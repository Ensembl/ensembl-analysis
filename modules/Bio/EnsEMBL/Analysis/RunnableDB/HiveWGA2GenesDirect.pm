
package Bio::EnsEMBL::Analysis::RunnableDB::HiveWGA2GenesDirect;

use warnings ;
use strict;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::WGA2GenesDirect;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::ClusterFilter;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils
    qw(replace_stops_with_introns);

use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable',
            'Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild');

sub new {
  my ($class,@args) = @_;


  my $self = $class->SUPER::new(@args);
  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => "ProjectTranscripts",
                                             -module => "HiveWGA2GenesDirect",
                                            );
  $self->analysis($analysis);
  $self->read_and_check_config($WGA2GENES_CONFIG_BY_LOGIC);

  return $self;
}

sub fetch_input {
  my( $self) = @_;

  unless($self->input_id =~ /\"iid\" \=\> \"([a-zA-Z\d]+)\"/) {
    throw("Could not match an Ensembl stable id from the input id:\n".
          $self->input_id
         );
  }

  my $input_id = $1;

  print "YOUR INPUT ID:",$input_id,"\n";

  my $source_dna_db = $self->db('source_dna_db',$self->get_dba($self->param('source_dna_db')));
  my $target_dna_db = $self->db('target_dna_db',$self->get_dba($self->param('target_dna_db')));

  my $source_transcript_db = $self->db('source_transcript_db',$self->get_dba($self->param('source_transcript_db'),undef,'source_dna_db'));

  my $target_transcript_db = $self->db('target_transcript_db',$self->get_dba($self->param('target_transcript_db'),undef,'target_dna_db'));
  my $compara_db = $self->db('compara_db',$self->get_dba($self->param('compara_db'),'compara'));

  my $query_species = 
      $source_transcript_db->get_MetaContainerAdaptor->get_Species->binomial;
  my $target_species = 
      $target_transcript_db->get_MetaContainerAdaptor->get_Species->binomial;

  my $gdb_adap = $compara_db->get_GenomeDBAdaptor;
  my $q_gdb = $gdb_adap->fetch_by_name_assembly($query_species);
  my $t_gdb = $gdb_adap->fetch_by_name_assembly($target_species);

  ########
  # check that the default assembly for the query and target agrees
  # with that for the method_link_species_set GenomeDBs
  ########

  my ($q_assembly_version, $t_assembly_version);
  eval {
    $q_assembly_version = $source_transcript_db->get_CoordSystemAdaptor->
        fetch_by_name('toplevel',
                      $q_gdb->assembly);

    $t_assembly_version = $target_transcript_db->get_CoordSystemAdaptor->
        fetch_by_name('toplevel',
                      $t_gdb->assembly);
  };
  $@ and do {
    throw("Had trouble fetching coord systems for ".
          $q_gdb->assembly . " and " .
          $t_gdb->assembly . " from core dbs: $@");
  };


  ########
  # fetch the genes; need to work in the coordinate space of the
  # top-level slice to be consistent with compara
  ########
  my $gene = $source_transcript_db->get_GeneAdaptor->fetch_by_stable_id($input_id);
  $self->_check_gene($gene);
  $self->gene($gene);


  #########
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag,
  # and all GenomicAlignBlocks
  #########
  my $mlss = $compara_db->get_MethodLinkSpeciesSetAdaptor
      ->fetch_by_method_link_type_GenomeDBs($self->INPUT_METHOD_LINK_TYPE,
                                            [$q_gdb, $t_gdb]);
  throw("No MethodLinkSpeciesSet for :\n" .
        $self->INPUT_METHOD_LINK_TYPE . "\n" .
        $query_species . "\n" .
        $target_species)
      if not $mlss;

  my $dnafrag = $compara_db->get_DnaFragAdaptor->
      fetch_by_GenomeDB_and_name($q_gdb,
                                 $self->gene->slice->seq_region_name);

  my $gaba = $compara_db->get_GenomicAlignBlockAdaptor;

  my $gen_al_blocks = 
      $gaba->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,
                                                       $dnafrag,
                                                       $gene->start,
                                                       $gene->end);

  my (%chains, @chains);
  foreach my $block (@$gen_al_blocks) {
    my $qga = $block->reference_genomic_align;    
    my ($tga) = @{$block->get_all_non_reference_genomic_aligns};

    ###########################################################
    ###### INVESTIGATE: do we want to just use level 1 chains?
    ###########################################################
    #next if $qga->level_id != 1 or $tga->level_id != 1;
    
    # fetch the target slice for later reference
    if (not exists $self->target_slices->{$tga->dnafrag->name}) {

      $self->target_slices->{$tga->dnafrag->name} = 
          $target_transcript_db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                    $tga->dnafrag->name);
    }

    if ($block->reference_genomic_align->dnafrag_strand < 0) {
      $block->reverse_complement;
    }

    push @{$chains{$block->group_id}}, $block;
  }

  foreach my $chain_id (keys %chains) {
    push @chains, [
                   sort {
                     $a->reference_genomic_align->dnafrag_start <=> 
                     $b->reference_genomic_align->dnafrag_start;
                   } @{$chains{$chain_id}}
                   ];
  }
  $self->genomic_align_block_chains(\@chains);

}

sub run {
  my ($self) = @_;


  my @res_tran;
  my $tran_stable_id;
  my @final_tran;
  print scalar(@{$self->genomic_align_block_chains}),"\n";
  foreach my $chain (@{$self->genomic_align_block_chains}) {

    my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold->new(
                                                                                    -genomic_align_blocks => $chain,
                                                                                    -from_slice    => $self->gene->slice,
                                                                                    -to_slices     => $self->target_slices,
                                                                                    -transcripts   => $self->good_transcripts,
                                                                                    -max_readthrough_dist => $self->MAX_EXON_READTHROUGH_DIST,
                                                                                    -direct_target_coords => 1,
                                                                                    -ignore_translation => 1,
                                                                                    );
    foreach my $tran (@{$self->good_transcripts}) {
      $tran_stable_id = $tran->stable_id;
      my $proj_trans = $gene_scaffold->place_transcript($tran);

      if ($proj_trans) {
        push @res_tran, $proj_trans;
      }
    }
  }

  if ($self->TRANSCRIPT_FILTER){
    @res_tran = @{$self->filter->filter_results(\@res_tran)};
  }

  foreach my $res_tran (@res_tran){ 
    $res_tran = $self->process_transcript($res_tran,
                                          $tran_stable_id);
    push @final_tran, $res_tran;
  }
  # create new gene object for each transcript

  print "At the end of RUN, we had ", scalar(@final_tran), " transcripts\n";
  $self->output(\@final_tran);
}


sub write_output {
  my ($self) = @_;

  my $trans_count = 0;

  my $target_transcript_db = $self->db('target_transcript_db');
  my $t_gene_adaptor = $target_transcript_db->get_GeneAdaptor();

 foreach my $t (@{$self->output}) {
    $t->analysis($self->analysis);

#    unless(scalar(@{$self->output}) == 1) {
#      throw("FM2 exception WGA write: gene contains more than one transcript!!!!!");
#    }
    print "FM2 transcript stable id: ".$t->stable_id()."\n";
    my $gene = Bio::EnsEMBL::Gene->new( -analysis => $self->analysis,
                                        -biotype  => 'protein_coding',
                                        -stable_id => $t->stable_id(),
                                      );

    foreach my $tsf ( @{ $t->get_all_supporting_features }){
      $tsf->analysis($self->analysis);
    }

    foreach my $exon (@{$t->get_all_Exons()}){
      $exon->analysis($self->analysis);

      foreach my $esf (@{$exon->get_all_supporting_features()}){
        $esf->analysis($self->analysis);
      }
    }

    $gene->add_Transcript($t);
    empty_Gene($gene);
    $t_gene_adaptor->store($gene);
    $trans_count++;

    print "TRANSCRIPT:\n";
    foreach my $e (@{$t->get_all_Exons}) {
      printf("%s\tEnsembl\tExon\t%d\t%d\t%d\t%d\t%d\n", $e->slice->seq_region_name, $e->start, $e->end, $e->strand, $e->phase, $e->end_phase);
    }
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -fh => \*STDOUT);
    $seqio->write_seq($t->translate);
  }

  print "For gene " . $self->input_id . " you stored ", $trans_count, " transcripts\n";
  # to do: write gene back to core target database
}


sub db {
  my ($self,$adaptor_name,$value) = @_;

  if($value){
    $self->{$adaptor_name} = $value;
  }
  return $self->{$adaptor_name};
}

sub get_dba {
   my ($self,$connection_info, $non_standard_db_adaptor, $dna_db_name) = @_;
   my $dba;

   if(defined $non_standard_db_adaptor) {
     if($non_standard_db_adaptor eq 'compara') {
       $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
                                                            %$connection_info
                                                          );
       print "Not attaching a dna db to: ".$dba->dbname."\n";
     }
   }

   else {
          $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                     %$connection_info
                                                   );

          if($dna_db_name) {
            my $dnadb = $self->db($dna_db_name);

            # try to get default asm+ species name for OTHER db - does not work
            # for comapra database
            my $core_db_asm = $dba->get_MetaContainer->get_default_assembly();
            my $core_db_species =
            $dba->get_MetaContainer->get_common_name();

            # get the same for dna-db
            my $dna_db_asm = $dnadb->get_MetaContainer->get_default_assembly();
            my $dna_db_species =
            $dnadb->get_MetaContainer()->get_common_name();

            my $dna_db_and_core_db_are_compatible = 1;

            unless ( $core_db_asm eq $dna_db_asm ) {
            warning( "You try to add  a DNA_DB with assembly $dna_db_asm to "
              . "a core/cdna/otherfeatures DB with assembly $core_db_asm ...\n\t"
              . "that's incompatbile. I will not add dna_database "
              . $dnadb->dbname . " to core " . $dba->dbname . "\n" );

              $dna_db_and_core_db_are_compatible = 0;
            }

            unless ( $core_db_species eq $dna_db_species ) {
            warning( "You try to add a DNA_DB with species ".$dna_db_species." to "
                . "a core database with species: '" .$core_db_species . "' - this does not work. \n"
                . "Check that you are using the correct DNA_DB and that the species.common_name values in the meta tables match\n"
            );
            $dna_db_and_core_db_are_compatible = 0;

            }

            if ($dna_db_and_core_db_are_compatible) {
            $dba->dnadb($dnadb);
            print "\nAttaching DNA_DB "
              . $dnadb->dbname . " to "
              . $dba->dbname . "\n";
            }
          }

          else {
            print "Not attaching a dna db to: ".$dba->dbname."\n";
          }
   }

  $dba->dbc->disconnect_when_inactive(1) ;
  return $dba;

 }

######################################
# internal methods
#####################################

sub _check_gene {
  my ($self, $gene) = @_; 

  my @good_transcripts;  

  foreach my $t (@{$gene->get_all_Transcripts}){
    #print "Transcript Start: ",$t->start, "  END: ", $t->end,"\n";
    #print "translation: ",$t->translateable_seq,"\n";

    print "FM2 in _check_gene transcript: ".$t->stable_id.":srname:".$t->seq_region_name.":srs:".$t->seq_region_start.":sre:".$t->seq_region_end."\n";
    if ((length($t->translateable_seq) % 3) == 0){

      push @good_transcripts, $t;

      }
    else{
      warn ("Gene ", $gene->display_id(), " contains no valid transcripts");
      #push @good_transcripts, $t;
      #print "Done\n";
      #exit(1);
    }
  }
  
  warn ("Gene ", $gene->display_id(), " contains no valid transcripts") if (scalar(@good_transcripts) == 0); 

  $self->good_transcripts(\@good_transcripts);
  
  # throw an exception here if:
  # 1. gene is not protein-coding (translations in all transcripts)
  # 2. gene contains at least one transcript with a coding region
  #    that is not multiple-of-three in length

}


#################################################################
# FUNCTION  : process_transcript
#
# Description:
#    Subjects the given transcript to a few tests, returning 
#    the transcript if they succeed, undef if not. If the
#    transcript contains less than $max_stops stops, these
#    are "spliced out"; otherwise the transcripts is rejected
#################################################################

sub process_transcript {
  my ($self, 
      $tran,
      $source_id
) = @_;
  
  return 0 if not $tran;

  my (@processed_transcripts);

  my @exons = @{$tran->get_all_Exons};
  my ($tsf) = @{$tran->get_all_supporting_features};
  my $pep = $tran->translate->seq;


  if (CORE::length($pep) == 0) {
    logger_info("Rejecting proj of $source_id because was just a stop codon");
    return 0;
  }

  my $num_stops = $pep =~ s/\*/\*/g;

  ##################
  # number of stops is non-zero but acceptable. Need to 
  # operate on the transcript to jump over the stops
  if($num_stops) {
    $tran = replace_stops_with_introns($tran,10);
  }

  return $tran;
}



###########################
# gets/sets
###########################

sub genomic_align_block_chains {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gen_al_chains} = $val;
  }

  return $self->{_gen_al_chains};
}


sub gene {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gene_recs} = $val;
  }

  return $self->{_gene_recs};
}

sub good_transcripts {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_good_transcripts} = $val;
  }

  return $self->{_good_transcripts};
}

sub target_slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_target_slices} = $val;
  }
  
  if (not exists $self->{_target_slices}) {
    $self->{_target_slices} = {};
  }

  return $self->{_target_slices};
}



####################################
# config variable holders
####################################

sub read_and_check_config {
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);
 
  my $logic = $self->analysis->logic_name;

  foreach my $var (qw(INPUT_METHOD_LINK_TYPE
                      QUERY_CORE_DB
                      TARGET_CORE_DB
                      COMPARA_DB)) {

    throw("You must define $var in config for logic '$logic'" . 
          " or in the DEFAULT entry")
        if not $self->$var;
  }

 # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
  if ($self->TRANSCRIPT_FILTER) {
    if (not ref($self->TRANSCRIPT_FILTER) eq "HASH" or
        not exists($self->TRANSCRIPT_FILTER->{OBJECT}) or
        not exists($self->TRANSCRIPT_FILTER->{PARAMETERS})) {

      throw("FILTER in config foR '$logic' must be a hash ref with elements:\n" .
            "  OBJECT : qualified name of the filter module;\n" .
            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
    } else {
      my $module = $self->TRANSCRIPT_FILTER->{OBJECT};
      my $pars   = $self->TRANSCRIPT_FILTER->{PARAMETERS};

      (my $class = $module) =~ s/::/\//g;
      eval{
        require "$class.pm";
      };
      throw("Couldn't require ".$class." Exonerate2Genes:require_module $@") if($@);

      $self->filter($module->new(%{$pars}));
    }
  }




}

#
# core options
#

sub INPUT_METHOD_LINK_TYPE {
  my ($self, $type) = @_;

  if (defined $type) {
    $self->{_input_method_link_type} = $type;
  }

  return $self->{_input_method_link_type};
}


sub COMPARA_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_compara_db} = $db;
  }

  return $self->{_compara_db};
}


sub QUERY_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_query_core_db} = $db;
  }

  return $self->{_query_core_db};
}


sub TARGET_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_target_core_db} = $db;
  }

  return $self->{_target_core_db};
}


#
# transcript editing and filtering
#


sub TRANSCRIPT_FILTER {
   my ($self, $val) = @_;

  if (defined $val) {
    $self->{_transcript_filter} = $val; 
  }

  return $self->{_transcript_filter};
}

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->{_filter} = $val;
  }
  return $self->{_filter};
}



sub MAX_EXON_READTHROUGH_DIST {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_max_ex_rt_dist} = $val; 
  }

  return $self->{_max_ex_rt_dist};
}



sub MIN_COVERAGE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_min_coverage} = $val;
  }

  return $self->{_min_coverage};
}



1;
