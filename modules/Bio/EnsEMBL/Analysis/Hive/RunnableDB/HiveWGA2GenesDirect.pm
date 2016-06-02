#!/usr/bin/env perl

=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::HiveWGA2GenesDirect -

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect
    ->new (-db         => $pipelinedb,
           -input_id   => $input_id
           -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to stdout


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect;

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

#use Bio::EnsEMBL::Analysis::Config::General;
#use Bio::EnsEMBL::Analysis::Config::WGA2GenesDirect;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::ClusterFilter;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils
    qw(replace_stops_with_introns);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my($self) = @_;

  $self->hive_set_config();

  my $input_id = $self->param('iid');

  my $max_internal_stops = $self->param('max_internal_stops');
  unless(defined($max_internal_stops)) {
    $self->warning("No max_internal_stops param found, defaulting to 1");
    $max_internal_stops = 1;
  }
  $self->max_internal_stops($max_internal_stops);

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->QUERY_CORE_DNA_DB);
  my $target_dna_dba = $self->hrdb_get_dba($self->TARGET_CORE_DNA_DB);
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');
  my $source_dna_dbc = $self->hrdb_get_con('source_dna_db');
  my $target_dna_dbc = $self->hrdb_get_con('target_dna_db');

  # Define the source transcript and target transcript dbs
  my $source_transcript_dba = $self->hrdb_get_dba($self->QUERY_CORE_DB,undef,'source_dna_db');
  my $target_transcript_dba = $self->hrdb_get_dba($self->TARGET_CORE_DB,undef,'target_dna_db');
  $self->hrdb_set_con($source_transcript_dba,'source_transcript_db');
  $self->hrdb_set_con($target_transcript_dba,'target_transcript_db');
  my $source_transcript_dbc = $self->hrdb_get_con('source_transcript_db');
  my $target_transcript_dbc = $self->hrdb_get_con('target_transcript_db');

  # Define the compara db
  my $compara_dba = $self->hrdb_get_dba($self->COMPARA_DB,'compara');
  $self->hrdb_set_con($compara_dba,'compara_db');
  my $compara_dbc = $self->hrdb_get_con('compara_db');

  # Get the genome db adpator
  my $gdb_adap = $compara_dbc->get_GenomeDBAdaptor;

  # Retrieve the production names for the query and target species
  my $query_species = $source_transcript_dbc->get_MetaContainerAdaptor->get_production_name();
  my $target_species = $target_transcript_dbc->get_MetaContainerAdaptor->get_production_name();


  my $q_gdb = $gdb_adap->fetch_by_core_DBAdaptor($source_transcript_dbc);
  my $t_gdb = $gdb_adap->fetch_by_core_DBAdaptor($target_transcript_dbc);


  ########
  # check that the default assembly for the query and target agrees
  # with that for the method_link_species_set GenomeDBs
  ########

  my $q_assembly = $q_gdb->assembly;
  my $t_assembly = $t_gdb->assembly;



  my ($q_assembly_version, $t_assembly_version);
  eval {
    $q_assembly_version = $source_transcript_dbc->get_CoordSystemAdaptor->fetch_by_name('toplevel',$q_gdb->assembly);
    $t_assembly_version = $target_transcript_dbc->get_CoordSystemAdaptor->fetch_by_name('toplevel',$t_gdb->assembly);
  };
  if($@) {
    throw("Had trouble fetching coord systems for ".$q_gdb->assembly . " and " .$t_gdb->assembly . " from core dbs:\n".$@);
  }


  ########
  # fetch the genes; need to work in the coordinate space of the
  # top-level slice to be consistent with compara
  ########
  my $transcript = $source_transcript_dbc->get_TranscriptAdaptor->fetch_by_dbID($input_id);
  $transcript->analysis($self->analysis);
  unless($transcript->biotype eq 'protein_coding') {
    $self->input_job->autoflow(0);
    $self->complete_early('Transcript does not have protein_coding biotype!');
  }

  my $gene = Bio::EnsEMBL::Gene->new();
  $gene->analysis($self->analysis);
  $gene->biotype('projection');
  $gene->add_Transcript($transcript);
  $gene->stable_id($transcript->stable_id);
  $gene->start($transcript->start);
  $gene->end($transcript->end);
  $gene->strand($transcript->strand);
  $gene->slice($transcript->slice);

  $self->_check_gene($gene);
  $self->gene($gene);

  #########
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag,
  # and all GenomicAlignBlocks
  #########
  my $mlss = $compara_dbc->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs($self->INPUT_METHOD_LINK_TYPE,
                                                                                                [$q_gdb, $t_gdb]);
  unless($mlss) {
    throw("No MethodLinkSpeciesSet for :\n" .$self->INPUT_METHOD_LINK_TYPE . "\n" .$query_species . "\n" .$target_species);
  }

  my $dnafrag = $compara_dbc->get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($q_gdb,$self->gene->slice->seq_region_name);

  my $gaba = $compara_dbc->get_GenomicAlignBlockAdaptor;


  my $gen_al_blocks = $gaba->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,$dnafrag,$gene->start,$gene->end);


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
      $self->target_slices->{$tga->dnafrag->name} = $target_transcript_dbc->get_SliceAdaptor->fetch_by_region('toplevel',$tga->dnafrag->name);
    }

    if ($block->reference_genomic_align->dnafrag_strand < 0) {
      $block->reverse_complement;
    }

    push @{$chains{$block->group_id}}, $block;
  }

  foreach my $chain_id (keys %chains) {
    push @chains, [
                   sort {
                     $a->reference_genomic_align->dnafrag_start <=> $b->reference_genomic_align->dnafrag_start;
                   } @{$chains{$chain_id}}
                  ];
  }
  $self->genomic_align_block_chains(\@chains);


}


sub run {
  my ($self) = @_;

  unless(scalar(@{$self->good_transcripts}) > 0) {
    warning("No transcripts in the good_transcripts list, so nothing to project");
    return 1;
  }

  my @res_tran;
  my @final_tran;
  print scalar(@{$self->genomic_align_block_chains}),"\n";
  foreach my $chain (@{$self->genomic_align_block_chains}) {

    #    my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffoldDirect->new(
    my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold->new(
                                                                                     -genomic_align_blocks => $chain,
                                                                                     -from_slice    => $self->gene->slice,
                                                                                     -to_slices     => $self->target_slices,
                                                                                     -transcripts   => $self->good_transcripts,
                                                                                     -max_readthrough_dist => $self->MAX_EXON_READTHROUGH_DIST,
                                                                                     -direct_target_coords => 1,
                                                                                   );

    foreach my $tran (@{$self->good_transcripts}) {
      # Remember to remove the 1 in place transcripts as right now this is only used for testing purposes
      #$gene_scaffold->place_transcript($tran,1);
      my $proj_trans;
      eval {
        $proj_trans  = $gene_scaffold->place_transcript($tran);
      };

      if($@) {
        $self->runnable_failed(1);
      }
      if ($proj_trans) {
        push @res_tran, $proj_trans;
      }
    }
  }


  if ($self->TRANSCRIPT_FILTER){
    @res_tran = @{$self->filter->filter_results(\@res_tran)};
  }

  foreach my $res_tran (@res_tran){ 
    $res_tran = $self->process_transcript($res_tran);

    unless($res_tran) {
      next;
    }

    if($self->param('calculate_coverage_and_pid') && $res_tran) {
      my $gene = $self->gene();
      my $transcripts = $gene->get_all_Transcripts;
      my $source_transcript = $$transcripts[0];
      if($res_tran->translation->length < 20000 && $source_transcript->translation->length < 20000) {
        $self->realign_translation($source_transcript,$res_tran);
      } else {
        $self->warning('Not realigning translation as translation length >= 20000');
      }
    }
    push @final_tran, $res_tran;
  }

  # create new gene object for each transcript
  print "At the end of RUN, we had ", scalar(@final_tran), " transcripts\n";
  $self->output(\@final_tran);
}


sub write_output {
  my ($self) = @_;

  my $trans_count = 0;
  my $target_transcript_dbc = $self->hrdb_get_con('target_transcript_db');
  my $t_gene_adaptor = $target_transcript_dbc->get_GeneAdaptor();
  my $failure_branch_code = -3;

  if($self->runnable_failed && scalar(@{$self->output}) == 0) {
    $self->warning("Issue with projection, will dataflow input id on branch -3");
    my $output_hash = {};
    $output_hash->{'iid'} = $self->param('iid');
    $self->dataflow_output_id($output_hash,$failure_branch_code);
  } else {
    foreach my $t (@{$self->output}) {
      $t->analysis($self->analysis);
      $t->biotype('projection');

      my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis,-biotype  => 'projection');

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

    print "For gene " . $self->param('iid'). " you stored ", $trans_count, " transcripts\n";
    # to do: write gene back to core target database
  }

  return 1;
}


sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  if($runnable_failed) {
    $self->param('_runnable_failed',$runnable_failed);
  }
  return($self->param('_runnable_failed'));
}


sub calculate_coverage_and_pid {
  my ($self, $value) = @_;
  if($value){
    $self->{_calculate_coverage_and_pid} = $value;
  }
  return $self->{_calculate_coverage_and_pid};
}

sub realign_translation {
  my ($self,$source_transcript,$projected_transcript) = @_;

  my $query_seq = $source_transcript->translate->seq();
#  my $transcripts = $gene->get_all_Transcripts();
#  my $transcript = ${$transcripts}[0];
  my $projected_seq = $projected_transcript->translate->seq();

  my $align_input_file = "/tmp/projected_align_".$$.".fa";
  my $align_output_file = "/tmp/projected_align_".$$.".aln";

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $query_seq;
  say INPUT ">target";
  say INPUT $projected_seq;
  close INPUT;

  my $align_program_path = 'muscle';

  my $cmd = $align_program_path." -in ".$align_input_file." -out ".$align_output_file;
  my $result = system($cmd);

  if($result) {
    throw("Got a non-zero exit code from alignment. Commandline used:\n".$cmd);
  }

  my $file = "";
  open(ALIGN,$align_output_file);
  while(<ALIGN>) {
    $file .= $_;
  }
  close ALIGN;

  unless($file =~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/) {
    throw("Could not parse the alignment file for the alignment sequences. Alignment file: ".$align_output_file);
  }

  my $aligned_query_seq = $1;
  my $aligned_projected_seq = $3;

  $aligned_query_seq =~ s/\n//g;
  $aligned_projected_seq =~ s/\n//g;

  `rm $align_input_file`;
  `rm $align_output_file`;

  # Work out coverage
  my $coverage;
  my $temp = $aligned_projected_seq;
  my $projected_gap_count = $temp =~ s/\-//g;
  my $ungapped_query_seq = $aligned_query_seq;
  $ungapped_query_seq  =~ s/\-//g;

  if(length($ungapped_query_seq) == 0) {
    $coverage = 0;
  } else {
    $coverage = 100 - (($projected_gap_count/length($ungapped_query_seq)) * 100);
  }

  # Work out precent identity
  my $match_count = 0;
  my $aligned_positions = 0;
  for(my $j=0; $j<length($aligned_query_seq); $j++) {
    my $char_query = substr($aligned_query_seq,$j,1);
    my $char_target = substr($aligned_projected_seq,$j,1);
    if($char_query eq '-' || $char_target  eq '-') {
      next;
    }
    if($char_query eq $char_target) {
      $match_count++;
    }
    $aligned_positions++;
  }

  unless($aligned_positions) {
    throw("Pairwise alignment between the query sequence and the translation shows zero aligned positions. Something has gone wrong");
  }

  my $percent_id = ($match_count / $aligned_positions) * 100;

  # Get all exons and transcript supporting features
  my $transcript_supporting_features = $projected_transcript->get_all_supporting_features();
  my $exons = $projected_transcript->get_all_Exons();

  # Now clean these out
  $projected_transcript->flush_Exons();
  $projected_transcript->flush_supporting_features();

  # Loop through the TSFs and add the coverage and pid, then add back into transcript
  foreach my $transcript_supporting_feature (@{$transcript_supporting_features}) {
    $transcript_supporting_feature->hcoverage($coverage);
    $transcript_supporting_feature->percent_id($percent_id);
    $transcript_supporting_feature->hseqname($source_transcript->stable_id);
    $projected_transcript->add_supporting_features($transcript_supporting_feature);
  }


  # Loop through exons, get supporting features for each, flush existing SF, add coverage and pid, add back to exon, add exon to transcript
  foreach my $exon (@{$exons}) {
    my $exon_supporting_features = $exon->get_all_supporting_features();
    $exon->flush_supporting_features();
    foreach my $exon_supporting_feature (@{$exon_supporting_features}) {
      $exon_supporting_feature->hcoverage($coverage);
      $exon_supporting_feature->percent_id($percent_id);
      $exon_supporting_feature->hseqname($source_transcript->stable_id);
      $exon->add_supporting_features($exon_supporting_feature);
    }
    $projected_transcript->add_Exon($exon);
  }

}

sub hive_set_config {
  my $self = shift;

  # Throw is these aren't present as they should both be defined
  unless($self->param_is_defined('logic_name') && $self->param_is_defined('module')) {
    throw("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, ".
          "even if they are already defined in the analysis hash itself. This is because the hive will not allow the runnableDB ".
          "to read values of the analysis hash unless they are in the parameters hash. However we need to have a logic name to ".
          "write the genes to and this should also include the module name even if it isn't strictly necessary"
         );
  }

  # Make an analysis object and set it, this will allow the module to write to the output db
  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                           );
  $self->analysis($analysis);

  # Now loop through all the keys in the parameters hash and set anything that can be set
  my $config_hash = $self->param('config_settings');
  foreach my $config_key (keys(%{$config_hash})) {
    if(defined &$config_key) {
      $self->$config_key($config_hash->{$config_key});
    } else {
      throw("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does ".
            "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending ".
            "key:\n".$config_key
           );
    }
  }

  my $logic = $self->analysis->logic_name;
  foreach my $var (qw(INPUT_METHOD_LINK_TYPE QUERY_CORE_DB TARGET_CORE_DB COMPARA_DB)) {
    unless($self->$var) {
      throw("You must define $var in config for logic '".$logic."' or in the DEFAULT entry");
    }
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

sub hrdb_set_con {
  my ($self,$dba,$dba_con_name) = @_;

  if($dba_con_name){
      $self->param('_'.$dba_con_name,$dba);
    } else {
      $self->param('_hrdbadaptor',$dba);
    }

}


sub hrdb_get_con {
  my ($self,$dba_con_name) = @_;

  if($dba_con_name) {
    return $self->param('_'.$dba_con_name);
  } else {
    return $self->param('_hrdbadaptor');
  }
}

sub hrdb_get_dba {
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

            my $dnadb = $self->hrdb_get_con($dna_db_name);

            # try to get default asm+ species name for OTHER db - does not work
            # for comapra database
            my $core_db_asm = $dba->get_MetaContainer->get_default_assembly();
            my $core_db_species = $dba->get_MetaContainer->get_common_name();

            # get the same for dna-db
            my $dna_db_asm = $dnadb->get_MetaContainer->get_default_assembly();
            my $dna_db_species = $dnadb->get_MetaContainer()->get_common_name();

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

    if ((length($t->translateable_seq) % 3) == 0){

      push(@good_transcripts,$t);

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
  my ($self,$tran) = @_;

  return 0 if not $tran;

  my (@processed_transcripts);

  my @exons = @{$tran->get_all_Exons};
  my ($tsf) = @{$tran->get_all_supporting_features};
  my $pep = $tran->translate->seq;


  if (CORE::length($pep) == 0) {
    logger_info("Rejecting proj of ".$tran->dbID()." because was just a stop codon");
    return 0;
  }

  my $num_stops = $pep =~ s/\*/\*/g;

  if($num_stops > $self->max_internal_stops()) {
    $self->warning("The internal stop count (".$num_stops.") is greater than the allowed value: ".$self->max_internal_stops());
    return 0;
  }

  ##################
  # number of stops is non-zero but acceptable. Need to 
  # operate on the transcript to jump over the stops
  if($num_stops) {
    $tran = replace_stops_with_introns($tran,$self->max_internal_stops());
  }

  return $tran;
}



###########################
# gets/sets
###########################

sub genomic_align_block_chains {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_gen_al_chains',$val);
  }

  return $self->param('_gen_al_chains');
}


sub gene {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_gene_recs',$val);
  }

  return $self->param('_gene_recs');
}

sub good_transcripts {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_good_transcripts',$val);
  }

  return $self->param('_good_transcripts');
}

sub target_slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_target_slices',$val);
  }

  unless($self->param_is_defined('_target_slices')) {
    $self->param('_target_slices',{});
  }

  return $self->param('_target_slices');
}

sub max_internal_stops {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_max_internal_stops',$val);
  }

  return $self->param('_max_internal_stops');
}

####################################
# config variable holders
####################################

#sub read_and_check_config {
#  my ($self, $hash) = @_;

#  $self->SUPER::read_and_check_config($hash);

#  my $logic = $self->analysis->logic_name;

#  foreach my $var (qw(INPUT_METHOD_LINK_TYPE
#                      QUERY_CORE_DB
#                      TARGET_CORE_DB
#                      COMPARA_DB)) {

#    throw("You must define $var in config for logic '$logic'" .
#          " or in the DEFAULT entry")
#        if not $self->$var;
#  }

 # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
#  if ($self->TRANSCRIPT_FILTER) {
#    if (not ref($self->TRANSCRIPT_FILTER) eq "HASH" or
#        not exists($self->TRANSCRIPT_FILTER->{OBJECT}) or
#        not exists($self->TRANSCRIPT_FILTER->{PARAMETERS})) {

#      throw("FILTER in config foR '$logic' must be a hash ref with elements:\n" .
#            "  OBJECT : qualified name of the filter module;\n" .
#            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
#    } else {
#      my $module = $self->TRANSCRIPT_FILTER->{OBJECT};
#      my $pars   = $self->TRANSCRIPT_FILTER->{PARAMETERS};

#      (my $class = $module) =~ s/::/\//g;
#      eval{
#        require "$class.pm";
#      };
#      throw("Couldn't require ".$class." Exonerate2Genes:require_module $@") if($@);

#      $self->filter($module->new(%{$pars}));
#    }
#  }




#}

#
# core options
#

sub INPUT_METHOD_LINK_TYPE {
  my ($self, $type) = @_;

  if (defined $type) {
    $self->param('_input_method_link_type',$type);
  }

  return $self->param('_input_method_link_type');
}


sub COMPARA_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->param('_compara_db',$db);
  }

  return $self->param('_compara_db');
}


sub QUERY_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->param('_query_core_db',$db);
  }

  return $self->param('_query_core_db');
}

sub QUERY_CORE_DNA_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->param('_query_core_dna_db',$db);
  }

  return $self->param('_query_core_dna_db');
}

sub TARGET_CORE_DB {
  my ($self,$db) = @_;

  if (defined $db) {
    $self->param('_target_core_db',$db);
  }

  return $self->param('_target_core_db');
}

sub TARGET_CORE_DNA_DB {
  my ($self,$db) = @_;

  if (defined $db) {
    $self->param('_target_core_dna_db',$db);
  }

  return $self->param('_target_core_dna_db');
}


#
# transcript editing and filtering
#


sub TRANSCRIPT_FILTER {
   my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_transcript_filter',$val);
  }

  return $self->param('_transcript_filter');
}


sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_filter',$val);
  }

  return $self->param('_filter');
}


sub MAX_EXON_READTHROUGH_DIST {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_max_ex_rt_dist',$val);
  }

  return $self->param('_max_ex_rt_dist');
}


sub MIN_COVERAGE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_min_coverage',$val);
  }

  return $self->param('_min_coverage');
}

1;
