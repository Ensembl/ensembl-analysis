=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion;

use strict;
use warnings;

use Bio::Seq;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw ( create_file_name write_seqfile);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my( $self) = @_;


  ##########################################
  # set up the target (genome)
  ##########################################

  $self->create_analysis;
  my ($slice_name, $accession) = $self->input_id =~ /^(.*):+([^:]+)$/;
  $self->query_acc($accession);
  #repeat masking logic names
  $self->throw("Repeat logic names are not in an array") if(!(ref($self->SOFT_MASKED_REPEATS) eq "ARRAY"));
  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $dba->dnadb($dna_dba);
  }
  $self->hrdb_set_con($dba,'target_db');

  foreach my $repeat_logic_name ( @{ $self->SOFT_MASKED_REPEATS } ) {
    my $repeat_analysis =
      $dba->get_AnalysisAdaptor->fetch_by_logic_name($repeat_logic_name);
    if ( !$repeat_analysis ) {
      $self->throw(   "Failed to find an analysis with the logic_name "
             . $repeat_logic_name . ". Cannot continue." );
    }
  }

  my $slice = $dba->get_SliceAdaptor->fetch_by_name($slice_name)->get_repeatmasked_seq($self->SOFT_MASKED_REPEATS,1);
  my $genomic_seq = Bio::Seq->new(-ID => $slice_name, -SEQ => $slice->seq);

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################
  my $query_file;
  # check if QUERYSEQ dir exists and file exists
  if($self->param_is_defined('iid_type') and $self->param('iid_type') eq 'db_seq') {
    $query_file = $self->output_query_file();
  }
  elsif(defined $self->QUERYSEQS){
    if(-d $self->QUERYSEQS && -e $self->QUERYSEQS.'/'.$self->query_acc){
      # read use this file as seq file
      print "Using existing file\n";
      $query_file = $self->QUERYSEQS.'/'.$self->query_acc;
    }
  }
  #retrieve and write seq
  else{
    my $seqfetcher = $self->seqfetcher;
    my $query_seq  = $seqfetcher->get_Seq_by_acc( $self->query_acc );
    if(!$query_seq){
      $self->throw("No entry in sequence index for ".$self->query_acc."\n");
    }
    $query_file    = create_file_name( "query_", "fa", $self->ANALYSIS_WORK_DIR );

    write_seqfile($query_seq, $query_file, 'fasta');
  }
  ##########################################
  # Annotation file with CDS positions
  ##########################################

  ##########################################
  # setup the runnables
  ##########################################

  my %parameters = %{$self->parameters_hash};
  if (not exists($parameters{-options}) and
      defined $self->OPTIONS) {
    $parameters{-options} = $self->OPTIONS
  }
  if (not exists($parameters{-coverage_by_aligned}) and
      defined $self->COVERAGE_BY_ALIGNED) {
    $parameters{-coverage_by_aligned} = $self->COVERAGE_BY_ALIGNED;
  }

  if (defined $self->PROGRAM && defined $self->analysis->program_file) {
    if ($self->PROGRAM ne $self->analysis->program_file) {
      $self->throw("CONFLICT: You have defined -program in your config file and ".
            "-program_file in your analysis table.");
    }
  }
  if ($parameters{-options} =~ /cdna2genome/ and !$self->QUERYANNOTATION) {
    $self->param('tempfile', File::Temp->new());
    my $tempfile = $self->param('tempfile');
    $self->QUERYANNOTATION($tempfile->filename);
    foreach my $seq (@$query_file) {
      my $best_start = 0;
      my $best_length = 0;
      my $bestM_start = 0;
      my $bestM_length = 0;
      my $start = 1;
      my $length = 0;
      foreach my $frame (0, 1, 2) {
        $start = $frame+1;
        foreach my $sequence (split('\*', $seq->translate(undef, undef, $frame)->seq)) {
          next unless ($sequence);
          $length = length($sequence);
          if ($best_length < $length) {
            $best_length = $length;
            $best_start = $start;
          }
          if ($length > $bestM_length) {
            if ($sequence =~ /M/) {
              my $pos = $-[0];
              if ($bestM_length < $length-$pos) {
                $bestM_length = $length-$pos;
                $bestM_start = $start+($pos)*3;
                if ($best_length < $bestM_length) {
                  $best_length = $bestM_length;
                  $best_start = $bestM_start;
                }
              }
            }
          }
          $start += ($length+1)*3;
        }
      }
      if ($bestM_length > $best_length*0.8) {
        print $tempfile $seq->id, ' + ', $bestM_start, ' ', (($bestM_length+1)*3), "\n";
      }
      else {
        print $tempfile $seq->id, ' + ', $best_start, ' ', (($best_length+1)*3), "\n";
      }
    }
  }

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
              -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_seqs => [$genomic_seq],
              -query_type     => $self->QUERYTYPE,
              -annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
              %parameters,
              );
    if ($self->debug) {
      $runnable->_verbose(1);
    }
    if (ref($query_file) eq 'ARRAY') {
      $runnable->query_seqs($query_file);
    }
    else {
      $runnable->query_file($query_file);
    }
    $self->runnable($runnable);

}

sub output_query_file {
  my ($self) = @_;

  my $accession_array = [$self->query_acc];

  my $table_adaptor = $self->db->get_NakedTableAdaptor();

  # table name here should probably be changed to something more general
  $table_adaptor->table_name('cdna_sequences');


#  my $output_dir = $self->param('query_seq_dir');
#
#  # Note as each accession will occur in only one file, there should be no problem using the first one
#  my $outfile_name = "exonerate_".${$accession_array}[0].".fasta";
#  my $outfile_path = $output_dir."/".$outfile_name;

  my $biotypes_hash = {};

#  unless(-e $output_dir) {
#    `mkdir $output_dir`;
#  }
#
#  if(-e $outfile_path) {
#    $self->warning("Found the query file in the query dir already. Overwriting. File path:\n".$outfile_path);
#  }

#  open(QUERY_OUT,">".$outfile_path);
  my @query_sequences;
  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw("Did not find an entry in the cdna_sequences table matching the accession. Accession:\n".$accession);
    }

    my $seq = $db_row->{'seq'};
    $biotypes_hash->{$accession} = $db_row->{'biotype'};

#    my $record = ">".$accession."\n".$seq;
    push(@query_sequences, Bio::Seq->new(-id => $accession, -seq => $seq));

#    say QUERY_OUT $record;
  }
#  close QUERY_OUT;

  #$self->files_to_delete($outfile_path);
#  $self->get_biotype($biotypes_hash);

#  return($outfile_path);
  return \@query_sequences;
}
############################################################

sub run {

  my ($self) = @_;
  my @results;

  $self->throw("Can't run - no runnable objects") unless ( $self->runnable );

  foreach my $runnable ( @{ $self->runnable } ) {

    $runnable->run;
    push( @results, @{ $runnable->output } );
  }

  if ( $self->filter ) {
    my $filtered_transcripts = $self->filter->filter_results( \@results );
    @results = @$filtered_transcripts;
  }

  my @genes = $self->make_genes(@results);

  $self->output( \@genes );
}


############################################################

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb        = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $outdb->get_GeneAdaptor;

  my $fails = 0;
  my $total = 0;
  foreach my $gene (@{$self->output}) {

    eval { $gene_adaptor->store($gene); };
    if ($@) {
      $self->warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ( $fails > 0 ) {
    $self->throw(   "Not all genes could be written successfully "
           . "($fails fails out of $total)" );
  }
  $self->input_job->autoflow(0) unless ($total);
} ## end sub write_output

############################################################

sub make_genes {
  my ( $self, @transcripts ) = @_;

  my (@genes);

  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;

  my %genome_slices;

  foreach my $tran (@transcripts) {
    $tran->analysis( $self->analysis );
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis( $self->analysis );
    $gene->biotype( $self->analysis->logic_name );

    ############################################################
    # put a slice on the transcript

    my $slice_id = $tran->start_Exon->seqname;
    if ( not exists $genome_slices{$slice_id} ) {
      # assumes genome seqs were named in the Ensembl API Slice naming
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};

    foreach my $exon ( @{ $tran->get_all_Exons } ) {
      $exon->slice($slice);
      foreach my $evi ( @{ $exon->get_all_supporting_features } ) {
        $evi->slice($slice);
        $evi->analysis( $self->analysis );
      }
    }
    foreach my $evi ( @{ $tran->get_all_supporting_features } ) {
      $evi->slice($slice);
      $evi->analysis( $self->analysis );
    }

    if ( !$slice ) {
      my ($sf);

      if ( @{ $tran->get_all_supporting_features } ) {
        ($sf) = @{ $tran->get_all_supporting_features };
      } else {
        my @exons = @{ $tran->get_all_Exons };
        ($sf) = @{ $exons[0]->get_all_supporting_features };
      }
      print $sf->hseqname . "\t$slice_id\n";
    }

    $self->throw("Have no slice") if ( !$slice );
    $tran->slice($slice);
    $tran->biotype($gene->biotype);
    $gene->add_Transcript($tran);
    push( @genes, $gene );
  } ## end foreach my $tran (@transcripts)
  return @genes;
} ## end sub make_genes

############################################################
############################################################

sub seqfetcher {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->throw(   "RunnableDB::Exonerate2GenesRegion " . $arg
           . " must have a method get_Seq_by_acc" )
      unless ( $arg->can("get_Seq_by_acc") );
    $self->param('seqfetcher', $arg);
  }
  if ( !$self->param_is_defined('seqfetcher') ) {
    $self->require_module( $self->SEQFETCHER_OBJECT );
    my %params = %{ $self->SEQFETCHER_PARAMS };
    print $params{-db}->[0], "\n";
    $self->param('seqfetcher', $self->SEQFETCHER_OBJECT->new( %params, ));
  }
  return $self->param('seqfetcher');
}


sub get_output_db {
  my ($self) = @_;

  my $outdb;

  if ($self->OUTDB) {
    if ( ref($self->OUTDB)=~m/HASH/) {

      $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$self->OUTDB},
                                                -dnadb => $self->hrdb_get_con('dna_db'));
    }else{
      $outdb = $self->get_dbadaptor($self->OUTDB);
    }
  } else {
    $outdb = $self->hrdb_get_con('target_db');
  }
  return $outdb;
}



############################################################
#
# get/set methods
#
############################################################


############################################################

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_transcript_filter',$val);
  }
  elsif (!$self->param_is_defined('_transcript_filter')
    and $self->param_is_defined('FILTER')
    and exists $self->param('FILTER')->{OBJECT}) {
    my $module = $self->require_module($self->param('FILTER')->{OBJECT});
    $self->param('_transcript_filter', $module->new(%{$self->param('FILTER')->{PARAMETERS}}));
  }
  if ($self->param_is_defined('_transcript_filter')) {
    return $self->param('_transcript_filter');
  }
  else {
    return;
  }
}


#############################################################
# Declare and set up config variables
#############################################################

sub QUERYSEQS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('QUERYSEQS', $value);
  }

  if ($self->param_is_defined('QUERYSEQS')) {
    return $self->param('QUERYSEQS');
  } else {
    return;
  }
}

sub QUERYTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('QUERYTYPE', $value);
  }

  if ($self->param_is_defined('QUERYTYPE')) {
    return $self->param('QUERYTYPE');
  } else {
    return;
  }
}


sub QUERYANNOTATION {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('annotation_file', $value);
  }

  if ($self->param_is_defined('annotation_file')) {
    return $self->param('annotation_file');
  } else {
    return;
  }
}



sub GENOMICSEQS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('GENOMICSEQS', $value);
  }

  if ($self->param_is_defined('GENOMICSEQS')) {
    return $self->param('GENOMICSEQS');
  } else {
    return;
  }
}

sub query_acc {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->param('query_acc', $value);
  }

  if ($self->param_is_defined('query_acc')) {
    return $self->param('query_acc');
  } else {
    return;
  }
}


sub IIDREGEXP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('IIDREGEXP', $value);
  }

  if ($self->param_is_defined('IIDREGEXP')) {
    return $self->param('IIDREGEXP');
  } else {
    return;
  }
}

sub OUTDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('target_db', $value);
  }

  if ($self->param_is_defined('target_db')) {
    return $self->param('target_db');

  } else {
    return;
  }
}

sub COVERAGE_BY_ALIGNED {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('COVERAGE', $value);
  }

  if ($self->param_is_defined('COVERAGE')) {
    return $self->param('COVERAGE');
  } else {
    return;
  }
}


sub FILTER {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('FILTER', $value);
  }

  if ($self->param_is_defined('FILTER')) {
    return $self->param('FILTER');
  } else {
    return;
  }
}

sub OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('OPTIONS', $value);
  }

  if ($self->param_is_defined('OPTIONS')) {
    return $self->param('OPTIONS');
  } else {
    return;
  }
}

sub NONREF_REGIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('NONREF_REGIONS', $value);
  }

  if ($self->param_is_defined('NONREF_REGIONS')) {
    return $self->param('NONREF_REGIONS');
  } else {
    return;
  }
}

sub PROGRAM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('PROGRAM', $value);
  }

  if ($self->param_is_defined('PROGRAM')) {
    return $self->param('PROGRAM');
  } else {
    return;
  }
}

sub USE_KILL_LIST {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('USE_KILL_LIST', $value);
  }

  if ($self->param_is_defined('USE_KILL_LIST')) {
    return $self->param('USE_KILL_LIST');
  } else {
    return;
  }
}

sub KILL_TYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('KILL_TYPE', $value);
  }

  if ($self->param_is_defined('KILL_TYPE')) {
    return $self->param('KILL_TYPE');
  } else {
    return;
  }
}

sub SEQFETCHER_OBJECT {
  my ($self, $value) = @_;
  if($value){
    $self->param('SEQFETCHER_OBJECT', $value);
  }
  return $self->param('SEQFETCHER_OBJECT');
}

sub SEQFETCHER_PARAMS {
  my ($self, $value) = @_;
  if($value){
    $self->param('SEQFETCHER_PARAMS', $value);
  }
  return $self->param('SEQFETCHER_PARAMS');
}

sub SOFT_MASKED_REPEATS {
  my ($self, $value) = @_;
  if($value){
    $self->param('SOFT_MASKED_REPEATS', $value);
  }
  return $self->param('SOFT_MASKED_REPEATS');
}

sub ANALYSIS_WORK_DIR {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_ana_dir', $val);
  }
  return $self->param('_ana_dir');
}


###############################################
###     end of config
###############################################


1;
