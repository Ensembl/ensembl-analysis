=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise - 

=head1 SYNOPSIS

my $obj = Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise->
    new(
	-db        => $db,
        -input_id  => $id,
        -analysis  => $analysis, 
        );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

This module runs the runnable BlastMiniGenewise to use protein alignments to seed
genewise runs and then filter and write back the predictions to an ensembl core 
database

It is configured by the code

Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise and
Bio::EnsEMBL::Analysis::Config::Databases

and it takes 3 input id formates

The first is the standard slice name

coord_system:version:name:start:end:strand

The second includes a logic name and a protein id

coord_system:version:name:start:end:strand:logic_name:protein_id

Note this will only work if there is a feature of that type in the databae currently

Lastly it expects a couple of numbers at the end of the input id to allow it to take
a certain number of ids from the full list in a consistent manner

coord_system:version:name:start:end:strand:id_pool_bin:id_pool_index

The first indicates the bin size the ids were split up into
The second the index at which ids should be take

e.g 5:2

This would mean the 2nd of each group of 5 ids should be taken

=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenewise;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(convert_to_genes Transcript_info set_start_codon set_stop_codon list_evidence attach_Slice_to_Transcript);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info print_Gene attach_Analysis_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::KillList::KillList;
use Bio::Seq;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Function  : This processes the input id, fetches both the sequence and the
  align features to run with and filters the align features where appropriate on
  the basis of a kill list and comparision to existing gene models
  Returntype: 1
  Exceptions: 
  Example   : 

=cut


sub fetch_input{
  my ($self) = @_;

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($dba,'target_db');

  $self->create_analysis;

  my $iid_type = $self->param('iid_type');
  my $calculate_coverage_and_pid = $self->param('calculate_coverage_and_pid');

  my $slice;
  my $accession;
  my $peptide_seq;
  my $transcript_features;

 if($iid_type eq 'feature_region') {
    my $feature_region_id = $self->param('iid');
    $self->parse_feature_region_id($feature_region_id);
    $peptide_seq = $self->get_peptide_seq($self->accession());
 } elsif($iid_type eq 'feature_id') {
    my $feature_type = $self->param('feature_type');
    if($feature_type eq 'transcript') {
      my $transcript_id = $self->param('iid');
      my $transcript_dba = $self->hrdb_get_dba($self->param('transcript_db'));
      if($dna_dba) {
        $transcript_dba->dnadb($dna_dba);
      }

      $self->hrdb_set_con($transcript_dba,'transcript_db');
      $transcript_features = $self->get_transcript_features($transcript_id);
#      $transcript_features = $self->correct_features($transcript_features);

      $peptide_seq = $self->get_peptide_seq($self->accession());
    } else {
      $self->throw("The feature_type you passed in is not supported! Type:\n".$feature_type);
    }
  } else {
    $self->throw("The iid_type you passed in is not supported! Type:\n".$iid_type);
  }


  $accession = $self->accession();
  $slice = $self->query();

  my $use_killlist = $self->param('use_killlist');
  my $kill_query = 0;
  if($use_killlist) {
    # code that can set kill_query
  }

  # Some code for the PAF threshold if I want it (that would also set kill_query)

  # If the query is to be killed then return instead of making a runnable
  if($kill_query) {
    return;
  }

  my $peptide_obj = Bio::Seq->new( -display_id => $self->accession(),
                                   -seq => $peptide_seq);

  my %params;# = %{{                   ENDBIAS => 1,
             #                       MATRIX => 'BLOSUM80.bla',
             #                        GAP => 20,
             #                        EXTENSION => 8,
             #                        SPLICE_MODEL => 0
             #   }};


  if($iid_type eq 'feature_region') {
    my %params = %{$self->genewise_options} if($self->genewise_options);
    my $runnable = new Bio::EnsEMBL::Analysis::Runnable::Genewise
       (
         -query   => $slice,
         -protein  => $peptide_obj,
         -reverse  => 0,
         -analysis => $self->analysis,
         %params
     );
    $self->runnable($runnable);
  } elsif($iid_type eq 'feature_id') {
    my %params = %{$self->MINIGENEWISE_PARAMETERS} if($self->MINIGENEWISE_PARAMETERS);
    my $runnable  = Bio::EnsEMBL::Analysis::Runnable::MiniGenewise->new
       (
         -query            => $slice,
         -protein          => $peptide_obj,
         -features         => $transcript_features,
         -genewise_options => $self->GENEWISE_PARAMETERS,
         -analysis         => $self->analysis,
         -calculate_coverage_and_pid => $calculate_coverage_and_pid,
         -best_in_genome   => $self->best_in_genome_transcript(),
         %params,
       );
    $self->runnable($runnable);
  }
  return 1;
}



sub get_transcript_features {
  my ($self,$transcript_id) = @_;

  my $transcript_dba = $self->hrdb_get_con('transcript_db');
  my $transcript = $transcript_dba->get_TranscriptAdaptor()->fetch_by_dbID($transcript_id);
  my $tsf = $transcript->get_all_supporting_features();

  my $slice = $transcript->slice();
  $self->query($slice);

  my $feature_pair = ${$tsf}[0];
  my $accession = $feature_pair->hseqname();
  $self->accession($accession);

  if($self->param('use_genblast_best_in_genome')) {
    my $logic_name = $transcript->analysis->logic_name();
    if($logic_name =~ /_not_best$/) {
      $self->best_in_genome_transcript(0);
    } else {
      $self->best_in_genome_transcript(1);
    }
  }

  return($tsf);
}


sub best_in_genome_transcript {
   my ($self,$val) = @_;

   if(defined($val) && $val==1) {
     $self->param('_best_in_genome_transcript', 1);
   } elsif(defined($val) && $val==0) {
     $self->param('_best_in_genome_transcript', 0);
   }

   return($self->param('_best_in_genome_transcript'));
}


=head2 generate_ids_to_ignore

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Features
  Function  : to produce a list of ids which shouldnt be considered as their features
  overlap with existing genes
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut



sub generate_ids_to_ignore{
  my ($self, $features) = @_;
  my $exon_mask_list = $self->create_mask_list;
  my @exonmask_regions = @{$exon_mask_list};
  my %features;
  foreach my $f (@$features) {
    push @{$features{$f->hseqname}}, $f;
  }
  my %ids_to_ignore;
 SEQID:foreach my $sid (keys %features) {
    my $ex_idx = 0;
  FEAT: foreach my $f (sort {$a->start <=> $b->start} 
                       @{$features{$sid}}) {
      for( ; $ex_idx < @exonmask_regions; ) {
        my $mask_exon = $exonmask_regions[$ex_idx];
        if ($mask_exon->{'start'} > $f->end) {
          # no exons will overlap this feature
          next FEAT;
        } elsif ( $mask_exon->{'end'} >= $f->start) {
          # overlap
          $ids_to_ignore{$f->hseqname} = 1;
          next SEQID;
        }  else {
          $ex_idx++;
        }
      }
    }
  }
  return \%ids_to_ignore;
}



sub parse_feature_region_id {
  my ($self,$feature_region_id) = @_;

  my $dba = $self->hrdb_get_con('target_db');
  my $sa = $dba->get_SliceAdaptor();

  unless($feature_region_id =~ s/\:([^\:]+)$//) {
    $self->throw("Could not parse the accession from the feature region id. Expecting a normal slice id, with an extra colon ".
                 "followed by the accession. Offending feature_region_id:\n".$feature_region_id);
  }

  my $slice_name = $feature_region_id;
  my $accession = $1;

  my $slice = $sa->fetch_by_name($slice_name);
  $self->query($slice);
  $self->accession($accession);
}



sub get_peptide_seq {
  my ($self,$accession) = @_;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('uniprot_sequences');

  my $output_dir = $self->param('query_seq_dir');

  my $biotypes_hash = {};

  my $db_row = $table_adaptor->fetch_by_dbID($accession);
  unless($db_row) {
    $self->throw("Did not find an entry in the uniprot_sequences table matching the accession. Accession:\n".$accession);
  }

  my $peptide_seq = $db_row->{'seq'};
  $biotypes_hash->{$accession} = $db_row->{'biotype'};

  $self->get_biotype($biotypes_hash);

  return($peptide_seq);
}

sub get_biotype {
  my ($self,$biotype_hash) = @_;
  if($biotype_hash) {
    $self->param('_biotype_hash',$biotype_hash);
  }
  return($self->param('_biotype_hash'));
}



=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Function  : This runs the created runnables, attached unmasked slices to the
  transcripts, converts the transcripts into genes, sets start and stop codons if 
  possible then filters the results
  Returntype: arrayref of Bio::EnsEMBL::Gene objects
  Exceptions: throws if the runnable fails to run
  Example   : 

=cut



sub run{
  my ($self) = @_;
  my @genes;

  my $runnable = shift @{$self->runnable};
  my $output;
  eval{
    $runnable->run;
  };
  if($@){
    my $except = $@;
    $self->runnable_failed(1);
    $self->warning("Issue with running genewise, will dataflow input id on branch -3. Exception:\n".$except);
  } else {
    $output =  $runnable->output;
    my $transcript = ${$output}[0];

    unless($self->check_coverage($transcript)) {
      next;
    }

    unless($self->check_pid($transcript)) {
      next;
    }

    $transcript = $self->add_slice_to_features($transcript);
    $transcript = $self->process_transcript($transcript);
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);

    my $biotypes = $self->get_biotype();
    my $accession = $self->accession();
    my $transcript_biotype = $biotypes->{$accession};
    $transcript->analysis($self->analysis);
    $transcript->biotype($transcript_biotype);
    $gene->add_Transcript($transcript);
    push(@genes,$gene);
    $self->param('output_genes',\@genes);
  }

#  my $complete_transcripts = $self->process_transcripts(\@transcripts);
#  my @genes = @{convert_to_genes($complete_transcripts, $self->analysis, 
#                                 $self->OUTPUT_BIOTYPE)};
  #my $hash = $self->group_genes_by_id(\@genes);
  #foreach my $name(keys(%{$hash})){
  #  print "FILTER ".$name." has ".$hash->{$name}." genes before filter\n";
  #}
#  my @masked_genes;
#  logger_info("HAVE ".@genes." genes before masking") if($self->POST_GENEWISE_MASK);
#  if($self->POST_GENEWISE_MASK){
#    @masked_genes = @{$self->post_mask(\@genes)};
#  }else{
#    @masked_genes = @genes;
#  }
#  logger_info("HAVE ".@masked_genes." genesafter masking") if($self->POST_GENEWISE_MASK);
#  $self->filter_genes(\@masked_genes);
  #$hash = $self->group_genes_by_id($self->output);
  #foreach my $name(keys(%{$hash})){
  #  print "FILTER ".$name." has ".$hash->{$name}." genes after filter\n";
  #}
#  logger_info("HAVE ".@{$self->output}." genes to return");

  return 1;
}

sub check_coverage {
  my ($self,$transcript) = @_;

  my $transcript_supporting_features = $transcript->get_all_supporting_features();
  my $transcript_supporting_feature = $$transcript_supporting_features[0];
  my $coverage = $transcript_supporting_feature->hcoverage();
  my $coverage_cutoff = $self->param('genewise_cov');
  if(defined($coverage_cutoff) && $coverage < $coverage_cutoff) {
    return 0;
  } else {
    return 1;
  }

}

sub check_pid {
  my ($self,$transcript) = @_;

  my $transcript_supporting_features = $transcript->get_all_supporting_features();
  my $transcript_supporting_feature = $$transcript_supporting_features[0];
  my $pid = $transcript_supporting_feature->percent_id();
  my $pid_cutoff = $self->param('genewise_pid');
  if(defined($pid_cutoff) && $pid < $pid_cutoff) {
    return 0;
  } else {
    return 1;
  }

}


sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  if($runnable_failed) {
    $self->param('_runnable_failed',$runnable_failed);
  }
  return($self->param('_runnable_failed'));
}


sub add_slice_to_features {
  my ($self,$transcript) = @_;

    my $slice = $self->query();
    $transcript->slice($slice);

    my @exons = @{$transcript->get_all_Exons()};
    my @transcript_supporting_features = @{$transcript->get_all_supporting_features()};

    $transcript->flush_Exons;
    $transcript->flush_supporting_features();

    my $i=0;
    for($i=0; $i<scalar @transcript_supporting_features; $i++) {
      $transcript_supporting_features[$i]->slice($slice);
      $transcript->add_supporting_features($transcript_supporting_features[$i]);
    }

    for($i=0; $i<scalar(@exons); $i++) {
      $exons[$i]->slice($slice);
      my $exon_supporting_features = $exons[$i]->get_all_supporting_features;
      my @esf_with_slice;
      $exons[$i]->flush_supporting_features();
      foreach my $esf (@{$exon_supporting_features}) {
        $esf->slice($slice);
        push(@esf_with_slice,$esf);
      }
      $exons[$i]->add_supporting_features(@esf_with_slice);
      $transcript->add_Exon($exons[$i]);
    }

  return($transcript);

}


#sub correct_features {
#  my ($self,$features) = @_;

#  my $feature_string = "";
 # for(my $i=0; $i<scalar(@{$features})) {
#    $feature_string .= $$features[$i]->start.'..'.$$features[$i]->end.":";
#  }
#  chop $feature_string;
#  say "Starting FS: ".$feature_string;


#  my $feature_string = "1..10:13..100:101..150:500..550:552..554:";

#  my $feature_string_temp =  $feature_string;
#  while($feature_string_temp =~ s/\.\.(\d+)\:(\d+)\.\./\.\./) {
#    my $d1 = $1;
#    my $d2 = $2;
#    if($d2 - $d1 + 1 < 10) {
#      say "Removing frameshift from feature";
#      $feature_string =~ s/\.\.$d1\:$d2\.\./\.\./;
#    }
#    say "TS: ".$feature_string_temp;
#    say "FS: ".$feature_string;
#  }

  
#}

=head2 filter_genes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : arrayref of Genes
  Function  : filters the gene sets using the defined filter
  Returntype: n/a, stores the genes in the output array
  Exceptions: 
  Example   : 

=cut


sub filter_genes{
  my ($self, $genes) = @_;
  if($self->filter_object){
    my ($filtered_set, $rejected_set) = 
      $self->filter_object->filter_genes($genes);
    $self->output($filtered_set);
    $self->rejected_set($rejected_set);
  }else{
   $self->output($genes);
  }
}

=head2 process_transcripts

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Transcripts
  Function  : attaches unmasked slices, sets start and stop codons on transcripts
  Returntype: arrayref of Bio::EnsEMBL::Transcripts
  Exceptions: 
  Example   : 

=cut



sub process_transcript {
  my ($self, $transcript) = @_;
  my $start_t = set_start_codon($transcript);
  my $end_t = set_stop_codon($start_t);
  return $end_t;
}


=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Function  : writes the accepted and if desired rejected genes to the specified
  database
  Returntype: 1
  Exceptions: throws if failed to store any genes, otherwise warnings are issued
  Example   : 

=cut



sub write_output{
  my ($self) = @_;

  my $failure_branch_code = -3;
  if($self->runnable_failed == 1) {
    my $output_hash = {};
    $output_hash->{'iid'} = $self->param('iid');
    $self->dataflow_output_id($output_hash,$failure_branch_code);
  } else {
    my $ga = $self->get_adaptor;
    my @output = @{$self->param('output_genes')};
    foreach my $gene (@output) {
      empty_Gene($gene);
      eval {
        $ga->store($gene);
      };
      if($@){
        $self->throw("Failed to write gene: ".$@);
      }
    }
  }

  return 1;
}



=head2 exon_mask_list

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : various, String, Int, Arrayref
  Function  : various accessor methods for holding variables during the
  run. Some will create and return the approprite object if not already holding
  it
  Returntype: as given
  Exceptions: some thrown if given the wrong type
  Example   : 

=cut


sub exon_mask_list{
  my ($self, $arg) = @_;
  if(!$self->param_is_defined('exon_mask_list')){
    $self->param('exon_mask_list', []);
  }
  if($arg){
    $self->param('exon_mask_list', $arg);
  }
  return $self->param('exon_mask_list');
}

=head2 gene_mask_list

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : various, String, Int, Arrayref
  Function  : various accessor methods for holding variables during the
  run. Some will create and return the approprite object if not already holding
  it
  Returntype: as given
  Exceptions: some thrown if given the wrong type
  Example   : 

=cut
sub gene_mask_list{
  my ($self, $arg) = @_;
  if(!$self->param_is_defined('gene_mask_list')){
    $self->param('gene_mask_list', []);
  }
  if($arg){
    $self->param('gene_mask_list', $arg);
  }
  return $self->param('gene_mask_list');
}


sub id_pool_bin{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('id_pool_bin', $arg);
  }
  return $self->param('id_pool_bin');
}

sub id_pool_index{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('id_pool_index', $arg);
  }
  return $self->param('id_pool_index');
}

sub use_id{
  my ($self, $id) = @_;
  if($id){
    $self->param('use_id', $id);
  }
  return $self->param('use_id');
}


sub paf_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->param('paf_source_db', $db);
  }
  if(!$self->param_is_defined('paf_source_db')){
    my $db = $self->get_dbadaptor($self->PAF_SOURCE_DB); 
    $self->param('paf_source_db', $db);
  }
  return $self->param('paf_source_db');
}

sub gene_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->param('gene_source_db', $db);
  }
  if(!$self->param_is_defined('gene_source_db')){
    if (@{$self->BIOTYPES_TO_MASK}) {
      my $db = $self->get_dbadaptor($self->GENE_SOURCE_DB);
      $self->param('gene_source_db', $db);
    }
  }
  return $self->param('gene_source_db');
}

sub output_db{
  my ($self, $db) = @_;
  if($db){
    $self->param('_output_db',$db);
  }
  if(!$self->param('_output_db')){
    my $db = $self->hrdb_get_con('target_db');
    $self->param('_output_db',$db);
  }
  return $self->param('_output_db');
}

sub get_adaptor{
  my ($self) = @_;
  return $self->output_db->get_GeneAdaptor;
}

sub paf_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->param('paf_slice', $slice);
  }
  if(!$self->param_is_defined('paf_slice')){
    my $slice = $self->fetch_sequence($self->input_id, $self->paf_source_db,
                                      $self->REPEATMASKING);
    $self->param('paf_slice', $slice);
  }
  return $self->param('paf_slice');
}

sub gene_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->param('gene_slice', $slice);
  }
  if(!$self->param_is_defined('gene_slice')){
    my $slice = $self->fetch_sequence($self->input_id, $self->gene_source_db,
                                      $self->REPEATMASKING);
    $self->param('gene_slice', $slice);
  }
  return $self->param('gene_slice');
}

sub output_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->param('output_slice', $slice);
  }
  if(!$self->param_is_defined('output_slice')){
    my $slice = $self->fetch_sequence($self->input_id, $self->output_db);
    $self->param('output_slice', $slice);
  }
  return $self->param('output_slice');
}

sub kill_list{
  my ($self, $arg) = @_;

  if($arg){
    $self->param('kill_list', $arg);
  }
  if(!$self->param_is_defined('kill_list')){
    my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => 'protein');
    $self->param('kill_list', $kill_list_object)->get_kill_list;
  }
  return $self->param('kill_list');
}

sub rejected_set{
  my ($self, $arg) = @_;
  $self->param('rejected_set', []) if(!$self->param_is_defined('rejected_set'));
  if($arg){
    $self->param('rejected_set', $arg);
  }
  return $self->param('rejected_set');
}

sub filter_object{
  my ($self, $arg) = @_;
  if($arg){
    $self->throw("RunnableDB::BlastMiniGenewise ".
          $arg." must have a method filter_genes") 
      unless($arg->can("filter_genes"));
    $self->param('filter_object', $arg);
  }
  if(!$self->param_is_defined('filter_object') && $self->FILTER_OBJECT){
    $self->require_module($self->FILTER_OBJECT);
    my %params = %{$self->FILTER_PARAMS};
    $arg = $self->FILTER_OBJECT->new(
                                     -slice => $self->query,
                                     -seqfetcher => $self->seqfetcher,
                                     %params,
                                    );
    $self->param('filter_object', $arg);
  }
  return $self->param('filter_object');
}

sub seqfetcher{
  my ($self, $arg) = @_;

  if($arg){
    $self->throw("RunnableDB::BlastMiniGenewise ".
          $arg." must have a method get_Seq_by_acc") 
      unless($arg->can("get_Seq_by_acc"));
    $self->param('seqfetcher', $arg);
  }
  if(!$self->param_is_defined('seqfetcher')){
    $self->require_module($self->SEQFETCHER_OBJECT);
    my %params = %{$self->SEQFETCHER_PARAMS};
    print $params{-db}->[0], "\n";
    $arg = $self->SEQFETCHER_OBJECT->new(
                                    %params,
                                   );
    $self->param('seqfetcher', $arg);
  }
  return $self->param('seqfetcher');
}

=head2 require_module

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Arg [2]   : string, module path
  Function  : uses perls require to use the past in module
  Returntype: returns module name with / replaced by ::
  Exceptions: throws if require fails
  Example   : my $parser = $self->require('Bio/EnsEMBL/Analysis/Tools/BPliteWrapper');

=cut

sub require_module{
  my ($self, $module) = @_;
  my $class;
  ($class = $module) =~ s/::/\//g;
  eval{
    require "$class.pm";
  };
  $self->throw("Couldn't require ".$class." Blast:require_module $@") if($@);
  return $module;
}



=head2 overlaps_fiveprime_end_of_slice

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Bio::EnsEMBL::Feature
  Arg [3]   : Bio::EnsEMBL::Slice
  Function  : returns 1 if the features starts of the 5 end of the slice
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut



sub overlaps_fiveprime_end_of_slice{
  my ($self, $feature, $slice) = @_;
  if($feature->start < 1){
    return 1;
  }
  return 0;
}

#CONFIG METHODS


=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : hashref from config file
  Function  : call the superclass method to set all the varibles and carry
  out some sanity checking
  Returntype: N/A
  Exceptions: throws if certain variables arent set properlu
  Example   : 

=cut

sub read_and_check_config{
  my ($self, $hash) = @_;
  $self->SUPER::read_and_check_config($hash);

  #######
  #CHECKS
  #######
  foreach my $var(qw(PAF_LOGICNAMES
                     PAF_SOURCE_DB
                     GENE_SOURCE_DB
                     OUTPUT_DB
                     SEQFETCHER_OBJECT)){
    $self->throw("RunnableDB::BlastMiniGenewise $var config variable is not defined")
      unless($self->$var);
  }
  $self->OUTPUT_BIOTYPE($self->analysis->logic_name) if(!$self->OUTPUT_BIOTYPE);
  $self->throw("PAF_LOGICNAMES must contain at least one logic name")
    if(@{$self->PAF_LOGICNAMES} == 0);
}


=head2 PAF_LOGICNAMES

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise

sub accession {
  my ($self, $val) = @_;
  if($val){
    $self->param('_accession',$val);
  }
  return $self->param('_accession');
}


sub PAF_LOGICNAMES{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PAF_LOGICNAMES', $arg);
  }
  return $self->param('PAF_LOGICNAMES');
}

sub PAF_MIN_SCORE_THRESHOLD {
  my ($self, $arg) = @_;
  if( defined($arg) ) {
    $self->param('PAF_MIN_SCORE_THRESHOLD', $arg);
  }
  return $self->param('PAF_MIN_SCORE_THRESHOLD');
}

sub PAF_UPPER_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PAF_UPPER_SCORE_THRESHOLD', $arg);
  }
  return $self->param('PAF_UPPER_SCORE_THRESHOLD');
}

sub PAF_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PAF_SOURCE_DB', $arg);
  }
  return $self->param('PAF_SOURCE_DB');
}

sub GENE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('GENE_SOURCE_DB', $arg);
  }
  return $self->param('GENE_SOURCE_DB');
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OUTPUT_DB', $arg);
  }
  return $self->param('OUTPUT_DB');
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OUTPUT_BIOTYPE', $arg);
  }
  return $self->param('OUTPUT_BIOTYPE');
}

sub GENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('GENEWISE_PARAMETERS',$arg);
  }
  return $self->param('GENEWISE_PARAMETERS');
}

sub MINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MINIGENEWISE_PARAMETERS',$arg);
  }
  return $self->param('MINIGENEWISE_PARAMETERS');
}

sub MULTIMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MULTIMINIGENEWISE_PARAMETERS', $arg);
  }
  return $self->param('MULTIMINIGENEWISE_PARAMETERS');
}

sub BLASTMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('BLASTMINIGENEWISE_PARAMETERS', $arg);
  }
  return $self->param('BLASTMINIGENEWISE_PARAMETERS');
}


sub EXONERATE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('EXONERATE_PARAMETERS', $arg);
  }
  return $self->param('EXONERATE_PARAMETERS');
}


sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FILTER_PARAMETERS',$arg);
  }
  return $self->param('FILTER_PARAMETERS');
}

sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FILTER_OBJECT', $arg);
  }
  return $self->param('FILTER_OBJECT');
}

sub BIOTYPES_TO_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('BIOTYPES_TO_MASK', $arg);
  }
  return $self->param('BIOTYPES_TO_MASK');
}

sub EXON_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('EXON_BASED_MASKING', $arg);
  }
  return $self->param('EXON_BASED_MASKING');
}

sub GENE_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('GENE_BASED_MASKING', $arg);
  }
  return $self->param('GENE_BASED_MASKING');
}


sub POST_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('POST_GENEWISE_MASK', $arg);
  }
  return $self->param('POST_GENEWISE_MASK');
}

sub PRE_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PRE_GENEWISE_MASK', $arg);
  }
  return $self->param('PRE_GENEWISE_MASK');
}

sub REPEATMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('REPEATMASKING', $arg);
  }
  return $self->param('REPEATMASKING');
}

sub SEQFETCHER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('SEQFETCHER_OBJECT', $arg);
  }
  return $self->param('SEQFETCHER_OBJECT');
}

sub SEQFETCHER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('SEQFETCHER_PARAMS', $arg);
  }
  return $self->param('SEQFETCHER_PARAMS');
}

sub USE_KILL_LIST{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('USE_KILL_LIST', $arg);
  }
  return $self->param('USE_KILL_LIST');
}

sub LIMIT_TO_FEATURE_RANGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('LIMIT_TO_FEATURE_RANGE',$arg);
  }
  return $self->param('LIMIT_TO_FEATURE_RANGE');
}


sub FEATURE_RANGE_PADDING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FEATURE_RANGE_PADDING',$arg);
  }
  return $self->param('FEATURE_RANGE_PADDING');
}

sub WRITE_REJECTED{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->param('WRITE_REJECTED', $arg);
  }
  return $self->param('WRITE_REJECTED');
}

sub REJECTED_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('REJECTED_BIOTYPE', $arg);
  }
  return $self->param('REJECTED_BIOTYPE');
}

sub SOFTMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('SOFTMASKING', $arg);
  }
  return $self->param('SOFTMASKING');
}


sub MAKE_SIMGW_INPUT_ID_PARMAMS {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MAKE_SIMGW_INPUT_ID_PARMAMS', $arg);
  }
  return $self->param('MAKE_SIMGW_INPUT_ID_PARMAMS');
}

=head2 group_genes_by_id

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Genes
  Function  : groups the genes on the basis of transcript supporting features
  Returntype: hashref
  Exceptions: throws if a transcript has no supporting features
  Example   : 

=cut


sub group_genes_by_id{
  my ($self, $genes) = @_;
  my %hash;
  foreach my $gene(@$genes){
    my $transcripts = $gene->get_all_Transcripts;
    my $transcript = $transcripts->[0];
    my $sfs = $transcript->get_all_supporting_features;
    $self->throw(id($transcript)." appears to have no supporting features") if(!@$sfs);
    my $id = $sfs->[0]->hseqname;
    if(!$hash{$id}){
      $hash{$id} = 1;
    }else{
      $hash{$id}++;
    }
  }
  return \%hash;
}

sub OPTIMAL_LENGTH {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OPTIMAL_LENGTH', $arg);
  }
  return $self->param('OPTIMAL_LENGTH');
}

1;
