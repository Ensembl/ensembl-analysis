=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::miRNA - 

=head1 SYNOPSIS

     my $runnableDB = Bio::EnsEMBL::Analysis::RunnableDB::miRNA->new(
            	-db       => $db_adaptor,
		-input_id => 'analysis logic name',
		-analysis => $analysis,
	        );
    $runnabledb->fetch_input();
    $runnabledb->run();
    $runnabledb->write_output();

=head1 DESCRIPTION

RunnableDB to provide database access for miRNA detection.
Runs as an accumulator job on miRNA blast hits found by Bio::EnsEMBL::RunnableDB::BlastmiRNA
Takes an analysis logic name as an input id and uses it to fetch all dna align features associated with 
that analysis.
It then groups the dna align features by miRNA families and ignores families with > 50
members as there is a high probability that these are hitting repetitive sequences.
Creates and runs the miRNA runnable.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivemiRNA;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::miRNA;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

  Title      : fetch_input
  Usage      : $miRNA->fetch_input();
  Function   : opens and stores connections to databases in Bio::EnsEMBL::Analysis::Config::Databases
             : fetches all dna align features by analysis logic name specified in the 
             : input id 
  Returns    : Hash reference
  Exceptions : throws if the analysis object is not found or no dna align features are retrieved
  Args       : None

=cut

sub fetch_input{
  my ($self) = @_;

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                             -program_file => $self->param('cmsearch_exe_path'),
                                             -db_file => $self->param('blast_db_dir_path'),
                                           );
  $self->analysis($analysis);

  # The output db should be the one that the dafs to check have been written to
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));


  if($dna_dba) {
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');


  my $daf_adaptor = $output_dba->get_DnaAlignFeatureAdaptor;
  my $slice_adaptor = $output_dba->get_SliceAdaptor;


  my $daf_ids = $self->param('iid');
  my $dafs = [];
  foreach my $db_id (@{$daf_ids}) {
    push(@{$dafs},$daf_adaptor->fetch_by_dbID($db_id));
  }


  my %families = %{$self->family($dafs)};
  undef($dafs);

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::miRNA->new
    (
     -queries => \%families,
     -analysis => $self->analysis,
    );
  $self->runnable($runnable);
}


=head2 family

  Title      : family
  Usage      : my %families = %{$self->family(\@dafs)};
  Function   : order dafs by family, removes families with blast hits to repeats
  Returns    : Hash reference
  Exceptions : None
  Args       : Array ref of Bio::EnsEMBL::DnaDnaAlignFeatures

=cut

sub family{
  my ($self,$dafs_ref) = @_;
  my %families;
  foreach my $daf (@$dafs_ref){
    push @{$families{$daf->hseqname}},$daf;
  }
  my %filtered_fam;
  foreach my $key (keys %families){
  if (scalar @{$families{$key}} <= 50){
    $filtered_fam{$key} = $families{$key};
    } else {
      # take top scoring 50 hits
      my @array = sort {$a->p_value <=> $b->p_value} @{$families{$key}};
      my @filtered_array =  splice(@array,0,50);
      $filtered_fam{$key} = \@filtered_array;
    }
 }
  return \%filtered_fam;
}

=head2 write_output

  Args       : none
  Description: Writes the single exon miRNA genes into the final genebuild database,
             : also stores attributes associated with the transcript
  Exceptions : Throws if gene or transcript attribute fail to write to the database
  Returntype : scalar

=cut

sub write_output{
  my ($self) = @_;

  my $db_adaptor = $self->hrdb_get_con('output_db');
  my $gene_adaptor = $db_adaptor->get_GeneAdaptor;
  my $transcript_adaptor = $db_adaptor->get_TranscriptAdaptor;
  my $attribute_adaptor = $db_adaptor->get_AttributeAdaptor;
  my @attributes;
  my $xref;

  foreach my $gene_hash (@{$self->output}){
    my $gene = $gene_hash->{'gene'};
    @attributes = @{$gene_hash->{'attrib'}};
    $gene->analysis($self->analysis);
    foreach my $trans (@{$gene->get_all_Transcripts}){
      $trans->analysis($self->analysis);
    }
    $gene->slice($self->query) if(!$gene->slice);
    $self->feature_factory->validate($gene);
    eval{
      $gene_adaptor->store($gene);
    };
    if($@){
      $self->throw("miRNA:store failed, failed to write ".$gene." to ".
		   "the database $@");
    }
    foreach my $trans (@{$gene->get_all_Transcripts}){
      eval{
	$attribute_adaptor->store_on_Transcript($trans->dbID,\@attributes);
	$transcript_adaptor->update($trans);
	$gene_adaptor->update($gene);
      };
      if($@){
	$self->throw("miRNA:store failed, failed to write ".@attributes." on transcript ".
		     $trans." in the database $@");
      }
    }
  }
  return 1;
}

1;
