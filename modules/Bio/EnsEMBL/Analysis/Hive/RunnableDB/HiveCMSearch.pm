=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Infernal -

=head1 SYNOPSIS

    my $runnableDB = Bio::EnsEMBL::Analysis::RunnableDB::Infernal->new(
              -db       => $db_adaptor,
    -input_id => $flag_id,
    -analysis => $analysis,
          );
    $runnabledb->fetch_input();
    $runnabledb->run();
    $runnabledb->write_output();


=head1 DESCRIPTION

RunnableDB to wrap cmsearch - part of the Infernal suite of programs by Sean Eddy.
Uses RFAM blast hits to identify regions where cmsearch should be run. Blast hits
are run using the RfamBlast module and then are filtered on a familly by familly basis.
Blast hits that look promising are flagged by predict_ncRNAs.pl and the flags are used
as input_ids.
Uses the Bio::EnsEMBL::Analysis::Config::Databases config file. Writes non coding genes to
the $GB_FINALDB database and gets DNA from $GB_ database. The blast hits from BlastRfam are
stored in the pipeline database;
Creates single exon non-coding genes with gene descriptions obtained from Rfam.descriptions
file. The dna align feature representing the initial blast hit is added as a supporting feature
and the RNA secondary structure predicted by Infernal is added as a transcript attribute
in a length encoded string form to take up less room in the db.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCMSearch;

use Bio::Seq;
use Bio::EnsEMBL::Analysis::Runnable::CMSearch;

use strict;
use warnings;
use feature 'say';

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

  Title   :   fetch_input
  Usage   :   $self->fetch_input
  Function:   Opens connections to the final gene build database to store the genes in.
          :   Also opens connection to the genebuild DNA database to fetch DNA from.
          :   Parses flags use as input_ids and fetches dna_align_features to run cmsearch
          :   on.
  Returns :   none
  Args    :   none

=cut

sub fetch_input{
  my ($self)=@_;

  $self->setup_fasta_db;
  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                             -program_file => $self->param('cmsearch_exe_path'),
                                             -gff_feature => $self->param('rfam_seeds'),
                                             -db_file => $self->param('rfam_cm'),
                                           );
  $self->analysis($analysis);

  # The output db should be the one that the dafs to check have been written to
  my $repeat_dba = $self->get_database_by_name('output_db');
  my $output_dba = $self->get_database_by_name('output_db');
  $self->hrdb_set_con($repeat_dba,'repeat_db');
  $self->hrdb_set_con($output_dba,'output_db');

  my $repeat_logic_names = $self->param('repeat_logic_names');
  
  my $slice_adaptor = $output_dba->get_SliceAdaptor();

  my $input_ids = $self->param('iid');

  foreach my $input_id (@{$input_ids}) {
    my $slice = $self->fetch_sequence($input_id,$repeat_dba,$repeat_logic_names);
    
     unless ($slice->seq =~ /[CATG]{3}/) {
      say "The following slice is > 3bp after applying repeatmasking, will skip:";
       say $slice->name;
       next;
     }
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::CMSearch->new
                   (
                    -queries  => $slice, 
                    -datadir => $self->param('output_dir'),
                    -analysis => $self->analysis,
                    -program => $self->analysis->program_file
                   );
    $self->runnable($runnable);
  }

  # Make  the runnable
   if ($self->param('disconnect_jobs')) {
     $repeat_dba->dbc->disconnect_when_inactive(1);
  }
  return 1;
}

=head2 run

  Args       : none
  Description: Runs the runnable.
  Exceptions : Throws if the the runnable module is not set
  Returntype : scalar

=cut

sub run{
  my ($self) = @_;
  $self->dbc->disconnect_when_inactive(1) if ($self->param('disconnect_jobs'));
  foreach my $runnable (@{$self->runnable}) {
    $self->throw("Runnable module not set") unless ($runnable->isa("Bio::EnsEMBL::Analysis::Runnable"));
    $runnable->run();
    $self->output($runnable->output);
  }
  $self->dbc->disconnect_when_inactive(0);
  return 1;
}


=head2 write_output

  Args       : none
  Description: Writes the single exon ncRNA genes into the final genebuild databse,
             : also stores the structure attribute associated with the transcript
  Exceptions : Throws if gene or transcript attribute fail to write to the database
  Returntype : scalar

=cut

sub write_output{
  my ($self) = @_;

  my $adaptor = $self->hrdb_get_con('output_db');
  $adaptor->dbc->disconnect_when_inactive(0);
  my $aa = $adaptor->get_AttributeAdaptor;
  my $dbea = $adaptor->get_DBEntryAdaptor;
  my $gene_adaptor = $adaptor->get_GeneAdaptor;
  my $transcript_adaptor = $adaptor->get_TranscriptAdaptor;
  my @attributes;
  my $xref;

  say "In write output, total genes for output: ".scalar(@{$self->output});
  foreach my $gene_hash (@{$self->output}){
    my $gene = $gene_hash->{'gene'};
    @attributes = @{$gene_hash->{'attrib'}};
    $xref = $gene_hash->{'xref'};
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
      $self->throw("Infernal:store failed, failed to write ".$gene." to ".
       "the database $@");
    }
    foreach my $trans (@{$gene->get_all_Transcripts}){
      eval{
  $aa->store_on_Transcript($trans->dbID,\@attributes);
  $trans->display_xref($xref);
  $gene->display_xref($xref);
  $transcript_adaptor->update($trans);
  $gene_adaptor->update($gene);
      };
      if($@){
  $self->throw("Infernal:store failed, failed to write xrefs or attributes  on transcript ".
         $trans." in the database $@");
      }
    }
  }
  return 1;
}


1;
