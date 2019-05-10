=head1 LICENSE

 Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Pmatch -

=head1 SYNOPSIS


=head1 DESCRIPTION


Pmatch is a fast alignment program written by Richard Durbin we used
to align species specific proteins to the genome, (We also use it to
align very closely related species proteins sets e.g Mouse to Rat or
Fugu to Tetraodon).

The pmatch source code is available from the
sanger cvs respository in module rd-utils,
(http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/rd-utils/?root=ensembl).
The code to run this process in the ensembl code
base can be found in 2 RunnableDBs and a config
file. Bio::EnsEMBL::Analysis::RunnableDB::Pmatch,
Bio::EnsEMBL::Analysis::RunnableDB::BestPmatch and
Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut



package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePmatch;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::Runnable::Pmatch;
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDBSeqFiles');


sub fetch_input{
  my ($self) = @_;

  $self->create_analysis;
  my $dna_db = $self->get_database_by_name('dna_db');
  $self->hrdb_set_con($self->get_database_by_name('target_db', $dna_db), 'target_db');
  my $slice = $self->fetch_sequence($self->input_id, $self->hrdb_get_con('target_db'),
                                    $self->REPEAT_MASKING);
  $self->query($slice);
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  $program = $self->BINARY_LOCATION if(!$program);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Pmatch->new(
          -query => $self->query,
          -program => $program,
          -protein_file => $self->PROTEIN_FILE,
          -analysis => $self->analysis,
          -max_intron_length => $self->MAX_INTRON_LENGTH,
          -min_coverage => $self->MIN_COVERAGE,
          -options => $self->OPTIONS,
         );
#  my $proteins;
#  if ($self->param_required('iid_type') eq 'db_seq') {
#    my $accs = $self->input_id;
#    if (ref($accs) ne 'ARRAY') {
#      $accs = [$accs];
#    }
#    my $filename = $runnable->write_seq_file($self->get_query_seqs($accs));
#    $self->PROTEIN_FILE($filename);
#  }
#  $runnable->protein_file($self->PROTEIN_FILE);
  $self->runnable($runnable);
}


sub PROTEIN_FILE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PROTEIN_FILE', $arg);
  }
  return $self->param('PROTEIN_FILE');
}


sub MIN_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MIN_COVERAGE', $arg);
  }
  return $self->param('MIN_COVERAGE');
}

sub BINARY_LOCATION{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('BINARY_LOCATION', $arg);
  }
  return $self->param('BINARY_LOCATION');
}


sub MAX_INTRON_LENGTH{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MAX_INTRON_SIZE', $arg);
  }
  return $self->param('MAX_INTRON_SIZE');
}


sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('target_db', $arg);
  }
  return $self->param('target_db');
}

sub REPEAT_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    throw("Runnable::Pmatch ".$arg." must be an array ref of logic names not .".$arg)
      unless(ref($arg) eq 'ARRAY');
    $self->param('REPEAT_MASKING', $arg);
  }
  return $self->param('REPEAT_MASKING');
}

sub OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('_CONFIG_OPTIONS', $value);
  }

  if ($self->param_is_defined('_CONFIG_OPTIONS')) {
    return $self->param('_CONFIG_OPTIONS');
  } else {
    return undef;
  }
}





sub get_adaptor{
  my ($self) = @_;
  my $output_db = $self->hrdb_get_con('target_db');
  return $output_db->get_ProteinAlignFeatureAdaptor;
}


sub write_output{
  my ($self) = @_;

  my $adaptor = $self->get_adaptor;
  my %unique;
  foreach my $feature (@{$self->output}) {
    $feature->analysis($self->analysis);
    $feature->slice($self->query) unless ($feature->slice);
    my $unique_string = $feature->start.' '.$feature->end.' '.$feature->score.' '.$feature->hseqname;
    next if ( exists $unique{$unique_string});
    $unique{$unique_string} = 1;
    $self->feature_factory->validate($feature);
    $adaptor->store($feature);
  }
}

1;
