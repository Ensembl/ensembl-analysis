=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadGenomeFlatfiles;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Bio::EnsEMBL::Hive::Utils qw(destringify);
use POSIX;
use File::Spec::Functions;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  my $assembly_registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(%{$self->param_required('assembly_registry_db')});
  $self->hrdb_set_con($assembly_registry_dba,'assembly_registry_db');
}


sub run {
  my $self = shift;

  my $gca = $self->param_required('iid');
  my ($chain,$version) = $self->split_gca($gca);
  $chain =~ /^GCA\_(\d{3})(\d{3})(\d{3})/;
  my $gca_p1 = $1;
  my $gca_p2 = $2;
  my $gca_p3 = $3;

  my $ncbi_base_ftp = $self->param_required('ncbi_base_ftp');

  my $assembly_registry_dba = $self->hrdb_get_con('assembly_registry_db');
  my $assembly_name = $assembly_registry_dba->fetch_assembly_name_by_gca($gca);
  $assembly_name =~ s/[\- \/]/\_/g;

  my $assembly_dir_name = $gca."_".$assembly_name;
  my $fasta_file_name = $assembly_dir_name.'_genomic.fna.gz';
  my $full_ftp_path = $ncbi_base_ftp."/".catfile('GCA',$gca_p1,$gca_p2,$gca_p3,$assembly_dir_name,$fasta_file_name);
  my $output_dir = catfile($self->param_required('base_output_path'),$gca);
  unless(-e $output_dir) {
    my $make_output_dir_result = system('mkdir -p '.$output_dir);
    if($make_output_dir_result) {
      $self->throw("Failed to create output dir on the following path: ".$output_dir);
    }

    my $stripe_dir_result = system('lfs setstripe -c -1 '.$output_dir);
    if($stripe_dir_result) {
      $self->throw("Failed to stripe the output dir with: lfs setstripe -1 -c ".$output_dir);
    }
  }

  my $command = "wget -O ".catfile($output_dir,'genomic.fna.gz')." ".$full_ftp_path;
  my $wget_result = system($command);
  if($wget_result) {
    $self->throw("Failed to download genomic fna file. Commandline used:\n".$command);
  }

  my $gunzip_result = system('gunzip '.$output_dir.'/genomic.fna.gz');
  if($gunzip_result) {
    $self->throw("Failed to unzip genomic fna file. Commandline used:\ngunzip ".$output_dir."/genomic.fna.gz");
  }

  my $sed_result = system('sed -i "s/\([a-z]\+\)/\U\1\E/g" '.$output_dir.'/genomic.fna');
  if($sed_result) {
    $self->throw("Failed to uppercase genomic fna file. Commandline used:\nsed -i sed -i 's/\([a-z]\+\)/\U\1\E/g' ".$output_dir."/genomic.fna");
  }

  $self->path_to_genomic_fasta($output_dir);
}


sub write_output {
  my $self = shift;

  my $path_to_genomic_fasta = $self->path_to_genomic_fasta;
  my $job_params = destringify($self->input_job->input_id);
  $job_params->{'path_to_genomic_fasta'} = $path_to_genomic_fasta;
  $self->input_job->input_id($job_params);

}


# Should consider removing redundancy with the same method in the assembly registry adaptor
sub split_gca {
  my ($self,$chain_version) = @_;

  unless($chain_version =~ /^(GCA\_\d{9})\.(\d+)$/) {
    $self->throw("Could not parse versioned GCA. GCA used: ".$chain_version);
  }

  my $chain = $1;
  my $version = $2;

  return($chain,$version);
}

sub path_to_genomic_fasta {
  my ($self,$path) = @_;
  if($path) {
    $self->param('_path_to_genomic_fasta',$path);
  }

  return($self->param('_path_to_genomic_fasta'));
}

1;
