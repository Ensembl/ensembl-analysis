# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 NAME 

HiveListUnprocessed.pm

=head1 DESCRIPTION

This module creates a file which contains the list of genes which have not been processed during the merge. It requires merge.pl output file in 'output_dir'.

=head1 OPTIONS

-processed_genes_filename   File name of the merge process output file containing the list of processed genes (usually 'havana_merge_list_processed_genes.ids'). 
-output_file                File name of the file containing the list of gene identifiers which were not processed during the merge analysis.
-output_dir                 Directory where 'output_file' will be written.
-host_secondary             Secondary database host.
-port_secondary             Secondary database port.
-user_secondary             Secondary database user.
-password_secondary         Secondary database password.
-database_secondary         Secondary database username.
-secondary_include          Secondary database list of gene logic names to include. 
-secondary_exclude          Secondary database list of gene logic names to exclude (if 'secondary_include' is used, 'secondary_exclude' will not have any effect).

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveListUnprocessed -processed_genes_filename -output_dir OUTPUT_DIR -output_file -host_secondary ENSEMBL_DB_HOST -user_secondary => READ_ONLY_USER -password_secondary READ_ONLY_PASS -database_secondary => ENSEMBL_DB_NAME

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveListUnprocessed;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      processed_genes_filename => undef,
      output_file => 'genes_to_copy.ids',
      output_dir => undef,
      host_secondary => undef,
      port_secondary => undef,
      user_secondary => undef,
      password_secondary => undef,
      database_secondary => undef,
      secondary_include => '',
      secondary_exclude => '',
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  
  $self->param_required('processed_genes_filename');
  $self->param_required('output_file');
  $self->param_required('output_dir');
  $self->param_required('host_secondary');
  $self->param_required('port_secondary');
  $self->param_required('user_secondary');
  $self->param_required('password_secondary');
  $self->param_required('database_secondary');

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_dir') =~ /\/$/)) {
    $self->param('output_dir',$self->param('output_dir')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_dir')) {
    run_command("mkdir -p ".$self->param('output_dir'),"Create output path.");
  }

  my @gene_ids_to_copy = $self->get_unprocessed_genes($self->param('processed_genes_filename'),
                                               $self->param('output_dir'),
                                               $self->param('host_secondary'),
                                               $self->param('port_secondary'),
                                               $self->param('user_secondary'),
                                               $self->param('password_secondary'),
                                               $self->param('database_secondary'),
                                               $self->param('secondary_include'),
                                               $self->param('secondary_exclude'));

  my $filepath = $self->param('output_dir').$self->param('output_file');
  open(my $fh, ">$filepath") or die $!;
  foreach my $gene_id (@gene_ids_to_copy) {
    printf $fh $gene_id."\n";
  }

  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub get_unprocessed_genes() {
  my ($self, $filename,$output_dir,$dbhost,$dbport,$user,$pass,$dbname,$include,$exclude) = @_;
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbname => $dbname,
                                              -host   => $dbhost,
                                              -user   => $user,
                                              -pass   => $pass,
                                              -port   => $dbport);
  my $ga = $db->get_GeneAdaptor();
  my %genes_to_copy = ();
  
  if ($include) {
    # include the genes whose logic name matches
    foreach my $gene (@{$ga->fetch_all_by_logic_name($include)}) {
      $genes_to_copy{$gene->dbID()} = 1;
    }
  } elsif ($exclude) {
    # include the genes whose logic name does not match
    foreach my $gene (@{$ga->fetch_all()}) {
      my $gene_analysis = $gene->analysis();
      if ($gene_analysis->logic_name() ne $exclude) {
        $genes_to_copy{$gene->dbID()} = 1;
      }
    }
  } else {
    # include all genes
    %genes_to_copy = map { $_ => 1 } @{$ga->list_dbIDs()};
  }
  
  # only copy the genes which have not been processed by the merge code
  if (-e "$output_dir/$filename") {
      open PROCCESED_GENES, $output_dir."/".$filename or die $!;
      while (my $proccesed_gene_id = <PROCCESED_GENES>) {
        chomp($proccesed_gene_id);
        delete $genes_to_copy{$proccesed_gene_id};
      }
  }
  else {
      $self->warning("$output_dir/$filename does not exists, check that it's OK!\n");
  }
  
  return keys(%genes_to_copy);
}

1;
