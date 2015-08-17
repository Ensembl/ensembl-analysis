# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

This module creates a file which contains the list of genes which have not been processed during the merge. It requires merge.pl output file 'genes-processed.txt' in 'output_dir'.

=head1 OPTIONS

-dbhost         database host name

-dbport         database port (default 3306)

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveListUnprocessed 

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveListUnprocessed;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      processed_genes_filename => '',
      output_file => 'genes_to_copy.ids',
      output_dir => '',
      host_secondary => '',
      user_secondary => '',
      password_secondary => '',
      database_secondary => '',
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
  
  if (not $self->param('output_file') or not $self->param('output_dir') or not $self->param('host_secondary') or not $self->param('user_secondary') or not $self->param('database_secondary') or not $self->param('database_secondary') ) {
    throw("Parameters missing");
  }

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_dir') =~ /\/$/)) {
    $self->param('output_dir',$self->param('output_dir')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_dir')) {
    run_command("mkdir -p ".$self->param('output_dir'),"Create output path.");
  }

  my @gene_ids_to_copy = get_unprocessed_genes($self->param('processed_genes_filename'),
                                               $self->param('output_dir'),
                                               $self->param('host_secondary'),
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
  my ($filename,$output_dir,$dbhost,$user,$pass,$dbname,$include,$exclude) = @_;  
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbname => $dbname,
                                              -host   => $dbhost,
                                              -user   => $user,
                                              -pass   => $pass,
                                              -port   => '3306');
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
  open PROCCESED_GENES, $output_dir."/".$filename or die $!;
  while (my $proccesed_gene_id = <PROCCESED_GENES>) {
    chomp($proccesed_gene_id);
    delete $genes_to_copy{$proccesed_gene_id};
  }
  
  return keys(%genes_to_copy);
}

1;
