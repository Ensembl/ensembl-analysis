#!/usr/bin/env perl

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

use strict;
use warnings;
use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::SeqIO;


my @dbID;
my $count =0;
my $num=0;
my $start;
my $logic_name;
my $SE_logic_name;
my @input_ids;
my @multiexon_files;
my $ref_db = 'REFERENCE_DB';

GetOptions( '-pseudogene_logic_name:s' => \$logic_name,
            '-SE_logic_name:s'         => \$SE_logic_name,
            '-ref_db:s'                => \$ref_db );

my $config = read_config("Bio::EnsEMBL::Analysis::Config::Pseudogene");
my $pseudo_config;

foreach my $key ( keys %{$config->{"PSEUDOGENE_CONFIG_BY_LOGIC"}->{"DEFAULT"} } ) {
  if ( $config->{"PSEUDOGENE_CONFIG_BY_LOGIC"}->{$logic_name}->{$key} ) {
    $pseudo_config->{$key} =  $config->{"PSEUDOGENE_CONFIG_BY_LOGIC"}->{$logic_name}->{$key};
  } else {
    $pseudo_config->{$key} =  $config->{"PSEUDOGENE_CONFIG_BY_LOGIC"}->{"DEFAULT"}->{$key};
  }
}

my $usage = "prepare_SplicedElsewhere.pl 
-pseudogene_logic_name < logic name used to run the pseudogene analysis > 
-SE_logic_name <dummy analysis which will run spliced elsewhere>
-ref_db  $ref_db < hash key in databases.pm of the reference database
Loads input ids into pseudo_db and gets output files from pseudogene and makes them into a blast db";


die $usage unless($logic_name && $SE_logic_name );

my $db = get_db_adaptor_by_string($ref_db);
# we need a pipeline adaptor
$ref_db =
  Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( -dbname => $db->dbc->dbname,
                                                 -dbport => $db->dbc->port,
                                                 -host   => $db->dbc->host,
                                                 -user   => 'ensadmin',
                                                 -pass => $db->dbc->password,
  );

if ( $pseudo_config->{"SINGLE_EXON"} ) {
  print "Making input ids for single exon genes\n";

  my $fa       = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($ref_db);
  my $aa       = $ref_db->get_AnalysisAdaptor;
  my $analysis = $aa->fetch_by_logic_name( $pseudo_config->{"SINGLE_EXON"} );
  my $multifile =
    $pseudo_config->{"PS_MULTI_EXON_DIR"} . "/all_multi_exon_genes.fasta";
  die( "analysis object not found " . $pseudo_config->{"SINGLE_EXON"} . "\n" )
    unless ($analysis);
  my @ids = @{ $fa->fetch_by_analysis($analysis) };
  @ids = sort { $a->dbID <=> $b->dbID } @ids;

  foreach my $id (@ids) {
    if ( $count == 0 ) {
      $start = $id->dbID;
    }
    $count++;
    if ( $count == $pseudo_config->{"PS_CHUNK"} ) {
      $num++;
      push @input_ids, "$start:" . $id->dbID;
      $count = 0;
    }
  }
  if ( $count > 0 ) {
    push @input_ids, "$start:" . $ids[$#ids]->dbID;
  }

  my $inputIDFactory =
    new Bio::EnsEMBL::Pipeline::Utils::InputIDFactory(
                                                -db         => $ref_db,
                                                -top_level  => 'top_level',
                                                -logic_name => $SE_logic_name,
                                                -file       => 1 );
  $inputIDFactory->input_ids( \@input_ids );
  $inputIDFactory->store_input_ids;

  print STDERR "Pooling multiexon genes into single blastDB .\n";

  my $db_output = Bio::SeqIO->new( -file   => ">$multifile",
                                   -format => 'fasta' );

  unless ( opendir( DIR, $pseudo_config->{"PS_MULTI_EXON_DIR"} ) ) {
    closedir(DIR);
    die "cannot read files from "
      . $pseudo_config->{"PS_MULTI_EXON_DIR"} . "\n";
  }
  foreach ( readdir(DIR) ) {
    my $file = "$_";
    if ( $file =~ m/^multi_exon_seq.*\.fasta$/ ) {
      my $bioseq = Bio::SeqIO->new(
                 -file => $pseudo_config->{"PS_MULTI_EXON_DIR"} . "/" . $file,
                 -format => 'fasta' );
      while ( my $seq = $bioseq->next_seq ) {
        $db_output->write_seq($seq);
      }
    }
  }
  #   system ("rm $file");
  system("xdformat -n $multifile");
} ## end if ( $pseudo_config->{...

print "Finished\n";

exit 0;
