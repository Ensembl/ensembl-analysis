#!//usr/local/ensembl/bin/perl

#this script takes database options and a file of gene ids and copies them between two 
#databases. It can if asked split multi transcript genes into single genes
#
#an example commandline might look like
#
#perl copy_genes.pl -in_config_name COALESCER_DB -out_config_name UTR_DB gene_ids_to_copy
#
# or
#
#perl copy_genes.pl -dbhost host -dbuser ensro -dbname est_db -outhost host1 -outuser ensadmin -outpass ensembl -outdbname utr_db gene_ids_to_copy
#
#When using the in and out_config_name options it reads the equivalent database from the
#Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases file

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

my $host      = '';
my $user      = '';
my $pass      = undef;
my $dbname    = '';
my $port      = 3306;

my $outhost   = '';
my $outuser   = '';
my $outpass   = '';
my $outdbname = undef;
my $outport   = 3306;

my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnapass   = '';
my $dnadbname = undef;
my $dnaport   = 3306;
my $in_config_name;
my $out_config_name;

my $split = 0;

&GetOptions(
            'dbhost:s'        => \$host,
            'dbuser:s'        => \$user,
            'dbpass:s'        => \$pass,
            'dbname:s'      => \$dbname,
            'dbport:n'        => \$port,
            'in_config_name:s' => \$in_config_name,
            'out_config_name:s' => \$out_config_name,
            'outhost:s'     => \$outhost,
            'outuser:s'     => \$outuser,
            'outpass:s'     => \$outpass,
            'outdbname:s'   => \$outdbname,
            'outport:n'     => \$outport,
            'split!' => \$split,
           );


#print "Connecting to ".$dbname." at ".$host." ".$user."\n";
my $db;
if($in_config_name){
  $db = get_db_adaptor_by_string($in_config_name);
}else{
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host   => $host,
                                           -user   => $user,
                                           -pass   => $pass,
                                           -port   => $port,
                                           -dbname => $dbname
                                          );
}





my $ga = $db->get_GeneAdaptor;

my @genes;
my @copy_genes;

while(<>){
  chomp;
  my $gene_id= $_;
  my $gene = $ga->fetch_by_dbID($gene_id);
  empty_Gene($gene);
  push(@copy_genes, $gene);
}

my $outdb;
if($out_config_name){
  $outdb = get_db_adaptor_by_string($out_config_name);
}else{
  $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                              -host   => $outhost,
                                              -user   => $outuser,
                                              -pass   => $outpass,
                                              -port   => $outport,
                                              -dbname => $outdbname
                                             );
}

if($split){
  foreach my $gene(@copy_genes){
    push(@genes, @{convert_to_genes($gene->get_all_Transcripts, $gene->analysis, $gene->biotype)});
  }
}else{
  @genes = @copy_genes;
}

my $outga = $outdb->get_GeneAdaptor;

foreach my $gene(@genes){
  eval{
    fully_load_Gene($gene);
    empty_Gene($gene);
    $outga->store($gene);
  }
}
