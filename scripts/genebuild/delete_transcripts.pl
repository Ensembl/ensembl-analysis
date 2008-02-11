#!/usr/local/bin/perl -w
=head1 NAME

  delete_transcripts.pl

=head1 SYNOPSIS
 
  delete_transcripts.pl
  deletes transcripts from given database whose ids are passed in through STDIN

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    -help      summary of options


=head2 EXAMPLES

./delete_transcripts.pl -dbhost ecs2b -dbuser ensadmin -dbpass **** -dbname rat_Jun03_mk2 transcripts_to_delete

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

my $host;
my $port;
my $dbname;
my $user;
my $pass;
my $config_dbname;

&GetOptions( 
            'dbhost:s'      => \$host,
            'dbport:n'      => \$port,
            'dbname:s'    => \$dbname,
            'dbuser:s'    => \$user,
            'dbpass:s'      => \$pass,
            'config_dbname:s' => \$config_dbname,
           );


my $db;

if($config_dbname){
  $db = get_db_adaptor_by_string($config_dbname);
}elsif($dbname && $host){
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host   => $host,
                                            -user   => $user,
                                            -port   => $port,
                                            -dbname => $dbname,
                                            -pass => $pass,
                                           );
}else{
  throw("Need to pass either -dbhost $host and -dbname $dbname or ".
        "-config_dbname $config_dbname for the script to work");
}

my $transcript_adaptor = $db->get_TranscriptAdaptor;

while(<>){
  chomp;
  my $transcript_id= $_;
  
  #my $sth = $db->prepare("delete from transcript where transcript_id = $transcript_id");
  #$sth->execute;
  eval{
    my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    $transcript_adaptor->remove($transcript);
    print STDERR "Deleted $transcript_id\n";
  };
  if($@){
    print "Couldn't remove transcript $transcript_id ($@)\n";
  }
}

