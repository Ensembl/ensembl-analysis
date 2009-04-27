#!/usr/local/ensembl/bin/perl -w

=head1 NAME

ensembl-analysis/scripts/sequence_dump.pl

=head1 SYNOPSIS

This script will dump the sequence of all the seq_regions in a particular
coordinate system 

=head1 DESCRIPTION

The script can dump into individual fasta files for each seq_region or one
file containing all the sequences. The format is fasta as standard but others
can be specified. The sequence can also be masked for repeats in either
normal uppercase Ns or with softmasking

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (user= in locator)

    -dbpass    For RDBs, what password to use (pass= in locator)

    -species   Species name/alias.  Only required for multispecies dna DBs

    -multi_species Boolean for multi-species DBs
    
    -species_id Species ID in the database. Only required for multispecies DBs

    -coord_system_name the name of the coordinate system you want to dump

    -output_dir the directory to dump the files too

    -toplevel to indicate you want to dump the top level seq regions

    -sequence_level to indicate you want to dump the sequence level seq
                   regions
    -nonref to indicate you want available non-reference regions

    -format the format you want your sequence dumped in. As standard this is
            fasta format and this can be any format Bio::SeqIO can dump

    -extension the file extention you want to give the dumped files, by 
               default this is fa

    -mask   to indicate you want the sequence repeatmasked

    -mask_repeat, which logic name of repeats you want masked. This can 
            appear on the command line several times. If it doesnt but -mask
            does all repeats will be masked

    -softmask to indicate you want to softmask the repeats (ie lower case
              rather than upper case Ns)

    -onefile to indicate you want all the sequences in one file

    -help this will print out the docs
=head1 EXAMPLES

perl sequence_dump.pl -dbhost myhost -dbuser myuser -dbpass mypass -dbname
  mydatabase -dbport 3306 -coord_system_name chromosome 
  -output_dir /path/to/output

this will dump each chromosome in a separate file in the output directory

perl sequence_dump.pl -dbhost myhost -dbuser myuser -dbpass mypass -dbname
  mydatabase -dbport 3306 -coord_system_name contig -onefile
  -output_dir /path/to/output

this will dump the sequence of all the contigs into one file in the output
directory

perl sequence_dump.pl -dbhost myhost -dbuser myuser -dbpass mypass -dbname
  mydatabase -dbport 3306 -coord_system_name chromosome 
  -output_dir /path/to/output -mask -mask_repeat RepeatMask -mask_repeat
  Dust -softmask

this will dump each chromosome sequence softmasked for RepeatMasker and
Dust repeats into to separate files


=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my ($species, $filename, $multi_species);
my $species_id = 1;
my $coord_system_name;
my $output_dir;
my $onefile;
my $top_level;
my $seq_level;
my $non_ref;
my $format = 'fasta';
my @logic_names;
my $mask;
my $softmask = 0;
my $extension = 'fa';
my $help;
&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'species=s'  => \$species,
            'multi_species'  => \$multi_species,
            'species_id=i' => \$species_id,
            'coord_system_name:s' => \$coord_system_name,
            'output_dir:s' => \$output_dir,
            'extension:s' => \$extension,
            'toplevel!' => \$top_level,
            'seqlevel!' => \$seq_level,
             'nonref!' => \$non_ref,
            'format!' => \$format,
            'mask!' => \$mask,
            'mask_repeat:s@' => \@logic_names,
            'softmask!' => \$softmask,
            'onefile!' => \$onefile,
            'filename=s' => \$filename,
            'help!' => \$help,
           ) or ($help = 1);



if ($help) {
    exec('perldoc', $0);
}

if(!$host || !$dbname || !$dbuser){
  throw("Need -dbhost $host -dbuser $dbuser and -dbname $dbname to run ".
        " use -help for docs");
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -dbname => $dbname,
   -host   => $host,
   -user   => $dbuser,
   -port   => $port,
   -pass   => $dbpass,
   -species => $species,
   -multispecies_db => $multi_species,
   -species_id => $species_id
);


if(!$coord_system_name && !$top_level && !$seq_level){
  throw("Must specify at least one of -coord_system_name -toplevel or ".
        "-seqlevel to run using -help for docs");
}

if($top_level && $seq_level){
  throw("Can't specify both -toplevel and -seqlevel must be one or the ".
        "other see -help for docs"); 
}

if(! $filename && (!$output_dir || ! -e $output_dir)){
  throw("Can't dump sequence into ".$output_dir." it doesn't exist ".
        "use -help for docs");
}

if($top_level){
  $coord_system_name = 'toplevel';
}
if($seq_level){
  $coord_system_name = 'seqlevel';
}


my $oneout;
if($onefile || $filename){
  $filename ||= $output_dir."/".$coord_system_name.".".$extension;  
  $oneout = Bio::SeqIO->new(
                            -file => ">".$filename,
                            -format => $format,
                           );
}




my $sa = $db->get_SliceAdaptor;

my $slices = $sa->fetch_all($coord_system_name, undef, $non_ref);

foreach my $slice(@$slices){
  my $seq;

  if($mask){
    $seq = $slice->get_repeatmasked_seq(\@logic_names, $softmask);
  }else{
    $seq = $slice;
  }

  my $seqout;
  if($filename){
    $seqout = $oneout;
  }else{

        my $filename = $output_dir."/".$slice->seq_region_name.".".$extension;
    $seqout = Bio::SeqIO->new(
                              -file => ">".$filename,
                              -format => $format,
                             );
  }
  $seqout->write_seq($seq);
}



