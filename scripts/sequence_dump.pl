#!/usr/local/ensembl/bin/perl -w
use Data::Dumper;

#print Dumper(\@ARGV); die;
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

    -coord_system_version the version of the coordinate system you want to dump

    -output_dir the directory to dump the files too

    -toplevel to indicate you want to dump the top level seq regions

    -sequence_level to indicate you want to dump the sequence level seq
                   regions
    -nonref to indicate you want available non-reference regions

    -include_duplicates Returns duplicate regions. In order to get non-PAR regions of
                        chrY padded with N's it needs to be turned off (default)

    -padded_human_Chr_Y Returns full length human chrY with non-PAR regions padded with N's.
                        Needed for the FuncGen pipeline. Only works if -include_duplicates is disabled

    -human_female Creates a second output file that ommits chrY. Needs -onefile or -filename

    -format Deprecated. Following suggestion of the Ensembl core team, Bio::SeqIO has been replaced with
                        Bio::EnsEMBL::Utils::IO::FASTASerializer. Therefore output is limited to FASTA.

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
use Bio::EnsEMBL::PaddedSlice;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;

my $serializer = 'Bio::EnsEMBL::Utils::IO::FASTASerializer';

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my ($species, $filename, $multi_species);
my $species_id = 1;
my $coord_system_name;
my $coord_system_version;
my $output_dir;
my $single_file;
my $top_level;
my $seq_level;
my $non_ref;
my $padded_human_Chr_Y;
my $human_female;
my $include_duplicates;
my $format = 'fasta';
my @logic_names;
my $mask;
my $softmask = 0;
my $extension = 'fa';
my $help;

&GetOptions(
            'dbhost:s'                => \$host,
            'dbport:n'                => \$port,
            'dbname:s'                => \$dbname,
            'dbuser:s'                => \$dbuser,
            'dbpass:s'                => \$dbpass,
            'species=s'               => \$species,
            'multi_species'           => \$multi_species,
            'species_id=i'            => \$species_id,
            'coord_system_name:s'     => \$coord_system_name,
            'coord_system_version:s'  => \$coord_system_version,
            'output_dir:s'            => \$output_dir,
            'extension:s'             => \$extension,
            'toplevel!'               => \$top_level,
            'seqlevel!'               => \$seq_level,
            'nonref!'                 => \$non_ref,
            'padded_human_Chr_Y!'     => \$padded_human_Chr_Y,
            'human_female!'           => \$human_female,
            'include_duplicates!'     => \$include_duplicates,
            'format!'                 => \$format,
            'mask!'                   => \$mask,
            'mask_repeat:s@'          => \@logic_names,
            'softmask!'               => \$softmask,
            'onefile!'                => \$single_file,
            'filename=s'              => \$filename,
            'help!'                   => \$help,
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

########################## Sanity tests ####################################
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

if($include_duplicates && $padded_human_Chr_Y){
  my $message =
    'Retrieving padded human ChrY only works if duplicates are excluded. '.
    'Run again either not using -include_duplicates or not using '.
    '-padded_human_Chr_Y';
  throw($message);
}

if( ($mask || scalar(@logic_names > 0) ) && $padded_human_Chr_Y){
  my $message =
    'Masked padded chrY not implemented yet.';
  throw($message);
}

if($human_female && ( (!$filename) && (!$single_file) ) ){
  my $message =
    'When using -human_female, you need to pass a filename or use -onefile';
  throw($message);
}

if($human_female && $species ne 'homo_sapiens' ){
  my $message =
    "-human_female only works for homo_sapiens, not '$species'";
  throw($message);
}

if($format !~ /fasta/i ){
  my $message =
    'Following advice from the core team, the output method has been changed '.
    'from Bio::SeqIO to Bio::EnsEMBL::Utils::IO::FASTASerializer. ' .
    'Output is therefore in FASTA format only at the moment. ' .
    'FASTASerializer handles Slices more efficiently memory-wise and also ' .
    'allows changing the FASTA-header. ';
  throw($message);
}

if($top_level){
  $coord_system_name = 'toplevel';
}
if($seq_level){
  $coord_system_name = 'seqlevel';
}

###############################################################################
# All sequences go into a single file
# If a filename was passed, use that one, otherwise build it
###############################################################################
my ($singleSerializer, $fh_singleFile) = (undef, undef);
if($single_file || $filename){
  $filename ||= $output_dir.'/'.$coord_system_name.'.'.$extension;
  open($fh_singleFile, '>', $filename) or die
    "Cant open stream to $filename: $!";
  $singleSerializer = $serializer->new($fh_singleFile);
}
###############################################################################

###############################################################################
# If using a single file, removing chrY can be very memory intensive as the
# resulting file can be 3GB+ (release 66). An option has been added that
# simultaneously writtes a female FASTA file by simply ommiting chrY.
# If the filename is given and contains male, it is replaced by female
# homo_sapiens_male_GRCh37_66_37_unmasked.fasta ->
# homo_sapiens_female_GRCh37_66_37_unmasked.fasta
###############################################################################
my ($singleSerializer_female, $fh_singleFile_female) = (undef, undef);
if($human_female) {

  my $tmpName;
  if($filename){
    $tmpName = $filename;
    if($tmpName =~ /male/){$tmpName =~ s/male/female/}
    else{$tmpName .= '.female'}
  }
  else{
    $tmpName= $output_dir.'/'.$coord_system_name.'.female.'.$extension;
  }

  open($fh_singleFile_female, '>', $tmpName) or die
    "Cant open stream to $tmpName: $!";
  $singleSerializer_female = $serializer->new($fh_singleFile_female);
}
###############################################################################

###############################################################################
# get slice adaptor
# fetch slices where args are eg:
# '('toplevel','GRCh37',1,undef,undef)
# will toplevel chromosomes and the unique regions of haplotypes and Y
# but no Locus Reference Genomic (LRG) sequences
###############################################################################
my $sa = $db->get_SliceAdaptor;
my $slices =
    $sa->fetch_all(
      $coord_system_name,
      $coord_system_version,
      $non_ref,
      $include_duplicates,
      undef
      );
###############################################################################

###############################################################################
# Process slices
###############################################################################
SLICE:
foreach my $slice(@$slices){
  # Remove after testing:
#  next unless ( $slice->name() =~ /^chromosome:GRCh\d\d:Y/);

  # For FuncGen pipeline
  my ($y_header, $y_slice) = (undef, undef);

  if($padded_human_Chr_Y && $slice->name() =~ /^chromosome:GRCh\d\d:Y/ ){
    ($y_header, $y_slice) = _build_complete_PAR_padded_chrY($slice);
    next SLICE if (!$y_header);
  }

  if($mask){
    $slice = $slice->get_repeatmasked_seq(\@logic_names, $softmask);
  }

  # printing output
  # An existing filename at this stage implies that output goes
  # into a single file
  if ($filename) {
    if($y_header && $y_slice){
      _print_padded_y($singleSerializer, $y_header, $y_slice);
    }
    else{
      $singleSerializer->print_Seq($slice);
      # Write female file, if demanded
      if( $human_female && $slice->name() !~ /^chromosome:GRCh\d\d:Y/){
        $singleSerializer_female->print_Seq($slice);
      }
    }
  }
  # Write into seperate files
  else {
    my $name = $output_dir.'/'.$slice->seq_region_name.'.'.$extension;
    open(my $fh, '>', $name) or die "Cant open stream to $name: $!";
      my $multiSerial = $serializer->new($fh);
      if($y_header && $y_slice){
        _print_padded_y($multiSerial, $y_header, $y_slice);
      }
      else{
        $multiSerial->print_Seq($slice);
      }
    close($fh);
  }
}
close($fh_singleFile) if($fh_singleFile);
print "Finished\n";

=head2 _print_padded_y

  Arg [1]    : Bio::EnsEMBL::Utils::IO::FASTASerializer
  Arg [2]    : Modified header line
  Arg [3]    : Bio::EnsEMBL::Slice
  Example    : _print_padded_y($singleSerializer, $y_header, $y_slice)
  Description: Replaces the original FASTA-header from the slice with
               modified one. After writing the sequence, the original
               header is restored.
  Returntype : none
  Exceptions : none
  Caller     : main method
  Status     : at risk

=cut
sub _print_padded_y {
  my ($out, $y_header, $y_slice) = @_;

  my $original_header_function = $out->header_function();
  $out->header_function($y_header);
  $out->print_Seq($y_slice);
  $out->header_function($original_header_function);
}

=head2 _build_complete_PAR_padded_chrY

  Arg [1]    : Bio::EnsEMBL::Slice
  Example    : ($y_header, $y_slice) = _build_complete_PAR_padded_chrY($slice);
  Description: This method was solely written to create human chrY version that
               has non-PAR regions padded with N's. Therefore the script has to
               be run using
  Returntype : none
  Exceptions : none
  Caller     : main method
  Status     : at risk

=cut
sub _build_complete_PAR_padded_chrY {
  my ($slice) = @_;

  if ($species ne 'homo_sapiens'){
    my $message = "Only tested for homo_sapiens, not for '$species'";
    throw($message);
  }


  # As this method is highly specific, test that we get what we expect
  # The 1st 10,000bp are an accepted guess for the telomeric region of chrY
  # The second slice contains non-PAR regions
  if    ($slice->name() eq 'chromosome:GRCh37:Y:1:10000:1'){return 0}
  elsif ($slice->name() eq 'chromosome:GRCh37:Y:2649521:59034049:1'){
    print STDERR "Chromosome Y will have padded PAR regions\n";
    my $y_slice = Bio::EnsEMBL::PaddedSlice->new($slice);
    my $y_header = sub {
      my ($slice) = @_;
      my $original = $slice->name();
      my @header = split(/:/, $original);
      $header[3] = 1;
      $header[4] = $slice->length();
      my $tmp = join(q{:}, @header);
      my $newHeader = "Y dna:chromosome $tmp";
      return ($newHeader);
    };
    return($y_header, $y_slice);
  }
  else {
    my $message =
      'This method has been specifically written for the Ensembl FuncGen '.
      'pipeline. The slice it expect are specific for GRCh37. Compliance '.
      'with any other assembly has not been tested. The header found is '.
      "not expected: '".$slice->name()."'";
    throw($message);
  }

}
