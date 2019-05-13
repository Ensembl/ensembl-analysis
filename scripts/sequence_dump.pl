#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

    -padded_nonref Returns all non-reference sequences padded with N's to match the corresponding full length
                   of the reference chromosome for the given non-reference sequence. Only works if -include_duplicates is disabled.

    -padded_human_Chr_Y Returns full length human chrY with non-PAR regions padded with N's.
                        Needed for the FuncGen pipeline. Only works if -include_duplicates is disabled

    -human_female Creates a second output file that ommits chrY. Needs -onefile or -filename

    -format Deprecated. Following suggestion of the Ensembl core team, Bio::SeqIO has been replaced with
                        Bio::EnsEMBL::Utils::IO::FASTASerializer. Therefore the output is limited to FASTA.

    -header [default | basic | funcgen | rnaseq]
             default: chromosome:GRCh37:11:1:135006516:1 chromosome 11
             basic:   chromosome:GRCh37:11:1:135006516:1
             funcgen: 18 dna:chromosome chromosome:GRCh37:18:1:78077248:1
             rnaseq:  18


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

use warnings ;
use strict;
use feature 'say';
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use File::Spec;
use Cwd;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::PaddedSlice;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;


my $host   = '';
my $port   = '3306';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';

my $serializer = 'Bio::EnsEMBL::Utils::IO::FASTASerializer';
my $species_id  = 1;
my $header      = 'default';
my $format      = 'fasta';
my $softmask    = 0;
my $extension   = 'fa';
my $output_dir  = getcwd();

my $species;
my $filename;
my $multi_species;
my $coord_system_name;
my $coord_system_version;
my $single_file;
my $top_level;
my $seq_level;
my $non_ref;
my $padded_human_Chr_Y;
my $human_female;
my $include_duplicates;
my @logic_names;
my $mask;
my $help;
my $padded_nonref;
my $patch_only = 0;
my $alt_as_scaffolds = 0;

GetOptions( 'dbhost|host|h:s'               => \$host,
            'dbport|port|P:n'               => \$port,
            'dbname|db|D:s'               => \$dbname,
            'dbuser|user|u:s'               => \$dbuser,
            'dbpass|pass|p:s'               => \$dbpass,
            'species=s'              => \$species,
            'multi_species'          => \$multi_species,
            'species_id=i'           => \$species_id,
            'header:s'               => \$header,
            'coord_system_name|cs_name:s'    => \$coord_system_name,
            'coord_system_version|cs_version:s' => \$coord_system_version,
            'output_dir:s'           => \$output_dir,
            'extension:s'            => \$extension,
            'toplevel!'              => \$top_level,
            'seqlevel!'              => \$seq_level,
            'nonref!'                => \$non_ref,
            'padded_nonref!'         => \$padded_nonref,
            'padded_human_Chr_Y!'    => \$padded_human_Chr_Y,
            'human_female!'          => \$human_female,
            'include_duplicates!'    => \$include_duplicates,
            'format!'                => \$format,
            'mask!'                  => \$mask,
            'mask_repeat:s@'         => \@logic_names,
            'softmask!'              => \$softmask,
            'onefile!'               => \$single_file,
            'filename=s'             => \$filename,
            'patch_only!'            => \$patch_only,
            'alt_as_scaffolds!'      => \$alt_as_scaffolds,
            'help!'                  => \$help,
) or ( $help = 1 );


if ($help) {
    exec('perldoc', $0);
}



############################# Sanity tests ################################
if(!$host || !$dbname || !$dbuser){
  my $message =
    "Need -dbhost '$host' -dbuser '$dbuser' and -dbname '$dbname' to run. ".
    'Use -help for more detailed documentation.';
  throw($message);
}
if(!$coord_system_name && !$top_level && !$seq_level){
  my $message =
  'Must specify either -coord_system_name, -toplevel, or -seqlevel to run. '.
  'Use -help for more detailed documentation.';
  throw($message);
}

if($top_level && $seq_level){
  my $message =
  'Cannot  specify both -toplevel and -seqlevel must be one or the other. '.
  'Use -help for more detailed documentation.';
  throw($message);
}

if(! $filename && (!$output_dir || ! -e $output_dir)){
  my $message =
  "Cannot dump sequence into '$output_dir' it does not exist ".
  'Use -help for more detailed documentation.';
  throw($message);
}

if($include_duplicates && $padded_human_Chr_Y){
  my $message =
    'Retrieving padded human ChrY only works if duplicates are excluded. '.
    'Run again either not using -include_duplicates or not using '.
    '-padded_human_Chr_Y';
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

if($format !~ /^fasta$/i ){
  my $message =
    'Following advice from the core team, the output method has been changed '.
    'from Bio::SeqIO to Bio::EnsEMBL::Utils::IO::FASTASerializer. ' .
    'Output is therefore in FASTA format only at the moment. ' .
    'FASTASerializer handles Slices more efficiently memory-wise and also ' .
    'allows changing the FASTA-header. ';
  throw($message);
}
############################# Global settings ###############################
if($top_level){
  $coord_system_name = 'toplevel';
}
if($seq_level){
  $coord_system_name = 'seqlevel';
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -dbname          => $dbname,
   -host            => $host,
   -user            => $dbuser,
   -port            => $port,
   -pass            => $dbpass,
   -species         => $species,
   -multispecies_db => $multi_species,
   -species_id      => $species_id
);

################################################################################
# All sequences go into a single file
# If a filename was passed, use that one, otherwise build it
################################################################################
my ($singleSerializer, $fh_singleFile) = (undef, undef);
if($single_file || $filename){
  $filename ||= File::Spec->catfile($output_dir , "$coord_system_name.$extension");
  open($fh_singleFile, '>', $filename) or die
    "Cannot open stream to $filename: $!";
  $singleSerializer = $serializer->new($fh_singleFile);
}
################################################################################

################################################################################
# If using a single file, removing chrY can be very memory intensive as the
# resulting file can be 3GB+ (release 66). An option has been added that
# simultaneously writtes a female FASTA file by simply ommiting chrY.
# If the filename is given and contains male, it is replaced by female
# homo_sapiens_male_GRCh37_66_37_unmasked.fasta ->
# homo_sapiens_female_GRCh37_66_37_unmasked.fasta
################################################################################
my ($singleSerializer_female, $fh_singleFile_female) = (undef, undef);
if($human_female) {

  my $tmpName;
  if($filename){
    $tmpName = $filename;
    if($tmpName =~ /male/){$tmpName =~ s/male/female/}
    else{$tmpName .= '.female'}
  }
  else{
    $tmpName = File::Spec->catfile($output_dir, "$coord_system_name.female.$extension");
  }

  open($fh_singleFile_female, '>', $tmpName) or die
    "Cant open stream to $tmpName: $!";
  $singleSerializer_female = $serializer->new($fh_singleFile_female);
}
################################################################################

################################################################################
# get slice adaptor
# fetch slices where args are eg:
# '('toplevel','GRCh37',1,undef,undef)
# will toplevel chromosomes and the unique regions of haplotypes and Y
# but no Locus Reference Genomic (LRG) sequences
###############################################################################
my $sa = $db->get_SliceAdaptor;
my $slices;
if ($alt_as_scaffolds) {
  $slices = fetch_alts_as_scaffolds($db,$sa,$coord_system_name,$coord_system_version);
} else {
  $slices = $sa->fetch_all(
                            $coord_system_name,
                            $coord_system_version,
                            $non_ref,
                            $include_duplicates,
                            undef,
                          );
}
################################################################################

################################################################################
# Dispatch table
# Used for creating project-specific headers
# -default:
#   chromosome:GRCh37:11:1:135006516:1 chromosome 11
# -basic:
#   chromosome:GRCh37:11:1:135006516:1
# -funcgen
#   18 dna:chromosome chromosome:GRCh37:18:1:78077248:1
# -rnaseq
#   18
################################################################################
my $dispatch = {
  default => sub {
    my ($slice) = @_;
    return sprintf(
        '%s %s %s',
        $slice->name(),
        $slice->coord_system_name(),
        $slice->seq_region_name(),
        );
  },

  basic   => sub {
    my ($slice) = @_;
    return $slice->name();
  },

  funcgen => sub {
    my ($slice) = @_;
    return sprintf (
      '%s %s %s',
      $slice->seq_region_name(),
      $slice->moltype().':'.$slice->coord_system_name(),
      $slice->name(),
      );
  },
  rnaseq => sub {
    my ($slice) = @_;
    return sprintf (
      '%s',
      $slice->seq_region_name(),
      );
  },
};
# Sanity check for header passed
if(not exists $dispatch->{$header}){
  my $message =
    "'$header' not defined. Please adjust dispatch table.\nHeaders currently ".
    "available: " . join (', ', keys (%{$dispatch}));
  throw($message);
}
################################################################################

################################################################################
# Process slices
################################################################################
SLICE:
foreach my $slice(@$slices){

  if($patch_only) {
    my $is_patch = 0;
    my @pt = ('patch_novel','patch_fix');
    foreach my $type (@pt) {
      my @slice_attributes = @{$slice->get_all_Attributes($type)};
      if (scalar(@slice_attributes) > 0) {
        $is_patch = 1;
        last;
      }
    }

    unless($is_patch) {
      next;
    }
  }

  # Compliance with header format used in previous version of this script
  $singleSerializer->header_function($dispatch->{$header}) if($filename);

  # For FuncGen pipeline, print a PAR-padded chrY
  my ($padded_header, $padded_slice) = (undef, undef);
  if($padded_human_Chr_Y && $slice->name() =~ /^chromosome:GRCh\d\d:Y/ ){
    ($padded_header, $padded_slice) = _build_complete_PAR_padded_chrY($slice);
    next SLICE if (!$padded_header);
  }

  if($mask){
    $slice = $slice->get_repeatmasked_seq(\@logic_names, $softmask);
  }

  if ($padded_nonref and (!($slice->is_reference()))) {
    $padded_slice = Bio::EnsEMBL::PaddedSlice->new($slice);
    $padded_header = sub {
      my ($padded_slice) = @_;
      my $original = $padded_slice->name();
      my @header = split(/:/, $original);
      $header[3] = 1;
      $header[4] = $padded_slice->length();
      my $tmp = join(q{:}, @header);
      my $newHeader = "$tmp";
      return ($newHeader);
    };
  }

  # printing output
  # An existing filename at this stage implies that output goes
  # into a single file
  if ($filename) {
    if ($padded_header && $padded_slice) {
      _print_padded($singleSerializer,$padded_header,$padded_slice);
    } else {
      $singleSerializer->print_Seq($slice);
      # Write female file, if demanded
      if( $human_female && $slice->name() !~ /^chromosome:GRCh\d\d:Y/){
        $singleSerializer_female->print_Seq($slice);
      }
    }
  }
  # Write into separate files
  else {
    my $name = File::Spec->catfile($output_dir, $slice->seq_region_name.'.'.$extension);
    print "Multi: $name\n";
    open(my $fh, '>', $name) or die "Cant open stream to $name: $!";
      my $multiSerializer = $serializer->new($fh);
      if ($padded_header && $padded_slice) {
        _print_padded($multiSerializer,$padded_header,$padded_slice);
      } else {
        $multiSerializer->header_function($dispatch->{$header});
        $multiSerializer->print_Seq($slice);
      }
    close($fh);
  }
}
close($fh_singleFile) if($fh_singleFile);
print "Finished\n";


=head2 _print_padded

  Arg [1]    : Bio::EnsEMBL::Utils::IO::FASTASerializer
  Arg [2]    : Modified header line
  Arg [3]    : Bio::EnsEMBL::Slice
  Example    : _print_padded($singleSerializer, $y_header, $y_slice)
  Description: Replaces the original FASTA-header from the slice with
               modified one. After writing the sequence, the original
               header is restored.
  Returntype : none
  Exceptions : none
  Caller     : main method
  Status     : at risk

=cut
sub _print_padded {
  my ($out, $padded_header, $padded_slice) = @_;

  my $original_header_function = $out->header_function();
  $out->header_function($padded_header);
  $out->print_Seq($padded_slice);
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
  if    ( ($slice->name() eq 'chromosome:GRCh38:Y:1:10000:1') or ($slice->name() eq 'chromosome:GRCh38:Y:57217416:57227415:1') ){return 0}
  elsif ($slice->name() eq 'chromosome:GRCh38:Y:2781480:56887902:1'){
    print STDERR "Chromosome Y will have padded PAR regions\n";
    if($mask){
    $slice = $slice->get_repeatmasked_seq(\@logic_names, $softmask);
    }
    my $y_slice = Bio::EnsEMBL::PaddedSlice->new($slice);
    my $y_header = sub {
      my ($slice) = @_;
      my $original = $slice->name();
      my @header = split(/:/, $original);
      $header[3] = 1;
      $header[4] = $slice->length();
      my $tmp = join(q{:}, @header);
      #my $newHeader = "Y dna:chromosome $tmp";
      my $newHeader = "$tmp";
      return ($newHeader);
    };
    return($y_header, $y_slice);
  }
  else {
    my $message =
      'This method has been specifically written for the Ensembl FuncGen '.
      'pipeline. The slice it expect are specific for GRCh38. Compliance '.
      'with any other assembly has not been tested. The header found is '.
      "not expected: '".$slice->name()."'";
    throw($message);
  }

}


=head2 fetch_alts_as_scaffolds

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [2]    : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Arg [3]    : Coord system name
  Arg [4]    : Coord system version
  Example    : $slices = fetch_alts_as_scaffolds($slice_adaptor,$coord_system_name,$coord_system_version);
  Description: We are likely moving to a system where we no longer have the patches represented as pseudochromosomes and instead
               they'll just be scaffolds. This would be good for many reasons. Until then this method will allow for dumping of
               the sequence in a way that will mimic our future plans by replacing toplevel alts with the underlying scaffold
  Returntype : Listref of Bio::EnsEMBL::Slice
  Exceptions : die if the slice is not found
  Caller     : main method
  Status     : at risk

=cut
sub fetch_alts_as_scaffolds {
  my ($db,$sa,$coord_system_name,$coord_system_version) = @_;

  my  $slices = $sa->fetch_all(
                                $coord_system_name,
                                $coord_system_version,
                                0,
                              );

  # This query is really only designed for human, but it will probably work on other thing. This will homefully be removed when we switch to representing alts as scaffolds
  my $query = "select distinct(seq_region_id) from seq_region where coord_system_id=(select coord_system_id from coord_system where name='scaffold' and version='".$coord_system_version.
              "') and seq_region_id in (select cmp_seq_region_id from assembly where asm_seq_region_id in (select seq_region_id from assembly_exception where exc_type != 'PAR'));";

  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  while (my $seq_region_id = $sth->fetchrow_array) {
    my $slice = $sa->fetch_by_seq_region_id($seq_region_id);
    unless($slice) {
      die "Could not find slice for the following seq_region_id: ".$seq_region_id;
    }
    push(@{$slices},$slice);
  }
  return($slices);
}
