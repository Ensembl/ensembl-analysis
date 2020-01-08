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

#
# Adds GRC genomic mapping to the dna align feature table
# in the HAP pipeline with their cigar strings. Note that
# in some cases there is more than one alignment per patch.
#
# Example:
#
# perl add_GRC_align_features.pl -dbhost genebuildn \
#      -dbname homo_sapiens_core_nn_nn -dbuser user -dbpass pass \
#      -gca_patch_release GCA_000001405.20_GRCh38.p5 -verbose

use strict;
use warnings;

use Getopt::Long;
use Net::FTP;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis;

$| = 1;

my $dbname         = '';
my $dbhost         = '';
my $dbuser         = 'ensro';
my $dbpass         = '';
my $dbport         = '3306';
my $gca_patch_release  = '';
my $store          = 0;
my $external_db_id = '50692'; # GRC_alignment_import
my $syn_external_db_id = '50710'; # seq_region_synonym slice type - i.e. INSDC

my $verbose        = 0; 

my @patch_types = ('PATCH_FIX','PATCH_NOVEL','HAP');
my @ftpdir_types = ('ALT_REF_LOCI_','PATCHES'); # all ALT_REF_LOCI_* will be included
my @dna_align_features = ();

&GetOptions( 
  'dbhost:s'                 => \$dbhost,
  'dbuser:s'                 => \$dbuser,
  'dbpass:s'                 => \$dbpass,
  'dbname:s'                 => \$dbname,
  'dbport:n'                 => \$dbport,
  'gca_patch_release:s'      => \$gca_patch_release,
  'external_db_id:n'         => \$external_db_id,
  'write!'                   => \$store,
  'verbose!'                 => \$verbose,
);

if(!$gca_patch_release){
  throw ("Need to specify gca accession and assembly version with -gca_patch_release.\n");
}

# get alt_scaffold_placement.txt to generate filename that we will need
# to retrieve files from the NCBI ftp site. Populates...

my %ftp_filename=();                 

my ($content, $remote_file_handle) = "";
open($remote_file_handle, '>', \$content);

my $ftp = Net::FTP->new('ftp.ncbi.nlm.nih.gov', Debug => 0)
  or die "Can't connect to NCBI FTP: $@";

$ftp->login('anonymous', '-anonymous@')
  or die 'Cannot login ', $ftp->message;

chomp $gca_patch_release;

my $ncbi_patch_release_wd = "/genomes/genbank/vertebrate_mammalian/Mus_musculus/all_assembly_versions/".$gca_patch_release."/".$gca_patch_release."_assembly_structure";
my $ncbi_wd = $ncbi_patch_release_wd;
$ftp->cwd($ncbi_wd);

# get the list of ftp dirs which contain patches and haplotypes data
my @patches_ftpdirs = ();
my @ftpdirs = $ftp->ls();
foreach my $ftpdir (@ftpdirs) {
  foreach my $ftpdir_type (@ftpdir_types) {
    if ($ftpdir =~ m/$ftpdir_type/) {
      push(@patches_ftpdirs,$ftpdir);
    }
  }
}

my %align_str = ();
foreach my $patches_ftpdir (@patches_ftpdirs) {

  %ftp_filename = ();

  print "---Processing directory $patches_ftpdir\n";

  my $ncbi_wd = $ncbi_patch_release_wd."/".$patches_ftpdir."/alt_scaffolds/";

  $ftp->cwd($ncbi_wd)
    or die 'Cannot change working directory ', $ftp->message;

  close $remote_file_handle;
  open($remote_file_handle, '>', \$content);
  $ftp->get('alt_scaffold_placement.txt', $remote_file_handle)
    or die "get failed ", $ftp->message;

  my @asp_lines = split /\n/, $content;

  foreach my $asp_line (@asp_lines) {
    next if $asp_line =~ /^#/;
    my @elem = split /\t/, $asp_line;
    my $patch_name = $elem[2];
    my $file_name = $elem[3]."_".$elem[6].".gff";
    print "Filename:  $file_name\t\tPatchname:  $patch_name\n" if $verbose;
    $ftp_filename{$patch_name} = $file_name;
  }

  # change directory to where the GRC alignments are kept:

  $ncbi_wd = "alignments"; 
  $ftp->cwd($ncbi_wd) or die 'Cannot change working directory ', $ftp->message;

  # hash of arrays - there way be more than one alignment per file if they
  # have been manually annotated, However the GRC may change all to one
  # line in the near future, in the meantime, we need to deal with them.

  foreach my $patch (keys %ftp_filename) {
    close $remote_file_handle;
    open($remote_file_handle, '>', \$content);
    $ftp->get($ftp_filename{$patch}, $remote_file_handle)
      or die "get failed ", $ftp->message;

    my @lines = split "\n", $content;
    foreach my $line (@lines) {
      next if $line =~ /^\#/;
      # We'll parse the data later because we need most of it.
      push @{$align_str{$patch}},$line;

      # In GRCh37, the HAP names were shortened like HSCHR17_1 instead of HSCHR17_1_CTG5
      # In GRCh38 and GRCm38, the HAP and PATCHES names were extended like CHR_HSCHR17_1 instead of HSCHR17_1
      # so I'll add a 'duplicated' line associated to the new name too
      # so that when the HAP and PATCHES names are fetched from our DB, there can be a match
      my $new_hap_name = "CHR_".$patch;
      push @{$align_str{$new_hap_name}},$line;
    }
  }
} # endif patches_ftpdir


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dbhost,
                                             -user   => $dbuser,
                                             -pass   => $dbpass,
                                             -port   => $dbport,
                                             -dbname => $dbname );

my $sa = $db->get_SliceAdaptor();

my $analysis = new Bio::EnsEMBL::Analysis( -logic_name => "grc_alignment_import",
                                           -db_version => $gca_patch_release);



# TODO - leave $asm_exc_adaptor in - that way we can compare
# information from file with what we already know as a sanity check.


# Now get the patches, they come in pairs, the assembly exception and the reference
print "Getting patches...\n" if $verbose;
my $asm_exc_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();
my @exceptions = @{$asm_exc_adaptor->fetch_all()};
my @patches;
EXC: foreach my $exc (@exceptions){
  foreach my $type (@patch_types){
    if($exc->type() =~ m/$type/){
      push(@patches, $exc);
      next EXC;
    }
  }
}
# Assuming that AssemblyExceptionFeatureAdaptor's fetch_all will always
# return 2 entries for each patch and that those two entries are adjacent
my $num_patches = scalar(@patches)/2;

print "Have ".$num_patches." patches.\n";

# for each patch
for (my $i = 0; $i < $num_patches; $i++) {
  # get the two slices
  my $ref_slice;
  my $patch_slice;
 
  for(my $j = 0; $j < 2; $j++) {
    my $exc = pop(@patches);
    # if this is the ref version
    if($exc->type =~ m/REF/){
      # alt is only the patch slice
      $patch_slice = $exc->alternate_slice();
    }
    else{
      # alt is replaced region of ref
      $ref_slice = $exc->alternate_slice();
    }
  }
  if(!($patch_slice and $ref_slice)){
    throw("Something is wrong, the patch and ref slices were not set correctly.\n");
  }

  my @patch_vals = split /:/, $patch_slice->display_id; 
  my $patch_name = $patch_vals[2]; 
  foreach my $string ( @{ $align_str{$patch_name}}) {
    my @el = split /\t/, $string;
    
    my $num = $#el;
    throw ("Incorrect number of elements in gtf file: $num") unless $num == 8;
 
    my ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attr) = split /\t/, $string;
    
    $strand = fix_strand($strand);

    my %attribute = ();
    foreach my $kvp (split ';', $attr) {
      my ($key, $value) = split '=', $kvp;
      $attribute{$key} = $value;
    }

    my $target = $attribute{"Target"};
    my ($hseqname, $hstart, $hend, $hstrand ) = split " ", $target;
    
    $hstrand = fix_strand($hstrand); 
    
    my $length = ($hend - $hstart) + 1;    

    my $cigar_line;    
    $cigar_line = $attribute{"Gap"};
    if (defined $cigar_line) {
      sanity_check_cigar_line($cigar_line, $length);
      $cigar_line = reformat_cigar_line($cigar_line);
    } else {
      $cigar_line = $length."M";
    }
print $cigar_line."\n";
                      
    # need the seq_region_id from seq_region_synonym
    my @synonyms = @{$ref_slice->get_all_synonyms()};
    my $seq_region_id = '';
    foreach my $syn (@synonyms) {
      if ($syn->external_db_id() == $syn_external_db_id) {
        $seq_region_id = $syn->seq_region_id();
        last();
      }
    }
print "about to print seq_region_id\n";
print "seq_region_id is $seq_region_id\n";
    # ...to obtain the slice:

    my $slice = $sa->fetch_by_seq_region_id($seq_region_id);


    my $daf = new Bio::EnsEMBL::DnaDnaAlignFeature(
      -slice          => $slice,
      -start          => $start,
      -end            => $end,
      -strand         => $strand,
      -analysis       => $analysis,                        
      -score          => $score,
      -hstart         => $hstart,
      -hend           => $hend,
      -hstrand        => $hstrand,
      -hseqname       => $hseqname,
      -hcoverage      => $attribute{"pct_coverage"},
      -percent_id     => $attribute{"pct_identity_ungap"}, 
      -external_db_id => $external_db_id,
      -cigar_string   => $cigar_line,
    );

    push @dna_align_features, $daf;
  }
}


# now store all the dna_align features
if (scalar(@dna_align_features) > 0) {
  write_dna_align_features_to_db($db,\@dna_align_features,$store)
}
print "There are ".scalar (@dna_align_features)." new dna_align_features.\n";


sub write_dna_align_features_to_db {
  my ($db,$dna_align_features,$store) = @_;

  DAF: foreach my $dna_align_feat (@$dna_align_features) {
    if ($store) {
      $db->get_DnaAlignFeatureAdaptor->store($dna_align_feat);
      if ($@) {
        throw("ERROR: Can't write dna_align_feat ".$dna_align_feat->hseqname." [$@]");
      }  else {
        print "Written ".$dna_align_feat->hseqname." on chr ".$dna_align_feat->slice->name
              ." strand ".$dna_align_feat->hstrand." with start ".$dna_align_feat->start
              ." end ".$dna_align_feat->end."\n" if $verbose;
      }
    } else {
      print "Not storing ".$dna_align_feat->hseqname."\n" if $verbose;
    }
  } # DAF 
  return 1;  
}

sub fix_strand { 
  my $strand = shift;
  if ($strand eq '+') {
    $strand = 1;
  } elsif ($strand eq '-') {
    $strand = -1;
  } else {
    throw("Strand problem :".$strand);
  }
  return $strand;
}
 
sub sanity_check_cigar_line {
  # ok, it sanity checks the GRCs idea of a cigar line which is close to GFF3 format
  my ($line, $len) = @_;
  my $cl_length = 0;
  throw("Can only sanity check cigar lines with whitespace") unless $line =~ /\s/;
  my @elements = split /\s/, $line;
  foreach my $el (@elements) {
    my ($operator, $num) = ($1, $2) if $el =~ /^(\w{1})(\d+)$/;
    if ($operator =~ /[MI]/) {
      $cl_length += $num;
    } elsif ($operator eq 'D') {
      # nothing to do
    } else {
      throw("Unknown alignment operator: $operator acting on $num");
    }
  }
  if ($cl_length != $len) {
    warn("Cigar_line length: $cl_length does not match length: $len for this line:\n$line\n\n");
  }
}

sub reformat_cigar_line {
  my $line = shift;
  # hack
  # the GRC cigar line format turns out to be back to front - fix it
  # this is a retrospective hack, with hindsight the logic of the script
  # would be different and probably incorporated into the sub above.
  my @elements = split /\s/, $line;  
  $line = ''; 
  foreach my $el (@elements) { 
    my ($operator, $num) = ($1, $2) if $el =~ /^(\w{1})(\d+)$/;   
    $line .= $num.$operator;
  }
  return $line;
}


exit;
