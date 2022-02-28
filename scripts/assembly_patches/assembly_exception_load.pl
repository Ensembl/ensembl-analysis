#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Getopt::Long;

use strict;
use warnings;

my $assembly_acc;
my $assembly_name;
my $pass;
my $mapping_file   = "./data/alt.scaf.agp";
my $txt_file       = "./data/alt_scaffold_placement.txt";
my $patchtype_file = "./data/patch_type";
my $sql_file       = "./sql.txt";
my $dbname;
my $host;
my $user;
my $port;
my $scaf_syn_ext_id; # INSDC = 50710 but it's check by an SQL query
my $ensembl_syn_ext_id; # ensembl_internal_synonym = 50803 but it's check by an SQL query
my $central_coord_system;
my $chromosome_prefix = 'CHR_';

GetOptions(
            'pass=s'            => \$pass,
            'mapping_file=s'    => \$mapping_file,
            'txt_file=s'        => \$txt_file,
            'patchtype_file=s'  => \$patchtype_file,
            'sql_file=s'        => \$sql_file,
            'host=s'            => \$host,
            'dbname=s'          => \$dbname,
            'user=s'            => \$user,
            'port=n'            => \$port,
            'scaf_syn_ext_id=n' => \$scaf_syn_ext_id,
            'ensembl_syn_ext_id=n' => \$ensembl_syn_ext_id,
            'assembly_name=s'   => \$assembly_name,
            'assembly_acc=s'    => \$assembly_acc,
            'coord_system=s'    => \$central_coord_system
           );

#needed to update the meta table
if (!defined $assembly_name)  {
    throw("Please enter -assembly_name eg. GRCh37.p4");
} 

if (!defined $assembly_acc)  {
    throw("Please enter -assembly_acc eg. GCA_000001405.4");
}

if (!defined $central_coord_system)  {
    throw("Please enter -coord_system to specify the central coord system eg. supercontig");
}


#connect to the database

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    '-host'    => $host,
    '-port'    => $port,
    '-user'    => $user,
    '-pass'    => $pass,
    '-dbname'  => $dbname,
    '-species' => "load"
    );


my $sth = $dba->dbc->prepare("select max(seq_region_id) from seq_region")
  || die "Could not get max seq_region_id";

$sth->execute || die "problem executing get max seq_region_id";
my $max_seq_region_id;
$sth->bind_columns(\$max_seq_region_id) || die "problem binding";
$sth->fetch() || die "problem fetching";
$sth->finish;

if (! $scaf_syn_ext_id) {
    $sth = $dba->dbc->prepare("SELECT external_db_id FROM external_db WHERE db_name = 'INSDC'")
        || die "Could not get external_db_id for INSDC";
    $sth->execute || die "problem executing INSDC query";
    $sth->bind_columns(\$scaf_syn_ext_id) || die "problem binding INSDC query";
    $sth->fetch() || die "problem fetching INSDC query";
    $sth->finish;
}

if (! $ensembl_syn_ext_id) {
    $sth = $dba->dbc->prepare("SELECT external_db_id FROM external_db WHERE db_name = 'ensembl_internal_synonym'")
        || die "Could not get external_db_id for INSDC";
    $sth->execute || die "problem executing INSDC query";
    $sth->bind_columns(\$ensembl_syn_ext_id) || die "problem binding INSDC query";
    $sth->fetch() || die "problem fetching INSDC query";
    $sth->finish;
}

$sth = $dba->dbc->prepare("select count(assembly_exception_id) from assembly_exception")
  || die "Could not get number of rows in assembly_exception";
my $count_assembly_exception_id;
$sth->execute || die "problem executing";
$sth->bind_columns(\$count_assembly_exception_id) || die "problem binding";
$sth->fetch() || die "problem fetching";
$sth->finish;

my $max_assembly_exception_id;
if (defined $count_assembly_exception_id && $count_assembly_exception_id > 0) {
  $sth = $dba->dbc->prepare("select max(assembly_exception_id) from assembly_exception")
    || die "Could not get max assembly_exception_id";

  $sth->execute || die "problem executing get max assembly_exception_id";
  $sth->bind_columns(\$max_assembly_exception_id) || die "problem binding";
  $sth->fetch() || die "problem fetching";
  $sth->finish;
} else {
  warning("The assembly_exception table is empty. This is OK if it's the first ever patch release for a species, otherwise there is a problem.");
  $max_assembly_exception_id = 0;
}

print "starting new seq_region at seq_region_id of $max_seq_region_id\n";
print "\nTo reset\nDELETE FROM dna WHERE seq_region_id > $max_seq_region_id\nDELETE FROM seq_region WHERE seq_region_id > $max_seq_region_id\n";
print "DELETE FROM srs USING seq_region_synonym srs LEFT JOIN seq_region sr ON srs.seq_region_id = sr.seq_region_id WHERE sr.seq_region_id IS NULL\nDELETE FROM sra USING seq_region_attrib sra LEFT JOIN seq_region sr ON sra.seq_region_id = sr.seq_region_id WHERE sr.seq_region_id IS NULL\nDELETE FROM a USING assembly a LEFT JOIN seq_region sr ON a.asm_seq_region_id = sr.seq_region_id WHERE sr.seq_region_id IS NULL\nDELETE FROM a USING assembly a LEFT JOIN seq_region sr ON a.cmp_seq_region_id = sr.seq_region_id WHERE sr.seq_region_id IS NULL\n";
if (defined $count_assembly_exception_id && $count_assembly_exception_id > 0) {
  print "DELETE FROM assembly_exception WHERE assembly_exception_id > $max_assembly_exception_id\n\n";
}

$max_seq_region_id++;
$max_assembly_exception_id++;



my ($chr_coord_id, $contig_coord_id);

$sth = $dba->dbc->prepare('select coord_system_id from coord_system where name like ? and attrib like "%default_version%"')
  || die "Could not prepare coord_system sql";

$sth->execute("chromosome") || die "Could not execute coord_system sql";
$sth->bind_columns(\$chr_coord_id) || die "What a bind NOT";
$sth->fetch() || die "Could not fetch coord_system_id";

$sth->execute("contig") || die "Could not execute coord_system sql";
$sth->bind_columns(\$contig_coord_id) || die "What a bind NOT";
$sth->fetch() || die "Could not fetch coord_system_id";
$sth->finish;

my $get_code_sth = $dba->dbc->prepare("select attrib_type_id from attrib_type where code like ?");

my $get_seq_region_sth = $dba->dbc->prepare("select seq_region_id, length from seq_region where name like ?")
  || die "Could not prepare get_seq_region_sth";

my $get_seq_region_id_sth = $dba->dbc->prepare("select seq_region_id, length from seq_region where coord_system_id = ? and name like ?")
  || die "Could not prepare get_seq_region_id_sth";

my $get_seq_id_sth = $dba->dbc->prepare("select sr.seq_region_id from seq_region sr , coord_system cs where sr.name like ? and cs.coord_system_id = sr.coord_system_id and cs.attrib like '%sequence%'");

my $coord_sth = $dba->dbc->prepare("select coord_system_id from coord_system where name like ? order by rank limit 1");

if(!defined($chr_coord_id) or !defined($contig_coord_id)){
  die "No coord_system_id for chromosome ($chr_coord_id) or contig ($contig_coord_id)\n";
}

my $non_ref;
my $toplevel;
my $patch_novel;
my $patch_fix;

$get_code_sth->execute("toplevel");
$get_code_sth->bind_columns(\$toplevel);
$get_code_sth->fetch;

$get_code_sth->execute("non_ref");
$get_code_sth->bind_columns(\$non_ref);
$get_code_sth->fetch;

$get_code_sth->execute("patch_novel");
$get_code_sth->bind_columns(\$patch_novel);
$get_code_sth->fetch;

$get_code_sth->execute("patch_fix");
$get_code_sth->bind_columns(\$patch_fix);
$get_code_sth->fetch;

print "toplevel $toplevel, non_ref $non_ref, patch_novel $patch_novel, patch_fix $patch_fix\n";


#      Store contigs as seq_region + dna.

open(MAPPER,"<".$mapping_file)  || die "Could not open file $mapping_file";
open(TXT,"<".$txt_file)         || die "Could not open file $txt_file";
open(TYPE,"<".$patchtype_file)         || die "Could not open file $patchtype_file";

open(SQL,">".$sql_file) || die "Could not open $sql_file for writing\n";

my %acc_to_name;
my %txt_data; 
my %key_to_index;
my %name_to_seq_id;
my %name_to_sc_seq_id;
my %seq_id_to_start;
my %name_to_type;
my %seq_id_to_alt_scaf_start;
my %seq_id_to_alt_scaf_stop;
my %seq_id_to_stop;
my %seq_id_to_strand;
my %old_patch_names;

while (<TYPE>) {
  chomp;
  next if (/^#/);
  my ($alt_scaf_name,$alt_scaf_acc,$type) = split(/\t/,$_);
  if (!$alt_scaf_name || !$alt_scaf_acc || !$type) {
    throw("Unable to find name, accession or type"); 
  }
  $name_to_type{$alt_scaf_name} = $type;
}


#alt_asm_name	prim_asm_name	alt_scaf_name	alt_scaf_acc	parent_type	parent_name	parent_acc	region_name	ori	alt_scaf_start	alt_scaf_stop	parent_start	parent_stop	alt_start_tail	alt_stop_tail
#PATCHES	Primary Assembly	HG104_HG975_PATCH	GL383535.1	CHROMOSOME	8	CM000670.1	EPPK1_SPATC1	+	1	429806	144743526	145146062	0	0
#PATCHES	Primary Assembly	HG991_PATCH	GL383525.1	CHROMOSOME	3	CM000665.1	SLC25A26	+	27267	65063	66270271	66308065	27266	0
#PATCHES	Primary Assembly	HG987_PATCH	GL383561.1	CHROMOSOME	17	CM000679.1	REGION27	+	1	314281	21250948	21566608	0	92682
#PATCHES	Primary Assembly	HSCHR12_1_CTG2	GL383549.1	CHROMOSOME	12	CM000674.1	REGION21	-	1	120804	28148967	28263711	0	0

# See http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml

SCAF: while(<TXT>){
  if(/^#/){
    chomp;
    my @arr = split(/\t/,$_);
    # please note: this split was originally on whitespace and not tabs for human. 
    # the $i was then set to 1 instead of zero as a workaround because the human assembly name was 2 words
    # however this workaround broke for mouse where the assembly name is one word without whitespace
    my $i = 0; 
    foreach my $name (@arr){
      $key_to_index{$name} = $i;
      $i++;
    }
    foreach my $name (qw(alt_scaf_name alt_scaf_acc parent_type parent_name ori alt_scaf_start alt_scaf_stop parent_start parent_stop alt_start_tail alt_stop_tail)){
      if(!defined($key_to_index{$name})){
        print STDERR "PROBLEM could not find index for $name\n";
      }
      else{
        print $name."\t".$key_to_index{$name}."\n";
      }
    }
  }
  else{
    my @arr = split(/\t/,$_);
    my $alt_acc = $arr[$key_to_index{'alt_scaf_acc'}];
    my $alt_name = $arr[$key_to_index{'alt_scaf_name'}];

    my ($tmp_id,$tmp_len);
    $get_seq_region_sth->execute($alt_name)
      || die "Could not execute get seq_region for $alt_name";
    $get_seq_region_sth->bind_columns(\$tmp_id,\$tmp_len)
      || die "Could not bind get_seq_region";
    $get_seq_region_sth->fetch();
#     || die "Could not fetch seq_region id for name $name";;

    #go to next
    if(defined($tmp_id)){

      print $alt_name." already in db\n";
      $old_patch_names{$alt_name} = 1;
      next SCAF;
    }
    else{
      print $alt_name." is new\n";
    }

    my $coord_sys;
    $coord_sth->execute($arr[$key_to_index{'parent_type'}]);
    $coord_sth->bind_columns(\$coord_sys);
    $coord_sth->fetch();

    if(!defined($coord_sys)){
      die "Could not get coord_system_id for ".$arr[$key_to_index{'parent_type'}]." col index ".$key_to_index{'parent_type'}."\n";
    }

    my ($seq_region_id, $seq_length);
    $get_seq_region_id_sth->execute($coord_sys, $arr[$key_to_index{'parent_name'}]);
    $get_seq_region_id_sth->bind_columns(\$seq_region_id,\$seq_length);
    $get_seq_region_id_sth->fetch();
    if(!defined($seq_region_id)){
      die "Could not get seq_region_id for ".$arr[$key_to_index{'parent_name'}]. " on coord system $coord_sys";
    }

    $acc_to_name{$alt_acc} = $alt_name;
    $name_to_seq_id{$alt_name} = $max_seq_region_id;
    $name_to_seq_id{$alt_acc}  = $max_seq_region_id;
    $seq_id_to_start{$max_seq_region_id} = $arr[$key_to_index{'parent_start'}];

    if ($arr[$key_to_index{'ori'}] eq "+"){
      $seq_id_to_strand{$max_seq_region_id} = 1;
    }
    elsif($arr[$key_to_index{'ori'}] eq "-"){
      $seq_id_to_strand{$max_seq_region_id} = -1;
    }
    elsif($arr[$key_to_index{'ori'}] eq "b" ){
      $seq_id_to_strand{$max_seq_region_id} = 1;
      print "WARNING: Encountered ori = b for $alt_name - defaulting to ori = 1 (although this may be totally incorrect!)\n";
    }
    else{
      print "Problem with strand.\n";
    }

    my $exc_len    = $arr[$key_to_index{'parent_stop'}] - $arr[$key_to_index{'parent_start'}] + 1;
    my $new_length = $arr[$key_to_index{'alt_scaf_stop'}] - $arr[$key_to_index{'alt_scaf_start'}] + 1 +
                     $arr[$key_to_index{alt_start_tail}] + $arr[$key_to_index{alt_stop_tail}];

    my $length = $seq_length + ($new_length-$exc_len);
    $seq_id_to_stop{$max_seq_region_id} = $arr[$key_to_index{'parent_start'}] + $new_length -1;

    print SQL "insert into seq_region (seq_region_id, name, coord_system_id, length)\n";
    print SQL "\tvalues ($max_seq_region_id, '$chromosome_prefix$alt_name', $coord_sys, $length);\n\n";

    print SQL "insert into seq_region_synonym (seq_region_id, synonym, external_db_id)\n";
    print SQL "\tvalues ($max_seq_region_id, '$alt_name', $ensembl_syn_ext_id);\n\n";

    print SQL "insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values ($max_seq_region_id, $toplevel, 1);\n";
    print SQL "insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values ($max_seq_region_id, $non_ref, 1);\n";

    # is this a hap, a patch_novel or a patch_fix?
    my $hap_type;
    if (exists $name_to_type{$alt_name} && defined $name_to_type{$alt_name}) {
      if ($name_to_type{$alt_name} =~ /fix/i) {
        print SQL "insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values ($max_seq_region_id, $patch_fix, NOW());\n";
        $hap_type = "'PATCH_FIX'";
      } elsif ($name_to_type{$alt_name} =~ /novel/i) {
        print SQL "insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values ($max_seq_region_id, $patch_novel, NOW());\n";
        $hap_type = "'PATCH_NOVEL'";
      } elsif ($name_to_type{$alt_name} =~ /hap/i) {
		# note there is no need to insert anything into the seq_region_attrib because the only attrib for haps has already been inserted above (non_ref)
		$hap_type = "'HAP'"; # 'HAP' type does not appear in the original patch_type files since there is no such file for haps
		                     # but fake patch_type files are created by us to be able to use this script with haps too
      } else {
        throw("Patch type ".$name_to_type{$alt_name}." for $alt_name not recognised");
      }
    } else {
      throw("No alt_name $alt_name found in patchtypes_file");
    }

    print SQL "insert into assembly_exception (assembly_exception_id, seq_region_id, seq_region_start, seq_region_end, exc_type, exc_seq_region_id, exc_seq_region_start, exc_seq_region_end, ori)\n";
    print SQL "\tvalues($max_assembly_exception_id, $max_seq_region_id, ",
     $arr[$key_to_index{'parent_start'}], ",", 
     ($new_length+$arr[$key_to_index{'parent_start'}]-1), ",",
     "$hap_type, ",
     $seq_region_id,", ",
     $arr[$key_to_index{'parent_start'}], ", ",
     $arr[$key_to_index{'parent_stop'}], ",",
     " 1);\n\n";

    #add supercontig seq_region and chrom-supercontig mapping
    $max_seq_region_id++;

    my $central_coord_sys;
    $coord_sth->execute($central_coord_system);
    $coord_sth->bind_columns(\$central_coord_sys);
    $coord_sth->fetch();

    print SQL "insert into seq_region (seq_region_id, name, coord_system_id, length)\n";
    print SQL "\tvalues ($max_seq_region_id, '$alt_name', $central_coord_sys, $new_length);\n\n";

    #add GL accession as synonym for supercontig
    print SQL "insert into seq_region_synonym (seq_region_id, synonym, external_db_id)\n";
    print SQL "\tvalues ($max_seq_region_id, '$alt_acc', $scaf_syn_ext_id);\n\n";

    $name_to_sc_seq_id{$alt_name} = $max_seq_region_id;
    $name_to_sc_seq_id{$alt_acc}  = $max_seq_region_id;

    my $asm_start = $arr[$key_to_index{'parent_start'}];
    my $asm_end = ($arr[$key_to_index{'parent_start'}]-1) + $new_length;
    my $cmp_start = $arr[$key_to_index{'alt_scaf_start'}] - $arr[$key_to_index{alt_start_tail}];
    my $cmp_end = $arr[$key_to_index{'alt_scaf_stop'}] + $arr[$key_to_index{alt_stop_tail}];
    my $strand = $seq_id_to_strand{$name_to_seq_id{$alt_name}};
    print SQL "insert into assembly (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start , cmp_end, ori) \n";
    print SQL "               values(".$name_to_seq_id{$alt_name}.", ".$name_to_sc_seq_id{$alt_name}.", ".$asm_start.", ".$asm_end.", $cmp_start, $cmp_end, $strand);\n\n";


    $max_seq_region_id++;
    $max_assembly_exception_id++;
  }
}



# Create haplotype seq_region (calculate length from mapping data)
MAP: while(<MAPPER>){
  next if /^#/; # ignore comments
  chomp;
  my ($acc, $p_start, $p_end, $part, $type, $contig, $c_start, $c_end, $strand) = split(/\t/,$_);

  # check if this is one of the new patches
  if(!$name_to_seq_id{$acc}){
    next MAP;
  }
  else{
    print "About to process an agp line for ".$acc_to_name{$acc}."\n";
  }
  if($type =~ /^[ADFGOPW]$/){ 
    if($strand eq "+"){
      $strand = 1;
    }
    elsif($strand eq "-"){
      $strand = -1;
    }
    else{
      die "Strand is neither + or - ??\n";
    }
    my $cmp_seq_id = undef;

    if(defined($name_to_seq_id{$contig})){
      $cmp_seq_id = $name_to_seq_id{$contig};
    }
    else{
      $get_seq_id_sth->execute($contig);
      $get_seq_id_sth->bind_columns(\$cmp_seq_id);
      $get_seq_id_sth->fetch;

      if(!defined($cmp_seq_id)){
        print "Could not get seq_region_id for $contig trying eutils\n";
        #die "failed to fetch $contig" unless ( system("wget -o $contig.wget.log -O $contig.fa 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$contig&rettype=fasta&retmode=text'"));
        if(system("wget -o $contig.wget.log -O $contig.fa 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$contig&rettype=fasta&retmode=text'")){
          die "failed to fetch $contig";
        }
        my $contig_seq ="";
        open(FA,"<./$contig.fa") || die "Could not open file ./$contig.fa\n";
        while (<FA>){
          if(/^>/){
#            $contig_name =~ /^>(S+)/;
#            print $contig_name."\n"
          }
          else{
            chomp;
            $contig_seq .= $_;
          }
        }
        close FA;
        system("rm ./$contig.fa");
        load_seq_region($contig, \$contig_seq);
        $cmp_seq_id = $name_to_seq_id{$contig};
      }
      else{
        print "$contig already stored :-)\n";
      }
    }

    if(!defined($cmp_seq_id)){
      die "Could not get seq id for $contig\n\n";
    }
    #add supercontig-contig mapping
    my $sc_seq_id = $name_to_sc_seq_id{$acc};
    print SQL "insert into assembly (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start , cmp_end, ori) \n";
    print SQL "               values(".$sc_seq_id.", ".$cmp_seq_id.", ".$p_start.", ".$p_end.", $c_start, $c_end, $strand);\n\n";

    #add chrom-contig mapping
    my $seq_id = $name_to_seq_id{$acc};
    my $asm_start;
    my $asm_end;
    my $cmp_start;
    my $cmp_end;

    $asm_start = $p_start + $seq_id_to_start{$seq_id} - 1;
    $cmp_start = $c_start;
    $asm_end = $seq_id_to_start{$seq_id} + $p_end - 1;
    $cmp_end = $c_end;

    #for negative chromosome strand orientations
      #alter the chromosome level starts and ends if needed

    if ($strand == 1 and $seq_id_to_strand{$seq_id} == -1){

      $asm_start = $seq_id_to_stop{$seq_id} - ($p_end -1);
      $asm_end   = $seq_id_to_stop{$seq_id} - ($p_start -1);
      $strand    = -1;
    } elsif ($strand == -1 and $seq_id_to_strand{$seq_id} == -1){

      $asm_start = $seq_id_to_stop{$seq_id} - ($p_end -1);
      $asm_end   = $seq_id_to_stop{$seq_id} - ($p_start -1);
      $strand    = 1;
    } else {
      # Positive seq_id_to_strand - DO NOTHING
    }

    print SQL "insert into assembly (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start , cmp_end, ori) \n";
    print SQL "               values(".$seq_id.", ".$cmp_seq_id.", ".$asm_start.", ".$asm_end.", $cmp_start, $cmp_end, $strand);\n\n";
  }
  else{
    print "GAP of $contig length\n";
  }
}
print SQL "\nUPDATE meta SET meta_value = '$assembly_name' WHERE meta_key = 'assembly.name';\n";
print SQL "\nUPDATE meta SET meta_value = '$assembly_acc' WHERE meta_key = 'assembly.accession';\n";
close SQL;


sub load_seq_region {
  my $name = shift;
  my $seq =  shift;

  my ($tmp_id,$tmp_len);
  $get_seq_region_sth->execute($name)
    || die "Could not execute get seq_region for $name";
  $get_seq_region_sth->bind_columns(\$tmp_id,\$tmp_len)
    || die "Could not bind get_seq_region";
  $get_seq_region_sth->fetch();
#    || die "Could not fetch seq_region id for name $name";;

  if(defined($tmp_id)){
    if($tmp_len == length($$seq)){
      print "seq region $name already found and in database\n";
      $name_to_seq_id{$name} = $tmp_id;
      return;
    }
    die "Seq region found for $name but not the same length. In datatabse length = $tmp_len and in the new fasta file ".length($$seq);
  }
  print "adding new seq region $name length of sequence = ".length($$seq)."\n";

  print SQL "insert into seq_region (seq_region_id, name, coord_system_id, length) values ($max_seq_region_id, '$name', $contig_coord_id, ".length($$seq).");\n";
  print SQL "insert into dna (seq_region_id, sequence) values ($max_seq_region_id,'".$$seq."');\n\n";

  $name_to_seq_id{$name} = $max_seq_region_id;
  $max_seq_region_id++;

  return;
}


