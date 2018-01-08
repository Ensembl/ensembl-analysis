#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs;

use strict;
use warnings;
use feature 'say';


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
  unless($self->param('wgs_ids') && $self->param('output_path') && $self->param('contigs_source')) {
    $self->throw("Must pass in the following parameters:\n".
          "wgs_id e.g AAEX for Dog".
          "output_path e.g /path/to/work/dir\n".
          "contigs_source e.g. 'ENA' or 'NCBI'");
  }

}

sub run {
  my $self = shift;

  my $wgs_id = $self->param('wgs_id');
  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');
  my $contig_accession_path = $self->param('output_path')."/".$primary_assembly_dir_name."/AGP/contigs.txt";
  my $output_path = $self->param('output_path')."/".$primary_assembly_dir_name."/contigs";
  my $source = $self->param('contigs_source');

  $self->download_ftp_contigs($source,$wgs_id,$output_path);
  $self->unzip($output_path);
  $self->fix_contig_headers($source,$output_path);
  $self->find_missing_accessions($output_path,$contig_accession_path);

  say "Finished downloading contig files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}


sub download_ftp_contigs {
  my ($self,$source,$wgs_id,$output_path) = @_;

  say "The contigs will be downloaded from the ".$source." ftp site";

  # check if the output dir for contigs exists; otherwise, create it
  if (-e "$output_path") {
    say "Output path ".$output_path." found";
  } else {
    `mkdir -p $output_path`;
    if (-e "$output_path") {
      say "Output path $output_path not found.\n".$output_path." created successfully";
    } else {
      $self->throw("Cannot create output path for contigs ".$output_path);
    }
  }

  # It's only letters
  $wgs_id =~ s/\d+//; # remove any digit 0-9
  $wgs_id = uc($wgs_id);
  my @wgs_ids = split(',',$wgs_id); # multiple wgs ids are allowed (ie for human: AADC,AADD,ABBA)

  $source = lc($source);
  if($source eq 'ncbi') {
    my $base = 'wget -nv "ftp://ftp.ncbi.nlm.nih.gov/genbank/wgs';
    foreach my $a_wgs_id (@wgs_ids) {
      $a_wgs_id = uc($a_wgs_id);
      my $file = 'wgs.'.$a_wgs_id.'.*.fsa_nt.gz';
      my $wget = "$base/$file\" -P $output_path";
      my $return = system($wget);
      if($return) {
        $self->throw("wget failed on the following command line:\n".$wget);
      }
    }
  } #elsif($source eq 'ena') {
    #$wgs_id =~ /^(..)/;
    #my $prefix = $1;

   # my $base = 'wget -nv "ftp://ftp.ebi.ac.uk/pub/databases/embl/release/wgs/'.$prefix.'/';
   # foreach my $a_wgs_id (@wgs_ids) {
   #   my $file = 'wgs_'.$a_wgs_id.'_*.dat.gz';
   #   my $wget = "$base/$file\" -P $output_path";
   ##   my $return = system($wget);
   #   if($return) {
   #     $self->throw("wget failed on the following command line:\n".$wget);
   #   }
   # }
  #}
   else {
    $self->throw("You have specified an unknown source! Source must be NCBI or ENA! Source specified:\n".$source);
  }
}

sub unzip {
  my ($self,$output_path) = @_;
  say "Unzipping the compressed files...";
  $self->throw("gunzip operation failed. Please check your error log file.") if (system("gunzip -r $output_path") == 1);
  say "Unzipping finished!";
}

sub fix_contig_headers {
  my ($self,$source,$output_path) = @_;

  $source = lc($source);
  if($source eq 'ncbi') {
    my $contigs_unfixed = $output_path.'/contigs_unfixed_header.fa';
    my $contigs_fixed = $output_path.'/contigs.fa';

    my $cat_files = 'cat '.$output_path.'/*.fsa_nt > '.$contigs_unfixed;
    my $return = system($cat_files);
    if($return) {
      $self->throw("Problem concatenating the contig files. Commandline used:\n".$cat_files);
    }
    open(IN,$contigs_unfixed);
    open(OUT,">$contigs_fixed");
    while(<IN>) {
      my $line = $_;
      if($line =~ /^>.*gb\|([^\|]+\.\d+)\|/) {
        say OUT '>'.$1;
      } elsif($line =~ /^>/) {
        $self->throw("Found a header line that could not be parsed for the unversioned accession. Header:\n".$line);
      } else {
        print OUT $line; # print instead of say since the newline will be there already
      }
    }
    close OUT;
    close IN;

    my $contig_count1 = int(`grep -c '>' $contigs_unfixed`);
    my $contig_count2 = int(`grep -c '>' $contigs_fixed`);

    unless($contig_count1 == $contig_count2) {
      $self->throw("The contig count in contigs_unfixed_header.fa (".$contig_count1.") did not match the count in contigs.fa (".
                   $contig_count2."). They should match");
    }

  } elsif($source eq 'ena') {

  }
}

sub find_missing_accessions {
  my ($self,$output_path,$contig_accession_path) = @_;

  my $contig_file = $output_path.'/contigs.fa';
  my $contig_headers_file = $output_path.'/contig_headers.txt';

  my $cmd = "grep '^>' ".$contig_file." > ".$contig_headers_file;
  my $result = system($cmd);
  if($result) {
    $self->throw("Issue parsing headers from contigs.fa into contig_headers.txt. Commandline used:\n".$cmd);
  }

  open(IN,$contig_accession_path);
  my @all_accessions = <IN>;
  close IN;

  my $missing_accessions = [];
  foreach my $accession (@all_accessions) {
    chomp $accession;
    $cmd = "grep -c '^>".$accession."\$' ".$contig_headers_file;
    my $found_accession = int(`$cmd`);
    unless($found_accession) {
      $self->warning("No match in initial wgs download for accession: ".$accession);
      push(@{$missing_accessions},$accession);
    }
  }

  if(scalar(@{$missing_accessions}) > 1000) {
    $self->throw("Found a large amount of missing accessions (".scalar(@{$missing_accessions})."), something might be wrong");
  } elsif(scalar(@{$missing_accessions})) {
    $self->recover_missing_accessions($output_path,$missing_accessions);
  } else {
    say "All accession accounted for in the initial wgs download";
  }

}

sub recover_missing_accessions {
  my ($self,$output_path,$missing_accessions) = @_;

    say "Using efetch for missing accessions";

    my $fetchbase = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=text&rettype=fasta';
    my $searchbase = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=';
    my $fasta = '';
    my $ids = '';
    my $max_batch_size = 100;
    my $id_string_array = [];
    my $count = 0;
    foreach my $accession (@{$missing_accessions}) {
      my $search_cmd = 'wget -q -O - "'.$searchbase.$accession.'"';
      open(HTTP, $search_cmd.' | ') or $self->throw("Could not run $search_cmd");
      while (my $line = <HTTP>) {
        $ids .= '&id='.$1 if ($line =~ /<Id>(\d+)<\/Id>/);
        $count++;
        if($count==$max_batch_size) {
           push(@{$id_string_array},$ids);
          $ids = "";
          $count=0;
        }
      }
      close(HTTP) || $self->throw("Could not close $search_cmd");
    }
  if($ids) {
     push(@{$id_string_array},$ids);
  }

  unless(scalar(@{$id_string_array})) {
    $self->throw("In recovery mode but did not find any ids from efetch");
  }

  $count = 0;
  foreach my $id_string (@{$id_string_array}) {
    my $efetch_path = $output_path.'/efetched.'.$count.'.fa';
    my $fetch_cmd = 'wget -q -O '.$efetch_path.' "'.$fetchbase.$id_string.'"';
    $self->throw($fetch_cmd.' did not run successfully') if (system($fetch_cmd));

    open(IN,$efetch_path);
    open(OUT,">".$efetch_path.".fixed");
    while(<IN>) {
        my $line = $_;
        if($line =~ /^>.*gb\|([^\|]+\.\d+)\|/) {
          say OUT '>'.$1;
        } elsif($line =~ /^>/) {
          $self->throw("Found a header line that could not be parsed for the unversioned accession. Header:\n".$line);
        } else {
          print OUT $line; # print instead of say since the newline will be there already
        }
      }
    close OUT;
    close IN;

    $count++;
  }

   my $concat_efetch_path = $output_path.'/all_efetched.fa';
   my $contigs_path = $output_path.'/contigs.fa';
   my $cmd = 'cat '.$output_path.'/efetched.*.fa.fixed > '.$concat_efetch_path;
   my $return = system($cmd);
   if($return) {
     $self->throw("Problem concatenating the efetched files together. Commandline used:\n".$cmd);
   }

  $cmd = 'cp '.$contigs_path.' '.$output_path.'/contigs_wgs_only.fa';
  $return = system($cmd);
  if($return) {
    $self->throw("Problem backing up original contig file before inserting recovered sequences. Commandline used:\n".$cmd);
  }

  $cmd = 'cat '.$contigs_path.' '.$concat_efetch_path.' > '.$output_path.'/contigs_temp.fa';
  $return = system($cmd);
  if($return) {
    $self->throw("Problem concatenating the recovered accessions with the original contig file. Commandline used:\n".$cmd);
  }

  $cmd = 'mv '.$output_path.'/contigs_temp.fa '.$contigs_path;
  $return = system($cmd);
  if($return) {
    $self->throw("Problem overwriting the original contig file with the concatenated file containing recovered seqs. Commandline used:\n".$cmd);
  }

  system('rm '.$output_path.'/efetched.*.fa');
  system('rm '.$output_path.'/efetched.*.fa.fixed');

  foreach my $accession (@{$missing_accessions}) {
    $cmd = "grep -c '^>".$accession."\$' ".$contigs_path;
    my $found_accession = int(`$cmd`);
    unless($found_accession) {
      $self->throw("The following accession was missing from the contigs.fa file: ".$accession);
    }
  }

}

#sub recover_missing_accessions_new {
#  my ($self,$output_path,$missing_accessions) = @_;

#  foreach my $accession (@{$missing_accessions}) {
#      my $accession_file = $output_path.'/'.$accession.'.fa';
#      my $cmd = 'wget -q -O '.$accession_file.' "http://www.ebi.ac.uk/ena/data/view/'.$accession.'&display=fasta&download=gzip"';
#      my $return;
##      open(IN,$output_path.'/'.$accession);
#      open(OUT,">".$efetch_path.".fixed");
#      while(<IN>) {
#        my $line = $_;
#        if($line =~ /^>.*gb\|([^\|]+\.\d+)\|/) {
#          say OUT '>'.$1;
#        } elsif($line =~ /^>/) {
#          $self->throw("Found a header line that could not be parsed for the unversioned accession. Header:\n".$line);
#        } else {
#          print OUT $line; # print instead of say since the newline will be there already
#        }
#      }
#      close OUT;
#      close IN;
#      $cmd = 'cat '.$output_path.'/contigs.fa '.$output_path.'/efetched.fa > contig_temp.fa';
#      $return = system($cmd);
#      if($return) {
#        $self->throw("Problem concatenating the recovered accessions with the original contig file. Commandline used:\n".$cmd);
#      }

#      $cmd = 'mv '.$output_path.'/contig_temp.fa '.$output_path.'/contigs.fa';
#      $return = system($cmd);
#      if($return) {
#        $self->throw("Problem overwriting the original contig file with the concatenated file containing recovered seqs. Commandline used:\n".$cmd);
#      }
#    }

#}

1;
