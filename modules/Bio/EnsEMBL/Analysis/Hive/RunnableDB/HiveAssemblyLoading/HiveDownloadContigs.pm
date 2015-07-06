#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs;

use strict;
use warnings;
use feature 'say';


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
  unless($self->param('wgs_id') && $self->param('output_path') && $self->param('contigs_source')) {
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
  my $output_path = $self->param('output_path')."/".$self->param('species_name')."/".$primary_assembly_dir_name."/contigs";
  my $source = $self->param('contigs_source');

  $self->download_ftp_contigs($source,$wgs_id,$output_path);
  $self->unzip($output_path);
  $self->fix_contig_headers($source,$output_path);

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
      my $file = 'wgs.'.$a_wgs_id.'.*.fsa_nt.gz';
      my $wget = "$base/$file\" -P $output_path";
      my $return = system($wget);
      if($return) {
        $self->throw("wget failed on the following command line:\n".$wget);
      }
    }
  } elsif($source eq 'ena') {
    $wgs_id =~ /^(..)/;
    my $prefix = $1;

    my $base = 'wget -nv "ftp://ftp.ebi.ac.uk/pub/databases/embl/release/wgs/'.$prefix.'/';
    foreach my $a_wgs_id (@wgs_ids) {
      my $file = 'wgs_'.$a_wgs_id.'_*.dat.gz';
      my $wget = "$base/$file\" -P $output_path";
      my $return = system($wget);
      if($return) {
        $self->throw("wget failed on the following command line:\n".$wget);
      }
    }
  } else {
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
      if($line =~ /^>[^\|]+\|[^\|]+\|[^\|]+\|([^\|]+\.\d+)\|/) {
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

1;
