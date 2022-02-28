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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions;
use Bio::EnsEMBL::IO::Parser::Fasta;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  unless($self->param('wgs_id') && $self->param('output_path') && $self->param('contigs_source')) {
    $self->throw("Must pass in the following parameters:\n".
          "wgs_id e.g AAEX for Dog".
          "output_path e.g /path/to/work/dir\n".
          "contigs_source e.g. 'ENA' or 'NCBI'");
  }
  return 1;
}

sub run {
  my $self = shift;

  my $wgs_id = $self->param('wgs_id');
  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');
  my $contig_accession_path = catfile($self->param('output_path'), $primary_assembly_dir_name, 'AGP', 'contigs.txt');
  my $output_path = catdir($self->param('output_path'), $primary_assembly_dir_name, 'contigs');
  my $source = $self->param('contigs_source');

  my $single_level = 0;
  unless(-e $contig_accession_path) {
    $self->warning("No contig accession file found. Assuming assembly is single level");
    $single_level = 1;
    system('mkdir -p '.$output_path);
  }

  say "Downloading contig files";
  $self->download_ftp_contigs($source,$wgs_id,$output_path);
  say "Unzipping contig files";
  $self->unzip($output_path);
  say "Fixing contig headers";
  $self->fix_contig_headers($source,$output_path);

  unless($single_level) {
    say "Checking for accessions from nuclear AGP that are missing in the contig accessions";
    $self->find_missing_accessions($output_path,$contig_accession_path);
    say "Comparing contigs to AGP file, all contigs must match nuclear AGP accession. Non-nuclear AGP accession will have contigs removed";
    $self->compare_to_agp($output_path,$contig_accession_path);
  }

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
  if (-e $output_path) {
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
  $wgs_id = uc($wgs_id);
  my @wgs_ids = split(',',$wgs_id); # multiple wgs ids are allowed (ie for human: AADC,AADD,ABBA)

  $source = lc($source);
  if($source eq 'ncbi') {
    my $base = 'rsync -av ftp.ncbi.nlm.nih.gov::genbank/wgs/' ;  # wgs.AALT.*.fsa_nt.gz .  wget -nv "ftp://ftp.ncbi.nlm.nih.gov/genbank/wgs';
    foreach my $a_wgs_id (@wgs_ids) {
      $a_wgs_id =~ s/\d+//; # remove any digit 0-9
      $a_wgs_id = uc($a_wgs_id);
      my $wgs_sub_dir = substr($a_wgs_id,0,1);
      my $file = 'wgs.'.$a_wgs_id.'.*.fsa_nt.gz';
      my $wget = "$base/$wgs_sub_dir/$file   $output_path";
      my $return = system($wget);
      if($return) {
        $self->throw("wget/rsync failed on the following command line:\n".$wget);
      }
    }
  } elsif($source eq 'ena') {
    $wgs_id =~ /^(..)/;
    my $prefix = $1;
    $prefix = lc($prefix);

    my $base = 'rsync -av "rsync://ftp.ebi.ac.uk/pub/databases/ena/wgs_fasta/'.$prefix.'/';
    foreach my $a_wgs_id (@wgs_ids) {
      my $file = $a_wgs_id.'*.fasta.gz';
      print $file ."\n";
      my $rsync = "$base/$file\" $output_path";
      my $return = system($rsync);
      if($return) {
        $self->throw("rsync failed on the following command line:\n".$rsync);
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

  my $contigs_unfixed = catfile($output_path, 'contigs_unfixed_header.fa');
  my $contigs_fixed = catfile($output_path, 'contigs.fa');

  my $cat_files;
  my $return = 1;
  if($source eq 'ncbi') {
    $cat_files = 'cat '.$output_path.'/*.fsa_nt > '.$contigs_unfixed;
    $return = system($cat_files);
  } else {
    $cat_files = 'cat '.$output_path.'/*.fasta > '.$contigs_unfixed;
    $return = system($cat_files);
  }

  if($return) {
    $self->throw("Problem concatenating the contig files. Commandline used:\n".$cat_files);
  }

  open(IN,$contigs_unfixed);
  open(OUT,">$contigs_fixed");
  while(<IN>) {
    my $line = $_;
    if($source eq 'ncbi') {
      # fix headers for ncbi
      if($line =~ /^>.*gb\|([^\|]+\.\d+)\|/) {
        say OUT '>'.$1;
      } elsif($line =~ /^>gi\|[^\|]+\|[^\|]+\|([^\|]+)\|/) {
        say OUT '>'.$1;
      } elsif($line =~ /^\>([A-Z\d\.]+\.\d+) /) {
        say OUT '>'.$1;
      } elsif($line =~ /^>/) {
        $self->throw("Found a header line that could not be parsed for the unversioned accession. Header:\n".$line);
      } else {
        print OUT $line; # print instead of say since the newline will be there already
      }
    } else {
      # fix headers for ena
      if($line =~ /^>ENA\|[^|]+\|([A-Z\d]+\.\d+) /) {
        say OUT '>'.$1;
      } elsif($line =~ /^>/) {
        $self->throw("Found a header line that could not be parsed for the unversioned accession. Header:\n".$line);
      } else {
        print OUT $line; # print instead of say since the newline will be there already
      }
    } # end else
  }
  close OUT;
  close IN;

  my $contig_count1 = int(`grep -c '>' $contigs_unfixed`);
  my $contig_count2 = int(`grep -c '>' $contigs_fixed`);

  unless($contig_count1 == $contig_count2) {
    $self->throw("The contig count in contigs_unfixed_header.fa (".$contig_count1.") did not match the count in contigs.fa (".
                 $contig_count2."). They should match");
  }

}


=head2 find_missing_accessions

  Arg [1]    : string $output_path
  Arg [2]    : string $contig_accession_path
  Example    : $self->find_missing_accessions();
  Description: Compares a list of contig accessions from the assembly AGP files to the accession the the contigs.fa file
               and determines if any contigs listed in the AGP file are missing from the contigs.fa file. If there are less missing
               than the value of $max_allowed_missing, the missing contigs are retrieved by calling recover_missing_accessions
               If more are missing than the limit, the subroutine will throw. If there are duplicate contigs in the contigs.fa file
               the subroutine will also throw. If this subroutine finishes running successfully, all contigs that are listed in the assembly
               should be present in contigs.fa. Note that at that point there could be extra contigs that are not defined in the nuclear AGP
               files (for example contigs from the mitochondrion). These are checked later by the compare_to_agp subroutine
  Returntype : none
  Exceptions : Duplicate headers in contig file
               Large number of missing contigs
  Caller     : run
  Status     : Stable

=cut
sub find_missing_accessions {
  my ($self,$output_path,$contig_accession_path) = @_;

  my $contig_file = $output_path.'/contigs.fa';
  my $max_allowed_missing = 3000;
  my $agp_accession_hash = {};
  my $fasta_header;
  my $fasta_accession_hash = {};

  # Load the agp accessions into a hash
  open(IN,$contig_accession_path);
  while(<IN>) {
    my $agp_accession = $_;
    chomp $agp_accession;
    $agp_accession_hash->{$agp_accession} = 1;
  }
  close IN;

  # Load the fasta headers into a hash
  # Note: I have tested doing the contig accession parsing with EnsEMBL::IO and it is much much slower
  #       A grep on the command line might be slightly quicker, but the code below is cleaner
  open(IN,$contig_file);
  while(<IN>) {
    my $line = $_;
    unless($line =~ /^\>(.+)\n/) {
      next;
    }

    my $fasta_accession = $1;
    unless($fasta_accession_hash->{$fasta_accession}) {
      $fasta_accession_hash->{$fasta_accession} = 1;
    } else {
      $self->throw("There appears to be a duplicate header in ".$contig_file."\nHeader: ".$fasta_header);
    }

  }
  close IN;

  # Store these to use in the reverse comparison that's done in compare_to_agp
  $self->agp_accessions($agp_accession_hash);
  $self->fasta_accessions($fasta_accession_hash);

  my $missing_accessions = [];
  foreach my $agp_accession (keys(%$agp_accession_hash)) {
    unless($fasta_accession_hash->{$agp_accession}) {
      $self->warning("No match in initial wgs download for agp accession: ".$agp_accession);
      push(@{$missing_accessions},$agp_accession);
    }
  }

  if(scalar(@{$missing_accessions}) > $max_allowed_missing) {
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
      say "Attempting to fetch missing accession: ".$accession;
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
        } elsif($line =~ /\>([^ ]+)/) {
          say OUT '>'.$1;
	} elsif($line =~ /^>(\w+)/) {
          my $tmp_1 = $1 . '.1';
          say OUT '>'.$tmp_1;
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


=head2 compare_to_agp

  Arg [1]    : string $output_path
  Arg [2]    : string $contig_accession_path
  Example    : $self->compare_to_agp($output_path,contig_accession_path);
  Description: Compares the contig accessions in the fasta files to the AGP file contig accession. This is the reverse comparison of the
               find_missing_accessions subroutine in that it looks to accessions present in the contig files that are not listed in the
               AGP files describing the nuclear assembly. Compares the unidentified contigs to the non-nuclear AGP files and removes them
               (with a warning) as we currently do not load MTs from anywhere other than RefSeq. If there are contigs that are not described
               by any AGP file then the subroutine will throw
  Returntype : none
  Exceptions : Throw on finding contigs in the fasta files that are not described in either the nuclear or non-nuclear AGP files
               Throw on failing to remove non-nuclear contigs from the contig fasta file
               Warn on finding non-nuclear contigs
  Caller     : find_missing_accessions (set), compare_to_agp (get)
  Status     : Likely to change to allow loading of the assembly mitochondrion data in future

=cut
sub compare_to_agp {
  my ($self,$output_path,$contig_accession_path) = @_;

  my $contig_headers_file = $output_path.'/contig_headers.txt';
  my $agp_accessions = $self->agp_accessions;
  my $fasta_accessions = $self->fasta_accessions;
  my $extra_accessions = {};

  foreach my $fasta_accession (keys(%$fasta_accessions)) {
    unless($agp_accessions->{$fasta_accession}) {
      $self->warning("Found a contig accession in the contig.fa file that is not in any nuclear AGP file. Accession: ".$fasta_accession);
      $extra_accessions->{$fasta_accession} = 1;
    }
  }

  my @extra_accessions = keys(%$extra_accessions);
  if(scalar(@extra_accessions)) {
    my $non_nuclear_agp_path = $output_path."/../../non_nuclear_agp/concat_non_nuclear.agp";

#    unless(-e $non_nuclear_agp_path) {
#      $self->throw("Could not find a non-nuclear AGP file, but there appears to be extra contigs present that aren't nuclear. Offending accessions:\n".@extra_accessions);
#    }

    my $non_nuclear_accessions = {};
    if(-e $non_nuclear_agp_path) {
      open(IN,$non_nuclear_agp_path);
      while(<IN>) {
        my $line = $_;
        my @cols = split(/\t/,$line);
        if($cols[5]) {
          $non_nuclear_accessions->{$cols[5]} = 1;
        }
      }
      close IN;
    }

    my $unknown_contig_count = 0;
    foreach my $accession (@extra_accessions) {
      unless($non_nuclear_accessions->{$accession}) {
        $self->warning("An extra contig is present is present that is not present in the nuclear or non-nuclear AGP files. Offending accession: ".$accession);
        $unknown_contig_count++;
      }
    }

    if($unknown_contig_count) {
      $self->warning("Found ".$unknown_contig_count." contigs that were not present in the nuclear or non-nuclear AGP files");
    }

    say "Will delete extra non-nuclear or unknown contigs from contigs.fa";

    my $cmd = 'cp '.$output_path.'/contigs.fa'.' '.$output_path.'/contigs_with_extra_accessions.fa';
    my $return = system($cmd);
    if($return) {
      $self->throw("Problem backing up original contig file before inserting recovered sequences. Commandline used:\n".$cmd);
    }

    $cmd = 'rm '.$output_path.'/contigs.fa';
    $return = system($cmd);
    if($return) {
      $self->throw("Problem deleting the contigs.fa file. Commandline used:\n".$cmd);
    }

    open(OUT,">".$output_path.'/contigs.fa');
    my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($output_path.'/contigs_with_extra_accessions.fa');
    while($parser->next()) {
      my $seq = $parser->getSequence();
      my $header = $parser->getHeader();
      unless($extra_accessions->{$header}) {
        say OUT ">".$header;
        say OUT $seq;
      }
    }
   close OUT;

#    open(OUT,">".$output_path.'/contigs.fa');
#    open(IN,$output_path.'/contigs_with_extra_accessions.fa');
#    my $header = <IN>;
#    chomp $header;
#    $header =~ s/^\>//;
#    my $seq = "";
#    while(<IN>) {
#      my $line = $_;
#      if($line =~ /^\>/) {
#        unless($extra_accessions->{$header}) {
#          say OUT ">".$header;
#          print OUT $seq;
#        }
#        chomp $line;
#        $line =~ s/^\>//;
#        $header = $line;
#        $seq = "";
#      } else {
#        $seq .= $line;
#      }
#    }

#    unless($extra_accessions->{$header}) {
#      say OUT ">".$header;
#      print OUT $seq;
#    }
#    close IN;
#    close OUT;

    $cmd = 'rm '.$output_path.'/contigs_with_extra_accessions.fa';
    $return = system($cmd);
    if($return) {
      $self->throw("Problem deleting the contigs_with_extra_accessions.fa file. Commandline used:\n".$cmd);
    }
  }
}


=head2 agp_accessions

  Arg [1]    : hashref $val
  Example    : $self->agp_accessions($agp_accesion_hash);
  Description: Getter/setter for storing the contig accessions listed in the AGP files
  Returntype : none
  Exceptions : none
  Caller     : find_missing_accessions (set), compare_to_agp (get)
  Status     : Stable

=cut
sub agp_accessions {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_agp_accessions',$val);
  }

  return($self->param('_agp_accessions'));
}


=head2 fasta_accessions

  Arg [1]    : hashref $val
  Example    : $self->fasta_accessions($fasta_accesion_hash);
  Description: Getter/setter for storing the contig accessions listed in the fasta files
  Returntype : none
  Exceptions : none
  Caller     : find_missing_accessions (set), compare_to_agp (get)
  Status     : Stable

=cut
sub fasta_accessions {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_fasta_accessions',$val);
  }

  return($self->param('_fasta_accessions'));
}

1;
