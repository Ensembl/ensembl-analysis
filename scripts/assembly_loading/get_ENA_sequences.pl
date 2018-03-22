#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Find;

my $format = '.gz';
my $taxon;
my $wgs_id;
my $outdir;
my $fastadir;
my $agpdir;
my $component;
my $batch_size = 100;
my $ENA_ARCHIVE_BATCH_SIZE = 100;
my $help;
my $no_desc;
my $verbose;
my $use_ena = 0;
my $remove_non_atg = 0;

GetOptions (
            'tax=s'         => \$taxon,
            'id=s'          => \$wgs_id,
            'outdir=s'      => \$outdir,
            'fastadir=s'    => \$fastadir,
            'agp=s'         => \$agpdir,
            'component=s'   => \$component,
            'batch=s'       => \$batch_size,
            'verbose!'      => \$verbose,
            'no_desc!'      => \$no_desc,
            'ena!'          => \$use_ena,
            'n_correction!' => \$remove_non_atg,
            'help!'         => \$help,
        );
if ($help) {
    &Usage;
    exit(0);
}

warn('You should use bioperl >= 1.5.2!', "\n", '-help for more options', "\n\n");

if (!$taxon and !$wgs_id and !$outdir) {
    &Usage;
    die("You're missing some parameters!\n");
}

# It's only letters
$wgs_id =~ s/\d+//; # remove any digit 0-9
$wgs_id = uc($wgs_id);
my @wgs_ids = split(',',$wgs_id); # multiple wgs ids are allowed (ie for human: AADC,AADD,ABBA)

if (!$agpdir and !$component) {
    &Usage;
    die("You're missing some parameters!\n");
}

if ($component and !-e $component) {
    &Usage;
    die('Your component file '.$component." does not exist!\n");
}
#if ($agpdir and !-e $agpdir) {
#    &Usage;
#    die("-agp should be an existing directory/file!\n");
#}

if (!$fastadir) {
    $fastadir = $outdir;
}
my $contigfile = $fastadir.'/contig.fa';

die('Missing directory '.$outdir."\n") unless (-d $outdir);
die('Missing directory '.$fastadir."\n") unless (-d $fastadir);

my %contigs_id;
my @agpfiles;
my @clones;
print STDOUT "Getting contigs names...\n" if $verbose;
if ($component) {
    open(RF, $component) || die('Could not open '.$component."\n");
    while (<RF>) {
        my $line = $_;
        next if ($line =~ /^#/);
        my ($id, $wgs_id) = $line =~ /^\S+\s+(([A-Z]{4})[A-Z0-9\.]+)/;
        $contigs_id{$id} = $id if (grep {/$wgs_id/} @wgs_ids);
    }
    close(RF);
}
elsif ($agpdir) {
    if (-d $agpdir) {
        find(\&get_agp_files, $agpdir);
    }
    else {
        push(@agpfiles, split(/,/, $agpdir));
    }
    my %clones_id;
    foreach my $agpfile (@agpfiles) {
        print STDOUT $agpfile, "\n" if $verbose;
        open(RF, $agpfile) || die('Could not open '.$agpfile."\n");
        while(<RF>) {
            my $line = $_;
            next if ($line =~ /^#/);
            my ($component_type, $id) = $line =~ /^\S+\s+\d+\s+\d+\s+\d+\s+(\w)\s+(\S+)/;
            next if ($component_type eq 'N' or $component_type eq 'U');
            if ($component_type eq 'W') {
                $contigs_id{$id} = $id;
            }
            else {
                $clones_id{$id} = $id;
            }
        }
        close(RF);
    }
    @clones = keys %clones_id;
}

print STDOUT "Mfetching non-wgs sequences...\n" if ($verbose and scalar(@clones) > 0);

my $mfetch_cmd = 'mfetch -d embl -v fasta "';
my $clones_file = $outdir.'/'.$wgs_id.'_clones_seq.fa';
system('rm '.$clones_file) if (-e $clones_file);
my $index = 0;
my $cids = '';
foreach my $cid (@clones) {
    if ($index == $batch_size) {
        system($mfetch_cmd.$cids.'" >> '.$clones_file);
        $index = 0;
        $cids = '';
    }
    $cids .= $cid.' ';
    ++$index;
}
system($mfetch_cmd.$cids.'" >> '.$clones_file) if ($cids);


my $seqio_format = 'fasta';
my @files;
if (scalar(keys %contigs_id) < 100) {
    print STDERR "Using efetch for less than 100 ids\n";
    my $fetchbase = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=text&rettype=fasta';
    my $searchbase = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=';
    my $fasta = '';
    my $ids = '';
    foreach my $accession (keys %contigs_id) {
        my $search_cmd = 'wget -q -O - "'.$searchbase.$accession.'"';
        open(HTTP, $search_cmd.' | ') or die("Could not run $search_cmd");
        while (my $line = <HTTP>) {
            $ids .= '&id='.$1 if ($line =~ /<Id>(\d+)<\/Id>/);
        }
        close(HTTP) || die("Could not close $search_cmd");;
    }
    if ($ids) {
        my $fetch_cmd = 'wget -q -O '.$outdir.'/'.$wgs_id.'_efetched.fa "'.$fetchbase.$ids.'"';
        throw($fetch_cmd.' did not run successfully') if (system($fetch_cmd));
        push(@files, $wgs_id.'_efetched.fa');
    }
}
else {
    if ($use_ena) {
        $seqio_format = 'embl';
        print STDOUT "Getting the contigs from the ENA...\n" if $verbose;
        foreach my $a_wgs_id (@wgs_ids) {
          my $base_file = 'wgs_'.$a_wgs_id.'*_'.$taxon;
          my $file = $base_file.'.dat';
          my $base = 'wget -nv -N "ftp://ftp.ebi.ac.uk/pub/databases/embl';
          my $wget = "$base/wgs/$file$format\"";
          system("$wget -O $outdir/$file$format");
          if ($? != 0) {
              $wget = "$base/release/wgs/$file$format\"";
              system("$wget -O $outdir/$file$format");
              if ($? != 0) {
                  $wget = "$base/release/wgs/$base_file*$format\"";
                  system($wget.' -N -P '.$outdir);
                  die("wget got a problem!\n") if ($?);
                  opendir(DR, $outdir) || die('Could not open '.$outdir."\n");
                  @files = grep {/$base_file\_.*$format$/} readdir(DR);
                  closedir(DR);
              }
              else {
                  @files = ($file);
              }
          }
          else {
              @files = ($file);
          }
        }
    }
    else {
        print STDOUT "Getting the contigs from the NCBI...\n" if $verbose;
        my $base = 'wget -nv "ftp://ftp.ncbi.nlm.nih.gov/genbank/wgs';
        foreach my $a_wgs_id (@wgs_ids) {
          my $first_letter_uc = uc(substr($a_wgs_id,0,1));
          my $file = 'wgs.'.$a_wgs_id.'.*.fsa_nt'.$format;
          my $wget = "$base/$first_letter_uc/$file\" -P $outdir";
          die("wget got a problem!\n$wget\n") if (system($wget));
          opendir(DR, $outdir) || die("Could not open directory $outdir");
          push(@files, grep { /$file/ } readdir(DR));
          close(DR);
        }
    }
}
print STDOUT "Unzipping ".scalar(@files)." flatfile(s)...\n" if $verbose;
my $error = 0;
for (my $i = 0; $i < scalar(@files); $i++) {
    if ($files[$i] =~ /gz$/) {
        system('gunzip '.$outdir.'/'.$files[$i]);
        $error += $?;
        $files[$i] =~ s/$format$//;
    }
}
die("gunzip got a problem!\n") if ($error);

my $file = $wgs_id.'_contigs';
print STDOUT "Converting to fasta...\n" if $verbose;
my $fh = new Bio::SeqIO(-format => 'fasta', -file => '>'.$fastadir.'/'.$file.'.fa');
my %written_ids;
foreach my $file (@files) {
    my $eh = new Bio::SeqIO(-format => $seqio_format, -file => $outdir.'/'.$file);
    while (my $seq = $eh->next_seq) {
        my ($seq_id, $wgs_id) = $seq->id =~ /(([A-Z]{4})[0-9\.]+)/;
        if (grep {/$wgs_id/} @wgs_ids) {
          $seq->id($seq_id);
          $seq->desc(' ') if ($no_desc);
          if (exists $contigs_id{$seq->id} or !$agpdir) {
              if (!(exists($written_ids{$seq->id}))) {
                  $fh->write_seq($seq);
                  $written_ids{$seq_id} = 1;
              }
          }
          else {
              print STDERR 'REJECTED: ', $seq->id, " Parsed id: ", $seq_id,"\n";
          }
        }
    }
}

if (-e $clones_file) {
    print STDOUT "Checking the mfetched contigs number...\n" if $verbose;
    open(RF, 'grep \> '.$clones_file.' | ') || die('Could not grep '.$clones_file."\n");
    my @res = <RF>;
    close(RF);
    if (scalar(@clones) != scalar(@res)) {
        my %infile;
        my @missing_ids;
        my $missing = $outdir.'/'.$wgs_id.'_missing.lst';
        foreach my $res (@res) {
            $res =~ s/>(\S+).*/$1/s;
            print STDOUT $res, "\n" if ($verbose);
            $infile{$res} = 1;
        }
        open(WF, '>'.$missing) || die('Could not open '.$missing."\n");
        foreach my $clone (@clones) {
            if (!(exists $infile{$clone})) {
              push(@missing_ids,$clone);
              print WF $clone, "\n";
            }
        }
        close(WF);
        warning("You had ".scalar(@clones)." agp sequence entries and only ".scalar(@res)." mfetched sequences, have a look at $missing\n");

        # try to get the mfetch missing sequences from the ENA SVA (mfetch doesn't cope with old versions but ENA SVA does)
        my $ena_archive_file = $outdir.'/'.$wgs_id.'_ena_archive_seq.fa';
        my $ena_archive_file_retry = $outdir.'/'.$wgs_id.'_ena_archive_seq_retry.fa';
        get_ENA_archive_sequences($ena_archive_file,$ENA_ARCHIVE_BATCH_SIZE,@missing_ids);

        print STDOUT "Checking the ENA SVArchive fetched sequences number...\n" if $verbose;

        open(RFILE, 'grep \> '.$ena_archive_file.' | ') || die('Could not grep '.$ena_archive_file."\n");
        my @res_ena = <RFILE>;
        close(RFILE);
        if (scalar(@missing_ids) != scalar(@res_ena)) {
          my %infile_ena;
          my $missing_ena = $outdir.'/'.$wgs_id.'_missing_ena.lst';
          my @missing_ids_ena;
          foreach my $res_ena (@res_ena) {
            $res_ena =~ s/>\S+\|\S+\|(\S+).*/$1/s;
            print STDOUT $res_ena, "\n" if ($verbose);
            $infile_ena{$res_ena} = 1;
          }
          open(WFILE, '>'.$missing_ena) || die('Could not open '.$missing_ena."\n");
          foreach my $clone_ena (@missing_ids) {
            if (!(exists $infile_ena{$clone_ena})) {
              print WFILE $clone_ena, "\n";
              push(@missing_ids_ena,$clone_ena);
            }
          }
          close(WFILE);
          warning("You had ".scalar(@missing_ids)." missing sequences from mfetch and only ".scalar(@res_ena)." were fetched from ENA SVArchive, have a look at $missing_ena\n");
          
          # try to get the ENA archive missing sequences from the ENA SVA again with
          # a batch size of 1
          
          get_ENA_archive_sequences($ena_archive_file_retry,1,@missing_ids_ena);
          
          print STDOUT "Checking the ENA SVArchive retry fetched sequences number...\n" if $verbose;

          open(RRFILE, 'grep \> '.$ena_archive_file_retry.' | ') || die('Could not grep '.$ena_archive_file_retry."\n");
          my @res_ena_retry = <RRFILE>;
          close(RRFILE);
          if (scalar(@missing_ids_ena) != scalar(@res_ena_retry)) {
            my %infile_ena_retry;
            my $missing_ena_retry = $outdir.'/'.$wgs_id.'_missing_ena_retry.lst';
            foreach my $res_ena_retry (@res_ena_retry) {
              $res_ena_retry =~ s/>\S+\|\S+\|(\S+).*/$1/s;
              print STDOUT $res_ena_retry, "\n" if ($verbose);
              $infile_ena_retry{$res_ena_retry} = 1;
            }
            open(RWFILE, '>'.$missing_ena_retry) || die('Could not open '.$missing_ena_retry."\n");
            foreach my $clone_ena_retry (@missing_ids_ena) {
              if (!(exists $infile_ena_retry{$clone_ena_retry})) {
                print RWFILE $clone_ena_retry, "\n";
              }
            }
            close(RWFILE);
            throw("You had ".scalar(@missing_ids_ena)." missing sequences from the ENA SVArchive (with a batch size of $ENA_ARCHIVE_BATCH_SIZE) and only ".scalar(@res_ena_retry)." were fetched from ENA SVArchive retry with a batch size of 1, have a look at $missing_ena_retry\n");
          }
       } else {
         system('cat '.$ena_archive_file_retry.' '.$ena_archive_file.' '.$clones_file.' '.$fastadir.'/'.$file.'.fa  | sed \'/no match/d\' > '.$contigfile);
       } #end if number ENA SVArchive sequences check
    } else {
      system('cat '.$clones_file.' '.$fastadir.'/'.$file.'.fa > '.$contigfile);
    } # end if number mfetch sequences check
}
else {
    system('mv '.$fastadir.'/'.$file.'.fa '.$contigfile);
}

sub get_agp_files {
    push(@agpfiles, $File::Find::name) if ($_ =~ /\.comp\.agp$/ or $_ =~ /\.unlocalized\.scaf\.agp$/ or $_ =~ /unplaced\.scaf\.agp$/);
}
if ($remove_non_atg) {
    # I do know it's ugly but I have no other idea now...
    my $count = `grep -c \\> $contigfile`;
    throw("Could not remove the non ATGC letters from $contigfile") if (system("sed -i '/>/ ! y/BDEFHIJKLMOPQRSUVWXYZ/NNNNNNNNNNNNNNNNNNNNN/' $contigfile"));
    throw('Oooops!!! missing some stuff after removing non ATCG letters... Count shoud be '.$count) if (`grep -c \\> $contigfile` != $count);
}

sub get_ENA_archive_sequences {
  # wgets the ENA sequences corresponding to the accessions 'ids' in batches of 'batch_size'
  # from the ENA SVA and writes the results into 'output_file' (after deleting the file if existed)
  my ($output_file,$batch_size,@ids) = @_;

  print STDOUT "Fetching sequences from the ENA Sequence Version Archive...\n";

  my $ena_archive_cmd = 'wget -nv -O - "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=emblsva;format=fasta;style=raw;id=';
  system('rm '.$output_file) if (-e $output_file);

  my $index = 0;
  my $ids_str = '';
  foreach my $id (@ids) {
    if ($index == $batch_size) {
        system($ena_archive_cmd.$ids_str.'" >> '.$output_file);
        $index = 0;
        $ids_str = '';
    }
    $ids_str .= $id.' ';
    ++$index;
  }
  system($ena_archive_cmd.$ids_str.'" >> '.$output_file) if ($ids_str);
}

sub Usage {
    print <<EOF
  perl $0 -tax <taxonomy> -id <contig prefixe> -outdir <directory> -agp <directory/file> [-fastadir <directory>] [-component <file>] [-batch <integer>] [-help] [-verbose] [-n_correction]
    -tax          Choose between: mam -> mammals
                                  vrt -> vertebrates
                                  rod -> other rodents
                                  mus -> mouse
                                  hum -> human
                                  inv -> invertebrates
    -id           The prefixe for the assembly, found on the master record, 4 letters
    -outdir       The directory where the downloaded archive will be stored
    -agp          The agp file with all the contigs (chr.comp.agp/unplaced.scaffold.agp) or a directory with agp files
    -component    The component file, use all the lines so it's safer with -agp
    -batch        The batch size to use with mfetch, default 100
    -fastadir     The directory where the fasta file will be written. If not set, it will be written
                   in the outdir directory. Default file name is contig.fa
    -n_correction Change non ATGC letters to N
    -verbose
    -help

    !!!!
    !! You need to add this hacked version of bioperl to your PERL5LIB to make it work:
    !! /software/ensembl/genebuild/lib/bioperl
    !!  OR
    !! bioperl >= 1.5.2
    !!!!

EOF
;
}
