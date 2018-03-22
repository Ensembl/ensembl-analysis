=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

Modified blast runnable for specific use with RFAMSEQ.
Use for running BLASTN of genomic vs RFAMSEQ prior to 
ncRNA analysis using Infernal.
Slice size seems best around 200k

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreatencRNADBs;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast;
use Bio::EnsEMBL::Analysis::Runnable::BlastRfam;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : fetch sequence out of database, instantiate the filter, 
  parser and finally the blast runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub run {
  my ($self) = @_;

  my $output_path = $self->param('output_path');

  unless(-e $output_path) {
    system("mkdir -p ".$output_path);
  }

  $self->download_files($output_path);
  $self->process_files($output_path);

  return 1;
}

sub write_output {
  my ($self) = @_;

  return 1;
}

sub download_files {
  my ($self,$output_path) = @_;
  # create Rfam.seed file
  my $exit;
  say "Updating RFAM descriptions file ... using Sanger FTP site (consider using EBI mirror)";
  system ("mkdir -p ".$output_path) unless -e $output_path;
  $exit =  system ("wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.cm.gz  -O ".$output_path."/Rfam.cm.gz");
  $self->throw("Error with obtaining Rfam covariance model file from ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.cm.gz\n") if $exit > 0;
  $exit =   system ("gunzip -f ".$output_path."/Rfam.cm.gz");
  $self->throw("Error decompressing Rfam.tar.gz\n") if $exit > 0;
  $exit =  system ("wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.seed.gz  -O ".$output_path."/Rfam.seed.gz");
  $self->throw("Error with obtaining Rfam.seed file from ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.seed.gz\n") if $exit > 0;
# Commented out code if for making the fasta file for Rfam 12 and up
#  $exit =   system ("wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/fasta_files/RF*.fa.gz -O ".$output_path."/Rfam.fasta.gz");
#  $self->throw("Error with obtaining Rfam fasta files from ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/fasta_files/RF*.fa.gz\n") if $exit > 0;
  $exit =   system ("wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.fasta.gz -O ".$output_path."/Rfam.fasta.gz");
  die ("Error with obtaining Rfam.fasta file from ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0//Rfam.fasta.gz\n") if $exit > 0;
  $exit =   system ("gunzip -f ".$output_path."/Rfam.seed.gz");
  $self->throw("Error decompressing Rfam.seed.gz\n") if $exit > 0;
  $exit =   system ("gunzip -f ".$output_path."/Rfam.fasta.gz");
  $self->throw("Error decompressing Rfam.fasta.gz\n") if $exit > 0;
  $exit =   system ("gunzip -f ".$output_path."/Rfam.thr.gz");
  say "done\nUpdating miRNA file...";
  $exit =   system ("wget ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz -O ".$output_path."/all_mirnas.embl.gz");
  $exit =   system ("gunzip -f ".$output_path."/all_mirnas.embl.gz");
  $self->throw("Error with obtaining miRNA.dat  file from ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat\n") if $exit > 0;
  $exit =   system ("wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz  -O ".$output_path."/all_mirnas.fa.gz");
  $exit =   system ("gunzip -f  ".$output_path."/all_mirnas.fa.gz");
  $self->throw("Error with obtaining hairpin.fa  file from ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa\n") if $exit > 0;
  say "done";

}

sub process_files {
  my ($self,$output_path) = @_;

  my %avoid = ( "RF00177" => 1, #ssu
                "RF00028" => 1, #GroupI Intron
                "RF00029" => 1, #Group II Intron
                "RF00005" => 1
              );
  my %families;

  print "Fetching all Rfam fasta sequences...";
  my $RFAM = Bio::SeqIO-> new
    (
     -file => $output_path."/Rfam.fasta",
     -format => "Fasta",
    );

  my $filtered_sequences = Bio::SeqIO-> new
    (
     -file => ">".$output_path."/filtered.fasta",
     -format => "Fasta",
    );

  $self->throw("Cannot open ".$output_path."/Rfam.fasta\n") unless $RFAM;
  $self->throw("Cannot open ".$output_path."/filtered.fasta\n") unless $filtered_sequences;
  say "Done";
  open (DESC,$output_path."/Rfam.seed") or die "Cannot open description file ".$output_path."/Rfam.seed\n";

  # determine if the ncRNA is a gene or cis acting etc
  # if it is a gene we want it, if not we dont!
  my $domain;
  while (<DESC>){
    chomp;
    $domain = $1 if ($_ =~ /^\#=GF AC   (RF.+)/);
    if ($_ =~ /^\#=GF TP   Gene;/){
      next if $_ =~ /miRNA/;
      next if $_ =~ /tRNA/;
      next if $_ =~ /antisense/;
      $families{$domain} = 1;
    }
  }
  close DESC;


  while (my $seq = $RFAM->next_seq){
    my $domain = $1 if $seq->display_id =~ /^(RF\d+);.+/;
    if ($families{$domain}){
      unless ($avoid{$domain}){
        $filtered_sequences->write_seq($seq);
      }
    }
  }

  say "Done\n Formatting blast databases....";
  system ("xdformat -n ".$output_path."/all_mirnas.fa");
  system ("xdformat -n ".$output_path."/filtered.fasta");
}

# This code is to deal with making the fasta file for Rfam 12 and up
sub process_files_new {
  my ($self,$output_path) = @_;

  my %avoid = ( "RF00177" => 1, #ssu
                "RF00028" => 1, #GroupI Intron
                "RF00029" => 1, #Group II Intron
                "RF00005" => 1
              );
  my %families;


  print "Fetching all Rfam fasta sequences...";
  my $RFAM = Bio::SeqIO-> new
    (
     -file => $output_path."/Rfam.fasta",
     -format => "Fasta",
    );

  my $filtered_sequences = Bio::SeqIO-> new
    (
     -file => ">".$output_path."/filtered.fasta",
     -format => "Fasta",
    );

  $self->throw("Cannot open ".$output_path."/Rfam.fasta\n") unless $RFAM;
  $self->throw("Cannot open ".$output_path."/filtered.fasta\n") unless $filtered_sequences;
  say "Done";
  open (DESC,$output_path."/Rfam.seed") or die "Cannot open description file ".$output_path."/Rfam.seed\n";

  # determine if the ncRNA is a gene or cis acting etc
  # if it is a gene we want it, if not we dont!
  my $domain;
  my %accessions;

  while (<DESC>){
    chomp;
    $domain = $1 if ($_ =~ /^\#=GF AC   (RF.+)/);
    if ($_ =~ /^\#=GF TP   Gene;/){
      next if $_ =~ /miRNA/;
      next if $_ =~ /tRNA/;
      next if $_ =~ /antisense/;
      $families{$domain} = 1;
    }

    if($_ =~ /^\#/) {
      next;
    }

    if($_ =~ /^([^\/]+)\//) {
      $accessions{$1} = $domain;
#      say "FM2 ".$1." ".$domain;
    }
  }
  close DESC;


  while (my $seq = $RFAM->next_seq){
    my $display_id = $seq->display_id;
    $display_id =~ /^([^\/]+)\/([0-9\-]+)/;
    my $fasta_accession = $1;
    my $range = $2;
    my $fasta_domain = $accessions{$fasta_accession};
    unless($fasta_domain) {
      next;
    }
    $seq->display_id($fasta_domain.";x;".$fasta_accession."/".$range);
    say "FM2 FASTA ACCESSION: ".$fasta_accession;
    say "FM2 FASTA RNAGE: ".$range;
    if ($families{$fasta_domain}){
        unless ($avoid{$fasta_domain}){
        $filtered_sequences->write_seq($seq);
      }
    }
  }

  say "Done\n Formatting blast databases....";
  system ("xdformat -n ".$output_path."/all_mirnas.fa");
  system ("xdformat -n ".$output_path."/filtered.fasta");

}
1;
