#!/usr/bin/env perl
#
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

use strict;
use warnings;

use Getopt::Long;
use POSIX qw(strftime);
use File::stat;
use File::Basename;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw info verbose);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

use Bio::EnsEMBL::IO::Parser::GTF;


my $host;
my $port = '3306';
my $user;
my $pass;
my $dbname;
my $dnahost;
my $dnaport = '3306';
my $dnauser;
my $dnapass;
my $dnadbname;
my $write;
my $gff_file;
my $logic_name = 'refseq_import';
my $cs = 'toplevel';
my $csv;
my $verbose = 0;

&GetOptions (
            'h|host|dbhost=s'   => \$host,
            'P|port|dbport=s'   => \$port,
            'u|user|dbuser=s'   => \$user,
            'p|pass|dbpass=s'   => \$pass,
            'd|dbname=s'        => \$dbname,
            'dnahost|ddbhost=s' => \$dnahost,
            'dnaport|ddbport=s' => \$dnaport,
            'dnauser|ddbuser=s' => \$dnauser,
            'dnapass|ddbpass=s' => \$dnapass,
            'dnadbname=s'       => \$dnadbname,
            'file=s'            => \$gff_file,
            'l|logic_name=s'    => \$logic_name,
            'write!'            => \$write,
            'verbose!'          => \$verbose,
        );

if ($verbose) {
  verbose('INFO')
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            '-host'   => $host,
            '-port'   => $port,
            '-user'   => $user,
            '-pass'   => $pass,
            '-dbname' => $dbname,
);

if ($dnahost and $dnaport and $dnauser and $dnadbname) {
  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
              '-host'   => $dnahost,
              '-port'   => $dnaport,
              '-user'   => $dnauser,
              '-pass'   => $dnapass,
              '-dbname' => $dnadbname,
  );
  $db->dnadb($dnadb);
}
else {
  warning('Check that your database has DNA or check your DNA DB connection as I could NOT create an adaptor');
}

my $sa = $db->get_SliceAdaptor;
my $ga = $db->get_GeneAdaptor;
my $aa = $db->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name($logic_name);
my $timestamp = strftime("%F %X", localtime(stat($gff_file)->mtime));
my $infile_name = basename($gff_file);
if ($analysis) {
  if ($infile_name ne $analysis->db_file) {
    warning('Your old analysis had '.$analysis->db_file." as db_file, it will be update to $infile_name\n");
    $analysis->db_file($infile_name);
  }

  warning('Old db_version: '.$analysis->db_version."\nNew db_version: $timestamp\n");
  $analysis->db_version($timestamp);
  $aa->update($analysis);
}
else {
    $analysis = Bio::EnsEMBL::Analysis->new(
        -logic_name => $logic_name,
        -db => 'RefSeq',
        -db_version => $timestamp,
        -db_file => $infile_name,
    );
}

my $codon_table = Bio::Tools::CodonTable->new;
my $gtf_parser = Bio::EnsEMBL::IO::Parser::GTF->open($gff_file, must_parse_metadata => 0);
my %sequences;
foreach my $slice (@{$sa->fetch_all('toplevel', undef, 1)}) {
  my $refseq_synonyms = $slice->get_all_synonyms();
  if ($slice->start != 1) {
    $slice = $sa->fetch_by_region($slice->coord_system->name, $slice->seq_region_name, 1, $slice->seq_region_length, 1, $slice->coord_system->version);
  }
  $sequences{$slice->seq_region_name} = $slice;
}
my %missing_sequences;
my $MT_acc;
my @par_regions;
my $par_srid;

foreach my $assemblyexception (@{$sa->db->get_AssemblyExceptionFeatureAdaptor->fetch_all}) {
  next unless ($assemblyexception->type eq 'PAR');
  push(@par_regions, [$assemblyexception->start, $assemblyexception->end]);
  $par_srid = $assemblyexception->slice->get_seq_region_id;
}

my @genes;
my %transcripts;
LINE: while ($gtf_parser->next) {
  my $seqname = $gtf_parser->get_seqname;
  next LINE if ($MT_acc && $seqname eq $MT_acc);
  next LINE if (exists $missing_sequences{$seqname});
  my $slice;
  if (exists $sequences{$seqname}) {
    $slice = $sequences{$seqname};
  }
  else {
    $slice = $sa->fetch_by_region($cs, $seqname, undef, undef, undef, $csv, undef);
# It's important to get this right: Ensembl and RefSeq have different styles of handling slices
# 1) Missing slices. Only report missing slices once. Missing slices are actually slices on
# alternative assemblies
    if ($slice) {
      $sequences{$seqname} = $slice;
    }
    else {
      warning('Slice not found '.$seqname);
      $missing_sequences{$seqname} = 1;
      next LINE;
    }
  }

  my $start = $gtf_parser->get_start;
  my $end = $gtf_parser->get_end;
# 2) Ignore annotations on pseudo-autosomal regions (currently only X/Y for human)
  if (@par_regions && $slice->get_seq_region_id() == $par_srid) {
    foreach my $aref (@par_regions) {
      if ( ($start >= $$aref[0] && $start <= $$aref[1]) || ($end >= $$aref[0] && $end <= $$aref[1]) ) {
        info( 'In PAR region, skip...');
        next LINE;
      }
    }
  }
  my $type = $gtf_parser->get_type;
  my $attributes = $gtf_parser->get_attributes;
  my $parent = $attributes->{transcript_id};
  if ($type eq 'exon') {
    my $exon = Bio::EnsEMBL::Exon->new();
    $exon->slice($slice);
    $exon->analysis($analysis);
    $exon->start($start);
    $exon->end($end);
    $exon->strand($gtf_parser->get_strand);
    $exon->stable_id($attributes->{exon_id});
    if (exists $transcripts{$parent}) {
      $transcripts{$parent}->add_Exon($exon);
    }
    else {
      throw('Missing transcript for '.$attributes->{exon_id}.' '.$start);
    }
  }
  elsif ($type eq 'CDS') {
    if (exists $transcripts{$parent}) {
      process_cds($gtf_parser, $transcripts{$parent}->{gene_object}, $transcripts{$parent}, $attributes, $start, $end);
    }
    else {
      throw('Missing transcript for '.$type.' '.$attributes->{exon_id}.' '.$start);
    }
  }
  elsif ($type eq 'start_codon' or $type eq 'stop_codon') {
    if (exists $transcripts{$parent}) {
      $transcripts{$parent}->{$type} = {
        start => $start,
        end => $end,
      };
    }
    else {
      throw('Missing transcript for '.$type.' '.$attributes->{exon_id}.' '.$start);
    }
  }
  elsif ($type eq 'transcript') {
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->slice($slice);
    $gene->analysis($analysis);
    $gene->stable_id($attributes->{gene_id});
    $gene->{__start} = $start;
    $gene->{__end} = $end;
    $gene->strand($gtf_parser->get_strand);
    $gene->biotype('non_coding');
    push(@genes, $gene);
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcripts{$attributes->{transcript_id}} = $transcript;
    $transcript->stable_id($attributes->{transcript_id});
    $transcript->slice($slice);
    $transcript->analysis($analysis);
    $transcript->biotype('non_coding');
    $transcript->{__start} = $start;
    $transcript->{__end} = $end;
    $transcript->strand($gtf_parser->get_strand);
    $transcript->{gene_object} = $gene;
    $gene->add_Transcript($transcript);
  }
  else {
    warning("$type is not supported. You need to add support for it if needed");
  }
}
$gtf_parser->close;
print "Finished processing GTF file\n";
my %stats;
GENE: foreach my $gene (@genes) {
  if ($gene->get_all_Transcripts) {
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      my $phase = 0;
      if ($transcript->translation) {
        $transcript->biotype('protein_coding');
        $gene->biotype('protein_coding');
        my $translation = $transcript->translation;
        my $genomic_start = $translation->start;
        my $genomic_end = $translation->end;
        if (exists $transcript->{start_codon}) {
          if ($transcript->strand == 1) {
            $genomic_start = $transcript->{start_codon}->{start} if ($genomic_start > $transcript->{start_codon}->{start});
          }
          else {
            $genomic_end = $transcript->{start_codon}->{end} if ($genomic_start < $transcript->{start_codon}->{end});
          }
          if ($genomic_start != $translation->start) {
            info("TRANSLATION START: $genomic_start ".$translation->start);
          }
        }
        if (exists $transcript->{stop_codon}) {
          if ($transcript->strand == 1) {
            $genomic_end = $transcript->{stop_codon}->{end} if ($genomic_end < $transcript->{stop_codon}->{end});
          }
          else {
            $genomic_start = $transcript->{stop_codon}->{start} if ($genomic_end > $transcript->{stop_codon}->{start});
          }
          if ($genomic_end != $translation->end) {
            info("TRANSLATION END: $genomic_end ".$translation->end);
          }
        }
        foreach my $exon (@{$transcript->get_all_Exons}) {
          if ($genomic_start >= $exon->seq_region_start and $genomic_start <= $exon->seq_region_end) {
            if ($exon->strand == 1) {
              $translation->start_Exon($exon);
              $translation->start($genomic_start-$exon->seq_region_start+1);
            }
            else {
              $translation->end_Exon($exon);
              $translation->end($exon->seq_region_end-$genomic_start+1);
            }
          }
          if ($genomic_end >= $exon->seq_region_start and $genomic_end <= $exon->seq_region_end) {
            if ($exon->strand == 1) {
              $translation->end_Exon($exon);
              $translation->end($genomic_end-$exon->seq_region_start+1);
            }
            else {
              $translation->start_Exon($exon);
              $translation->start($exon->seq_region_end-$genomic_end+1);
            }
          }
        }
        if (!$translation->start_Exon or !$translation->end_Exon) {
          warning('SKIPPING GENE because I could not find a start or end exon for a translation '.$gene->display_id.' '.$gene->biotype.' '.$gene->source);
          next GENE;
        }
        if (exists $transcript->translation->{phase}) {
          $phase = $transcript->translation->{phase};
        }
      }
      else {
        if ($gene->biotype ne 'protein_coding' and $gene->biotype ne $transcript->biotype) {
          info('Changing '.$transcript->stable_id.' '.$transcript->biotype.' to '.$gene->biotype.' '.$gene->stable_id);
          $transcript->biotype($gene->biotype);
        }
      }
      calculate_exon_phases($transcript, $phase);
      $stats{transcript}->{$transcript->biotype}++;
    }
    $gene->recalculate_coordinates; #Because I'm adding transcript first, I need to force recalculate
  }
  else {
    info('WARNING: Should not have happened '.$gene->stable_id);
    my $exon = Bio::EnsEMBL::Exon->new;
    $exon->start($gene->{__start});
    $exon->end($gene->{__end});
    $exon->strand($gene->strand);
    $exon->phase(-1);
    $exon->end_phase(-1);
    $exon->stable_id($gene->stable_id);
    $exon->slice($gene->slice);
    $exon->analysis($gene->analysis);
    my $transcript = Bio::EnsEMBL::Transcript->new;
    $transcript->{__start} = $gene->{__start};
    $transcript->{__end} = $gene->{__end};
    $transcript->add_Exon($exon);
    $transcript->biotype($gene->biotype);
    $transcript->stable_id($gene->stable_id);
    $transcript->description($gene->description);
    $transcript->analysis($gene->analysis);
    $stats{transcript}->{$transcript->biotype}++;
    $gene->add_Transcript($transcript);
    info('Was expecting gene '.$gene->stable_id.' '.$gene->biotype.' '.$gene->{__start}.' '.$gene->{__end}.' to have a transcript, I have added one...');
  }
# Another way to do it is to check that the canonical transcript is a NM if the gene has NM and XM
  $gene->canonical_transcript($gene->get_all_Transcripts->[0]) unless ($gene->canonical_transcript);
  $stats{gene}->{$gene->biotype}++;
  warning('Something is wrong with gene '.$gene->display_id) unless (check_gene($gene));
  print_gene($gene) if ($verbose);
  if ($write) {
    if ($gene->slice->is_toplevel) {
      $ga->store($gene);
    }
    else {
      my $toplevel_gene = $gene->transform('toplevel');
      if ($toplevel_gene) {
        $ga->store($toplevel_gene);
      }
      else {
        throw('Could not project '.$gene->display_id.' from '.$gene->slice->name.' to toplevel');
      }
    }
  }
}
foreach my $key (keys %stats) {
  foreach my $k (keys %{$stats{$key}}) {
    print $key, ' ', $k, ' ', $stats{$key}->{$k}, "\n";
  }
}

=head2 check_gene

 Arg [1]    : Bio::EnsEMBL::Gene
 Description: Check the current coordinates of the gene and transcripts against
              the values from the file. It also check for stop codon in the translation
              and for the absence of stop codon at the end of the CDS
 Returntype : Boolean, 1 if the gene is good, 0 otherwise
 Exceptions : None

=cut

sub check_gene {
  my ($gene) = @_;

  if ($gene->start == $gene->{__start} and $gene->end == $gene->{__end}) {
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      if ($transcript->start == $transcript->{__start} and $transcript->end == $transcript->{__end}) {
        if ($transcript->translation) {
          my $stop_codon = substr($transcript->translateable_seq,-3, 3);
          if ($codon_table->translate($stop_codon) ne '*') {
            info('MISSING STOP in '.$transcript->display_id." has $stop_codon instead");
          }
          if ($transcript->translation->seq =~ /\*/) {
            info('STOP in '.$transcript->translation->display_id."\n".$transcript->translation->seq);
          }
        }
      }
      else {
        info('Transcript coordinates '.$transcript->stable_id.' '.$transcript->source.' '.$transcript->biotype.' '.$transcript->start.' '.$transcript->end.' instead of '.$transcript->{__start}.' '.$transcript->{__end});
        return 0;
      }
    }
    return 1;
  }
  else {
    info('Gene coordinates '.$gene->stable_id.' '.$gene->source.' '.$gene->biotype.' '.$gene->{__start}.' '.$gene->{__end}.' instead of '.$gene->start.' '.$gene->end);
  }
  return 0;
}


=head2 process_cds

 Arg [1]    : Bio::EnsEMBL::IO::Parser::GFF3
 Arg [2]    : Bio::EnsEMBL::Gene
 Arg [3]    : Bio::EnsEMBL::Transcript
 Arg [4]    : Hash of String, the attributes for the current object
 Arg [5]    : Int start
 Arg [6]    : Int end
 Description: Process the CDS line of the file. It creates a
              Bio::EnsEMBL::Translation object and attach it
              to the transcript if it does not exists yet.
 Returntype : None
 Exceptions : None

=cut

sub process_cds {
  my ($gtf_parser, $gene, $transcript, $attributes, $start, $end) = @_;

  if ($transcript->translation) {
    $transcript->translation->start($start) if ($transcript->translation->start > $start);
    $transcript->translation->end($end) if ($transcript->translation->end < $end);
  }
  else {
    my $translation = Bio::EnsEMBL::Translation->new;
    $translation->start($start);
    $translation->end($end);
    $transcript->strand($gtf_parser->get_strand);
    $translation->stable_id($attributes->{transcript_id});
    $transcript->translation($translation);
  }
  if ($gtf_parser->get_phase) {
    if (($gtf_parser->get_strand == 1 and $transcript->translation->start >= $start)
        or ($gtf_parser->get_strand == -1 and $transcript->translation->end <= $end)) {
      if ($gtf_parser->get_phase == 1) {
        $transcript->translation->{phase} = 2;
      }
      elsif ($gtf_parser->get_phase == 2) {
        $transcript->translation->{phase} = 1;
      }
    }
  }
}


=head2 print_gene

 Arg [1]    : Bio::EnsEMBL::Gene
 Description: Print information abou the Arg[1], all the transcripts, translation and exons
 Returntype : None
 Exceptions : None

=cut

sub print_gene {
  my ($gene) = @_;

  print STDERR "GENE: ", join(' ', $gene->display_id, $gene->seq_region_name, $gene->seq_region_start, $gene->seq_region_end, $gene->strand, $gene->biotype), "\n";
  foreach my $transcript (@{$gene->get_all_Transcripts}) {
    print STDERR " TRANSCRIPT: ", join(' ', $transcript->display_id, $transcript->seq_region_name, $transcript->seq_region_start, $transcript->seq_region_end, $transcript->strand, $transcript->biotype), "\n";
    if ($transcript->translation) {
      print STDERR '   TRANSLATION: ', $transcript->translation->display_id, ' ', $transcript->translation->start, ' ', $transcript->translation->end, ' ', $transcript->translation->{phase} || 'NULL', "\n", $transcript->translation->seq, "\n";
      print STDERR "   CDNA: \n", $transcript->translateable_seq, "\n";
    }
    foreach my $exon (@{$transcript->get_all_Exons}) {
      print STDERR "  EXON: ", join(' ', $exon->display_id, $exon->seq_region_name, $exon->seq_region_start, $exon->seq_region_end, $exon->strand, $exon->phase, $exon->end_phase), "\n";
    }
  }
}
