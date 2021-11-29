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
my $gtf_file;
my $logic_name = 'braker';
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
            'file=s'            => \$gtf_file,
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
my $timestamp = strftime("%F %X", localtime(stat($gtf_file)->mtime));
my $infile_name = basename($gtf_file);
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
	#-db => $db,
	#-db_version => $timestamp,
	#-db_file => $infile_name,
	);
}

my $codon_table = Bio::Tools::CodonTable->new;
my $gtf_parser = Bio::EnsEMBL::IO::Parser::GTF->open($gtf_file, must_parse_metadata => 0);
my %sequences;
foreach my $slice (@{$sa->fetch_all('toplevel', undef, 1)}) {
  $sequences{$slice->seq_region_name} = $slice;
}

my %genes;
my %transcripts;
my $eid = 1;
while ($gtf_parser->next) {
  my $seqname = $gtf_parser->get_seqname;
  my $slice;
  if (exists $sequences{$seqname}) {
    $slice = $sequences{$seqname};
    my $start = $gtf_parser->get_start;
    my $end = $gtf_parser->get_end;
    my $type = $gtf_parser->get_type;
    my $strand = $gtf_parser->get_strand;
    my $source = retrieve_source($gtf_parser->get_source);
    if ($type eq 'exon') {
      my $exon = Bio::EnsEMBL::Exon->new();
      $exon->slice($slice);
      $exon->analysis($analysis);
      $exon->start($start);
      $exon->end($end);
      $exon->strand($strand);
      $exon->stable_id($eid++);
      my $attributes = $gtf_parser->get_attributes;
      if (!exists $transcripts{$attributes->{transcript_id}}) {
        my $transcript = Bio::EnsEMBL::Transcript->new();
        $transcripts{$attributes->{transcript_id}} = $transcript;
        $transcript->stable_id($attributes->{transcript_id});
        $transcript->slice($slice);
        $transcript->analysis($analysis);
        $transcript->biotype('non_coding');
        $transcript->strand($strand);
        $transcript->source($source);
      }
      $transcripts{$attributes->{transcript_id}}->add_Exon($exon);
      if (!exists $genes{$attributes->{gene_id}}) {
        my $gene = Bio::EnsEMBL::Gene->new();
        $genes{$attributes->{gene_id}} = $gene;
        $gene->slice($slice);
        $gene->analysis($analysis);
        $gene->stable_id($attributes->{gene_id});
        $gene->strand($strand);
        $gene->biotype('non_coding');
        $gene->source($transcripts{$attributes->{transcript_id}}->source);
        $genes{$attributes->{gene_id}} = $gene;
      }
      if (!exists $transcripts{$attributes->{transcript_id}}->{__gene}) {
        $transcripts{$attributes->{transcript_id}}->{__gene} = $genes{$attributes->{gene_id}};
        $genes{$attributes->{gene_id}}->add_Transcript($transcripts{$attributes->{transcript_id}});
      }
    }
    elsif ($type eq 'CDS') {
      my $attributes = $gtf_parser->get_attributes;
      if (!exists $transcripts{$attributes->{transcript_id}}) {
        my $transcript = Bio::EnsEMBL::Transcript->new();
        $transcripts{$attributes->{transcript_id}} = $transcript;
        $transcript->stable_id($attributes->{transcript_id});
        $transcript->slice($slice);
        $transcript->analysis($analysis);
        $transcript->biotype('non_coding');
        $transcript->strand($strand);
        $transcript->source($source);
      }
      process_cds($transcripts{$attributes->{transcript_id}}, $start, $end, $strand, $gtf_parser->get_phase);
    }
    elsif ($type eq 'intron') {
# just to remove warning for now
    }
    elsif ($type eq 'start_codon' or $type eq 'stop_codon') {
      my $attributes = $gtf_parser->get_attributes;
      if (!exists $transcripts{$attributes->{transcript_id}}) {
        my $transcript = Bio::EnsEMBL::Transcript->new();
        $transcripts{$attributes->{transcript_id}} = $transcript;
        $transcript->stable_id($attributes->{transcript_id});
        $transcript->slice($slice);
        $transcript->analysis($analysis);
        $transcript->biotype('non_coding');
        $transcript->strand($strand);
        $transcript->source($source);
      }
      $transcripts{$attributes->{transcript_id}}->{$type} = {
        start => $start,
        end => $end,
      };
    }
    elsif ($type eq 'transcript') {
      my $transcript_id = $gtf_parser->get_raw_attributes;
      if (!exists $transcripts{$transcript_id}) {
        my $transcript = Bio::EnsEMBL::Transcript->new();
        $transcripts{$transcript_id} = $transcript;
        $transcript->stable_id($transcript_id);
        $transcript->slice($slice);
        $transcript->analysis($analysis);
        $transcript->biotype('non_coding');
        $transcript->strand($strand);
        $transcript->source($source);
      }
    }
    elsif ($type eq 'gene') {
      my $gene_id = $gtf_parser->get_raw_attributes;
      if (!exists $genes{$gene_id}) {
        my $gene = Bio::EnsEMBL::Gene->new();
        $genes{$gene_id} = $gene;
        $gene->stable_id($gene_id);
        $gene->slice($slice);
        $gene->analysis($analysis);
        $gene->biotype('non_coding');
        $gene->strand($strand);
        $gene->source($source);
      }
    }
    else {
      warning("$type is not supported. You need to add support for it if needed");
    }
  }
  else {
    warning("Slice not found $seqname");
  }
}
$gtf_parser->close;
print "Finished processing GTF file\n";
my %stats;
GENE: foreach my $gene (sort {$a->slice->seq_region_name cmp $b->slice->seq_region_name} values %genes) {
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
    # Another way to do it is to check that the canonical transcript is a NM if the gene has NM and XM
    $gene->canonical_transcript($gene->get_all_Transcripts->[0]) unless ($gene->canonical_transcript);
    $stats{gene}->{$gene->biotype}++;
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
  else {
    warning('No transcripts for '.$gene->stable_id);
  }
}
foreach my $key (keys %stats) {
  foreach my $k (keys %{$stats{$key}}) {
    print $key, ' ', $k, ' ', $stats{$key}->{$k}, "\n";
  }
}


=head2 process_cds

 Arg [1]    : Bio::EnsEMBL::Transcript
 Arg [2]    : Int start
 Arg [3]    : Int end
 Arg [4]    : Int strand
 Arg [5]    : Int frame
 Description: Process the CDS line of the file. It creates a
              Bio::EnsEMBL::Translation object
 Returntype : None
 Exceptions : None

=cut

sub process_cds {
  my ($transcript, $start, $end, $strand, $phase) = @_;

  if ($transcript->translation) {
    $transcript->translation->start($start) if ($transcript->translation->start > $start);
    $transcript->translation->end($end) if ($transcript->translation->end < $end);
  }
  else {
    my $translation = Bio::EnsEMBL::Translation->new;
    $translation->start($start);
    $translation->end($end);
    $translation->stable_id($transcript->stable_id);
    $transcript->translation($translation);
  }
  if (($strand == 1 and $transcript->translation->start >= $start)
      or ($strand == -1 and $transcript->translation->end <= $end)) {
    if ($phase == 1) {
      $transcript->translation->{phase} = 2;
    }
    elsif ($phase == 2) {
      $transcript->translation->{phase} = 1;
    }
    else {
      $transcript->translation->{phase} = 0;
    }
  }
}


=head2 print_gene

 Arg [1]    : Bio::EnsEMBL::Gene
 Description: Print information about the Arg[1], all the transcripts, translation and exons
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


=head2 retrieve_source

 Arg [1]    : String source of the data
 Description: Return the source with 'braker_' prefixed
 Returntype : String
 Exceptions : None

=cut

sub retrieve_source {
  my ($source) = @_;

  my ($data) = split('\.', $source);
  return 'braker_'.lc($data);
}
