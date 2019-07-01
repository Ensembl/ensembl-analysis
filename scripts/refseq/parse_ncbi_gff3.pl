#!/usr/bin/env perl
#
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
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

use Bio::EnsEMBL::IO::Parser::GFF3;


# Connection to the target DB
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
my %dbconverter = (
  'GeneID'       => 'EntrezGene',
  'IMGT/GENE-DB' => 'IMGT/GENE_DB',
  'MIM'          => 'MIM_GENE',
);

my %nondisplaydb = (
  'MIM'          => 1,
);

my %unknowndbs = (
  HPRD => 1,
  Ensembl => 1,
);
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
my $gff_parser = Bio::EnsEMBL::IO::Parser::GFF3->open($gff_file, must_parse_metadata => 0);
my %sequences;
my $cs_rank2 = $db->get_CoordSystemAdaptor->fetch_by_rank(2);
foreach my $slice (@{$sa->fetch_all('toplevel', undef, 1)}) {
  my $refseq_synonyms = $slice->get_all_synonyms('RefSeq_genomic');
  if (!@$refseq_synonyms) {
    if ($cs_rank2) {
      my $projection = $slice->project($cs_rank2->name, $cs_rank2->version);
      if (@$projection and @$projection == 1) {
        $slice = $projection->[0]->to_Slice;
        $refseq_synonyms = $slice->get_all_synonyms('RefSeq_genomic');
        if (!@$refseq_synonyms) {
          warning('Could not find a RefSeq synonym for '.$slice->seq_region_name.' '.$cs_rank2->name);
          next;
        }
      }
      else {
        warning('Could not project '.$slice->seq_region_name.' to '.$cs_rank2->name);
        next;
      }
    }
    else {
      warning('Could not find a RefSeq synonym for '.$slice->seq_region_name.' on '.$slice->coord_system->name);
      next;
    }
  }
  $sequences{$refseq_synonyms->[0]->name} = $slice;
}
my %missing_sequences;
my $MT_acc;
my @par_regions;
my $par_srid;
my %xrefs;
my %objects_attributes = (
  cds_start_NF => Bio::EnsEMBL::Attribute->new(-code => 'cds_start_NF'),
  cds_end_NF => Bio::EnsEMBL::Attribute->new(-code => 'cds_end_NF'),
  initial_met => Bio::EnsEMBL::Attribute->new(-code => 'initial_met', -value => '1 1 M'),
);

foreach my $assemblyexception (@{$sa->db->get_AssemblyExceptionFeatureAdaptor->fetch_all}) {
  next unless ($assemblyexception->type eq 'PAR');
  push(@par_regions, [$assemblyexception->start, $assemblyexception->end]);
  $par_srid = $assemblyexception->slice->get_seq_region_id;
}

my %genes;
my %transcripts;
my %do_not_process;
LINE: while ($gff_parser->next) {
  my $seqname = $gff_parser->get_seqname;
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

  my $start = $gff_parser->get_start;
  my $end = $gff_parser->get_end;
# 2) Ignore annotations on pseudo-autosomal regions (currently only X/Y for human)
  if (@par_regions && $slice->get_seq_region_id() == $par_srid) {
    foreach my $aref (@par_regions) {
      if ( ($start >= $$aref[0] && $start <= $$aref[1]) || ($end >= $$aref[0] && $end <= $$aref[1]) ) {
        info( 'In PAR region, skip...');
        next LINE;
      }
    }
  }
  my $type = $gff_parser->get_type;
  my $attributes = $gff_parser->get_attributes;
  if (exists $attributes->{Parent}) {
    my $parent = $attributes->{Parent};
    if ($type eq 'exon') {
      next if (exists $do_not_process{$parent});
      my $exon = Bio::EnsEMBL::Exon->new();
      $exon->slice($slice);
      $exon->analysis($analysis);
      $exon->start($start);
      $exon->end($end);
      $exon->strand($gff_parser->get_strand);
      $exon->stable_id($attributes->{ID});
      if (exists $transcripts{$parent}) {
        $transcripts{$parent}->add_Exon($exon);
      }
      elsif (exists $genes{$parent}) {
        my $transcript = Bio::EnsEMBL::Transcript->new();
        $transcripts{$parent} = $transcript;
        $transcript->analysis($analysis);
        $transcript->strand($genes{$parent}->strand);
        $transcript->biotype($genes{$parent}->biotype);
        $transcript->stable_id($genes{$parent}->stable_id);
        $transcript->source($genes{$parent}->source);
        $transcript->description($genes{$parent}->description);
        $transcript->{__start} = $genes{$parent}->{__start};
        $transcript->{__end} = $genes{$parent}->{__end};
        foreach my $dbentry (@{$genes{$parent}->get_all_DBEntries}) {
          $transcript->add_DBEntry($dbentry);
        }
        $transcript->display_xref($genes{$parent}->display_xref);
        info('You should have 1 transcript in this gene '.$genes{$parent}->stable_id.' '.$genes{$parent}->biotype);
        $transcript->add_Exon($exon);
        $genes{$parent}->add_Transcript($transcript);
      }
      else {
        throw('Missing transcript for '.$attributes->{ID}.' '.$start);
      }
    }
    elsif ($type eq 'CDS') {
      if (exists $transcripts{$parent}) {
        process_cds($gff_parser, $transcripts{$parent}->{gene_object}, $transcripts{$parent}, $attributes, $start, $end);
      }
      elsif (exists $genes{$parent}) {
        if (!exists $transcripts{$attributes->{ID}}) {
          my $transcript = Bio::EnsEMBL::Transcript->new();
          $transcripts{$attributes->{ID}} = $transcript;
          $transcript->slice($slice);
          $transcript->start($start);
          $transcript->end($end);
          $transcript->strand($gff_parser->get_strand);
          $transcript->stable_id($attributes->{ID});
          $transcript->{gene_object} = $genes{$parent};
          $genes{$parent}->add_Transcript($transcript);
          info('You should have 1 transcript in this gene '.$genes{$parent}->stable_id.' '.$genes{$parent}->biotype. ' '.$genes{$parent}->source);
          process_cds($gff_parser, $genes{$parent}, $transcript, $attributes, $start, $end);
        }
      }
      else {
        throw('Missing transcript for '.$attributes->{ID}.' '.$start);
      }
    }
    elsif ($type eq 'miRNA') {
      $do_not_process{$attributes->{ID}} = 1;
    }
    else {
      if (!exists $transcripts{$attributes->{ID}}) {
        if (exists $genes{$parent}) {
          my $transcript = Bio::EnsEMBL::Transcript->new();
          $transcripts{$attributes->{ID}} = $transcript;
          $transcript->stable_id($attributes->{ID});
          $transcript->slice($slice);
          $transcript->analysis($analysis);
          $transcript->{__start} = $start;
          $transcript->{__end} = $end;
          $transcript->strand($gff_parser->get_strand);
          $attributes->{gbkey} =~ s/(\w_)(region|segment)/IG_$1gene/;
          $transcript->biotype($attributes->{gbkey});
          $transcript->external_name($attributes->{Name} || $attributes->{gene});
          $transcript->description($gff_parser->decode_string($attributes->{product} || $attributes->{standard_name} || $attributes->{gene}));
          $transcript->source(split(',', $gff_parser->get_source));
          if ($attributes->{Dbxref}) {
            add_xrefs(\%xrefs, $attributes->{Dbxref}, $transcript);
          }
          if (exists $attributes->{partial}) {
            if (exists $attributes->{start_range}) {
              $transcript->add_Attributes($objects_attributes{($transcript->strand == 1 ? 'cds_start_NF' : 'cds_end_NF')})
            }
            if (exists $attributes->{end_range}) {
              $transcript->add_Attributes($objects_attributes{($transcript->strand == 1 ? 'cds_end_NF' : 'cds_start_NF')})
            }
          }
          $transcript->{gene_object} = $genes{$parent};
          $genes{$parent}->add_Transcript($transcript);
          if (exists $attributes->{tag} and $attributes->{tag} eq 'refseq_select') {
            $genes{$parent}->canonical_transcript($transcript);
          }
        }
        else {
          throw('Could not find a gene for '.$attributes->{ID}.' '.$slice->seq_region_name.' '.$start.' '.$end.' '.$gff_parser->get_strand);
        }
      }
    }
  }
  elsif ($type eq 'gene' or $type eq 'pseudogene') {
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->slice($slice);
    $gene->analysis($analysis);
    $gene->stable_id($attributes->{ID});
    $gene->{__start} = $start;
    $gene->{__end} = $end;
    $gene->strand($gff_parser->get_strand);
    $attributes->{gene_biotype} =~ s/(\w_)(region|segment)/IG_$1gene/;
    $gene->biotype($attributes->{gene_biotype});
    $gene->external_name($attributes->{Name});
    $gene->description($gff_parser->decode_string($attributes->{description} || $attributes->{Name}));
    $gene->source(split(',', $gff_parser->get_source));
    if ($attributes->{Dbxref}) {
      add_xrefs(\%xrefs, $attributes->{Dbxref}, $gene);
    }
    $genes{$attributes->{ID}} = $gene;
  }
  elsif ($type eq 'region') {
    if (exists $attributes->{genome} and $attributes->{genome} eq 'mitochondrion') {
      $MT_acc = $seqname;
    }
  }
}
$gff_parser->close;
print "Finished processing GFF file\n";
my %stats;
GENE: foreach my $gene (values %genes) {
  if ($gene->get_all_Transcripts) {
    if ($gene->biotype =~ /^IG_/) {
      my $transcripts = $gene->get_all_Transcripts;
      if (@$transcripts == 1) {
        $transcripts->[0]->get_all_Exons;
      }
      else {
        $gene->flush_Transcripts;
        my @cds;
        foreach my $transcript (sort {my ($a1) = $a =~ /(\d+)/;my ($b1) = $b =~ /(\d+)/; return $a1 <=> $b1} @$transcripts) {
          eval {
            $transcript->get_all_Exons;
          };
          if ($@) {
            push(@cds, $transcript);
          }
          else {
            $gene->add_Transcript($transcript);
          }
        }
        foreach my $transcript (@{$gene->get_all_Transcripts}) {
          if (!$transcript->translation) {
            my $cds_transcript = shift @cds;
            if ($cds_transcript) {
              my $cds = $cds_transcript->translation;
              if ($transcript->start > $cds->start) {
                if ($transcript->start_Exon->seq_region_end > $cds->start) {
                  $transcript->start_Exon->move($cds->start, $transcript->start_Exon->end);
                  $transcript->recalculate_coordinates;
                }
                else {
                  $cds->start($transcript->start);
                }
              }
              if ($transcript->seq_region_end < $cds->end) {
                if ($transcript->end_Exon->seq_region_end < $cds->end) {
                  $transcript->end_Exon->move($transcript->end_Exon->start, $cds->end);
                  $transcript->recalculate_coordinates;
                }
                else {
                  $cds->end($transcript->end);
                }
              }
              $transcript->translation($cds);
            }
          }
        }
      }
    }
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      my $phase = 0;
      if (exists $transcript->{exception}) {
        if ($transcript->{exception}->{type} eq 'ribosomal slippage') {
          my ($direction, $length) = $transcript->{exception}->{note} =~ /([-+])(\d+)/;
          if ($direction eq '-') {
            for (my $index = 1; $index < @{$transcript->{exception}->{data}}-2; $index += 2) {
              if ($transcript->{exception}->{data}->[$index] == $transcript->{exception}->{data}->[$index+1]-$length+1) {
                my ($elm) = $transcript->genomic2cdna($transcript->{exception}->{data}->[$index+1], $transcript->{exception}->{data}->[$index], $transcript->strand);
                my $seq = $transcript->seq;
                $transcript->add_Attributes(Bio::EnsEMBL::SeqEdit->new(
                                                 -CODE    => '_rna_edit',
                                                 -NAME    => 'rna_edit',
                                                 -DESC    => 'RNA edit',
                                                 -START   => $elm->end,
                                                 -END     => $elm->end,
                                                 -ALT_SEQ => $seq->subseq($elm->start, $elm->end).$seq->subseq($elm->end, $elm->end),
                )->get_Attribute);
                $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(
                  -code => '_rib_frameshift',
                  -value => $direction.$length,
                ));
              }
            }
          }
          elsif ($direction eq '+') {
            my $tstrand = $transcript->strand;
            if ($tstrand == 1) {
              for (my $index = 1; $index < @{$transcript->{exception}->{data}}-2; $index += 2) {
                if ($tstrand == 1 and $transcript->{exception}->{data}->[$index]+$length+1 == $transcript->{exception}->{data}->[$index+1]) {
                  my ($elm) = $transcript->genomic2cdna($transcript->{exception}->{data}->[$index]+1, $transcript->{exception}->{data}->[$index+1]-1, $tstrand);
                  $transcript->add_Attributes(Bio::EnsEMBL::SeqEdit->new(
                                                   -CODE    => '_rna_edit',
                                                   -NAME    => 'rna_edit',
                                                   -DESC    => 'RNA edit',
                                                   -START   => $elm->start,
                                                   -END     => $elm->end,
                                                   -ALT_SEQ => '',
                  )->get_Attribute);
                  $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(
                    -code => '_rib_frameshift',
                    -value => $direction.$length,
                  ));
                }
              }
            }
            elsif ($tstrand == -1) {
              for (my $index = 0; $index < @{$transcript->{exception}->{data}}-2; $index += 2) {
                if ($tstrand == -1 and $transcript->{exception}->{data}->[$index] == $transcript->{exception}->{data}->[$index+3]+$length+1) {
                  my ($elm) = $transcript->genomic2cdna($transcript->{exception}->{data}->[$index+3]+1, $transcript->{exception}->{data}->[$index]-1, $tstrand);
                  $transcript->add_Attributes(Bio::EnsEMBL::SeqEdit->new(
                                                   -CODE    => '_rna_edit',
                                                   -NAME    => 'rna_edit',
                                                   -DESC    => 'RNA edit',
                                                   -START   => $elm->start,
                                                   -END     => $elm->end,
                                                   -ALT_SEQ => '',
                  )->get_Attribute);
                  $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(
                    -code => '_rib_frameshift',
                    -value => $direction.$length,
                  ));
                }
              }
            }
          }
          else {
            info("Cannot process this ribosomal slippage $direction$length");
          }
        }
      }
      if ($transcript->translation) {
        if ($transcript->biotype eq 'mRNA') {
          $transcript->biotype('protein_coding');
          if ($gene->biotype ne 'protein_coding') {
            info('Gene should be "protein_coding" but is '.$gene->biotype.', updating');
            $gene->biotype('protein_coding');
          }
        }
        else {
          info('Biotype is "'.$transcript->biotype.'" instead of "mRNA"');
        }
        my $translation = $transcript->translation;
        my $genomic_start = $translation->start;
        my $genomic_end = $translation->end;
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
      if ($transcript->translation) {
        foreach my $attribute (@{$transcript->translation->get_all_Attributes}) {
          if ($attribute->code eq '_selenocysteine' or $attribute->code eq 'amino_acid_sub') {
            my ($attribute_start, $attribute_end, $attribute_value) = $attribute->value =~ /(\d+) (\d+)(.*)/;
            my @coords = $transcript->genomic2pep($attribute_start, $attribute_end, $transcript->strand);
            if (@coords) {
              $attribute->value($coords[0]->start.' '.$coords[0]->end.$attribute_value);
              info('Attribute on split codon '.$transcript->display_id." $attribute_start $attribute_end$attribute_value")
                if (@coords > 1);
            }
            else {
              throw('There is a problem '.$transcript->display_id." $attribute_start $attribute_end$attribute_value");
            }
          }
        }
        foreach my $attribute (@{$transcript->get_all_Attributes('_rna_edit')}) {
          if ($attribute->value =~ /(\d+) (\d+) \*$/) {
            my ($attribute_start, $attribute_end) = ($1, $2);
            my @coords = $transcript->genomic2cdna($attribute_start, $attribute_end, $transcript->strand);
            my $sub_slice = $transcript->slice->sub_Slice($attribute_start, $attribute_end, $transcript->strand);
            my $codons = lc($sub_slice->seq);
            if (length($codons) < 3 and exists $transcript->{stop}) {
              $codons .= $transcript->{stop};
              if (@{$transcript->get_all_Attributes('_transl_end')}) {
                $transcript->get_all_Attributes('_transl_end')->[0]->value($transcript->get_all_Attributes('_transl_end')->[0]+length($transcript->{stop}));
              }
              else {
                my $attribute = Bio::EnsEMBL::Attribute->new(
                  -CODE        => '_transl_end',
                  -VALUE       => length($transcript->{stop})+$transcript->cdna_coding_end,
                );
                $transcript->add_Attributes($attribute);
              }
            }
            my @stops = $codon_table->revtranslate('*');
            my $best_codon;
            if (length($codons) == 3) {
              my $best_score = 4;
              foreach my $stop_codon (@stops) {
                my $score = 3;
                for (my $i = 0; $i < 3; $i++) {
                  --$score if (substr($stop_codon, $i, 1) eq substr($codons, $i, 1));
                }
                if ($score == 1) {
                  $best_codon = $stop_codon;
                  last;
                }
                elsif ($score < $best_score) {
                  $best_score = $score;
                  $best_codon = $stop_codon;
                }
              }
            }
            elsif (length($codons) == 2) {
              foreach my $stop_codon (@stops) {
                $best_codon = $stop_codon if ($stop_codon =~ /^$codons/);
              }
            }
            elsif (length($codons) == 1) {
              $best_codon = $stops[0];
            }
            $attribute->value($coords[0]->start.' '.$coords[0]->end.' '.uc($best_codon));
            info("Changing $codons with $best_codon at ".$coords[0]->start.' '.$coords[0]->end.' '.$attribute->value);
          }
        }
      }
      $stats{transcript}->{$transcript->biotype}++;
    }
    $gene->recalculate_coordinates; #Because I'm adding transcript first, I need to force recalculate
  }
  else {
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
    $transcript->external_name($gene->external_name);
    $transcript->analysis($gene->analysis);
    $stats{transcript}->{$transcript->biotype}++;
    foreach my $dbentry (@{$gene->get_all_DBEntries}) {
      $transcript->add_DBEntry($dbentry);
    }
    $transcript->display_xref($gene->display_xref);
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

sub check_gene {
  my ($gene) = @_;

  if ($gene->start == $gene->{__start} and $gene->end == $gene->{__end}) {
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      if ($transcript->start == $transcript->{__start} and $transcript->end == $transcript->{__end}) {
        if ($transcript->translation) {
          my $stop_codon = substr($transcript->translateable_seq,-3, 3);
          if ($codon_table->translate($stop_codon) eq '*') {
            info('STOP in '.$transcript->display_id.' but cds_end_NF is present')
              if (@{$transcript->get_all_Attributes('cds_end_NF')} != 0);
          }
          else {
            info('MISSING STOP in '.$transcript->display_id." has $stop_codon instead")
              if (@{$transcript->get_all_Attributes('cds_end_NF')} == 0);
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

=head2 add_xrefs

 Arg [1]    : Hashref of Bio::EnsEMBL::DBEntry, contains all the known Xrefs of the file
 Arg [2]    : String, the Dbxref value for the current line
 Arg [3]    : Bio::EnsEMBL::Feature, the object to add the Bio::EnsEMBL::DBEntry object on
 Description: Add a Bio::EnsEMBL::DBEntry object to Arg[3] unless the database is part of
              unknown databases (%unknowndbs). It sets the Bio::EnsEMBL::DBEntry as the
              display_xref unless the database is in %nondisplaydb. It also converts the
              name of the database to en Ensembl friendly name
 Returntype : None
 Exceptions : None

=cut

sub add_xrefs {
  my ($xrefs_hash, $xref_line, $object) = @_;

  foreach my $xref (split(',', $xref_line)) {
    if (!exists $xrefs_hash->{$xref}) {
      # HGNC ids have a : in their colon, which makes it a triplet
      my ($dbname, $symbol, $third) = split(':', $xref);
      if ($third) {
        info("External database $dbname has ':' in its id: $symbol, taken from $xref")
          unless ($symbol eq 'HGNC' or $symbol eq 'MGI');
        $symbol .= ':'.$third;
      }
      next if (exists $unknowndbs{$dbname});
      $xrefs_hash->{$xref} = Bio::EnsEMBL::DBEntry->new(
                                 -primary_id => $symbol,
                                 -display_id => $symbol,
                                 -analysis => $analysis,
                                 -dbname => exists $dbconverter{$dbname} ? $dbconverter{$dbname} : $dbname,
                               );
      $xrefs_hash->{$xref}->description($object->description) if ($object->can('description'));
      $xrefs_hash->{$xref}->display_id($object->external_name) if ($object->isa('Bio::EnsEMBL::Gene'));
    }
    $object->add_DBEntry($xrefs_hash->{$xref});
    if ($object->can('display_xref')) {
      if ($xrefs_hash->{$xref}->dbname eq 'Genbank') {
        $object->display_xref($xrefs_hash->{$xref});
      }
      elsif ($object->isa('Bio::EnsEMBL::Gene') and $xrefs_hash->{$xref}->dbname eq 'EntrezGene') {
        $object->display_xref($xrefs_hash->{$xref});
        # This should not be hard coded but the API does not allow us to fetch the db display name
        $object->description($object->description.' [Source:NCBI gene;Acc:'.$xrefs_hash->{$xref}->primary_id.']');
      }
    }
  }
}


=head2 process_cds

 Arg [1]    : Bio::EnsEMBL::IO::Parser::GFF3
 Arg [2]    : Bio::EnsEMBL::Gene
 Arg [3]    : Bio::EnsEMBL::Transcript
 Arg [4]    : Hash of String, the attributes for the current object
 Arg [5]    : Int start
 Arg [6]    : Int end
 Description: Process the CDS line of the file. Add attributes to the correct
              object for selenocysteine, non Meth start and get the correct
              information for processing ribosomal slippage.
              It creates a Bio::EnsEMBL::Translation object and attach it to the
              transcript if it does not exists yet.
 Returntype : None
 Exceptions : None

=cut

sub process_cds {
  my ($gff_parser, $gene, $transcript, $attributes, $start, $end) = @_;

  if ($transcript->translation) {
    $transcript->translation->start($start) if ($transcript->translation->start > $start);
    $transcript->translation->end($end) if ($transcript->translation->end < $end);
  }
  else {
    my $translation = Bio::EnsEMBL::Translation->new;
    $translation->start($start);
    $translation->end($end);
    $transcript->strand($gff_parser->get_strand);
    if ($attributes->{Dbxref}) {
      add_xrefs(\%xrefs, $attributes->{Dbxref}, $translation);
    }
    if (exists $attributes->{transl_except}) {
      foreach my $attribute_string (split(',', $attributes->{transl_except})) {
        my ($attribute_start, $attribute_end, $type) = $attribute_string =~ /(\d+)..(\d+)\)?[^:]+:(\w+)/;
        my $attribute;
        if ($type eq 'Sec') {
          $attribute = Bio::EnsEMBL::Attribute->new(
            -CODE        => '_selenocysteine',
            -VALUE       => "$attribute_start $attribute_end U",
          );
        }
        elsif ($type eq 'Met' and
            exists $attributes->{Note} and $attributes->{Note} =~ /non-AUG\s+\(\w{3}\)/) {
          $attribute = $objects_attributes{initial_met};
        }
        elsif ($type eq 'TERM') {
          $attribute = Bio::EnsEMBL::Attribute->new(
            -CODE        => '_rna_edit',
            -VALUE       => "$attribute_start $attribute_end *",
          );
          $transcript->add_Attributes($attribute);
          next;
        }
        else {
          my $string_seq = get_one_letter_code($type);
          $attribute = Bio::EnsEMBL::Attribute->new(
            -CODE        => 'amino_acid_sub',
            -VALUE       => "$attribute_start $attribute_end ".($string_seq ? $string_seq : ''),
          );
        }
        $translation->add_Attributes($attribute);
      }
    }
    if (exists $attributes->{Note}) {
      my $note = $attributes->{Note};
      if ($note =~ /([-+]\d+).*ribosomal frameshift/) {
        $transcript->{exception}->{note} = $1;
      }
      elsif ($note =~ /stop codon completed by the addition of 3' ([ATGC]+)/) {
        $transcript->{stop} = $1;
      }
    }
    $translation->stable_id($attributes->{ID});
    $transcript->translation($translation);
  }
  if (exists $attributes->{exception}) {
    if ($attributes->{exception} =~ /ribosomal slippage/) {
      push(@{$transcript->{exception}->{data}}, $start, $end);
      $transcript->{exception}->{type} = 'ribosomal slippage';
    }
    $gene->{exception} = $attributes->{exception};
  }
  if ($gff_parser->get_phase) {
    if (($gff_parser->get_strand == 1 and $transcript->translation->start >= $start)
        or ($gff_parser->get_strand == -1 and $transcript->translation->end <= $end)) {
      if ($gff_parser->get_phase == 1) {
        $transcript->translation->{phase} = 2;
      }
      elsif ($gff_parser->get_phase == 2) {
        $transcript->translation->{phase} = 1;
      }
    }
  }
}


sub get_one_letter_code {
  my ($value) = @_;

  my @codons = $codon_table->revtranslate($value);
  if (@codons) {
    return $codon_table->translate($codons[0]);
  }
  return;
}

sub print_gene {
  my ($gene) = @_;

  print STDERR "GENE: ", join(' ', $gene->display_id, $gene->seq_region_name, $gene->start, $gene->end, $gene->strand, $gene->biotype), "\n";
  foreach my $dbe (@{$gene->get_all_DBEntries}) {
    print STDERR "  DBE: ", $dbe->dbname, ' ', $dbe->description, "\n";
  }
  foreach my $attribute (@{$gene->get_all_Attributes}) {
    print STDERR '  ATTRIBUTE: ', $attribute->code, ' ', $attribute->value, "\n";
  }
  foreach my $transcript (@{$gene->get_all_Transcripts}) {
    print STDERR " TRANSCRIPT: ", join(' ', $transcript->display_id, $transcript->seq_region_name, $transcript->start, $transcript->end, $transcript->strand, $transcript->biotype), "\n";
    foreach my $dbe (@{$transcript->get_all_DBEntries}) {
      print STDERR "   DBE: ", $dbe->dbname, ' ', $dbe->description, "\n";
    }
    foreach my $attribute (@{$transcript->get_all_Attributes}) {
      print STDERR '   ATTRIBUTE: ', $attribute->code, ' ', $attribute->value || 'NULL', "\n";
    }
    if ($transcript->translation) {
      foreach my $attribute (@{$transcript->translation->get_all_Attributes}) {
        print STDERR '    ATTRIBUTE: ', $attribute->code, ' ', $attribute->value || 'NULL', "\n";
      }
      print STDERR '   TRANSLATION: ', $transcript->translation->display_id, ' ', $transcript->translation->start, ' ', $transcript->translation->end, ' ', $transcript->translation->{phase} || 'NULL', "\n", $transcript->translation->seq, "\n";
      print STDERR "   CDNA: \n", $transcript->translateable_seq, "\n";
    }
    foreach my $exon (@{$transcript->get_all_Exons}) {
      print STDERR "  EXON: ", join(' ', $exon->display_id, $exon->seq_region_name, $exon->start, $exon->end, $exon->strand, $exon->phase, $exon->end_phase), "\n";
    }
  }
}
