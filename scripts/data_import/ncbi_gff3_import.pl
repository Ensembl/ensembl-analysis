#!/usr/env perl
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

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use POSIX qw(strftime);
use File::stat;
use File::Basename;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
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
# Connection to the DNA DB
my $dnahost;
my $dnaport = '3306';
my $dnauser;
my $dnapass;
my $dnadbname;
# Parameters
my $write;
my $gff_file;
my $logic_name = 'refseq_import';
my $cs = 'toplevel';
my $csv;
# Xrefs processing
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
        );

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
  warning('Make sure that your database has DNA otherwise some features will not work');
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
}
else {
    $analysis = Bio::EnsEMBL::Analysis->new(
        -logic_name => $logic_name,
        -db => 'RefSeq',
        -db_version => $timestamp,
        -db_file => $infile_name,
    );
}

my $MT_acc;
my %sequences;
my %missing_sequences;
my %to_avoid;
my %xrefs;
my %objects_attributes = (
  cds_start_NF => Bio::EnsEMBL::Attribute->new(-code => 'cds_start_NF'),
  cds_end_NF => Bio::EnsEMBL::Attribute->new(-code => 'cds_end_NF'),
  initial_met => Bio::EnsEMBL::Attribute->new(-code => 'initial_met', -value => '1 1 M'),
);

my @par_regions;
my $par_srid;
foreach my $assemblyexception (@{$sa->db->get_AssemblyExceptionFeatureAdaptor->fetch_all}) {
  next unless ($assemblyexception->type eq 'PAR');
  push(@par_regions, [$assemblyexception->start, $assemblyexception->end]);
  $par_srid = $assemblyexception->slice->get_seq_region_id;
}

my %transcripts_to_fill;
my @genes;

my $gff_parser = Bio::EnsEMBL::IO::Parser::GFF3->open($gff_file);
while ($gff_parser->next) {
  my $seqname = $gff_parser->get_seqname;
  next if ($MT_acc && $seqname eq $MT_acc);
  next if (exists $missing_sequences{$seqname});
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
      next;
    }
  }

  my $start = $gff_parser->get_start;
  my $end = $gff_parser->get_end;
# 2) Ignore annotations on pseudo-autosomal regions (currently only X/Y for human)
  if (@par_regions && $slice->get_seq_region_id() == $par_srid) {
    foreach my $aref (@par_regions) {
      if ( ($start >= $$aref[0] && $start <= $$aref[1]) || ($end >= $$aref[0] && $end <= $$aref[1]) ) {
        print "In PAR region, skip...\n";
        next;
      }
    }
  }
  my $type = $gff_parser->get_type;
  my $attributes = $gff_parser->get_attributes;
  if (exists $attributes->{Parent}) {
    my $parent = $attributes->{Parent};
    if ($type eq 'exon' and !exists $to_avoid{$parent}) {
      my $exon = Bio::EnsEMBL::Exon->new();
      $exon->slice($slice);
      $exon->analysis($analysis);
      $exon->start($start);
      $exon->end($end);
      $exon->strand($gff_parser->get_strand);
      $exon->stable_id($attributes->{ID});
      for (my $index = @genes-1; $index > -1; $index--) {
        if ($genes[$index]->get_all_Transcripts) {
          foreach my $transcript (reverse @{$genes[$index]->get_all_Transcripts}) {
            if ($transcript->stable_id eq $parent) {
              $transcript->add_Exon($exon);
              $index = -1; #This is just to avoid GOTO, not sure it's better
              last;
            }
          }
        }
        elsif ($genes[$index]->stable_id eq $parent) {
          my $transcript;
          if (exists $transcripts_to_fill{$attributes->{Parent}}) {
            $transcript = $transcripts_to_fill{$attributes->{Parent}};
            delete $transcripts_to_fill{$attributes->{Parent}};
          }
          else {
            $transcript = Bio::EnsEMBL::Transcript->new();
          }
          $transcript->slice($slice);
          $transcript->analysis($analysis);
          $transcript->strand($genes[$index]->strand);
          $transcript->biotype($genes[$index]->biotype);
          $transcript->stable_id($genes[$index]->stable_id);
          $transcript->source($genes[$index]->source);
          $transcript->description($genes[$index]->description);
          $transcript->{__start} = $genes[$index]->{__start};
          $transcript->{__end} = $genes[$index]->{__end};
          foreach my $dbentry (@{$genes[$index]->get_all_DBEntries}) {
            $transcript->add_DBEntry($dbentry);
          }
          $transcript->display_xref($genes[$index]->display_xref);
          $transcript->add_Exon($exon);
          $genes[$index]->add_Transcript($transcript);
        }
      }
    }
    elsif ($type eq 'CDS') {
      for (my $index = @genes-1; $index > -1; $index--) {
        print STDERR "DEBUG ", $genes[$index]->stable_id, "\n";
        if ($genes[$index]->get_all_Transcripts) {
          foreach my $transcript (reverse @{$genes[$index]->get_all_Transcripts}) {
            if ($transcript->stable_id eq $parent) {
              print STDERR 'NORMAL ', $attributes->{ID}, "\n";
              process_cds($gff_parser, $genes[$index], $transcript, $attributes, $start, $end);
              $index = -1; #This is just to avoid GOTO, not sure it's better
              last;
            }
          }
          if ($genes[$index]->stable_id eq $parent) {
            print STDERR 'GOT gene ID', "\n";
              print STDERR $attributes->{ID}, "\n";
            if (@{$genes[$index]->get_all_Transcripts} == 1) {
              print STDERR $attributes->{ID}, "\n";
            print STDERR 'GOT 1 T ', $genes[$index]->get_all_Transcripts->[0]->stable_id, "\n";
              process_cds($gff_parser, $genes[$index], $genes[$index]->get_all_Transcripts->[0], $attributes, $start, $end);
            }
            else {
              if ($gff_parser->get_source eq 'Curated Genomic') {
                my $transcript = Bio::EnsEMBL::Transcript->new();
                $transcript->slice($slice);
                $transcript->start($start);
                $transcript->end($end);
                process_cds($gff_parser, $genes[$index], $transcript, $attributes, $start, $end);
                $transcripts_to_fill{$genes[$index]} = $transcript;
              }
            print STDERR 'GOT No T', "\n";
              warning('You should have 1 transcript in this gene '.$genes[$index]->stable_id);
            }
          }
        }
      }
    }
    elsif ($type eq 'miRNA') {
# We have no way of displaying mature miRNA, so we don't process them
      $to_avoid{$attributes->{ID}} = 1;
    }
    else {
      my $transcript;
      if (exists $transcripts_to_fill{$attributes->{Parent}}) {
        $transcript = $transcripts_to_fill{$attributes->{Parent}};
        delete $transcripts_to_fill{$attributes->{Parent}};
      }
      else {
        $transcript = Bio::EnsEMBL::Transcript->new();
      }
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
      if ($attributes->{Dbxref}) {
        add_xrefs(\%xrefs, $attributes->{Dbxref}, $transcript);
      }
      if (exists $attributes->{partial}) {
        if (exists $attributes->{start_range}) {
          $transcript->add_Attributes($objects_attributes{cds_start_NF})
        }
        if (exists $attributes->{end_range}) {
          $transcript->add_Attributes($objects_attributes{cds_end_NF})
        }
      }
      for (my $index = @genes-1; $index > -1; $index--) {
        if ($genes[$index]->stable_id eq $parent) {
          $genes[$index]->add_Transcript($transcript);
          last;
        }
      }
    }
  }
  elsif ($type eq 'region') {
    if (exists $attributes->{genome} and $attributes->{genome} =~ 'mitochondrion') {
      $MT_acc = $seqname;
    }
    if ($start != 1) {
      $slice = $slice->sub_Slice($start, $end, 1);
    }
  }
  else {
    if ($type eq 'gene') {
      push(@genes, Bio::EnsEMBL::Gene->new());
      $genes[-1]->slice($slice);
      $genes[-1]->analysis($analysis);
      $genes[-1]->stable_id($attributes->{ID});
      $genes[-1]->{__start} = $start;
      $genes[-1]->{__end} = $end;
      $genes[-1]->strand($gff_parser->get_strand);
      $attributes->{gene_biotype} =~ s/(\w_)(region|segment)/IG_$1gene/;
      $genes[-1]->biotype($attributes->{gene_biotype});
      $genes[-1]->external_name($attributes->{Name});
      $genes[-1]->description($gff_parser->decode_string($attributes->{description} || $attributes->{Name}));
      $genes[-1]->source(split(',', $gff_parser->get_source));
      if ($attributes->{Dbxref}) {
        add_xrefs(\%xrefs, $attributes->{Dbxref}, $genes[-1]);
      }
    }
  }
}
$gff_parser->close;
print "Finished processing GFF file\n";
my %stats;
foreach my $gene (@genes) {
  if ($gene->get_all_Transcripts) {
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
            warning("Cannot process this ribosomal slippage $direction$length");
          }
        }
      }
      if ($transcript->translation) {
        if ($transcript->biotype eq 'mRNA') {
          $transcript->biotype('protein_coding');
          if ($gene->biotype ne 'protein_coding') {
            warning('Gene should be "protein_coding" but is '.$gene->biotype.', updating');
            $gene->biotype('protein_coding');
          }
        }
        else {
          warning('Biotype is "'.$transcript->biotype.'" instead of "mRNA"');
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
        $phase = $transcript->translation->{phase} if (exists $transcript->translation->{phase});
      }
      else {
        if ($gene->biotype ne 'protein_coding' and $gene->biotype ne $transcript->biotype) {
          warning('Changing '.$transcript->stable_id.' '.$transcript->biotype.' to '.$gene->biotype.' '.$gene->stable_id);
          $transcript->biotype($gene->biotype);
        }
      }
      calculate_exon_phases($transcript, $phase);
      if ($transcript->{__start} != $transcript->start and $transcript->{__end} != $transcript->end) {
        warning('Was expecting transcript '.$transcript->stable_id.' '.$transcript->stable_id.' '.$transcript->biotype.' '.$transcript->{__start}.' '.$transcript->{__end}.' but have '.$transcript->start.' '.$transcript->end);
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
    warning('Was expecting gene '.$gene->stable_id.' '.$gene->biotype.' '.$gene->{__start}.' '.$gene->{__end}.' to have a transcript, I\'ve added one...')
      unless ($gene->biotype =~ /pseudo/ or $gene->source eq 'Curated Genomic');
  }
  if ($gene->{__start} != $gene->start and $gene->{__end} != $gene->end) {
    warning('Was expecting gene '.$gene->stable_id.' '.$gene->stable_id.' '.$gene->biotype.' '.$gene->{__start}.' '.$gene->{__end}.' but have '.$gene->start.' '.$gene->end)
      unless ($gene->biotype =~ /^IG_/ and exists $gene->{exception});
  }
  $stats{gene}->{$gene->biotype}++;
  print STDERR $gene->stable_id, ' ', $gene->biotype, ' ', $gene->display_id,  ' ', $gene->start, ' ', $gene->end, ' ', $gene->strand, ' ', $gene->description, "\n";
  foreach my $dbe (@{$gene->get_all_DBEntries}) {
    print STDERR '     DBE ', $dbe->primary_id, ' ', $dbe->display_id, ' ', $dbe->dbname, "\n";
  }
  foreach my $t (@{$gene->get_all_Transcripts}) {
    print STDERR '  ', $t->stable_id, ' ', $t->biotype, ' ', $t->display_id,  ' ', $t->start, ' ', $t->end, ' ', $t->strand, ' ', $t->description, "\n";
    foreach my $dbe (@{$t->get_all_DBEntries}) {
      print STDERR '     DBE ', $dbe->primary_id, ' ', $dbe->display_id, ' ', $dbe->dbname, "\n";
    }
    foreach my $e (@{$t->get_all_Exons}) {
      print STDERR '    ', $e->stable_id, ' ', $e->display_id,  ' ', $e->start, ' ', $e->end, ' ', $e->strand, ' ', $e->phase, ' ', $e->end_phase, "\n";
    }
    if ($t->translation) {
      print STDERR ' P ', $t->translation->stable_id, ' ', $t->translation->start, ' ', $t->translation->end, ' ', $t->translation->start_Exon->start, ' ', $t->translation->start_Exon->end, ' ', $t->translation->end_Exon->start, ' ', $t->translation->end_Exon->end, "\n";
      foreach my $dbe (@{$t->translation->get_all_DBEntries}) {
        print STDERR '     DBE ', $dbe->primary_id, ' ', $dbe->display_id, ' ', $dbe->dbname, "\n";
      }
      print STDERR ' P ', $t->translation->seq, "\n";
      warning("STOP") if ($t->translation->seq =~ /\*/);
    }
  }
  if ($write) {
    if ($gene->slice->asseembly_exception_type eq 'REF') {
      $ga->store($gene);
    }
    else {
      my $toplevel_gene = $gene->transform('toplevel');
      $ga->store($toplevel_gene);
    }
  }
}

# Just printing some stats about the number of genes and transcripts and their biotype
foreach my $key (keys %stats) {
  foreach my $k (keys %{$stats{$key}}) {
    print $key, ' ', $k, ' ', $stats{$key}->{$k}, "\n";
  }
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
      my ($dbname, $symbol) = split(':', $xref);
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
    if ($object->can('display_xref') and !exists $nondisplaydb{$dbname}) {
      $object->display_xref($xrefs_hash->{$xref});
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
    if ($attributes->{Dbxref}) {
      add_xrefs(\%xrefs, $attributes->{Dbxref}, $translation);
    }
    if (exists $attributes->{transl_except}) {
      if ($attributes->{transl_except} =~ /Sec/) {
        if (!exists $transcript->{__selenocysteine}) {
          foreach my $attribute (grep {$_ =~ /Sec/} split(',', $attributes->{transl_except})) {
            $attribute =~ /(\d+)\.\./;
            my $pos = $transcript->cdna2genomic($1, $1);
            my $seqedit = Bio::EnsEMBL::SeqEdit->new(
                                   -CODE    => '_selenocysteine',
                                   -NAME    => 'Selenocysteine',
                                   -DESC    => 'Selenocysteine',
                                   -START   => $pos,#$total_residues_before_stop_codon+1,
                                   -END     => $pos,#$total_residues_before_stop_codon+1,
                                   -ALT_SEQ => 'U'
                                   );
            $translation->add_Attributes($seqedit->get_Attribute);
          }
        }
      }
    }
    if (exists $attributes->{Note}) {
      if ($attributes->{Note} =~ /non-AUG\s+\(\w{3}\)/) {
        $translation->add_Attributes($objects_attributes{initial_met});
      }
      if ($attributes->{Note} =~ /([-+]\d+).*ribosomal frameshift/) {
        $transcript->{exception}->{note} = $1;
      }
    }
    $translation->stable_id($attributes->{ID});
    $transcript->translation($translation);
  }
  if (exists $attributes->{partial}) {
    if ($gff_parser->get_strand eq '+' and $transcript->translation->start >= $start
        or $gff_parser->get_strand eq '-' and $transcript->translation->end >= $start) {
      if ($gff_parser->get_phase == 1) {
        $transcript->translation->{phase} = 2;
      }
      elsif ($gff_parser->get_phase == 2) {
        $transcript->translation->{phase} = 1;
      }
    }
  }
  if (exists $attributes->{exception}) {
    if ($attributes->{exception} =~ /ribosomal slippage/) {
      push(@{$transcript->{exception}->{data}}, $start, $end);
      $transcript->{exception}->{type} = 'ribosomal slippage';
    }
    $gene->{exception} = $attributes->{exception};
  }
}
