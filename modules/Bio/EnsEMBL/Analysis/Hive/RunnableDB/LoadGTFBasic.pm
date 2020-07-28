=head1 LICENSE

 Copyright [2020] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::LoadGTFBasic

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::LoadGTFBasic->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module is for loading a simple GTF file. It works off the idea that the transcripts
should be loaded as single genes and that their cds should be calculated. It can work with
either an arrayref of gtf files (where it will load each file), or with a range within a
particular GTF file (used for parallel loading). The mode is dictated by the loading_type
param which should be either 'file' or 'range'. Defaults on file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::LoadGTFBasic;

use warnings;
use strict;
use feature 'say';
use File::Basename;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(contains_internal_stops compute_translation);

use Data::Dumper;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    loading_type => 'file',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch input expects either an arrayref of gtf files or an arrayref
              of [file_path,start_transcript_index,end_transcript_index] when doing
              range based loading
 Returntype : None
 Exceptions : Throws if it doesn't find input files

=cut

sub fetch_input {
  my ($self) = @_;

  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));
  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $target_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($target_dba,'target_db');


  my $slice_adaptor = $target_dba->get_SliceAdaptor();

  my $genome_file = $self->param_required('genome_file');
  unless(-e $genome_file) {
    $self->throw("Could not find the genome file. Path used:\n".$genome_file);
  }


  my $loading_type = $self->param_required('loading_type');
  if($loading_type eq 'range') {
    my $proto_transcripts = $self->load_range($slice_adaptor);
    unless(scalar(@$proto_transcripts)) {
      $self->throw("Loaded no transcripts. This module expects either an iid with an arrayref of gtf files or an iid with a filename and a start and end index");
    }
    $self->param('proto_transcripts',$proto_transcripts);
  } elsif($loading_type eq 'file') {
   #push(@$initial_transcripts,$self->load_full_gtf());
  } else {
    $self->throw("Unrecognised loading type selected. Refer to module for existing loading types");
  }


}


=head2 run

 Arg [1]    : None
 Description: Run will go through the gtf files, count the genes and then make batches
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $output_genes = [];
  if($self->param('loading_type') eq 'range') {
    my $proto_transcripts = $self->param('proto_transcripts');
    foreach my $proto_transcript (@$proto_transcripts) {
      push(@$output_genes,$self->build_gene($proto_transcript));
    }
  }

  unless(scalar(@$output_genes)) {
    $self->throw("No output genes created. Something went wrong");
  }

  $self->output($output_genes);
}


=head2 write_output

 Arg [1]    : None
 Description: Writes the output ids on branch 2
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_dba = $self->hrdb_get_con("target_db");
  my $gene_adaptor = $target_dba->get_GeneAdaptor();

  my $genes = $self->output();
  say "Total genes to output: ".scalar(@$genes);
  foreach my $gene (@$genes) {
    $gene_adaptor->store($gene);
  }
}


sub load_range {
  my ($self,$slice_adaptor) = @_;

  my $slice_hash = {};

  my $proto_transcripts = [];
  my $input_id = $self->param_required('gtf_array');
  unless($input_id) {
    $self->throw("No input id found via the iid parameter");
  }

  my ($file,$start_index,$end_index) = @{$input_id};
  say "Loading the following range: [".$start_index.",".$end_index."], ".$file;
  unless(-e $file) {
    $self->throw("The GTF file specified in the input id does not exist. Path specified:\n".$file);
  }

  my $source_name = basename($file);
  $source_name =~ s/\.gtf//;
  $self->param('source_name',$source_name);

  unless($start_index >= 0) {
    $self->throw("The start index is negative, this should not be. Should start at zero or higher. Start index: ".$start_index);
  }

#1       StringTie       transcript      9006    26621   1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "4.052028"; FPKM "0.425020"; TPM "0.272291";
#1       StringTie       exon    9006    10763   1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "4.108931";
#1       StringTie       exon    13491   13554   1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "4.421875";
#1       StringTie       exon    26570   26621   1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; cov "1.673077";

  my $proto_exons = [];
  my $current_index = 0;
  my $first_transcript = 1;
  open(IN,$file);
  while(<IN>) {
    my $line = $_;

    unless($current_index <= $end_index) {
      last;
    }

    if($line =~ /^#/) {
      next;
    }

    my @columns = split("\t",$line);
    my $type = $columns[2];

    if($current_index < $start_index and $type ne 'transcript') {
      next;
    }

    if($type eq 'transcript' and !$first_transcript) {
      if($current_index >= $start_index) {
        push(@$proto_transcripts,$proto_exons);
        $proto_exons = [];
      }
      $current_index++;
    } elsif($type eq 'transcript') {
      $first_transcript = 0;
    } elsif($type eq 'exon') {
      push(@$proto_exons,$line);
      my $seq_region_name = $columns[0];

      unless($slice_hash->{$seq_region_name}) {
        if($seq_region_name =~ /^[^\:]+\:[^\:]+\:([^\:]+)\:[^\:]+\:[^\:]+\:[1:\-1]$/) {
          $self->warning("Slice style name detected in first column of GTF. Switching to seq region name");
          $seq_region_name = $1;
        }
        my $slice = $slice_adaptor->fetch_by_region('toplevel',$seq_region_name);
        $slice_hash->{$seq_region_name} = $slice;
      }
    }
  } # End while
  close IN;

  if(scalar(@$proto_exons)) {
    push(@$proto_transcripts,$proto_exons);
  }

  foreach my $proto_transcript (@$proto_transcripts) {
    say "FERGAL DUMPER PT1: ".Dumper($proto_transcript);
  }

  $self->param('slice_hash',$slice_hash);
  return($proto_transcripts);
}


sub build_gene {
  my ($self,$proto_transcript) = @_;

  my $slice_hash = $self->param('slice_hash');
  my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
  my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => "gtf_import",
                                            -module     => "LoadGTFBasic");


  my $gene_id;
  my $transcript_id;
  my $strand;
  my $seq_region_name;
  my $source_name = $self->param('source_name');
  my $logic_name;
  if($source_name) {
    $logic_name = $source_name."_rnaseq_gene";
  }

  my $strand_undef = 0;
  my $exons = [];

  foreach my $proto_exon (@$proto_transcript) {
    my @columns = split("\t",$proto_exon);
    $seq_region_name = $columns[0];
    if($seq_region_name =~ /^[^\:]+\:[^\:]+\:([^\:]+)\:[^\:]+\:[^\:]+\:[1:\-1]$/) {
      $self->warning("Slice style name detected in first column of GTF. Switching to seq region name");
      $seq_region_name = $1;
    }

    my $start = $columns[3];
    my $end = $columns[4];
    $strand = $strand_conversion{$columns[6]};

    unless($strand) {
      $strand = 1;
      $strand_undef = 1;
    }

    my $phase = ($columns[7] =~ /\./) ? undef : $columns[7];
    my $attributes = $self->set_attributes($columns[8]);

    my $gene_id = $attributes->{'gene_id'};
    my $transcript_id = $attributes->{'transcript_id'};

    my $exon_number = $attributes->{'exon_number'};
    $exon_number--;

    my $exon = Bio::EnsEMBL::Exon->new(
                                      -START     => $start,
                                      -END       => $end,
                                      -STRAND    => $strand,
                                      -SLICE     => $slice_hash->{$seq_region_name},
                                      -PHASE     => -1,
                                      -END_PHASE => -1);
    $$exons[$exon_number] = $exon;
  } # End foreach my $proto_exon

  if($strand == 1) {
    $exons = [sort { $a->start <=> $b->start } @{$exons}];
  } else {
    $exons = [sort { $b->start <=> $a->start } @{$exons}];
  }

  if($logic_name) {
    $analysis->logic_name($logic_name);
  }


  my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
  $transcript->stable_id($transcript_id);
  $transcript->slice($slice_hash->{$seq_region_name});
  $transcript->source($source_name);
  $transcript->analysis($analysis);
  compute_translation($transcript);
  say "Created transcript ".$transcript_id." with ".scalar(@$exons)." exons";


  # This is a bit redundant but will do for now
  # Basically for single exon RNA-seq where the data are not stranded
  if(scalar(@$exons) == 1 and $strand_undef) {
    my $reverse_exon = Bio::EnsEMBL::Exon->new(
                                      -START     => $$exons[0]->start,
                                      -END       => $$exons[0]->end,
                                      -STRAND    => -1,
                                      -SLICE     => $slice_hash->{$seq_region_name},
                                      -PHASE     => -1,
                                      -END_PHASE => -1);

    my $reverse_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => [$reverse_exon]);
    $reverse_transcript->stable_id($transcript_id);
    $reverse_transcript->slice($slice_hash->{$seq_region_name});
    $reverse_transcript->source($source_name);
    $reverse_transcript->analysis($analysis);
    compute_translation($reverse_transcript);

    my $forward_translation = $transcript->translation;
    my $reverse_translation = $reverse_transcript->translation;
    my $ftl = 0;
    my $rtl = 0;
    if($forward_translation) {
      $ftl = length($forward_translation->seq);
    }

    if($reverse_translation) {
      $rtl = length($reverse_translation->seq);
    }

    if(length($reverse_translation) > length($forward_translation)) {
      $transcript = $reverse_transcript;
    }
  }


  my $gene = Bio::EnsEMBL::Gene->new();
  $gene->stable_id($gene_id);
  $gene->slice($slice_hash->{$seq_region_name});
  $gene->source($source_name);
  $gene->analysis($analysis);
  $gene->add_Transcript($transcript);

  return($gene);
}


sub set_attributes {
  my ($self,$attribute_string) = @_;
  my $attribute_pairs = {};

  my @attribute_array = split(";",$attribute_string);
  foreach my $attribute (@attribute_array) {
    my @pairs = split(" ",$attribute);
    if(scalar(@pairs) == 2) {
      $pairs[1] =~ s/\"//g;
      $attribute_pairs->{$pairs[0]} = $pairs[1];
    }
  }

  return($attribute_pairs);
}


1;
