=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
use File::Basename;

use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_best_translation);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
                loading_type => 'range', # you can use range or file
                use_transcript_id => 1, # set to 1 if you want single transcript genes
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    loading_type => 'range',
    use_transcript_id => 1,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch input expects either an arrayref of gtf files or an arrayref
              of [file_path,start_transcript_index,end_transcript_index] when doing
              range based loading
 Returntype : None
 Exceptions : Throws if it doesn't find input files
              Throws if 'use_genome_flatfile' is set and no 'genome_file' exists
              Throws if there is no transcripts to work on
              Throws if 'loading_type' is set to an unsupported type

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
  }
  elsif ($loading_type eq 'file' or $loading_type eq 'sorted_file') {
    my $logic_name;
    if ($self->param_is_defined('logic_name')) {
      $logic_name = $self->param('logic_name');
    }
    else {
      $logic_name = basename($self->param_required('filename'), '.gtf');
      my $scientific_name = $target_dba->get_MetaContainerAdaptor->get_scientific_name;
      $scientific_name =~ s/ /_/g;
      $logic_name = lc($scientific_name)."_${logic_name}_rnaseq_gene";
    }
    my $slice_adaptor = $target_dba->get_SliceAdaptor;
    my %slice_cache;
    foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
      $slice_cache{$slice->name} = $slice;
      $slice_cache{$slice->seq_region_name} = $slice;
    }
    $self->param('slice_cache', \%slice_cache);
    my $analysis = $target_dba->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    if (!$analysis) {
      $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name);
    }
    $self->analysis($analysis);
  }
  else {
    $self->throw("Unrecognised loading type selected. Refer to module for existing loading types");
  }
}


=head2 run

 Arg [1]    : None
 Description: Run will go through the gtf files and create all objects to be stored in the database
 Returntype : None
 Exceptions : Throws if no object has been generated

=cut

sub run {
  my ($self) = @_;

  if($self->param('loading_type') eq 'range') {
    my $proto_transcripts = $self->param('proto_transcripts');
    foreach my $proto_transcript (@$proto_transcripts) {
      $self->output([$self->build_gene($proto_transcript)]);
    }
  }
  elsif ($self->param('loading_type') eq 'file') {
    my $use_transcript_id = $self->param('use_transcript_id');
    my $slice_cache = $self->param('slice_cache');
    my @genes;
    my $analysis = $self->analysis;
    my $current_transcript;
    my $current_strand;
    open(FH, $self->param('filename')) || $self->throw('Could not open '.$self->param('filename'));
    while (my $line = <FH>) {
      next if (index($line, '#') == 0);
      my @data = split("\t", $line);
      my $slice = $slice_cache->{$data[0]};
      $self->throw('Unknown region '.$data[0]) unless ($slice);
      my $attributes = $self->set_attributes($data[8]);

      my $gene_id = $use_transcript_id ? $attributes->{'transcript_id'} : $attributes->{'gene_id'};
      my $transcript_id = $attributes->{'transcript_id'};

      if ($data[2] eq 'exon') {
        if ($current_transcript) {
          $current_transcript->add_Exon(Bio::EnsEMBL::Exon->new(
           -START     => $data[3],
           -END       => $data[4],
           -STRAND    => $current_strand,
           -SLICE     => $slice,
           -PHASE     => -1,
           -END_PHASE => -1
          ));
          if ($current_transcript->start == $current_transcript->{_gb_start} and $current_transcript->end == $current_transcript->{_gb_end}) {
            compute_best_translation($current_transcript);
            if (@{$current_transcript->get_all_Exons} == 1 and $current_transcript->strand == 1) {
              my %tmp_hash = %{$current_transcript->get_all_Exons->[0]};
              my $tmp_exon = Bio::EnsEMBL::Exon->new_fast(\%tmp_hash);
              $tmp_exon->strand(-1);
              my $tmp_transcript = Bio::EnsEMBL::Transcript->new;
              $tmp_transcript->add_Exon($tmp_exon);
              compute_best_translation($tmp_transcript);
              my $tmp_translation = $tmp_transcript->translation;
              if ($tmp_translation) {
                my $translation = $current_transcript->translation;
                # Because we have a single exon model, no need to look for the real sequence
                if (($translation and ($tmp_translation->end-$tmp_translation->start) > ($translation->end-$translation->start)) or !$translation) {
                  $tmp_transcript->source($current_transcript->source);
                  $tmp_transcript->stable_id($current_transcript->stable_id);
                  $tmp_transcript->analysis($genes[-1]->analysis);
                  $genes[-1]->flush_Transcripts;
                  $genes[-1]->add_Transcript($tmp_transcript);
                }
              }
            }
          }
        }
        else {
          $self->throw("No current transcript, something went wrong for '$line'");
        }
      }
      elsif ($data[2] eq 'transcript') {
        $current_transcript = Bio::EnsEMBL::Transcript->new(-stable_id => $transcript_id, -source => $data[1]);
        $current_transcript->analysis($analysis);
        $current_transcript->slice($slice); # I need to add the slice here as add_Exon which will add a slice will happen later but it will be needed for recalculate_coordinates
        $current_strand = $data[6] eq '-' ? -1 : 1;
        $current_transcript->strand($current_strand); # I need to add the strand here as add_Exon which will add a strand will happen later but it will be needed for recalculate_coordinates
        $current_transcript->{_gb_start} = $data[3];
        $current_transcript->{_gb_end} = $data[4];
        my $gene = Bio::EnsEMBL::Gene->new();
        $gene->stable_id($gene_id);
        $gene->source($data[1]);
        $gene->analysis($analysis);
        $gene->add_Transcript($current_transcript);
        push(@genes, $gene);
      }
    }
    close(FH) || $self->throw('Could not close '.$self->param('filename'));
    $self->output(\@genes);
  }
  elsif ($self->param('loading_type') eq 'unsorted_file') {
    my $use_transcript_id = $self->param('use_transcript_id');
    my $slice_cache = $self->param('slice_cache');
    my %genes;
    my %transcripts;
    my $analysis = $self->analysis;
    open(FH, $self->param('filename')) || $self->throw('Could not open '.$self->param('filename'));
    while (my $line = <FH>) {
      next if (index($line, '#') == 0);
      my @data = split("\t", $line);
      my $slice = $slice_cache->{$data[0]};
      $self->throw('Unknown region '.$data[0]) unless ($slice);
      my $attributes = $self->set_attributes($data[8]);

      my $gene_id = $use_transcript_id ? $attributes->{'transcript_id'} : $attributes->{'gene_id'};
      my $transcript_id = $attributes->{'transcript_id'};

      if ($data[2] eq 'exon') {
        my $exon = Bio::EnsEMBL::Exon->new(
                                           -START     => $data[3],
                                           -END       => $data[4],
                                           -STRAND    => $data[6] eq '-' ? -1 : 1,
                                           -SLICE     => $slice,
                                           -PHASE     => -1,
                                           -END_PHASE => -1
                                          );

        if (!exists $transcripts{$transcript_id}) {
          $transcripts{$transcript_id} = Bio::EnsEMBL::Transcript->new(-stable_id => $transcript_id, -source => $data[1]);
          $transcripts{$transcript_id}->analysis($analysis);
          $transcripts{$transcript_id}->slice($slice); # I need to add the slice here as add_Exon which will add a slice will happen later but it will be needed for recalculate_coordinates
          if ($gene_id) {
            if (!exists $genes{$gene_id}) {
              $genes{$gene_id} = Bio::EnsEMBL::Gene->new();
              $genes{$gene_id}->stable_id($gene_id);
              $genes{$gene_id}->source($data[1]);
              $genes{$gene_id}->analysis($analysis);
            }
            $genes{$gene_id}->add_Transcript($transcripts{$transcript_id});
          }
        }
        $transcripts{$transcript_id}->add_Exon($exon);
      }
      elsif ($data[2] eq 'transcript') {
        if (!exists $transcripts{$transcript_id}) {
          $transcripts{$transcript_id} = Bio::EnsEMBL::Transcript->new(-stable_id => $transcript_id, -source => $data[1]);
          $transcripts{$transcript_id}->analysis($analysis);
          $transcripts{$transcript_id}->slice($slice); # I need to add the slice here as add_Exon which will add a slice will happen later but it will be needed for recalculate_coordinates
          $transcripts{$transcript_id}->strand($data[6] eq '-' ? -1 : 1); # I need to add the strand here as add_Exon which will add a strand will happen later but it will be needed for recalculate_coordinates
        }
        if (!exists $genes{$gene_id}) {
          $genes{$gene_id} = Bio::EnsEMBL::Gene->new();
          $genes{$gene_id}->stable_id($gene_id);
          $genes{$gene_id}->source($data[1]);
          $genes{$gene_id}->analysis($analysis);
        }
        $genes{$gene_id}->add_Transcript($transcripts{$transcript_id});
      }
    }
    close(FH) || $self->throw('Could not close '.$self->param('filename'));
    foreach my $gene (sort {$a->slice <=> $b->slice || $a->slice->start <=> $b->slice->start} values %genes) {
      my $transcripts = $gene->get_all_Transcripts;
      $gene->flush_Transcripts;
      foreach my $transcript (@$transcripts) {
        compute_best_translation($transcript);
        if (@{$transcript->get_all_Exons} == 1 and $transcript->strand == 1) {
          my %tmp_hash = %{$transcript->get_all_Exons->[0]};
          my $tmp_exon = Bio::EnsEMBL::Exon->new_fast(\%tmp_hash);
          $tmp_exon->strand(-1);
          my $tmp_transcript = Bio::EnsEMBL::Transcript->new;
          $tmp_transcript->add_Exon($tmp_exon);
          compute_best_translation($tmp_transcript);
          my $tmp_translation = $tmp_transcript->translation;
          if ($tmp_translation) {
            my $translation = $transcript->translation;
            # Because we have a single exon model, no need to look for the real sequence
            if (($translation and ($tmp_translation->end-$tmp_translation->start) > ($translation->end-$translation->start)) or !$translation) {
              $gene->add_Transcript($tmp_transcript);
            }
            else {
              $gene->add_Transcript($transcript);
            }
          }
          else {
            $gene->add_Transcript($transcript);
          }
        }
        else {
          $gene->add_Transcript($transcript);
        }
      }
      $self->output([$gene]);
    }
  }

  unless(scalar(@{$self->output})) {
    $self->throw("No output genes created. Something went wrong");
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Store the genes in the 'target_db' database
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_dba = $self->hrdb_get_con("target_db");
  my $gene_adaptor = $target_dba->get_GeneAdaptor();

  my $genes = $self->output();
  $self->say_with_header("Total genes to output: ".scalar(@$genes));
  foreach my $gene (@$genes) {
    $gene_adaptor->store($gene);
  }
}


=head2 load_range

 Arg [1]    : Bio::EnsEMBL::DBSQL::SliceAdaptor
 Description: Process the tanscript and exon lines to create gene models
 Returntype : Arrayref of String, the lines corresponding to exons
              It also store a hash of Bio::EnsEMBL::Slice representing the sequences
 Exceptions : Throws if 'gtf_array' is not populated
              Throws if the file provided in 'gtf_array' does not exist
              Throws if the start index provided by 'gtf_array' is negative
              Throws if it fails opening the GTF file
              Throws if it fails cloing the GTF file

=cut

sub load_range {
  my ($self,$slice_adaptor) = @_;

  my $slice_hash = {};

  my $proto_transcripts = [];
  my $input_id = $self->param_required('gtf_array');
  unless($input_id) {
    $self->throw("No input id found via the iid parameter");
  }

  my ($file,$start_index,$end_index) = @{$input_id};
  $self->say_with_header("Loading the following range: [".$start_index.",".$end_index."], ".$file);
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
  open(IN, $file) or $self->throw("Could not open $file");
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
  close(IN) or $self->throw("Could not close $file");

  if(scalar(@$proto_exons)) {
    push(@$proto_transcripts,$proto_exons);
  }

  $self->param('slice_hash',$slice_hash);
  return($proto_transcripts);
}


=head2 build_gene

 Arg [1]    : Array of String, exon GTF lines
 Description: Create the gene, transcript and translation based on the exon lines
 Returntype : Bio::EnsEMBL::Gene
 Exceptions : None

=cut

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

    my $exon = Bio::EnsEMBL::Exon->new(
                                      -START     => $start,
                                      -END       => $end,
                                      -STRAND    => $strand,
                                      -SLICE     => $slice_hash->{$seq_region_name},
                                      -PHASE     => -1,
                                      -END_PHASE => -1);
    push(@$exons, $exon);
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
  compute_best_translation($transcript);
  $self->say_with_header("Created transcript ".$transcript_id." with ".scalar(@$exons)." exons");


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
    compute_best_translation($reverse_transcript);

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


=head2 set_attributes

 Arg [1]    : String, the attribute column of a GTF file
 Description: Split the attribute column into a dictionnary of key for the attribute
              name and its value
 Returntype : Hashref, key is the atribute, value is the value of the attribute
 Exceptions : None

=cut

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
