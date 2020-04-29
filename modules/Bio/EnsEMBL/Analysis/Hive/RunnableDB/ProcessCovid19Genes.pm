=head1 LICENSE

 Copyright [2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::WriteCDSCoordsToFile

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::WriteCDSCoordsToFile->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module dumps cds coord data from Ensembl dbs to file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessCovid19Genes;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name is_canonical_splice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene clone_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript exon_overlap features_overlap overlap_length set_stop_codon);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _branch_to_flow_to => 1,
    use_generic_output_type => 0,
    generic_output_type => 'cds',
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for cds coord dumping
 Returntype : None

=cut

sub fetch_input {
  my ($self) = @_;

  say "Fetching input";
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

  my $gene_adaptor = $target_dba->get_GeneAdaptor();
  my $genes = $gene_adaptor->fetch_all();
  $self->param('input_genes',$genes);
  $self->param('gene_adaptor',$gene_adaptor);
  $self->hrdb_set_con($target_dba,'target_db');
}


sub run {
  my ($self) = @_;
  my $input_genes = $self->param('input_genes');
  my $genes = [];
  foreach my $input_gene(@$input_genes) {
    push(@$genes,clone_Gene($input_gene));
  }

  my $output_genes = [];
  foreach my $gene (@$genes) {
    say "Processing gene with id: ".$gene->dbID();
    my $updated_gene = $self->process_gene($gene);
    if($updated_gene) {
      push(@$output_genes,$updated_gene);
    } else {
      push(@$output_genes,$gene);
    }
  }

  $output_genes = $self->finalise_geneset($output_genes);

  foreach my $output_gene (@$output_genes) {
    my $transcripts = $output_gene->get_all_Transcripts();
    $output_gene->flush_Transcripts();
    foreach my $transcript (@$transcripts) {
      $transcript = set_stop_codon($transcript);
      $output_gene->add_Transcript($transcript);
    }
  }


  $self->output($output_genes);
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $gene_adaptor = $self->param('gene_adaptor');
  my $input_genes = $self->param('input_genes');
  my $output_genes = $self->output();

  foreach my $output_gene (@$output_genes) {
    empty_Gene($output_gene);
    $gene_adaptor->store($output_gene);
  }

  foreach my $input_gene (@$input_genes) {
    $gene_adaptor->remove($input_gene);
  }
}


sub finalise_geneset {
  my ($self,$genes) = @_;

  # This is going to be highly specific to the current annotation for now
  $genes = [sort { $a->start <=> $b->start } @{$genes}];

  my $gene1 = shift(@$genes);
  my $gene2 = shift(@$genes);
  my $transcript1 = ${$gene1->get_all_Transcripts}[0];
  my $transcript2 = ${$gene2->get_all_Transcripts}[0];

  unless(features_overlap($gene1,$gene2)) {
    $self->throw("Expected first two transcripts to overlap");
  }

  my $new_gene = Bio::EnsEMBL::Gene->new(-slice    => $gene1->slice,
                                         -analysis => $gene1->analysis,
                                         -biotype  => $gene1->biotype);

  $new_gene->add_Transcript($transcript1);
  $new_gene->add_Transcript($transcript2);
  push(@$genes,$new_gene);
  return($genes);
}


sub process_gene {
  my ($self,$gene) = @_;

  my $updated_gene = 0;
  my $updated_transcripts = [];
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    my $updated_transcript = 0;
    my $exons = $transcript->get_all_Exons();
    unless(scalar(@$exons) > 1) {
      push(@$updated_transcripts,$transcript);
      next;
    }

    say "Transcript with id ".$transcript->dbID." has ".scalar(@$exons)." exons. Will process further";
    my $transcript_attribs = [];
    my $updated_exons = [];
    my $current_seq_offset = 0;
    for(my $i=0; $i<scalar(@$exons)-1; $i++) {
      my $current_exon = $$exons[$i];
      my $next_exon = $$exons[$i+1];
      my $intron_length = ($next_exon->seq_region_start - 1) - ($current_exon->seq_region_end + 1) + 1;
      say "Current exon end: ".$current_exon->seq_region_end;
      say "Next exon start: ".$next_exon->seq_region_start;
      say "Intron length: ".$intron_length;
      unless($intron_length < 3) {
        push(@$updated_exons,$current_exon);
        $current_seq_offset += $current_exon->length;
        next;
      }

      # Note at the moment the seq edit stuff is only going to work for the two exon modification scenario
      # This should be fine based on the data so far
      say "Found an intron of length ".$intron_length.", will merge exons";
      my $merged_exon = $self->merge_exons($current_exon,$next_exon);
      $$exons[$i+1] = $merged_exon;
      push(@$updated_exons,$merged_exon);
      say "Merged exon: (".$merged_exon->seq_region_start."..".$merged_exon->seq_region_end.")";
      my $edit_seq;
      if($intron_length == 2) {
        $edit_seq = substr($transcript->seq->seq, -1);
        $edit_seq .= $edit_seq;
      } else {
        $edit_seq = substr($transcript->seq->seq, -2);
        $edit_seq .= $edit_seq;
      }

      my $edit_attribute_value = ($current_exon->length + $current_seq_offset)." ".($current_exon->length + $current_seq_offset)." ".$edit_seq;
      my $edit_attribute = Bio::EnsEMBL::Attribute->new(-CODE => '_rna_edit', -VALUE => $edit_attribute_value);
      say "RNA edit attribute generated. Value: ".$edit_attribute_value;
      push(@$transcript_attribs,$edit_attribute);
      $current_seq_offset += $merged_exon->length;
      $updated_transcript = 1;
      $updated_gene = 1;
    }

    unless($updated_transcript) {
      push(@$updated_transcripts,$transcript);
      next;
    }

    my $new_transcript = Bio::EnsEMBL::Transcript->new(-exons    => $updated_exons,
                                                       -slice    => $transcript->slice,
                                                       -analysis => $transcript->analysis);

    foreach my $attribute (@$transcript_attribs) {
      $new_transcript->add_Attributes($attribute);
    }

    say "FERGAL DEBUG 1";
    my $new_exons = $new_transcript->get_all_Exons;
    foreach my $new_exon (@$new_exons) {
      say "FERGAL EXON: ".$new_exon->start."..".$new_exon->end;
    }

#    say "FERGAL TRANSCRIPT SEQ:\n".$new_transcript->seq->seq;
#    compute_translation($new_transcript);
#    say "FERGAL DEBUG 1.5";

    $$updated_exons[0]->phase(0);
    my $translation = Bio::EnsEMBL::Translation->new(
    -START_EXON => $$updated_exons[0],
    -END_EXON   => $$updated_exons[$#$updated_exons],
    -SEQ_START  => 1,
    -SEQ_END    => $$updated_exons[$#$updated_exons]->length,
    );

    $new_transcript->translation($translation);

    push((@$updated_transcripts,$new_transcript));
    say "FERGAL DEBUG 2";
  }

  unless($updated_gene) {
    return(0);
  }

  my $new_gene = Bio::EnsEMBL::Gene->new(-slice    => $gene->slice,
                                         -analysis => $gene->analysis);

 foreach my $transcript (@$updated_transcripts) {
    $new_gene->add_Transcript($transcript);
  }

  return($new_gene);
}

sub merge_exons {
  my ($self,$left_exon, $right_exon) = @_;

  my $merged_exon = Bio::EnsEMBL::Exon->new(-start     => $left_exon->start,
                                            -end       => $right_exon->end,
                                            -strand    => $right_exon->strand,
                                            -phase     => $left_exon->phase,
                                            -end_phase => $right_exon->end_phase,
                                            -analysis  => $left_exon->analysis,
                                            -slice     => $left_exon->slice);

  return($merged_exon);
}


sub combine_files {
  my ($self) = @_;

  my $full_cds_hash = {};
  my $initial_cds_files = $self->param_required('cds_file');
  foreach my $cds_file (@$initial_cds_files) {
    unless(open(IN,$cds_file)) {
      $self->throw("Could not open cds file: ".$cds_file);
    }
    while(<IN>) {
      chomp $_;
      my @line = split(',',$_);
      my $header = $line[0];
      my ($seq_region,$type) = split(':',$header);
      if($self->param('use_generic_output_type')) {
        $type = $self->param('generic_output_type');
      }

      foreach(my $i=1; $i<scalar(@line); $i++) {
        my $cds_string = $line[$i];
        my ($start,$end,$strand,$length,$count) = split(':',$cds_string);
        my $cds = $start.':'.$end.':'.$strand.":".$length;
        if($full_cds_hash->{$seq_region}->{$type}->{$cds}) {
          $full_cds_hash->{$seq_region}->{$type}->{$cds} += $count;
        } else {
          $full_cds_hash->{$seq_region}->{$type}->{$cds} = $count;
	}
      }
    }
    close IN;
  }

  open(OUT,">".$self->param_required('combined_cds_file'));
  foreach my $seq_region (keys(%$full_cds_hash)) {
    foreach my $type (keys(%{$full_cds_hash->{$seq_region}})) {
      my $line = $seq_region.":".$type.",";
      foreach my $cds (keys(%{$full_cds_hash->{$seq_region}->{$type}})) {
        my $count = $full_cds_hash->{$seq_region}->{$type}->{$cds};
        $line .= $cds.":".$count.",";
      }
      say OUT $line;
    }
  }
  close OUT;
}


sub create_filename {
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


1;
