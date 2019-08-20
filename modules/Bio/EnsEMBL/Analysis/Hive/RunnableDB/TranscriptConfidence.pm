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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::SpliceSiteCorrection

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::SpliceSiteCorrection->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module takes a set of input genes and a file containing intron frequency counts and then
edits exons where the corresponding intron bounrdy is not frequenctly observed. Everything,
whether edited or unedited is by default written to an output db

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::TranscriptConfidence;

use warnings;
use strict;
use feature 'say';

use Data::Dumper;

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name is_canonical_splice);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    small_intron_size        => 75,
    low_frequency            => 5,
    low_frequency_cds        => 5,
    low_frequency_single_exon_cds        => 10,
    low_frequency_stop_start => 10,
    min_cds_exon_size        => 10,
    max_cds_exon_size        => 20000,
    min_intron_ratio         => 0.1,
    min_cds_length           => 300,
    write_unmodified         => 1,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for intron dumping
 Returntype : None

=cut

sub fetch_input {
  my ($self) = @_;

# Ideas:
# instead of passing in gene ids, could pass in a genic slice (stranded). Then cluster
# and for each cluster look at all orf start/ends. This would allow frequency data for
# start/ends within a cluster to be worked out and could therefore enable a low frequency
# ratio to be worked out. Could have various problems though
# Should also do something in regards to single versus multi exon ORFs (cds length?)
# Process cds exons? If a high frequency intron overlaps a cds exon, then mark as low
# confidence. Similarly if an exon is unusually long, should it be low confidence?
# Certainly if > 20kb then yes as no known cds exon is that big


  say "Fetching input";

  my $intron_file = $self->param_required('combined_intron_file');
  unless(open(IN,$intron_file)) {
    $self->throw("Could not open intron file on specified path. Path: ".$intron_file);
  }

  # This bit will process the hash of introns
  my $intron_hash = {};
  my $combined_intron_hash = {};
  while(<IN>) {
    my $line = $_;
    chomp $line;
    $self->process_introns($line,$intron_hash,$combined_intron_hash);
  }
  close IN;

  $self->param('intron_hash',$intron_hash);
  $self->param('combined_intron_hash',$combined_intron_hash);


  my $cds_file = $self->param_required('combined_cds_file');
  unless(open(IN,$cds_file)) {
    $self->throw("Could not open cds file on specified path. Path: ".$cds_file);
  }

  # This bit will process the hash of introns
  my $cds_hash = {};
  my $combined_cds_hash = {};
  my $combined_start_hash = {};
  my $combined_end_hash = {};
  while(<IN>) {
    my $line = $_;
    chomp $line;
    $self->process_cds_coords($line,$cds_hash,$combined_cds_hash,$combined_start_hash,$combined_end_hash);
  }
  close IN;

  $self->param('cds_hash',$cds_hash);
  $self->param('combined_cds_hash',$combined_cds_hash);
  $self->param('combined_start_hash',$combined_start_hash);
  $self->param('combined_end_hash',$combined_end_hash);

  my $input_dba = $self->hrdb_get_dba($self->param('input_db'));
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
    $input_dba->dnadb($dna_dba);
  }

  # If we are dumping introns from genes the get the genes
  my $gene_adaptor = $input_dba->get_GeneAdaptor();
  my $gene_ids = $self->input_id();
#  my $gene_ids = [173011];
#  my $gene_ids= [172985];
#  my $gene_ids = [43161];
  my $genes = [];

  # Note that the original script in Ensembl common has some flagging here to reduce the set. I have
  # left it out for the the moment so that every gene gets porcessed
  foreach my $gene_id (@$gene_ids) {
    push(@$genes,$gene_adaptor->fetch_by_dbID($gene_id));
  }

  say "Fetched ".scalar(@$genes)." input genes";
  $self->param('input_genes',$genes);


  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));
  if($dna_dba) {
    $target_dba->dnadb($dna_dba);
  }
  $self->hrdb_set_con($target_dba,'target_db');

  my $target_slice_adaptor = $target_dba->get_SliceAdaptor;
  $self->param('target_slice_adaptor',$target_slice_adaptor);
}


sub run {
  my ($self) = @_;

  my $genes = $self->param('input_genes');
  my $intron_hash = $self->param('intron_hash');
  my $combined_intron_hash = $self->param('combined_intron_hash');
  my $target_slice_adaptor = $self->param('target_slice_adaptor');
  my $processed_output_genes = [];
  foreach my $gene (@$genes) {
    $self->assess_confidence($gene);
#    if($self->assess_confidence($gene)) {
#      $gene->biotype($gene->biotype().'_high');
#    } else {
#      $gene->biotype($gene->biotype().'_low');
#    }
  }
#  $self->output($processed_output_genes);
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_db = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $target_db->get_GeneAdaptor();

  say "Finished processing. Preparing to write genes to db";
  foreach my $gene (@{$self->output}) {
    my $confidence = $gene->{'confidence'};
    unless($confidence) {
      $self->throw("Could not find confidence assignment for gene with dbID ".$gene->dbID.", this should not happen");
    }
    $gene->biotype($gene->biotype()."_".$confidence);
    $gene_adaptor->update($gene);
  }
}


sub process_introns {
  my ($self,$line,$intron_hash,$combined_intron_hash) = @_;

  my @line = split(',',$line);
  my $header = $line[0];
  my ($seq_region,$type) = split(':',$header);

  foreach(my $i=1; $i<scalar(@line); $i++) {
    my $intron_details = $line[$i];
    my ($start,$end,$strand,$count) = split(':',$intron_details);
    my $intron_string = $start.':'.$end.':'.$strand;

    if($combined_intron_hash->{$seq_region}->{$intron_string}) {
      $combined_intron_hash->{$seq_region}->{$intron_string} += $count;
    } else {
      $combined_intron_hash->{$seq_region}->{$intron_string} = $count;
    }

    if($intron_hash->{$seq_region}->{$type}->{$intron_string}) {
      $self->warning("An entry for the intron already existed, the file is expected to have unique entries at this point. Skipping this intron entry");
      next;
    } else {
      $intron_hash->{$seq_region}->{$type}->{$intron_string} = $count;
    }
  }
}


sub process_cds_coords {
  my ($self,$line,$cds_hash,$combined_cds_hash,$combined_start_hash,$combined_end_hash) = @_;

  # Might update this later to diverge from the process_introns code, which is why it's separate
  # A possible thing to include would be a count of the coding exons the cds spanned, to show
  # different cds isoforms
  my @line = split(',',$line);
  my $header = $line[0];
  my ($seq_region,$type) = split(':',$header);

  foreach(my $i=1; $i<scalar(@line); $i++) {
    my $cds_details = $line[$i];
    my ($start,$end,$strand,$length,$count) = split(':',$cds_details);
    my $cds_string = $start.':'.$end.':'.$strand.":".$length;
    my $start_string = $start.':'.$strand;
    my $end_string = $end.':'.$strand;

    if($combined_cds_hash->{$seq_region}->{$cds_string}) {
      $combined_cds_hash->{$seq_region}->{$cds_string} += $count;
    } else {
      $combined_cds_hash->{$seq_region}->{$cds_string} = $count;
    }

    if($combined_start_hash->{$seq_region}->{$start_string}) {
      $combined_start_hash->{$seq_region}->{$start_string} += $count;
    } else {
      $combined_start_hash->{$seq_region}->{$start_string} = $count;
    }

    if($combined_end_hash->{$seq_region}->{$end_string}) {
      $combined_end_hash->{$seq_region}->{$end_string} += $count;
    } else {
      $combined_end_hash->{$seq_region}->{$end_string} = $count;
    }

    if($cds_hash->{$seq_region}->{$type}->{$cds_string}) {
      $self->warning("An entry for the cds already existed, the file is expected to have unique entries at this point. Skipping this cds entry");
      next;
    } else {
      $cds_hash->{$seq_region}->{$type}->{$cds_string} = $count;
    }
  }
}


sub assess_confidence {
  my ($self,$gene) = @_;

  my $intron_confidence = $self->intron_confidence($gene);
  if($intron_confidence) {
    say "Confident on the introns for the cds";
  }

  my $cds_confidence = $self->cds_confidence($gene);
  if($cds_confidence) {
    say "Confident on the cds";
  }

  # Note that single exon genes have a separate cds confidence cut off that should be set higher for something
  # to be called high confidence
  if($intron_confidence && $cds_confidence || $cds_confidence && $gene->{'is_single'}) {
    $gene->{'confidence'} = 'high';
    $self->output([$gene]);
  } elsif($intron_confidence) {
    $gene->{'confidence'} = 'int';
    $self->output([$gene]);
  } elsif($cds_confidence) {
    $gene->{'confidence'} = 'cds';
    $self->output([$gene]);
  }

}

sub intron_confidence {
  my ($self,$gene) = @_;

  say "Assessing intron confidence for cds";
  my $slice_adaptor = $self->param('target_slice_adaptor');
  my $small_intron_size = $self->param_required('small_intron_size');
  my $low_frequency = $self->param_required('low_frequency');
  my $intron_hash = $self->param('intron_hash');
  my $combined_intron_hash = $self->param('combined_intron_hash');

  my $min_intron_ratio = $self->param_required('min_intron_ratio');
  my $seq_region = $gene->seq_region_name;
  my $hash_intron_types = $intron_hash->{$seq_region};

  my $priority = {'cdna'      => 5,
                  'rnaseq'    => 4,
                  'long_read' => 3,
                  'projecton' => 2,
                  'protein'   => 1};

  my $transcripts = $gene->get_all_Transcripts();

  # This will unset if any transcript is multi exon
  $gene->{'is_single'} = 1;
  foreach my $transcript (@$transcripts) {
    # creating this because we're really only interested in cds introns and the API call for getting all cds introns
    # is exteremly costly cos of joins
    my $cds_transcript = $self->create_cds_transcript($transcript);
    my $introns = $cds_transcript->get_all_Introns();
    my $high_confidence = 1;

    foreach my $intron (@$introns) {
     $gene->{'is_single'} = 0;
      my $intron_start = $intron->seq_region_start();
      my $intron_end = $intron->seq_region_end();
      my $intron_strand = $intron->strand;
      my $intron_string = $intron_start.":".$intron_end.":".$intron_strand;
      say "FERGAL INT STRING: ".$intron_string;
      my $intron_count = 1;
      if($combined_intron_hash->{$seq_region}->{$intron_string}) {
        $intron_count = $combined_intron_hash->{$seq_region}->{$intron_string};
        say "FOUND INT STRING: ".$intron_count;
      }

      my ($is_canonical,$donor,$acceptor) = is_canonical_splice($intron,$slice_adaptor,$gene->slice);

#      my $observations = 1;

#      foreach my $type_key (keys(%$hash_intron_types)) {
#        my $hash_introns = $hash_intron_types->{$type_key};
#        foreach my $hash_intron (keys(%$hash_introns)) {
#          if($intron_string eq $hash_intron) {
#            $observations++;
#          }
#        } # foreach my $hash_intron (keys(%$hash_introns))
#      } # foreach my $type_key

      # This is very strict at the moment cos if there's a single dodgy cds intron it will call the whole thing
      # low confidence. Should try and normalise for number of introns, maybe 95 percent, then one intron in 20
      # could be wrong. That would require that the average transcript was fully supported while allowing for
      # small issues in longer ones
      if($intron_count <= $low_frequency || !$is_canonical || $intron->length <= $small_intron_size) {
        return(0);
      }
    } # foreach my $intron (@$introns)
  } # end foreach my $transcript (@$transcripts)

  if($gene->{'is_single'}) {
    return(0);
  }

  return(1);
}


sub cds_confidence {
  my ($self,$gene) = @_;

  say "Assessing cds confidence";
  my $low_frequency_cds = $self->param_required('low_frequency_cds');
  my $low_frequency_single_exon_cds = $self->param_required('low_frequency_single_exon_cds');
  my $low_frequency_stop_start = $self->param_required('low_frequency_stop_start');
  my $min_cds_length = $self->param_required('min_cds_length');

  my $cds_hash = $self->param('cds_hash');
  my $combined_cds_hash = $self->param('combined_cds_hash');

  my $combined_start_hash = $self->param('combined_start_hash');
  my $combined_end_hash   = $self->param('combined_end_hash');

  my $cds_strand = $gene->strand();
  my $seq_region = $gene->seq_region_name();


  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    unless($transcript->translation) {
      say "Found no translation so can't examine a cds";
      return(0);
    }

    my $cds_start = $transcript->coding_region_start();
    my $cds_end = $transcript->coding_region_end();
    my $cds_length = length($transcript->translateable_seq);
    unless($cds_length >= $min_cds_length) {
      say "cds is shorted than the min length, length found was ".$cds_length.", length required ".$min_cds_length;
      return(0);
    }

    my $cds_string = $cds_start.":".$cds_end.":".$cds_strand.":".$cds_length;
    my $start_string = $cds_start.":".$cds_strand;
    my $end_string = $cds_end.":".$cds_strand;

    say "FERGAL cds string: ".$cds_string;
    my $cds_count = $combined_cds_hash->{$seq_region}->{$cds_string};
    unless($cds_count) {
      $self->warning("Couldn't fine the cds string in tbe combined cds hash, this is unusual unless the cds set does not include the gene being assessed");
      return(0);
    }

    if($gene->{'is_single'} && ($cds_count < $low_frequency_single_exon_cds)) {
      # A single exon cds would need to be well supported in terms of the cds
      # to have confidence in it
      return(0);
    } elsif(!$gene->{'is_single'} && $cds_count < $low_frequency_cds) {
      my $cds_start_count = $combined_start_hash->{$seq_region}->{$start_string};
      my $cds_end_count = $combined_end_hash->{$seq_region}->{$end_string};
      unless($cds_start_count >= $low_frequency_stop_start && $cds_end_count >= $low_frequency_stop_start) {
        return(0);
      }
    }
  }
  return(1);
}


sub exon_confidence {
  my ($self,$gene) = @_;

  say "Assessing exon confidence for cds";
  my $combined_intron_hash = $self->param('combined_intron_hash');

  my $seq_region = $gene->seq_region_name;
  my $max_cds_exon_size = $self->param_required('max_cds_exon_size');
  my $min_cds_exon_size = $self->param_required('min_cds_exon_size');

  my $overlapping_introns = $self->get_overlapping_introns($gene,$combined_intron_hash);
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    # creating this because we're really only interested in cds introns and the API call for getting all cds introns
    # is exteremly costly cos of joins
    my $cds_transcript = $self->create_cds_transcript($transcript);
    foreach my $exon (@{$cds_transcript->get_all_Exons()}) {
      if($exon->length > $max_cds_exon_size || $exon->length < $min_cds_exon_size) {
        return(0);
      }
      if($self->check_exon_overlap($exon,$overlapping_introns)) {
        return(0);
      }
    }
  }
  return(1)
}


sub get_overlapping_introns {
  my ($self,$gene,$combined_intron_hash) = @_;

  my $overlapping_introns = {};
  my $gene_seq_region = $gene->seq_region_name();
  my $gene_strand = $gene->strand();
  my $gene_start = $gene->seq_region_start();
  my $gene_end = $gene->seq_region_end();
  my $low_frequency = $self->param_required('low_frequency');
  foreach my $intron_string (keys(%{$combined_intron_hash->{$gene_seq_region}})) {
    my ($start,$end,$strand) = split(':',$intron_string);
    my $count = $combined_intron_hash->{$gene_seq_region}->{$intron_string};
    unless($strand == $gene_strand && $count > $low_frequency) {
      next;
    }

    if ($start <= $gene_end && $end >= $gene_start) {
      $overlapping_introns->{$intron_string} = $count;
    }
  }
  return($overlapping_introns);
}


sub check_exon_overlap {
  my ($self,$exon,$overlapping_introns) = @_;

  foreach my $intron_string (keys(%{$overlapping_introns})) {
    my ($start,$end,$strand) = split(':',$intron_string);
    if ($start <= $exon->seq_region_end && $end >= $exon->seq_region_start) {
      return(1);
    }
  }
  return(0);
}

sub create_cds_transcript {
  my ($self,$transcript) = @_;

  my $cds_transcript = Bio::EnsEMBL::Transcript->new();
  my $cds_exons = $transcript->get_all_translateable_Exons;
  foreach my $cds_exon (@$cds_exons) {
    $cds_transcript->add_Exon($cds_exon);
  }

  $cds_transcript->biotype($transcript->biotype);
  $cds_transcript->analysis($transcript->analysis);
  $cds_transcript->slice($transcript->slice);


  compute_translation($cds_transcript);
#  my $start_exon = ${$cds_exons}[0];
 # my $end_exon = ${$cds_exons}[scalar(@$cds_exons)-1];
#  my $translation = Bio::EnsEMBL::Translation->new();
#  $translation->start_Exon($start_exon);
#  $translation->start(1);
#  $translation->end_Exon($end_exon);
#  $translation->end($end_exon->length());
#  foreach my $seq_edit (@{$transcript->translation->get_all_SeqEdits}) {
#    $translation->add_Attributes($seq_edit->get_Attribute);
#  }
#  $cds_transcript->translation($translation);
#  calculate_exon_phases($cds_transcript, 0);

  $cds_transcript->{'5_prime_utr'} = 0;
  $cds_transcript->{'3_prime_utr'} = 0;

  return($cds_transcript);
}


sub print_exon {
  my ($self,$exon) = @_;

  my $exon_string = "(".$exon->seq_region_start."..".$exon->seq_region_end."):".$exon->strand;
  say $exon_string;
}


1;
