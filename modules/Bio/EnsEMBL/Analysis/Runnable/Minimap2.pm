=head1 LICENSE

 Copyright [2019-2020] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Minimap2

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Minimap2;

use warnings;
use strict;
use feature 'say';

use File::Spec;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_best_translation);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);
use Bio::EnsEMBL::Utils::Argument qw( rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [DECOMPRESS]           : String as a command like 'gzip -c -'
 Arg [EXPECTED_ATTRIBUTES]  : String specify the attribute expected for the output, see STAR manual
 Description                : Creates a  object to align reads to a genome using STAR
 Returntype                 : 
 Exceptions                 : Throws if WORKDIR does not exist
                              Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($genome_index, $input_file, $database_adaptor, $delete_input_file, $skip_introns_check, $skip_compute_translation, $sensitive, $bestn, $coverage, $perc_id, $max_intron_size) = rearrange([qw (GENOME_INDEX INPUT_FILE DATABASE_ADAPTOR DELETE_INPUT_FILE SKIP_INTRONS_CHECK SKIP_COMPUTE_TRANSLATION SENSITIVE BESTN COVERAGE PERC_ID MAX_INTRON_SIZE)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->database_adaptor($database_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->skip_introns_check($skip_introns_check);
  $self->skip_compute_translation($skip_compute_translation);
  $self->sensitive($sensitive);
  $self->secondary_alignments($bestn);
  $self->coverage_cutoff($coverage);
  $self->perc_id_cutoff($perc_id);
  $self->max_intron_size($max_intron_size);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run Star to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $sam_file = $self->create_filename(undef,'sam');
  $self->files_to_delete($sam_file);

  my $genome_index  = $self->genome_index;
  my $input_file    = $self->input_file;
  if($self->delete_input_file) {
    $self->files_to_delete($input_file);
  }

  my $options       = $self->options;

  my $max_intron_size = $self->max_intron_size();

  unless(defined($max_intron_size)) {
    $max_intron_size = 200000;
  }

  my $splice_type = "splice:hq";
  # run minimap2
  if($self->sensitive()) {
    $splice_type = "splice";
  }

  # By default we have secondary alignments off, but if we want them on then we just remove the flag to turn them off and put in the -N flag with the number of secondary alignments we want
  my $secondary_alignments = '--secondary=no';
  if($self->secondary_alignments()) {
    $secondary_alignments = '-Y -N '.$self->secondary_alignments();
  }

  my $minimap2_command = $self->program." --cs ".$secondary_alignments." -G ". $max_intron_size." -ax ".$splice_type." -u b ".$genome_index." ".$input_file." > ".$sam_file;
  warning("Command:\n$minimap2_command");
  execute_with_wait($minimap2_command);

  $self->output($self->parse_results($sam_file));

  # This is mostly a repeat of the above but on the reads that were filtered because they had a high non-canonical rate (but passed cov/identity)
  # These could be samples where the reads where accidently reversed as was seen in pig
  say "Found ".scalar(@{$self->leftover_genes})." leftover genes";
  if(scalar(@{$self->leftover_genes})) {
    my $leftover_input_file = $self->create_filename(undef,'lo');
    $self->create_leftover_input($leftover_input_file,$input_file);
    $self->files_to_delete($leftover_input_file);

    my $sam_lo_file = $self->create_filename(undef,'samlo');
    my $bed_lo_file = $self->create_filename(undef,'bedlo');
    $self->files_to_delete($sam_lo_file);
    $self->files_to_delete($bed_lo_file);


    my $minimap2_lo_command = $self->program." --cs ".$secondary_alignments." -G ". $max_intron_size." -ax ".$splice_type." -uf ".$genome_index." ".$leftover_input_file." > ".$sam_lo_file;
    warning("Leftover command:\n$minimap2_lo_command");
    execute_with_wait($minimap2_lo_command);

    $self->output($self->parse_results($sam_lo_file));
  }

}


=head2 parse_cigar_line

 Arg [1]    : String $hname, target name, the reference
 Arg [2]    : Int $hstart, the 5' position on the forward strand of the hit
 Arg [3]    : Int $qstrand, strand of the query sequence
 Arg [4]    : String $qname, query name, the read name
 Arg [5]    : Int $qlength, the length of the query sequence
 Arg [6]    : String $cigar_line, CIGAR string of the alignment
 Description: Parses the CIGAR string of the alignment to generate a set of ungapped
              objects. THe SAM CIGAR is switched to an Ensembl CIGAR.
              I operations become D operations
              D operations become introns
              S is used to set the start and end of the query sequence but is not kept
              = and X operations become M operations
              On the reverse strand the CIGAR string is reversed.
 Returntype : ArrayRef of Bio::EnsEMBL::DnaDnaAlignFeature
 Exceptions : None

=cut

sub parse_cigar_line {
  my ($hname, $hstart, $qstrand, $qname, $qlength, $cigar_line) = @_;

  my $qstart = 1;
  my $qend = 0;
  if ($qstrand == -1) {
    $qstart = $qlength;
    $qend = $qlength+1;
  }
  my $hend = $hstart-1;
  my @features;
  my @object_cigar;
  my $h_offset;
  while ($cigar_line =~ /(\d*)([MIDNSHP=X])/gc) {
    my $len = $1;
    if ($2 eq 'M') {
      push(@object_cigar, "$len$2");
      $qend += $len*$qstrand;
      $hend += $len;
    }
    elsif ($2 eq 'N' or $2 eq 'D') {
      if (@object_cigar and !(@object_cigar == 1 and index($object_cigar[0], 'D') > 0)) {
        push(@features, Bio::EnsEMBL::DnaDnaAlignFeature->new(
          -seqname    => $hname,
          -start      => $hstart,
          -end        => $hend,
          -strand     => $qstrand,
          -hseqname   => $qname,
          -hstart     => $qstrand == 1 ? $qstart : $qend,
          -hend       => $qstrand == 1 ? $qend : $qstart,
          -hstrand    => 1,
          -cigar_string => $qstrand == 1 ? join('', @object_cigar) : join('', reverse(@object_cigar)),
          -align_type => 'ensembl',
        ));
        @object_cigar = ();
        $qstart = $qend+(1*$qstrand);
        $hend += $len;
        $hstart = $hend+1;
      }
      else {
        $hend += $len;
        $hstart = $hend+1;
      }
    }
    elsif ($2 eq 'I') {
      $qend += $len*$qstrand;
      push(@object_cigar, "${len}D");
    }
    elsif ($2 eq 'S') {
      if ($qstrand == 1 and $qstart == 1) {
        $qstart += $len;
        $qend += $len;
      }
      elsif ($qstrand == -1 and $qend == $qlength+1) {
        $qstart -= $len;
        $qend -= $len;
      }
    }
    elsif ($2 eq 'H') {
      $h_offset = $len;
    }
    elsif ($2 eq '=' or $2 eq 'X') {
      push(@object_cigar, "$len$2");
      $qend += $len*$qstrand;
      $hend += $len;
    }
    elsif ($2 eq 'P') {

    }
  }
  if (@object_cigar) {
    push(@features, Bio::EnsEMBL::DnaDnaAlignFeature->new(
      -seqname    => $hname,
      -start      => $hstart,
      -end        => $hend,
      -strand     => $qstrand,
      -hseqname   => $qname,
      -hstart     => $qstrand == 1 ? $qstart : $qend,
      -hend       => $qstrand == 1 ? $qend : $qstart,
      -hstrand    => 1,
      -cigar_string => $qstrand == 1 ? join('', @object_cigar) : join('', reverse(@object_cigar)),
      -align_type => 'ensembl',
    ));
  }
  else {
    warning("Ending without a match after a gap $cigar_line");
  }
  if ($h_offset) {
    foreach my $feature (@features) {
      $feature->hstart($feature->hstart+$h_offset);
      $feature->hend($feature->hend+$h_offset);
    }
  }
  return \@features;
}


=head2 parse_minimap2_cs

 Arg [1]    : String $cs_line, the cs line generated by minimap2
 Arg [2]    : Int $seq_length, the length of the original query sequence
 Description: Parse the cs string to calculate the coverage and percentage of identity of the alignment
 Returntype : (Float, Float), $perc_identity, $coverage
 Exceptions : Throws if the values are not within 0 and 100

=cut

sub parse_minimap2_cs {
  my ($cs_line, $seq_length) = @_;

  my $missmatches = $cs_line =~ tr/*/*/;
  my $matches = 0;
  while ($cs_line =~ /:(\d+)/gc) {
    $matches += $1;
  }

  my $aligned_count = ($matches+$missmatches);
  my $percent_identity = sprintf("%.2f", (100*($matches/$aligned_count)));

  my $coverage = sprintf("%.2f", (100*($aligned_count/$seq_length)));
  unless(($percent_identity >= 0 && $percent_identity <= 100) &&
         ($coverage >= 0 && $coverage <= 100)) {
    throw("Issue with coverage/percent id calculation. Got values outside of expected range.".
                 "\nPercent id: ".$percent_identity."\nCoverage: ".$coverage);
  }
  return $percent_identity, $coverage;
}


=head2 parse_results

 Arg [1]    : String $output_file, path to the minimap2 output file
 Description: Parse the output file to generate gene models from the alignments
              When the alignments have too many non canonical alignments, they are
              realigned to the genome after being complemented
              Cut off for coverage and percentage of identity can be set, the defaults
              are 90%
              The original CIGAR string is added as a transcript supporting features
 Returntype : Arrayref of Bio::EnsEMBL::Gene
 Exceptions : Throw if it fails to open or close the result file

=cut

sub parse_results {
  my ($self, $output_file) = @_;

#208853  16      JAHAME010000038.1       505     60      1176M48I312M322N145M1273N97M5306N71M    *       0       0
#GTGATCTGCATGTGTGACACTGATTCTTTGGAAATAAAGAGTGGAAGCTGCAGGTGACACGTGAAGGGTTATTTATGGTTATGATGACCCTGTCCTGCAACGAGGGACTGGCAGCCACTACTGAGGAGGAGGGTCCCATCTCTCTCCTGTCGGCTTTCACCGAGGTCACAGCCAGACGTGGGGCAAAGGTGTTCCCTGTCCTACCCAGCCATTCCTGGGCCTGCCGCCTAGGGGCTCACAGGGCCCAGGAGTCCCCAGCTCACAGGCCAGGGCATCAGGCCAGGCGCGCTCGGTGCACACCGCACCTGGGAGGACCTGGGTACACTCAGGAGACCAAGAGCACTGGCGGGTCAGGATGGTTGGCGTTCAGCTCCTACGGGGTGGGGAGAAGTCTGTAGCCGAGAGCCCAGCCCCCTCCTGCCAGGTCTTCCCAGGTTCGAGAGAGGCTGGAACTCAATTTCTGCAGAAAATCTCCCAGTTTTTCCTGTTTGGTTAGTTTTTTTACAAAGACAGGATCTTGCTGTGTTGCCCAGGCTGGTCTTGAACTCCTGGCCTCAAGCAATCCTCCCACCTCGGCCTCCGAAAGTGCTGGGATTACTGGCATGAACCACTGCGCCCGGCTGGAGCTCCCGGTTTTTAAGCACTGCACGATACTAGAAGAGCTGACCTTTTTTCTGGCCTCACAGCTTATGCTGAAGCTGAGTGTGAGGAACAGAGAGACCTTTCTGTGACGACCGCTGGGGCAGAGTGGTCTATGCGCCGAGATCCTGGCATCAGCAAGGGAGGCGGGTCCTCGGGGAGGGGCAGCTTCCACAGTGTGGCTGCAGCGTGCACAGCCAGGTAGGCCCTGGATGTTCACCCCTCACTGCCCTTGGGGAAGCACCTGACCGCTGGGGATGTCCACCAGGGAGAGGACGCTGTGTCGGGGACAACATGCAGCATCAGCACCCACAAGGGCCCGGCCTGGCCCCGCCTCTCCACTCGCCCGAGGTCTTGCTGTGGCCCAAAGCAGGAGTCCAGGGCTGGCGAGACCTCCGGCTGCAGAAAGGCAGGCCCAGGCCTCACGATGGAGAAAGTCTGGATGTCCTGGTCTGGCCTGCTGTTTCATGCAGTGTGCAAGCAGCAGCATGAGCAGGCAATAGGCCAACTCGGCTGGAGGCACCAACTGAAGCCAAGGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCATGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACAGTGGGGCTGTCAGCTCTGGCTTGCAGGGTCACCGGCCTCACTGTCGCTGCTCTGAGACAGCTCTGTGGGGCTCTCAAGGCAGTGCAGCTCCAGGAAGCTGGTGATGAGCTTGGCGATGGCTTCCACATCACTGCGGTGGCCCATGACTGCACTCTGCGTGAGGGGGTGGGAGAAGCCCGTCGCAAGCCGGAGCACCACCATGTAGCCTTTCCCGAAGTACCGGACCTTCTCCTCCTCCACGCTCACATCACGGACATCATGGAGCAGGACCACCACCTGGTCGTGGCCAGCTCTGAAAAGAGTCAGCAGCTTCTTGTAGAGGCTGAACGTCTTCAAAACAACCTTCCCTGTGCTCTTGTCGAAGATGGCTTCCGAGAGCACAGCGGTCCCCGCGCCGCAGCAGAGAGACGCGAACCGGCGGGGGCGGGGCCGCGCGCACTTCC
#*       NM:i:52 ms:i:1727       AS:i:1685       nn:i:0  ts:A:+  tp:A:P  cm:i:317        s1:i:1716       s2:i:83 de:f:0.0028     cs:Z::161*cg:593*ct:71*tc:348+gccacgggtcctgtgggcggtgcacgggctcaggcacaccgggaatgt:9*tc:302~ct322ac:145~ct1273ac:97~ct5306ac:71     rl:i:0
  my $percent_id_cutoff = $self->perc_id_cutoff();
  my $coverage_cutoff = $self->coverage_cutoff();
  unless(defined($percent_id_cutoff)) {
    $percent_id_cutoff = 90;
  }

  unless(defined($coverage_cutoff)) {
    $coverage_cutoff = 90;
  }

  my $canonical_intron_cutoff = 0.8;

  my $genes = [];
  my @leftover_genes;

  my %seq_lengths;
  open(IN, $output_file) or throw("Could not open $output_file");
  while(my $line = <IN>) {
    next if (index($line, '@') == 0);
    my @results = split("\t",$line, 12);
    my $flags = $results[1];
    next if ($flags&0x4);
    my $query_name = $results[0];
    my $hit_name = $results[2];
    my $hit_start = $results[3];
    my $mapping_quality = $results[4];
    my $cigar_line = $results[5];
    my $second_read_hit_name = $results[6];
    my $second_read_start = $results[7];
    my $tlen = $results[8];
    my $query_seq = $results[9];
    my $query_quality = $results[10];
    my $attributes = $results[11];
    my $strand = 1;
    if ($flags&0x10) {
      $strand = -1;
    }
    $seq_lengths{$query_name} = length($query_seq) unless ($flags&0x100);
    my $percent_identity = 100;
    my $coverage = 100;
    if ($attributes) {
      my ($cs_line) = $attributes =~ /cs:Z:(\S+)/;
      if ($cs_line) {
        ($percent_identity, $coverage) = parse_minimap2_cs($cs_line, $seq_lengths{$query_name});
        if ($percent_identity < $percent_id_cutoff) {
          warning("Percent id for the hit fails the cutoff.\nHit name: $hit_name\nPercent id: $percent_identity\nCut-off: $percent_id_cutoff");
          next;
        }

        if ($coverage < $coverage_cutoff) {
          warning("Coverage for the hit fails the cutoff.\nHit name: $hit_name\nCoverage: $coverage\nCut-off: $coverage_cutoff");
          next;
        }
      }
    }
    my $cigar_objects = parse_cigar_line($hit_name, $hit_start, $strand, $query_name, $seq_lengths{$query_name}, $cigar_line);
    my $slice;
    if ($self->query) {
      $slice = $self->query;
    }
    else {
      $slice = $self->slice_cache($hit_name);
    }
    if (!ref($slice)) {
      $slice = $self->database_adaptor->get_SliceAdaptor->fetch_by_region('toplevel', $hit_name);
    }

    my @exons;
    foreach my $sf (@$cigar_objects) {
      push(@exons, $self->create_exon($slice, $sf->start, $sf->end, $strand));
      $sf->analysis($self->analysis);
      $sf->slice($slice);
      $sf->percent_id($percent_identity);
      $sf->hcoverage($coverage);
      $exons[-1]->add_supporting_features($sf);
    }

    if($strand == -1) {
      @exons = reverse(@exons);
    }

    my $gene = $self->create_gene(\@exons,$slice,$query_name);

    $gene->get_all_Transcripts->[0]->add_supporting_features(Bio::EnsEMBL::DnaDnaAlignFeature->new(
      -seqname    => $hit_name,
      -start      => $gene->start,
      -end        => $gene->end,
      -strand     => $gene->strand,
      -hseqname   => $gene->stable_id,
      -hstart     => 1,
      -hend       => length($query_seq),
      -hstrand    => 1,
      -cigar_string => $cigar_line,
      -align_type => 'cigar',
      -slice => $slice,
      -analysis => $self->analysis,
    ));
    # We aren't going to store a supporting feature, but we can store the coverage and percent id on the gene
    $coverage = int($coverage);
    $gene->version($coverage);
    $gene->description($percent_identity);

    unless($self->skip_introns_check()) {
      my $transcript = ${$gene->get_all_Transcripts}[0];

      my $introns = $transcript->get_all_Introns;
      my $intron_count = scalar(@$introns);
      if($intron_count) {
        my $canonical_count = 0;
        foreach my $intron (@$introns) {
          if($intron->is_splice_canonical) {
            $canonical_count++;
	        }
        }

        # If it fails the canonical cutoff then we skip this gene, but in case it's just a stranded issue we put onto the leftover pile
        unless($canonical_count/$intron_count >= $canonical_intron_cutoff) {
          say "Gene fails canonical splice site cut-off";
          push(@leftover_genes,$gene);
          next;
        }
      }
   } # End $self->skip_introns_check()
    push(@$genes,$gene);
  }
  close(IN) or throw("Could not close $output_file");

  say "Finished parsing output";
  if (@leftover_genes) {
    $self->leftover_genes(\@leftover_genes);
  }
  return($genes);
}



sub create_leftover_input {
  my ($self,$leftover_input_file,$original_input_file) = @_;

  my $index = Bio::DB::HTS::Faidx->new($original_input_file);

  my $read_ids = [];
  unless(open(OUTLEFTOVER,">".$leftover_input_file)) {
    throw("Failed to open output file to use as input for leftover run");
  }

  foreach my $gene (@{$self->leftover_genes}) {
    my $read_id = $gene->stable_id;
    my $length = $index->length($read_id);
    my $location  = $read_id.":1-".$length;
    my $seq = $index->get_sequence_no_length($location);
    $seq = reverse($seq);
    $seq =~ tr/atgcATGC/tacgTACG/;
    say OUTLEFTOVER ">".$read_id;
    say OUTLEFTOVER $seq;
  }
  close OUTLEFTOVER or throw("Could not close '$leftover_input_file'");

}


sub create_exon {
  my ($self,$slice,$exon_start,$exon_end,$strand) = @_;

  my $exon = Bio::EnsEMBL::Exon->new(-start     => $exon_start,
                                     -end       => $exon_end,
                                     -strand    => $strand,
                                     -phase     => -1,
                                     -end_phase => -1,
                                     -analysis  => $self->analysis,
                                     -slice     => $slice);

  if($exon_start > $exon_end) {
    throw("Created an exon where the start > than the end, this shouldn't be possible: ".$slice->name." ".$exon->start."..".$exon->end." ".$strand);
  }
  return($exon);
}


sub create_gene {
  my ($self,$exons,$slice,$hit_name) = @_;

  my $transcript = Bio::EnsEMBL::Transcript->new(-exons    => $exons,
                                                 -slice    => $slice,
                                                 -analysis => $self->analysis);

  unless($self->skip_compute_translation()) {
    compute_best_translation($transcript);
  }

  my $gene = Bio::EnsEMBL::Gene->new(-slice    => $slice,
                                     -analysis => $self->analysis);

  $transcript->biotype('cdna');
  $gene->biotype('cdna');
  $transcript->stable_id($hit_name);
  $gene->stable_id($hit_name);

  $gene->add_Transcript($transcript);

  return($gene);
}


sub genome_index {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genome_index} = $val;
  }

  return $self->{_genome_index};
}


sub input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_input_file} = $val;
  }

  return $self->{_input_file};
}


sub database_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_database_adaptor} = $val;
  }

  return $self->{_database_adaptor};
}


sub delete_input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_delete_input_file} = $val;
  }

  return $self->{_delete_input_file};
}


sub skip_introns_check {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_skip_introns_check} = $val;
  }

  return $self->{_skip_introns_check};
}


sub skip_compute_translation {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_skip_compute_translation} = $val;
  }

  return $self->{_skip_compute_translation};
}


sub sensitive {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_sensitive} = $val;
  }

  return $self->{_sensitive};
}


sub secondary_alignments {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_secondary_alignments} = $val;
  }

  return $self->{_secondary_alignments};
}


sub coverage_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_coverage_cutoff} = $val;
  }

  return $self->{_coverage_cutoff};
}


sub perc_id_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_perc_id_cutoff} = $val;
  }

  return $self->{_perc_id_cutoff};
}


sub max_intron_size {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_max_intron_size} = $val;
  }

  return $self->{_max_intron_size};
}


sub leftover_genes {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_leftover_genes} = $val;
  }

  return $self->{_leftover_genes} || [];
}


=head2 slice_cache

 Arg [1]    : Hashref $val, a cache of Bio::EnsEMBL::Slice
 Description: Cache of Slices to avoid too many call to the Mysql server
              WIP
 Returntype : Hashref of Bio::EnsEMBL::Slice
 Exceptions : None

=cut

sub slice_cache {
  my ($self, $val) = @_;
}

1;
