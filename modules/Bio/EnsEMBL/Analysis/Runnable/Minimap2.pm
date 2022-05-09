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

Bio::EnsEMBL::Analysis::Runnable::Star

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Star->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence. Star is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

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
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(is_canonical_splice align_proteins execute_with_wait);
use Bio::EnsEMBL::Utils::Argument qw( rearrange);
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

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
  my ($genome_index, $input_file, $paftools_path, $database_adaptor, $delete_input_file, $skip_introns_check, $add_offset, $skip_compute_translation, $sensitive, $bestn, $coverage, $perc_id, $max_intron_size) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH DATABASE_ADAPTOR DELETE_INPUT_FILE SKIP_INTRONS_CHECK ADD_OFFSET SKIP_COMPUTE_TRANSLATION SENSITIVE BESTN COVERAGE PERC_ID MAX_INTRON_SIZE)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->database_adaptor($database_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->skip_introns_check($skip_introns_check);
  $self->add_offset($add_offset);
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

  my $paftools_path = $self->paftools_path;
  my $options       = $self->options;

  my $max_intron_size = $self->max_intron_size();

  unless(defined($max_intron_size)) {
    $max_intron_size = 200000;
  }

  unless($paftools_path) {
    $self->throw("Paftools path was empty");
  }

  my $splice_type = "splice:hq";
  # run minimap2
  if($self->sensitive()) {
    $splice_type = "splice";
  }

  # By default we have secondary alignments off, but if we want them on then we just remove the flag to turn them off and put in the -N flag with the number of secondary alignments we want
  my $secondary_alignments = '--secondary=no';
  if($self->secondary_alignments()) {
    $secondary_alignments = '-N '.$self->secondary_alignments();
  }

  my $minimap2_command = $self->program." --cs ".$secondary_alignments." -G ". $max_intron_size." -ax ".$splice_type." -u b ".$genome_index." ".$input_file." > ".$sam_file;
  $self->warning("Command:\n".$minimap2_command."\n");
  execute_with_wait($minimap2_command);

  $self->output($self->parse_results($sam_file));

  # This is mostly a repeat of the above but on the reads that were filtered because they had a high non-canonical rate (but passed cov/identity)
  # These could be samples where the reads where accidently reversed as was seen in pig
  say "Found ".scalar(@{$self->leftover_genes})." leftover genes";
  if(scalar(@{$self->leftover_genes})) {
    my $leftover_input_file = $self->create_filename(undef,'lo');
    $self->create_leftover_input($self->leftover_genes,$leftover_input_file,$input_file);
    $self->files_to_delete($leftover_input_file);

    my $sam_lo_file = $self->create_filename(undef,'samlo');
    my $bed_lo_file = $self->create_filename(undef,'bedlo');
    $self->files_to_delete($sam_lo_file);
    $self->files_to_delete($bed_lo_file);


    my $minimap2_lo_command = $self->program." --cs ".$secondary_alignments." -G ". $max_intron_size." -ax ".$splice_type." -uf ".$genome_index." ".$leftover_input_file." > ".$sam_lo_file;
    $self->warning("Leftover command:\n$minimap2_lo_command\n");
    execute_with_wait($minimap2_lo_command);

    $self->output($self->parse_results($sam_lo_file));
  }

}


sub parse_sam {
  my ($self,$sam_file,$percent_id_hash,$coverage_hash) = @_;

  unless(open(IN,$sam_file)) {
    $self->throw("Could not open the sam file for processing. Path:\n".$sam_file);
  }

  while(<IN>) {
    my $line = $_;
    if($line =~ /^@/) {
      next;
    }

    my @results = split("\t",$line);
    my $num_cols = scalar(@results);

    # Column number is variable. Should probably just use one of the column identifiers that is always present
    unless($num_cols >= 20) {
      $self->warning("Unexpected number of result columns, skipping line. Line:\n".$line."Number of cols: ".$num_cols);
      next;
    }

    my $read_name = $results[0];
    # Sometimes this does not have any seq for some reason, so this value becomes 1 incorrectly. This
    # is dealt with in the coverage calc later on. Really should look up the faidx for the seq if this
    # col does not have it
    my $seq_length = length($results[9]);

    # The index of this varies because the columns vary with each result
    my $cs;
    for(my $i=0; $i<scalar(@results); $i++) {
      if($results[$i] =~ /^cs\:Z\:/) {
        $cs = $results[$i];
        last;
      }
    }

    unless($cs) {
      $self->throw("CS column not parsed successfully. Line contents:\n".$line);
    }

    my $mismatch_count = () = $cs =~ /\*/gi;
    my $match_count = 0;
    while($cs =~ s/\:(\d+)//) {
      $match_count += $1;
    }

    my $aligned_count = ($match_count + $mismatch_count);

    my $percent_identity = 100 * ($match_count / $aligned_count);
    $percent_identity = sprintf("%.2f",$percent_identity);

    # Sometimes the read isn't included in the output for some reason. We could look it up in the file, though this could
    # be a little slow if the file is big
    # Will probably add this is later
    if($aligned_count > $seq_length) {
      $self->warning("The number of aligned bases listed in the cs:Z is greater than calculated seq length. Likely that the sequence ".
                     "was not included in the sam (represented by *?). Will set to same value as aligned bases. This will make the coverage ".
                     "100 percent. Seq column entry:\n".$results[9]);
      $seq_length = $aligned_count;
    }

    my $coverage = 100 * ($aligned_count / $seq_length);
    $coverage = sprintf("%.2f",$coverage);
    unless(($percent_identity >= 0 && $percent_identity <= 100) &&
           ($coverage >= 0 && $coverage <= 100)) {
      $self->throw("Issue with coverage/percent id calculation. Got values outside of expected range.".
                   "\nPercent id: ".$percent_identity."\nCoverage: ".$coverage."\nRead id: ".$read_name);
    }

    unless(exists $percent_id_hash->{$read_name}) {
      $percent_id_hash->{$read_name} = $percent_identity;
      $coverage_hash->{$read_name} = $coverage;
    } else {
      $self->warning("Found two result lines for a read. Only calculating percent id for the first one. ID: ".$read_name);
    }
  }
  close IN;
}


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
  print __LINE__, " $qstart $qend $qstrand $hstart $hend\n";
  while ($cigar_line =~ /(\d*)([MIDNSHP=X])/gc) {
    my $len = $1;
    if ($2 eq 'M') {
      push(@object_cigar, "$len$2");
      $qend += $len*$qstrand;
      $hend += $len;
      print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
    }
    elsif ($2 eq 'N' or $2 eq 'D') {
      print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
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
#          -cigar_string => join('', @object_cigar),
          -cigar_string => $qstrand == 1 ? join('', @object_cigar) : join('', reverse(@object_cigar)),
          -align_type => 'ensembl',
        ));
        @object_cigar = ();
        $qstart = $qend+(1*$qstrand);
        $hend += $len;
        $hstart = $hend+1;
        print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
      }
      else {
        $hend += $len;
#        $features[-1]->end($hend);
        $hstart = $hend+1;
        print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
      }
    }
    elsif ($2 eq 'I') {
      $qend += $len*$qstrand;
      push(@object_cigar, "${len}D");
      print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
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
      print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
    }
    elsif ($2 eq 'H') {
    }
    elsif ($2 eq '=' or $2 eq 'X') {
      push(@object_cigar, "$len$2");
      $qend += $len*$qstrand;
      $hend += $len;
      print __LINE__, " $len$2 $qstart $qend $qstrand $hstart $hend\n";
    }
    elsif ($2 eq 'P') {

    }
  }
  if (@object_cigar) {
    print __LINE__, " $qstart $qend $qstrand $hstart $hend\n";
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
  return \@features;
}


sub parse_minimap2_cs {
  my ($cs_line, $seq_length) = @_;

  my $missmatches = $cs_line =~ tr/*/*/;
  my $matches = 0;
  while ($cs_line =~ /:(\d+)/gc) {
    $matches += $1;
  }

  my $aligned_count = ($matches+$missmatches);
  my $percent_identity = sprintf("%.2f", (100*($matches/$aligned_count)));
  print __LINE__, " $matches $missmatches $aligned_count $percent_identity\n";

  my $coverage = sprintf("%.2f", (100*($aligned_count/$seq_length)));
  print __LINE__, " $aligned_count $seq_length $coverage\n";
  unless(($percent_identity >= 0 && $percent_identity <= 100) &&
         ($coverage >= 0 && $coverage <= 100)) {
    throw("Issue with coverage/percent id calculation. Got values outside of expected range.".
                 "\nPercent id: ".$percent_identity."\nCoverage: ".$coverage);
  }
  return $percent_identity, $coverage;
}

sub parse_results {
  my ($self, $output_file) = @_;

#208853  16      JAHAME010000038.1       505     60      1176M48I312M322N145M1273N97M5306N71M    *       0       0
#GTGATCTGCATGTGTGACACTGATTCTTTGGAAATAAAGAGTGGAAGCTGCAGGTGACACGTGAAGGGTTATTTATGGTTATGATGACCCTGTCCTGCAACGAGGGACTGGCAGCCACTACTGAGGAGGAGGGTCCCATCTCTCTCCTGTCGGCTTTCACCGAGGTCACAGCCAGACGTGGGGCAAAGGTGTTCCCTGTCCTACCCAGCCATTCCTGGGCCTGCCGCCTAGGGGCTCACAGGGCCCAGGAGTCCCCAGCTCACAGGCCAGGGCATCAGGCCAGGCGCGCTCGGTGCACACCGCACCTGGGAGGACCTGGGTACACTCAGGAGACCAAGAGCACTGGCGGGTCAGGATGGTTGGCGTTCAGCTCCTACGGGGTGGGGAGAAGTCTGTAGCCGAGAGCCCAGCCCCCTCCTGCCAGGTCTTCCCAGGTTCGAGAGAGGCTGGAACTCAATTTCTGCAGAAAATCTCCCAGTTTTTCCTGTTTGGTTAGTTTTTTTACAAAGACAGGATCTTGCTGTGTTGCCCAGGCTGGTCTTGAACTCCTGGCCTCAAGCAATCCTCCCACCTCGGCCTCCGAAAGTGCTGGGATTACTGGCATGAACCACTGCGCCCGGCTGGAGCTCCCGGTTTTTAAGCACTGCACGATACTAGAAGAGCTGACCTTTTTTCTGGCCTCACAGCTTATGCTGAAGCTGAGTGTGAGGAACAGAGAGACCTTTCTGTGACGACCGCTGGGGCAGAGTGGTCTATGCGCCGAGATCCTGGCATCAGCAAGGGAGGCGGGTCCTCGGGGAGGGGCAGCTTCCACAGTGTGGCTGCAGCGTGCACAGCCAGGTAGGCCCTGGATGTTCACCCCTCACTGCCCTTGGGGAAGCACCTGACCGCTGGGGATGTCCACCAGGGAGAGGACGCTGTGTCGGGGACAACATGCAGCATCAGCACCCACAAGGGCCCGGCCTGGCCCCGCCTCTCCACTCGCCCGAGGTCTTGCTGTGGCCCAAAGCAGGAGTCCAGGGCTGGCGAGACCTCCGGCTGCAGAAAGGCAGGCCCAGGCCTCACGATGGAGAAAGTCTGGATGTCCTGGTCTGGCCTGCTGTTTCATGCAGTGTGCAAGCAGCAGCATGAGCAGGCAATAGGCCAACTCGGCTGGAGGCACCAACTGAAGCCAAGGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCACGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACACCGGGAATGTGCCATGGGTCCTGTGGGCGGTGCACGGGCTCAGGCACAGTGGGGCTGTCAGCTCTGGCTTGCAGGGTCACCGGCCTCACTGTCGCTGCTCTGAGACAGCTCTGTGGGGCTCTCAAGGCAGTGCAGCTCCAGGAAGCTGGTGATGAGCTTGGCGATGGCTTCCACATCACTGCGGTGGCCCATGACTGCACTCTGCGTGAGGGGGTGGGAGAAGCCCGTCGCAAGCCGGAGCACCACCATGTAGCCTTTCCCGAAGTACCGGACCTTCTCCTCCTCCACGCTCACATCACGGACATCATGGAGCAGGACCACCACCTGGTCGTGGCCAGCTCTGAAAAGAGTCAGCAGCTTCTTGTAGAGGCTGAACGTCTTCAAAACAACCTTCCCTGTGCTCTTGTCGAAGATGGCTTCCGAGAGCACAGCGGTCCCCGCGCCGCAGCAGAGAGACGCGAACCGGCGGGGGCGGGGCCGCGCGCACTTCC
#*       NM:i:52 ms:i:1727       AS:i:1685       nn:i:0  ts:A:+  tp:A:P  cm:i:317        s1:i:1716       s2:i:83 de:f:0.0028     cs:Z::161*cg:593*ct:71*tc:348+gccacgggtcctgtgggcggtgcacgggctcaggcacaccgggaatgt:9*tc:302~ct322ac:145~ct1273ac:97~ct5306ac:71     rl:i:0
# 13  0   84793   ENST00000380152.7   1000    +   0   84793   0,128,255   27  194,106,249,109,50,41,115,50,112,1116,4932,96,70,428,182,188,171,355,156,145,122,199,164,139,245,147,2105,  0,948,3603,9602,10627,10768,11025,13969,15445,16798,20791,29084,31353,39387,40954,42268,47049,47705,54928,55482,61196,63843,64276,64533,79215,81424,82688,
  my $percent_id_cutoff = $self->perc_id_cutoff();
  my $coverage_cutoff = $self->coverage_cutoff();
  unless(defined($percent_id_cutoff)) {
    $percent_id_cutoff = 90;
  }

  unless(defined($coverage_cutoff)) {
    $coverage_cutoff = 90;
  }

  my $canonical_intron_cutoff = 0.8;

  say "Parsing minimap2 sam output";
  my $genes = [];
  my @leftover_genes;

  open(IN, $output_file) or throw("Could not open $output_file");
  while(my $line = <IN>) {
    next if (index($line, '@') == 0);
    say "Output:\n".$line;
    my @results = split("\t",$line, 12);
    my $query_name = $results[0];
    my $flags = $results[1];
    next if ($flags&0x4);
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
    my $percent_identity = 100;
    my $coverage = 100;
    if ($attributes) {
      my ($cs_line) = $attributes =~ /cs:Z:(\S+)/;
      if ($cs_line) {
        ($percent_identity, $coverage) = parse_minimap2_cs($cs_line, length($query_seq));
        if ($percent_identity < $percent_id_cutoff) {
          $self->warning("Percent id for the hit fails the cutoff.\nHit name: $hit_name\nPercent id: $percent_identity\nCut-off: $percent_id_cutoff");
          next;
        }

        if ($coverage < $coverage_cutoff) {
          $self->warning("Coverage for the hit fails the cutoff.\nHit name: $hit_name\nCoverage: $coverage\nCut-off: $coverage_cutoff");
          next;
        }
      }
    }
    else {
      print __LINE__, "NO attributes found\n";
    }
    my $cigar_objects = parse_cigar_line($hit_name, $hit_start, $strand, $query_name, length($query_seq), $cigar_line);
    my $slice;
    if ($self->query) {
      $slice = $self->query;
      print __LINE__, ' ', $slice, "\n";
    }
    else {
      $slice = $self->slice_cache($hit_name);
      print __LINE__, ' ', $slice, "\n";
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
    # We aren't going to store a supporting feature, but we can store the coverage and percent id on the gene
    $coverage = int($coverage);
    $gene->version($coverage);
    $gene->description($percent_identity);

    my $transcript = ${$gene->get_all_Transcripts}[0];

    unless($self->skip_introns_check()) {
      my $introns = $transcript->get_all_Introns;
      my $intron_count = scalar(@$introns);
      if($intron_count) {
        my $canonical_count = 0;
        foreach my $intron (@$introns) {
          my ($is_canonical,$donor,$acceptor) = is_canonical_splice($intron,$self->database_adaptor->get_SliceAdaptor,$slice);
          if($is_canonical) {
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
   print join(' ', __LINE__, $gene->display_id, $gene->seq_region_start, $gene->seq_region_end, $gene->strand), "\n";
   foreach my $transcript (@{$gene->get_all_Transcripts}) {
     print join(' ', __LINE__, $transcript->display_id, $transcript->seq_region_start, $transcript->seq_region_end, $transcript->strand), "\n";
     foreach my $exon (@{$transcript->get_all_Exons}) {
       print join(' ', __LINE__, $exon->seq_region_start, $exon->seq_region_end, $exon->strand, $exon->length, $exon->start, $exon->end), "\n";
       foreach my $sf (@{$exon->get_all_supporting_features}) {
         print join(' ', __LINE__, $sf->start, $sf->end, $sf->strand, $sf->hstart, $sf->hend, $sf->hseqname), "\n";
       }
     }
   }
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
  my ($self,$genes,$leftover_input_file,$original_input_file) = @_;

  my $index = Bio::DB::HTS::Faidx->new($original_input_file);

  my $read_ids = [];
  foreach my $gene (@$genes) {
    push(@$read_ids,$gene->stable_id);
  }

  unless(open(OUTLEFTOVER,">".$leftover_input_file)) {
    $self->throw("Failed to open output file to use as input for leftover run");
  }

  foreach my $read_id (@$read_ids) {
    my $length = $index->length($read_id);
    my $location  = $read_id.":1-".$length;
    my $seq = $index->get_sequence_no_length($location);
    $seq = reverse($seq);
    $seq =~ tr/atgcATGC/tacgTACG/;
    say OUTLEFTOVER ">".$read_id;
    say OUTLEFTOVER $seq;
  }
  close OUTLEFTOVER;

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
    $self->throw("Created an exon where the start > than the end, this shouldn't be possible: ".$slice->name." ".$exon->start."..".$exon->end." ".$strand);
  }
#  say "Created exon: ".$slice->name." (".$exon_start."..".$exon_end.":".$strand.")";
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


sub paftools_path {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_paftools_path} = $val;
  }

  return $self->{_paftools_path};
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


sub add_offset {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_add_offset} = $val;
  }

  if($self->{_add_offset}) {
    return $self->{_add_offset};
  } else {
    return(0);
  }
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

sub slice_cache {
  my ($self, $val) = @_;
}

1;
