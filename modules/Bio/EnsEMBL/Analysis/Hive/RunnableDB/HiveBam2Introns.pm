=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns

=head1 SYNOPSIS

my $runnable =
  Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns->new(
);

$runnable->run;
my @results = $runnable->output;

=head1 DESCRIPTION

Takes an input id of a rough transcript and fetches reads associated with that model
from a BAM file, then runs a spliced exonerate alignment against genomic DNA or the
rough transcript. Writes output as SAM files.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns;

use warnings ;
use strict;

use File::Spec;
use File::Path qw(make_path);

use Bio::Seq;
use Bio::DB::HTS;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable::Bam2Introns;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(convert_to_ucsc_name);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: It will fetch the transcript and all the reads ovelapping the region. Then it filters
              the reads so we realign only the reads with mismatches on either the genomic span of
              the transcript or on the exonic regions of the transcript if the transcript is longer
              than 'max_transcript'.
              If there is more reads in the region than 'batch_size', it will create a new set of input ids
              in the format: <stable_id>:<start_exon>:<batch_size>:<offset>
              If your genome use UCSC style names (chr1,...), set 'wide_use_ucsc_naming' to 1
 Returntype : None
 Exceptions : Throws if the BAM file 'bam_file' does not exist

=cut

sub fetch_input {
  my ($self) = @_;
  # do all the normal stuff
  # then get all the transcripts and exons
  $self->throw('Your bam file "'.$self->param('bam_file').'" does not exist!') unless (-e $self->param('bam_file'));
  my $dna_db = $self->get_database_by_name('dna_db');
  $self->hrdb_set_con($dna_db, 'dna_db');
  my $gene_db = $self->get_database_by_name('input_db', $dna_db);
  my $gene_adaptor = $gene_db->get_GeneAdaptor;
  my $slice_adaptor = $dna_db->get_SliceAdaptor;
  my $counters;
  $counters->{'start'} = 0;
  $counters->{'offset'} = 0;
  $counters->{'start_exon'} = 0;
  my $stable_id = $self->input_id;
  # check for batch info in the input id
  if ( $self->input_id =~ /(\S+):(\d+):(\d+):(\d+)/ ) {
    $stable_id = $1;
    $counters->{'start_exon'} = $2;
    $counters->{'start'} = $3;
    $counters->{'offset'} = $4;
  }
  my $rough = $gene_adaptor->fetch_by_stable_id($stable_id);
  $self->throw("Gene $stable_id not found \n") unless $rough;
    print 'Found model '.join(' ', $rough->stable_id, $rough->seq_region_name, $rough->start, $rough->end, $rough->strand)."\n";
  $self->rough($rough->get_all_Transcripts->[0]);
  # number of missmatches needed before using a read from the bam file
  # is calculated as the Exonerate word length - number of matches you
  # might expect by chance ie 1 in 4
  my $missmatch = $self->param('word_length') - int($self->param('word_length') / 4);
  $self->param('min_missmatch', $missmatch);
  print "Ignoring reads with < $missmatch mismatches\n";
  my @reads;
  my $bam = Bio::DB::HTSfile->open($self->param('bam_file'));
  my $header = $bam->header_read();
  my $seq_region_name = $rough->seq_region_name;
  if ($self->param('wide_use_ucsc_naming')) {
      $seq_region_name = convert_to_ucsc_name($seq_region_name, $rough->slice);
  }
  my ($tid, $start, $end) = $header->parse_region($seq_region_name);
  my $bam_index = Bio::DB::HTSfile->index_load($bam);
  $self->throw('Bam file ' . $self->param('bam_file') . "  not found \n") unless $bam_index;
  $counters->{'count'} = 0;
  $counters->{'batch'} = 0;
  $counters->{'batch_size'} = $self->param('batch_size');
  my $exon = 0;
  my @exons = sort { $a->start <=> $b->start } @{$rough->get_all_Exons};
  my @iids;
  $counters->{'stop_loop'} = 0;

  $self->param('seq_hash', {});
  for ( my $i = $counters->{'start_exon'} ; $i <= $#exons ; $i++ ){
    my $exon = $exons[$i];
    my $exon_start = $exon->start;
    $counters->{'read_offset'} = 0;
    $counters->{'last_start'} = 0;
    my %split_points;
    $exon_start = $counters->{'start'} if ($counters->{'start'} && $counters->{'start_exon'} == $i);
    my @callback_data = ($self, $i, $exon_start, $rough->stable_id, \@reads, \@iids, $counters, $self->param('seq_hash'));
    $bam_index->fetch($bam, $tid, $exon_start-1, $exon->end-1, \&_process_reads, \@callback_data);
  }
  if (  scalar(@reads) == 0 ) {
    $self->input_job->autoflow(0);
    $self->complete_early('No read with more than '.$missmatch.' missmatches');
  }
  else {
  if ( scalar(@iids) > 0  && !$counters->{'start'} ) {
    # We want the original input_id job to finish before submitting the new input ids otherwise we may have problems
    print 'Making ',  scalar(@iids), ' New input ids from ', $counters->{'count'}, " reads\n";
    $self->param('iids', \@iids);
  }
  my $fullseq = $self->fullslice;
  my $queryseq;
  # get the fullseq if required
  my $slice = $slice_adaptor->fetch_by_region('toplevel',$rough->seq_region_name,$rough->start,$rough->end,$rough->strand);
  $self->param('query', $slice);
  if ( $fullseq && $rough->length < $self->param('max_transcript') ) {
    if ( $self->param('mask') ) {
      $queryseq = $slice->get_repeatmasked_seq(undef,1)->seq;
    }
    else {
      $queryseq = $slice->seq;
    }
  }
  elsif ( $fullseq )  {
    # we dont want to use fullseq in the case of a very long transcript
    $self->fullslice(0);
    foreach my $exon ( @{$rough->get_all_Transcripts->[0]->get_all_Exons}) {
      my $slice = $exon->feature_Slice;
      if ( $self->param('mask') ) {
        $queryseq .= $slice->get_repeatmasked_seq(undef,1)->seq;
      }
      else {
        $queryseq .= $slice->seq;
      }
    }
  }
  my $masked_count = $queryseq =~ tr/atcgn/atcgn/;
  if (($masked_count/length($queryseq)) > .99) {
    $self->input_job->autoflow(0);
    $self->complete_early(sprintf("Highly repetitive sequence: %.2f%% masked", (($masked_count*100)/length($queryseq))));
  }
  my $seqio = Bio::Seq->new( -seq => $queryseq,
                 -display_id => $rough->stable_id
               );
  # set uo the runnable
  my $program = $self->param('program_file');
  $program = 'exonerate' unless $program;
  $self->param('saturatethreshold', scalar(@reads)) unless ($self->param_is_defined('saturatethreshold'));

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Bam2Introns->new(
     -analysis     => $self->create_analysis,
     -program      => $program,
     -basic_options => $self->get_aligner_options,
     -target_seqs => [$seqio],
     -query_seqs => \@reads,
     -percent_id   => $self->param('percent_id'),
     -coverage     => $self->param('coverage'),
     -missmatch     => $self->param('missmatch'),
    );
  $self->runnable($runnable);
  }
}


=head2 get_aligner_options

 Arg [1]    : None
 Description: Return the exonerate options needed for the alignemnt.
 Returntype : String, exonerate options
 Exceptions : None

=cut

sub get_aligner_options {
    my $self = shift;

    my $options =  '--showsugar false --showvulgar false --showalignment false --ryo "RESULT: %S %pi %ql %tl %g %V\n" '.
                '--model est2genome --forwardcoordinates false '.
                '--softmasktarget '.($self->param('mask') ? 'true' : 'false').' --exhaustive false --percent 80 '.
                '--dnahspthreshold 70 --minintron 20 --dnawordlen '.$self->param('word_length').' -i -12 --bestn 1';
    $options .= ' --saturatethreshold '.$self->param('saturate_threshold') if ($self->param_is_defined('saturate_threshold'));
    $options .= ' --maxintron '.$self->param('maxintron') if ($self->param_is_defined('maxintron'));
    return $options;
}


=head2 run

 Arg [1]    : None
 Description: Run exonerate on the genomic sequence using only the reads with missmatches.
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;
  $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
  $self->throw("Can't run - no runnable objects") unless ( $self->runnable );
  my ($runnable) = @{$self->runnable};
  $runnable->run;
  my $features = $runnable->output;
# It is not in the Runnable because it needs to do DB calls to get the sequence
# Might be worth to get the sequence in fetch_input then do some Perl substr
  my $genomic_features = $self->process_features($features);
  $self->output($genomic_features);
}


=head2 write_output

  Arg [1]   : None
  Function  : Write the alignments in SAM format and dataflow the new input ids via 'iid' on branch 2
  Returntype: 1
  Exceptions: Throws if the feature cannot be stored

=cut

sub write_output {
  my $self = shift;

  my $output = $self->output;
  print "Got " .  scalar(@$output) ." genomic features \n";
  if (scalar(@$output)) {
      # write to file
      my $iid = $self->input_id;
      # remove any batching info at the end
      $iid =~ s/:\d+:\d+:\d+$//;

      my $path;
      my $filename;
      # figure out a directory structure based on the stable ids
      if ( $iid =~ /^\w+\d+(\d)(\d)(\d)(\d)(\d\d)$/ ) {
        # make the directory structure
        $path = File::Spec->catdir($self->param('wide_output_sam_dir'), $1, $2, $3, $4);
        make_path($path);
        $filename = File::Spec->catfile($path, $self->input_id.'.sam');
      }
      else {
        $self->throw("Input id $iid structure not recognised should be something like BAMG00000002548\n");
      }
      open ( SAM ,">$filename" ) or $self->throw("Cannot open file for writing $filename\n");
      my $seq_hash = $self->param('seq_hash');
      my $line_count = 0;
      foreach my $feature ( @$output ) {
          my $line = $self->convert_to_sam($feature, $seq_hash);
          if ($line) {
              print SAM $line;
              $line_count++;
          }
      }
      # write an end of file marker so we can test that the job didnt die when writing
      print SAM '@EOF';
      close SAM;
      if ($line_count) {
          $self->dataflow_output_id([{filename => $filename}], 1);
      }
      else {
          $self->input_job->autoflow(0);
          unlink $filename || $self->warning("Could not remove empty SAM file $filename");
      }
  }
  else {
      $self->input_job->autoflow(0);
  }
  if ($self->param_is_defined('iids')) {
      $self->dataflow_output_id($self->param('iids'), 2);
  }
}


=head2 convert_to_sam

 Arg [1]    : Bio::EnsEMBL::DnaDnaAlignFeature, represents the read alignment
 Arg [2]    : Hashref of Bio::Seq, the name of the read is the key
 Example    : $self->convert_to_sam($feature, $seq_hash);
 Description: Convert a DnaDnaAlignFeature to a SAM record
 Returntype : String, the SAM record for the alignment or empty if the length of the cigar line is different
              from the sequence length
 Exceptions : Throws if it cannot find the original sequence

=cut

sub convert_to_sam {
  my ( $self, $feature, $seq_hash) = @_;
  my $line;
  # I have a feeling the hit strands are all backwards in the db
  # I will have a go at reversing them
  #$feature->hstrand($feature->hstrand * -1);
  # strip off the :NC of the end of the seq names
  my $tmp = $feature->hseqname;
  $tmp   =~ s/:NC$//;
  $feature->hseqname($tmp);
  my $feature_seq = $seq_hash->{$tmp};
  $self->throw("cannot find sequence for " . $feature->hseqname ."\n")
    unless $feature_seq;
  my $seq = "*";
  my $flag = 1;
  my $paired = 1;
  # reverses cigar line if feature is reversed
  $feature = $self->check_cigar($feature);
  my $cigar = $feature->cigar_string;
  # N is used to represent splicing in sam format
#  $cigar =~ s/I/N/g;
  # But SAM does not like D in the CIGAR line as it's not counted as part of the sequence. We need to change it to I
  $cigar =~ tr/ID/NI/;
  # is the feature paired in the alignment
  $flag = 3;
  # strand wrt the reference
  # do we know the strand?
  my $strand = $feature->strand ;
  if ( $strand == -1 ) {
    $flag +=16;
  }
  # strand of the mate = -1 so we add nothing
  my $length = length($feature_seq->seq);
  # add soft clip info to the cigar
  # cigar reversed if the hit strand is -1;
  if ( ( $feature->hstrand == -1 &&  $feature->strand == -1) or ( $feature->hstrand == 1 &&  $feature->strand == 1 ) ) {
    $cigar = $cigar . ($feature->hstart - 1) ."S"	if  $feature->hstart > 1;
    $cigar = ( $length - $feature->hend ) ."S$cigar"	if  $feature->hend < $length;
  } else {
    $cigar = join('', reverse ( $cigar =~ /(\d+[A-Z])/g ));
    $cigar = ( $feature->hstart -1 ) ."S$cigar"	if  $feature->hstart > 1;
    $cigar = $cigar . ( $length - $feature->hend ) ."S"	if  $feature->hend < $length;
  }
  my $check = $self->check_cigar_length($cigar);
  if ( $check != $length ) {
    print STDERR "Losing " .  $feature->hseqname . "  I get cigar length $check rather than $length, $cigar\n";
    return ;
  }
  if ( ($feature->hstrand*$feature->strand) == 1 ) {
    # need the reverse compliment
    $seq = $feature_seq->revcom->seq;
  }else {
    $seq = $feature_seq->seq;
  }
  # if the feature is aligned to the -ve strand we want the complement of he sequence as well
  if ( $feature->strand == -1 ) {
    $seq =~ tr/ACGTacgt/TGCAtgca/;
  }

  # strand of the mate = -1 so we add nothing
  $line .= $feature->hseqname ."\t";
  $line .= "$flag\t";
  my $seq_region_name = $feature->seq_region_name;
  $seq_region_name = convert_to_ucsc_name($seq_region_name, $self->param('query')) if ($self->param('wide_use_ucsc_naming'));
  $line .= $seq_region_name ."\t";
  $line .=  $feature->seq_region_start ."\t";
  $line .= "0\t";
  $line .=  "$cigar\t";
  if ( $paired ) {
    # what is the position of the other pair?
    $line .= $feature->seq_region_name ."\t" . $feature->hstart ."\t" . $feature->hend . "\t$seq\t*\tRG:Z:". $feature_seq->desc ."\n";
  } else {
    $line .= "*\t0\t0\t$seq\t*\tRG:Z:". $feature_seq->desc ."\n";
  }
  return $line;
}


=head2 check_cigar

 Arg [1]    : Bio::EnsEMBL::DnaDnaAlignFeature, represents the read alignment
 Example    : $self->check_cigar($feat);
 Description: Check that the length of the cigar line match the length of the sequence.
              Reverse the cigar line if it's on the reverse strand. Store the cigar line.
 Returntype : Bio::EnsEMBL::DnaDnaAlignFeature
 Exceptions : None

=cut

sub check_cigar {
  my ($self,$feat) = @_;
  my @ugfs;
  my $string = $feat->cigar_string;
  my @pieces = ( $string =~ /(\d*[MDI])/g );
  # if  strand is -1 strand cigar line is reversed
  @pieces = reverse @pieces if     $feat->strand == -1 ;
  my $cigar;
  foreach my $piece (@pieces) {
    my ($length,$match) = ( $piece =~ /^(\d*)(M|I|D)/ );
    if( $length eq "" ) { $length = 1 }
    $cigar.= $length;
    $cigar.= $match;
  }
  $feat->cigar_string($cigar);
  return $feat;
}


=head2 check_cigar_length

 Arg [1]    : String, representing the cigar line
 Example    : if ($self->check_cigar_length('30M488N70M') == $feat->length) { return 1 };
 Description: Check that the length of the cigar line matches the length of the sequence
              according to the SAM secifications
 Returntype : Integer, length of the cigar line
 Exceptions : None

=cut

sub check_cigar_length{
  my ($self,$cigar) = @_;
  my @pieces = ( $cigar =~ /(\d*[SMDN])/g );
  my $total_length;
  foreach my $piece (@pieces) {
    my ($length,$match) = ( $piece =~ /^(\d*)(S|M|D|N)/ );
    if( $length eq "" ) { $length = 1 };
    $total_length += $length unless $match eq "N";
  }
  return $total_length;
}


=head2 process_features

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Uses cdna2genome to convert the ungapped align
              features into genomic coords
  Returntype: 1
  Exceptions: Throws if the transcript the read has aligned to is
              not stored in the transcript hash

=cut

sub process_features {
    my ( $self, $flist ) = @_;

# first do all the standard processing, adding a slice and analysis etc.
# unless we are projecting in which case we dont really nead a slice

    my @dafs;
    my $trans = $self->rough;
    my $transcript_strand = $trans->strand;
# I hate to do that but let's uses Perl's features...
    my $splice_in_splice_limit = scalar(@$flist);
    my $splice_in_splice_found = 0;
    FEATURE: foreach my $f (@$flist) {
        my @features;
# get as ungapped features
#        print STDOUT 'DEBUG FEATURE ', ref($f), "\n";
        foreach my $ugf ( $f->ungapped_features ) {
#            print STDOUT join(' ', 'DEBUG UGF', $ugf->start, $ugf->end, $ugf->strand, $ugf->hstart, $ugf->hend, $ugf->hstrand), "\n";
            if ($self->fullslice ) {
                $ugf->hstrand($f->hstrand * $transcript_strand);
                $ugf->strand($f->strand * $transcript_strand);
                $ugf->slice($self->param('query'));
                push @features,$ugf;
            }
            else {
# Project onto the genome if it was not run with fullseq
# We only want the object. If there is more than one object per ungapped feature
# we have a problem. That would mean that the splice sites we found are wrong.
# In this case we don't want to write down the alignment
            my @ugfs = grep { $_->isa('Bio::EnsEMBL::Mapper::Coordinate') } $trans->cdna2genomic($ugf->start, $ugf->end);
            if (scalar(@ugfs) == 1) {
              my $obj = $ugfs[0];
# make into feature pairs?
              my $fp = Bio::EnsEMBL::FeaturePair->new(
                -start    => $ugfs[0]->start,
                -end      => $ugfs[0]->end,
                -strand   => $ugf->strand*$transcript_strand,
                -slice    => $trans->slice,
                -hstart   => $ugf->hstart,
                -hend     => $ugf->hend,
                -hstrand  => $ugf->hstrand,
                -percent_id => $f->percent_id,
                -score    => $f->score,
                -hseqname => $f->hseqname,
                -hcoverage => $f->hcoverage,
                -p_value   => $f->p_value,
              );
              push @features, $fp;
            }
            else {
              print STDOUT 'It seems that ', $f->hseqname, ' is mapping on a splice site.',
              ' This is wrong as we add an extra padding. ',
              join(' ', 'More info:', $ugf->start, $ugf->end, $f->start, $f->end, $f->strand, $f->hstart, $f->hend, $f->hstrand), "\n";
              ++$splice_in_splice_found;
              next FEATURE;
            }
          }
        }
        push(@dafs, @{$self->build_dna_align_features($f, \@features)});
    }
    $self->warning('Too many reads have splice site in their exons, something in wrong as these exons should not be spliced '.$splice_in_splice_found.'/'.$splice_in_splice_limit)
      if ($splice_in_splice_limit/100 < $splice_in_splice_found);
    return \@dafs;
}


=head2 build_dna_align_features

 Arg [1]    : Bio::EnsEMBL::DnaDnaAlignFeature representing the ungapped alignment
 Arg [2]    : Bio::EnsEMBL::FeaturePair representing the alignment
 Description: Create a Bio::EnsEMBL::DnaDnaAlignFeature object based on the alignment from exonerate
 Returntype : Arrayref of Bio::EnsEMBL::DnaDnaAlignFeature
 Exceptions : None

=cut

sub build_dna_align_features {
    my ($self, $f, $features) = @_;

    my @dafs;
    my @features = sort { $a->start <=> $b->start } @$features;
    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);
# corect for hstart end bug
#    $feat->hstart($f->hstart);
#    $feat->hend($f->hend);
    $feat->analysis($self->analysis);
# transfer the original sequence of the read
    $feat->{"_feature_seq"} = $f->{"_feature_seq"};
# dont store the same feature twice because it aligns to a different transcript in the same gene.
#    my $unique_id = join(':', $feat->seq_region_name, $feat->start, $feat->end, $feat->strand, $feat->hseqname);
    my $unique_id = join(':', $feat->seq_region_name, $feat->seq_region_start, $feat->seq_region_end, $feat->strand, $feat->hseqname);
    unless ($self->stored_features($unique_id)){
        push @dafs,$feat;
# keep tabs on it so you don't store it again.
        $self->stored_features($unique_id,1);
    }
    return \@dafs;
}


###########################################################
# containers

=head2 rough

 Arg [1]    : (optional) Bio::EnsEMBL::Gene, the proto transcript
 Example    : my $gene = $self->rough;
 Description: Getter/Setter for the proto-transcript
 Returntype : Bio::EnsEMBL::Gene
 Exceptions : None

=cut

sub rough {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('_rough', $value);
  }

  if ($self->param_is_defined('_rough')) {
    return $self->param('_rough');
  } else {
    return;
  }
}


=head2 fullslice

 Arg [1]    : (optional) boolean 0
 Example    : if ($self->fullslice) {};
 Description: Getter/Setter for knowing if you want to align on the genomic span
              of the transcript
 Returntype : Boolean, 1 if you want the genomic span of the transcript
 Exceptions : None

=cut

sub fullslice {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('fullseq', $value);
  }

  if ($self->param_is_defined('fullseq')) {
    return $self->param('fullseq');
  } else {
    return;
  }
}


=head2 stored_features

 Arg [1]    : String, name of the read
 Arg [2]    : (optional) Bio::Seq, sequence of the read
 Example    : if ($self->stored_features($readname)) {};
 Description: Getter/Setter for the sequence of a read
 Returntype : Bio::Seq
 Exceptions : None

=cut

sub stored_features {
  my ($self,$key,$value) = @_;

  return unless defined ($key);

  if (!$self->param_is_defined('_stored_features')) {
    $self->param('_stored_features', {});
  }
  if (defined $value) {
    $self->param('_stored_features')->{$key} = $value;
  }

  if ($self->param_is_defined('_stored_features') and exists($self->param('_stored_features')->{$key})) {
    return $self->param('_stored_features')->{$key};
  } else {
    return;
  }
}

####################################################
# batching works by running initialy on the first
# 100,000 reads, once it gets past this is just counts
# how many reads there are and submits new input_ids
# of the type STABLEID00000001:2 for the second bacth
# of 100,000 reads etc.

=head2 _process_reads

 Arg [1]    : Arrayref of objects needed for processing the reads
 Description: This function is called for each read in the region. It checks the number of missmatches,
              the nnumber of reads already processed, update the name of the read based on the FIRST_MATE
              flag and store the sequence in a Bio::Seq object which will be put in an arrayref
 Returntype : None
 Exceptions : None

=cut

sub _process_reads {
    my ($read, $callbackdata) = @_;

    my ($self, $i, $exon_start, $stable_id, $reads, $iids, $counters, $seq_hash) = @$callbackdata;
    return if ($counters->{'stop_loop'});
    my $missmatch  = $read->get_tag_values('NM');
    return unless ($missmatch and $missmatch >= $self->param('min_missmatch'));
    # get rid of any reads that might start before our start point
    return if (!$read->start or ($i == $counters->{'start_exon'} and $read->start < $exon_start));
    $counters->{'read_count'} ++;
    # get rid of any reads overlapping out start point
    return if ($counters->{'read_offset'} and $counters->{'read_count'} <= $counters->{'read_offset'});
    # calculate the offset
    if ( $read->start == $counters->{'last_start'} ) {
        $counters->{'read_offset'} ++;
    } else {
        $counters->{'read_offset'} = 0;
        $counters->{'last_start'} = $read->start;
    }

    my $rg = $read->get_tag_values('RG') || '*';

    $counters->{'batch'} ++;
    $counters->{'count'} ++;

    if ($counters->{'count'} <= $counters->{'batch_size'}) {
        my $name = $read->query->name;
        # is it the 1st or 2nd read?
        if ( $read->get_tag_values('FIRST_MATE')) {
            $name .= '/1';
        }
        elsif ( $read->get_tag_values('SECOND_MATE')) {
            $name .= '/2';
        }
        # write the seq files - store the read group information in case it is needed later
        my $bioseq = Bio::Seq->new(
                -seq        => $read->query->dna,
                -display_id => $name,
                -desc       => $rg
                );
        push(@$reads, $bioseq);
        # want to store the read sequence for making it into a SAM file later
        $seq_hash->{$name} = $bioseq;
    }
    else {
        # if running a batch finish here
        if ($counters->{'start'}) {
            $counters->{'stop_loop'} = 1;
            return;
        }
        if ($counters->{'batch_size'}+1 == $counters->{'batch'} ) {
            $counters->{'batch'} = 1;
            push(@$iids, { iid => $stable_id .":$i:" . $read->start .':'.$counters->{'read_offset'}});
        }
        # otherwise figure out the ids for the rest of the batches
    }
}

1;
