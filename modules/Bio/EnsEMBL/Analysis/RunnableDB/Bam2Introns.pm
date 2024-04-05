=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


package Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns;

use warnings ;
use vars qw(@ISA);
use strict;
use Bio::SeqFeature::Lite;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns;
use Bio::EnsEMBL::Analysis::Runnable::Bam2Introns;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::DB::HTS;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);
$| = 1;
sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->db->disconnect_when_inactive(1);
  $self->read_and_check_config($BAM2INTRONS_CONFIG_BY_LOGIC);
  return $self;
}

# fetch input
# get the transcript in question
# get all the reads that overlap
# do we need to hold them in memory or can we stream them into the file?



sub fetch_input {
  my ($self) = @_;
  # do all the normal stuff
  # then get all the transcripts and exons
  my $gene_db = $self->get_dbadaptor($self->TRANSDB);
  my $gene_adaptor = $gene_db->get_GeneAdaptor;
  my $slice_adaptor =  $gene_db->get_SliceAdaptor;
  my $start = 0;
  my $start_exon;
  my $offset;
  my $stable_id = $self->input_id;
  # check for batch info in the input id
  if ( $self->input_id =~ /(\S+):(\d+):(\d+):(\d+)/ ) {
    $stable_id = $1;
    $start = $3;
    $start_exon = $2;
    $offset = $4;
  }
  my $rough = $gene_adaptor->fetch_by_stable_id($stable_id);
  $self->throw("Gene $stable_id not found \n") unless $rough;
    print "Found model " . $rough->stable_id . " ". 
    $rough->seq_region_name ." " . 
       $rough->start ." " . 
	 $rough->end ." " .
	   $rough->strand . "\n";
  $self->rough($rough->get_all_Transcripts->[0]);
  my $fullseq =  $self->FULLSEQ;
  # get the fullseq if required
  if ( $fullseq && $rough->length < $self->MAX_TRANSCRIPT ) {
    my $slice = $slice_adaptor->fetch_by_region('toplevel',$rough->seq_region_name,$rough->start,$rough->end,$rough->strand);
    if ( $self->MASK ) {
      $fullseq = $slice->get_repeatmasked_seq(undef,1)->seq;
    } else {
      $fullseq = $slice->seq;
    }
    $self->fullslice($slice);
  } elsif ( $fullseq )  {
    # we dont want to use fullseq in the case of a very long transcript
    $fullseq = undef;
    $self->FULLSEQ(0);
  }
  # set uo the runnable
  my $mask = "TRUE";
  $mask = "FALSE" unless $self->MASK;
  my $program = $self->analysis->program_file;
  $program = "/software/ensembl/genebuild/bin/exonerate64-0.9.0" unless $program;

  my $options =  "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" " .
                 "--model est2genome --forwardcoordinates false ".
                 "--softmasktarget $mask --exhaustive false --percent 80 ".
                 "--dnahspthreshold 70 --minintron 20 --dnawordlen " . $self->WORD_LENGTH . " -i -12 --bestn 1";
  $options .= " --saturatethreshold " .$self->SATURATE_THRESHOLD if $self->SATURATE_THRESHOLD ;
  # number of missmatches needed before using a read from the bam file
  # is calculated as the Exonerate word length - number of matches you 
  # might expect by chance ie 1 in 4 
  my $missmatch = $self->WORD_LENGTH - int($self->WORD_LENGTH / 4);
  print "Ignoring reads with < $missmatch mismatches\n";
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Bam2Introns->new
    (
     -analysis     => $self->analysis,
     -program      => $program,
     -basic_options => $options,
     -model        => $rough,
     -missmatch    => $missmatch,
     -mask         => $self->MASK,
     -out_sam_dir  => $self->OUT_SAM_DIR,
     -bam_file     => $self->BAM_FILE,
     -percent_id   => $self->PERCENT_ID,
     -coverage     => $self->COVERAGE,
     -fullseq      => $fullseq,
     -max_tran     => $self->MAX_TRANSCRIPT,
     -start        => $start,
     -startexon    => $start_exon,     
     -offset       => $offset,
     -batch_size   => $self->BATCH_SIZE,
    );
  $self->runnable($runnable);
}

sub run {
  my ($self) = @_;
  throw("Can't run - no runnable objects") unless ( $self->runnable );
  my ($runnable) = @{$self->runnable};
  $runnable->run;
  my $features = $runnable->output;
  if ($self->MISSMATCH  ) {
    $features = $self->filter_solexa($features);
  }
  my $genomic_features = $self->process_features($features);
  $self->output($genomic_features);
  my $iids =  $runnable->iids;
  if ( $iids ) {
    # we have to make some new iids as the number of reasds is > batch size
    my $analysis = $self->analysis;
    my  $hash = $self->database_hash;
    # need to delete this hash ref in order 
    # to call the pipeline version
    $hash->{REFERENCE_DB} = undef;
    my $pipelinea = $self->get_dbadaptor('REFERENCE_DB','pipeline');
    my $pipeline_analysis = $pipelinea->get_AnalysisAdaptor->fetch_by_logic_name
      ($self->analysis->logic_name);
    my $ra = $pipelinea->get_RuleAdaptor;
    my $rule = $ra->fetch_by_goal($pipeline_analysis);
    my $sic = $pipelinea->get_StateInfoContainer;
    # check the ids have not already been stored
    my $check ;
    foreach my $id ( @$iids ) {
      if ( my $analysis = $sic->fetch_analysis_by_input_id($id) ) {
	$check++ if  scalar(@$analysis) > 0;
      }
    }
    unless ( $check) {
      foreach my $condition ( @{$rule->list_conditions} ) {
        # want to store the input id
        my $inputIDFactory = new Bio::EnsEMBL::Pipeline::Utils::InputIDFactory
	  (
	   -db => $pipelinea,
	   -file => 1,
	   -logic_name => $condition,
	  );
        $inputIDFactory->input_ids(\@$iids);
        $inputIDFactory->store_input_ids;
      }
    } else {
      print STDERR "Already stored these input ids\n";
    }
  }
}


sub filter_solexa {
  my ($self,$features_ref) = @_;
  my @features = @$features_ref;
  my @filtered_features;
  # allow no more than MISSMATCH missmatches and
  foreach my $feat ( @features ) {
    # restrict just to splice models where read is within an exon
    next unless $feat->{"_intron"};
    # check missmatches
    if ( $self->MISSMATCH ) {
      my $aligned_length = abs($feat->hend - $feat->hstart) +1;
      my $matches = $aligned_length *  $feat->percent_id / 100;
      my $missmatches = ( $aligned_length - $matches) / $aligned_length * 100;
      next if $missmatches > $self->MISSMATCH;      
    }
    push @filtered_features, $feat;
  }
  return \@filtered_features;
}

=head2 write_output

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Overrides the write_output method from the superclass
              to allow Writing in SAM format
  Returntype: 1
  Exceptions: Throws if the feature cannot be stored

=cut

sub write_output {
  my ( $self, $loop ) = @_;
  if ( $loop ) {
    # keep track of how many times we have tried to write this file
    $self->write_attempts($loop);
    $self->throw("Failed at writing output file twice - giving up\n")
      if $self->write_attempts > 1;
    # otherwise lets give it another try
  }
  my @output = @{$self->output};
  print "Got " .  scalar(@output) ." genomic features \n";
  return unless  scalar(@output) > 0 ;
  # write to file
  my $iid = $self->input_id;
  # remove any batching info at the end
  if ( $self->input_id =~ /(\S+):(\d+):(\d+):(\d+)/ ) {
    $iid = $1;
  } 

  my $path;
  # figure out a directory structure based on the stable ids
  if ( $iid =~ /^\w+\d+(\d)(\d)(\d)(\d)(\d\d)$/ ) {
    # make the directory structure
    $path = $self->OUT_SAM_DIR. '/' . "$1";
    print "PATH $path\n";
    mkdir($path) unless -d ($path);
    $path = $path . '/' . "$2";
    mkdir($path) unless -d ($path);
    $path = $path . '/' . "$3"; 
    mkdir($path) unless -d ($path);
    $path = $path . '/' . "$4";
    mkdir($path) unless -d ($path);
  } else {
    $self->throw("Input id $iid structure not recognised should be something like BAMG00000002548\n");
  }
  open ( SAM ,">".$path."/" . $self->input_id . ".sam" ) or
    $self->throw("Cannot open file for writing " .$path."/" . $self->input_id . ". sam\n");
  eval { 
    foreach my $feature ( @output ) {
      my $line = $self->convert_to_sam($feature);
      print SAM $line if $line;
    }
    # write an end of file marker so we can test that the job didnt die when writing
    print SAM '@EOF';
    close SAM;
  }; if ( $@ ) {
    print " What happened? \n$@\n";
  }
  # verify file wrote correctly
  open ( CHECK , $path."/" . $self->input_id . ".sam" ) or
    $self->throw("Cannot open file for checking " .$path."/" . $self->input_id . ". sam\n");
  my $line;
  my $line_count;
  while (<CHECK>) {
    # 1st file copy the header all the others just copy the data
    chomp;
    $line = $_;
    $line_count++;
  }
  return if ( $line eq '@EOF' or $line_count == 0 );
  # file must have been corrupted try and write it again but just once otherwise just fail
  system("rm  " . $path."/" . $self->input_id . ".sam" );
  print "Failed to verify file - deleting output and writing again\n";
  $self->write_output(1);    
}



sub convert_to_sam {
  my ( $self,$feature ) = @_;
  my $line;
  # I have a feeling the hit strands are all backwards in the db
  # I will have a go at reversing them 
  #$feature->hstrand($feature->hstrand * -1);
  # strip off the :NC of the end of the seq names
  my $tmp = $feature->hseqname;
  $tmp   =~ s/:NC$//;
  $feature->hseqname($tmp);
  my $seq = "*";
  my $flag = 1;
  my $paired = 1;
  # reverses cigar line if feature is reversed
  $feature = $self->check_cigar($feature);
  my $cigar = $feature->cigar_string;
  # N isued to represent splicing in sam format
  $cigar =~ s/I/N/g;
  # is the feature paired in the alignment
  $flag =3 ;
  # strand wrt the reference
  # do we know the strand?
  my $strand = $feature->strand ;
  if ( $strand == -1 ) {
    $flag +=16;
  }
  # strand of the mate = -1 so we add nothing
  my $feature_seq = $feature->{"_feature_seq"};
  $self->throw("cannot find sequence for " . $feature->hseqname ."\n")
    unless $feature_seq;
  my $length = length($feature_seq->seq);
  # add soft clip info to the cigar
  # cigar reversed if the hit strand is -1;
  if ( ( $feature->hstrand == -1 &&  $feature->strand == -1) or ( $feature->hstrand == 1 &&  $feature->strand == 1 ) ) {
    $cigar = $cigar . ($feature->hstart - 1) ."S"	if  $feature->hstart > 1;
    $cigar = ( $length - $feature->hend ) ."S$cigar"	if  $feature->hend < $length;
  } else {
    $cigar = ( $feature->hstart -1 ) ."S$cigar"	if  $feature->hstart > 1;
    $cigar = $cigar . ( $length - $feature->hend ) ."S"	if  $feature->hend < $length;
  }
  my $check = $self->check_cigar_length($cigar);
  if ( $check != $length ) {
    print STDERR "Losing " .  $feature->hseqname . "  I get cigar length $check rather than $length, $cigar\n";
    return ;
  }
  if ( $feature->hstrand == 1 &&  $feature->strand == 1 ) {
    # need the reverse compliment
    $seq = $feature_seq->revcom->seq;
  } elsif( $feature->strand == -1 &&  $feature->hstrand == -1 )   {
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
 # $line .= $feature->hstrand . "\t" . $feature->strand . "\t";
  $line .= $feature->seq_region_name ."\t";
  $line .=  $feature->start ."\t";
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

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my @dafs;
  my $count = 0;
FEATURE:  foreach my $f (@$flist) {
    $count++;
    my $trans = $self->rough;
    my @mapper_objs;
    my @features;
    my $start = 1;
    my $end = $f->length;
    my $out_slice = $slice_adaptor->fetch_by_name($trans->slice->name);
    # get as ungapped features
    foreach my $ugf ( $f->ungapped_features ) {
      # Project onto the genome unless it was run with fullseq
      unless ($self->FULLSEQ ) {
	foreach my $obj ($trans->cdna2genomic($ugf->start, $ugf->end)){
	  if( $obj->isa("Bio::EnsEMBL::Mapper::Coordinate")){
	    # make into feature pairs?
	    my $qstrand = $f->strand  * $trans->strand;
	    my $hstrand = $f->hstrand * $trans->strand;
	    my $fp;
	    $fp = Bio::EnsEMBL::FeaturePair->new
	      (-start    => $obj->start,
	       -end      => $obj->end,
	       -strand   => $qstrand,
	       -slice    => $trans->slice,
	       -hstart   => 1,
	       -hend     => $obj->length,
	       -hstrand  => $hstrand,
	       -percent_id => $f->percent_id,
	       -score    => $f->score,
	       -hseqname => $f->hseqname,
	       -hcoverage => $f->hcoverage,
	       -p_value   => $f->p_value,
	      );
	    push @features, $fp->transfer($out_slice);
	  }
	}
      } else {
	$ugf->hstrand($f->hstrand * $trans->strand);
#	$ugf->strand($f->strand * $trans->strand);
	$ugf->slice($self->fullslice);
	push @features,$ugf->transfer($out_slice);
      }
    }
    @features = sort { $a->start <=> $b->start } @features;
    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);
    # corect for hstart end bug
    $feat->hstart($f->hstart);
    $feat->hend($f->hend);
    $feat->analysis($self->analysis);
    # transfer the original sequence of the read
    $feat->{"_feature_seq"} = $f->{"_feature_seq"};
    # dont store the same feature twice because it aligns to a different transcript in the same gene.
    my $unique_id = $feat->seq_region_name .":" .
      $feat->start .":" .
	$feat->end .":" .	
	  $feat->strand .":" .
	    $feat->hseqname;
    unless ($self->stored_features($unique_id)){
      push @dafs,$feat;
      # keep tabs on it so you don't store it again.
      $self->stored_features($unique_id,1);
    }
  } 
  return \@dafs;
}



###########################################################
# containers


sub rough {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_rough'} = $value;
  }

  if (exists($self->{'_rough'})) {
    return $self->{'_rough'};
  } else {
    return undef;
  }
}

sub fullslice {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_fullslice'} = $value;
  }

  if (exists($self->{'_fullslice'})) {
    return $self->{'_fullslice'};
  } else {
    return undef;
  }
}

sub write_attempts {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_writeattempts'} += $value;
  }

  if (exists($self->{'_writeattempts'})) {
    return $self->{'_writeattempts'};
  } else {
    return undef;
  }
}

sub stored_features {
  my ($self,$key,$value) = @_;

  return undef unless defined ($key);

  if (defined $value) {
    $self->{'_stored_features'}{$key} = $value;
  }

  if (exists($self->{'_stored_features'}{$key})) {
    return $self->{'_stored_features'}{$key};
  } else {
    return undef;
  }
}


sub TRANSDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_TRANSDB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_TRANSDB'})) {
    return $self->{'_CONFIG_TRANSDB'};
  } else {
    return undef;
  }
}

sub MISSMATCH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MISSMATCH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MISSMATCH'})) {
    return $self->{'_CONFIG_MISSMATCH'};
  } else {
    return undef;
  }
}


sub OUT_SAM_DIR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUT_SAM_DIR'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUT_SAM_DIR'})) {
    return $self->{'_CONFIG_OUT_SAM_DIR'};
  } else {
    return undef;
  }
}

sub BAM_FILE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BAM_FILE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BAM_FILE'})) {
    return $self->{'_CONFIG_BAM_FILE'};
  } else {
    return undef;
  }
}

sub COVERAGE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_COVERAGE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_COVERAGE'})) {
    return $self->{'_CONFIG_COVERAGE'};
  } else {
    return undef;
  }
}

sub PERCENT_ID {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PERCENT_ID'} = $value;
  }
  
  if (exists($self->{'_CONFIG_PERCENT_ID'})) {
    return $self->{'_CONFIG_PERCENT_ID'};
  } else {
    return undef;
  }
}

sub WORD_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_WORD_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_WORD_LENGTH'})) {
    return $self->{'_CONFIG_WORD_LENGTH'};
  } else {
    return undef;
  }
}


sub SATURATE_THRESHOLD {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SATURATE_THRESHOLD'} = $value;
  }
  
  if (exists($self->{'_CONFIG_SATURATE_THRESHOLD'})) {
    return $self->{'_CONFIG_SATURATE_THRESHOLD'};
  } else {
    return undef;
  }
}

sub MASK {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MASK'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MASK'})) {
    return $self->{'_CONFIG_MASK'};
  } else {
    return undef;
  }
}


sub FULLSEQ {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_FULLSEQ'} = $value;
  }
  
  if (exists($self->{'_CONFIG_FULLSEQ'})) {
    return $self->{'_CONFIG_FULLSEQ'};
  } else {
    return undef;
  }
}

sub MAX_TRANSCRIPT {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_TRANSCRIPT'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_TRANSCRIPT'})) {
    return $self->{'_CONFIG_MAX_TRANSCRIPT'};
  } else {
    return undef;
  }
}

sub BATCH_SIZE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BATCH_SIZE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BATCH_SIZE'})) {
    return $self->{'_CONFIG_BATCH_SIZE'};
  } else {
    return undef;
  }
}

1;
