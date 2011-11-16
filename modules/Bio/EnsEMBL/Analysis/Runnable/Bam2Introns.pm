
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Bam2Introns

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

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Bam2Introns;

use vars qw(@ISA);
use strict;

use Bio::SeqFeature::Lite;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::DB::Sam;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ( $model,$missmatch,$mask,$sam_dir,$bam_file,$read_length,$fullseq,$max_tran,$start,$batch_size ) =
    rearrange( [qw(MODEL MISSMATCH MASK SAM_DIR BAM_FILE  READ_LENGTH FULLSEQ MAX_TRAN START BATCH_SIZE)],@args );
  $self->model($model);
  $self->MASK($mask);
  $self->MISSMATCH($missmatch);
  $self->OUT_SAM_DIR($sam_dir);
  $self->BAM_FILE($bam_file);
  $self->READ_LENGTH($read_length);
  $self->FULLSEQ($fullseq);
  $self->MAX_TRANSCRIPT($max_tran);
  $self->start($start);
  $self->BATCH_SIZE($batch_size);
  return $self;
}

####################################################
# batching works by running initialy on the first 
# 100,000 reads, once it gets past this is just counts
# how many reads there are and submits new input_ids
# of the type STABLEID00000001:2 for the second bacth 
# of 100,000 reads etc.

sub run  {
  my ( $self) = @_;
  my $batch_size = $self->BATCH_SIZE;
  print "Using batch size of $batch_size\n";
  # set up the output files
  my $rough = $self->model;
  my $query_seq = $self->create_filename("B2I_reads","fa");
  my $genomic_seq = $self->create_filename("B2I_transcript","fa");
  $self->files_to_delete($query_seq);
  $self->files_to_delete($genomic_seq);
  $self->query_file($query_seq);
  $self->target_file($genomic_seq);
  # do we want our transcript sequences to be masked
  
  my $seq;
  if ( $self->FULLSEQ ) {
    $seq = $self->FULLSEQ;
  } else {
    foreach my $exon ( @{$rough->get_all_Transcripts->[0]->get_all_Exons}) {
      my $slice = $exon->feature_Slice;
      if ( $self->MASK ) {
	$seq .= $slice->get_repeatmasked_seq(undef,1)->seq;
      } else {
	$seq .= $slice->seq;
      }
    }
  }
  my $seqio = Bio::Seq->new( -seq => $seq,
			     -display_id => $rough->stable_id
			   );
  $self->write_seq_file($seqio,$genomic_seq);
  # get the reads
  my @reads;
  my $bam = Bio::DB::Sam->new(   -bam => $self->BAM_FILE,
				 -autoindex => 1,
                             );
  $self->throw("Bam file " . $self->BAM_FILE . "  not found \n") unless $bam; 
  my $count = 0;
  my $batch = 0;

 EXON: foreach my $exon ( @{$rough->get_all_Exons} ) {
    my $segment = $bam->segment($rough->seq_region_name,$exon->start,$exon->end);
    my $iterator = $segment->features(-iterator=>1);
  READ:  while (my $read = $iterator->next_seq) {
      # in case you have bam files with mixed length reads
      # you will need the reads in each batch to be the same
      # length for the exonerate config to work
      next unless $read->length == $self->READ_LENGTH;
      # dont want reads that align perfectly as they won't splice
      my $num_missmatches = $read->get_tag_values('NM') ;
      next READ  unless $num_missmatches >= $self->MISSMATCH;
      $batch++;
      next unless $batch > $batch_size * $self->start ;
      $count++;
      if (  $count <= $batch_size ) {
	#print ".";
	# write the seq files
	my $bioseq = Bio::Seq->new( 
				   -seq => $read->query->dna,
				   -display_id => $read->name
				  );
	push @reads, $bioseq;
      # want to store the read sequence for making it into a SAM file later
      $self->seq_hash($read->name,$bioseq);
      } else {
	last EXON if $self->start;
      }
      ###      print "BATCH $batch $count\n";
    }
    #print "\n" . $exon->start . " " , $exon->end . " = $count reads \n";
  }
  if (  scalar(@reads) == 0 ) {
    $self->delete_files();
    return 0;
  }
  if ( $batch > $batch_size && !$self->start ) {
    my $num = int($batch / $batch_size);
    print "Making $num New input ids from $batch reads\n";
    my @iids;
    for ( my $i = 1 ; $i <= $num ; $i++ ) {
      push @iids, $rough->stable_id .":$i";
      print $rough->stable_id .":$i\n";
    }
    $self->iids(\@iids);
  }
  $self->write_seq_file(\@reads,$query_seq);
  # now to run the runnable
  $self->SUPER::run();  # attach the read seq to the output features
  $self->process_features;
}

sub process_features {
  my ($self) = @_;
  my $features = $self->output;
  my @new_features;
  foreach my $feat ( @$features ) {
    my $feature_seq = $self->seq_hash($feat->hseqname);
    # print $feat->hseqname ." $feature_seq\n";
    $self->throw("cannot find sequence for " . $feat->hseqname ."\n")
    unless $feature_seq;
    # hide it in the hash
    $feat->{"_feature_seq"} = $feature_seq;
    push @new_features, $feat;
  }

  # clear the output arrray
  $self->{'output'} = [];
  $self->output(\@new_features);
}


###########################################################
# containers

sub seq_hash {
  my ($self,$key,$value) = @_;

  return undef unless defined ($key);

  if (defined $value) {
    $self->{'_seq_hash'}{$key} = $value;
  }

  if (exists($self->{'_seq_hash'}{$key})) {
    return $self->{'_seq_hash'}{$key};
  } else {
    return undef;
  }
}

sub iids { 
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_IIDS'} = $value;
  }
  
  if (exists($self->{'_IIDS'})) {
    return $self->{'_IIDS'};
  } else {
    return undef;
  }
}

sub model { 
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_MODEL'} = $value;
  }
  
  if (exists($self->{'_MODEL'})) {
    return $self->{'_MODEL'};
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

sub READ_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_READ_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_READ_LENGTH'})) {
    return $self->{'_CONFIG_READ_LENGTH'};
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

sub start {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_start'} = $value;
  }
  
  if (exists($self->{'_start'})) {
    return $self->{'_start'};
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
