=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

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

=head1 METHODS

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
  my ( $model,$missmatch,$mask,$sam_dir,$bam_file,$percent_id,$coverage,$fullseq,$max_tran,$start,$batch_size ) =
    rearrange( [qw(MODEL MISSMATCH MASK SAM_DIR BAM_FILE  PERCENT_ID COVERAGE FULLSEQ MAX_TRAN START BATCH_SIZE)],@args );
  $self->model($model);
  $self->MASK($mask);
  $self->MISSMATCH($missmatch);
  $self->OUT_SAM_DIR($sam_dir);
  $self->BAM_FILE($bam_file);
  $self->PERCENT_ID($percent_id);
  $self->COVERAGE($coverage);
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
      # dont want reads that align perfectly as they won't splice
      my $num_missmatches = $read->get_tag_values('NM') ;	
      my $rg = "*";
      if ( $read->get_tag_values('RG') ) {
        $rg = $read->get_tag_values('RG') ;
      }
      next READ  unless $num_missmatches >= $self->MISSMATCH;
      $batch++;
      next unless $batch > $batch_size * $self->start ;
      $count++;
      if (  $count <= $batch_size ) {
	my $suffix ="";
	# is it the 1st or 2nd read?
	if ( $read->get_tag_values('FLAGS') =~ /FIRST_MATE/ ) {
	  $suffix = "/1";
	}
	if ( $read->get_tag_values('FLAGS') =~ /SECOND_MATE/ ) {
	  $suffix = "/2";
	}
	my $name = $read->name.$suffix;
	# write the seq files - store the read group information in case it is needed later
	my $bioseq = Bio::Seq->new( 
				   -seq        => $read->query->dna,
				   -display_id => $name,
				   -desc       => $rg
				  );
	push @reads, $bioseq;
      # want to store the read sequence for making it into a SAM file later
      $self->seq_hash($name,$bioseq);
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
    # filter on coverage and percent id
    next unless $feat->percent_id >= $self->PERCENT_ID;
    next unless $feat->hcoverage >= $self->COVERAGE;
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
