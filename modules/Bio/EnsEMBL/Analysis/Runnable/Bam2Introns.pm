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

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::SeqFeature::Lite;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::DB::Sam;

$| = 1;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ( $model,$missmatch,$mask,$sam_dir,$bam_file,$percent_id,$coverage,$fullseq,$max_tran,$start,$start_exon, $offset, $batch_size ) =
    rearrange( [qw(MODEL MISSMATCH MASK SAM_DIR BAM_FILE  PERCENT_ID COVERAGE FULLSEQ MAX_TRAN START STARTEXON OFFSET BATCH_SIZE)],@args );
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
  $self->start_exon($start_exon);
  $self->offset($offset);
  $self->BATCH_SIZE($batch_size);
  return $self;
}

####################################################
# batching works by running initialy on the first 
# 100,000 reads, once it gets past this is just counts
# how many reads there are and submits new input_ids
# of the type STABLEID00000001:2 for the second bacth 
# of 100,000 reads etc.

sub _increase_read_offset {
    my $self = shift;

    $self->{'_read_offset'}++;
}

sub _reset_read_offset {
    my $self = shift;

    $self->{'_read_offset'} = 0;
}

sub read_offset {
    my $self = shift;

    return $self->{'_read_offset'};
}

sub _last_start {
    my ($self, $value) = @_;

    $self->{'_last_start'} = $value if (defined $value);
    return $self->{'_last_start'};
}

sub _increase_read_count {
    my $self = shift;

    $self->{'_read_count'}++;
}

sub _reset_read_count {
    my $self = shift;

    $self->{'_read_count'} = 0;
}

sub read_count {
    my $self = shift;

    return $self->{'_read_count'};
}

sub _increase_count {
    my $self = shift;

    $self->{'_count'}++;
}

sub _reset_count {
    my $self = shift;

    $self->{'_count'} = 0;
}

sub count {
    my $self = shift;

    return $self->{'_count'};
}

sub _reset_batch {
    my $self = shift;

    $self->{'_batch'} = 0;
}

sub _increase_batch {
    my $self = shift;

    $self->{'_batch'}++;
}

sub batch {
    my $self = shift;

    return $self->{'_batch'};
}

sub _process_reads {
    my ($read, $callbackdata) = @_;

    my ($self, $i, $exon_start, $stable_id, $reads, $iids) = @$callbackdata;
    return if ($self->kill_loop);
    my $missmatch  = $read->get_tag_values('NM');
    return unless ($missmatch and $missmatch >= $self->MISSMATCH);
#    unless ($missmatch and $missmatch >= $self->MISSMATCH) {
#        print STDERR 'TIBO missmatch: ', $read->query->name, ' : ', $missmatch, ' @ ', $read->get_tag_values('FLAGS'), ' ? ', $read->get_tag_values('MD'), "\n";
#        return;
#    }
    # get rid of any reads that might start before our start point
    return if (!$read->start or ($i == $self->start_exon and $read->start < $exon_start));
#    if (!$read->start or ($i == $self->start_exon and $read->start < $exon_start)) {
#        print STDERR 'TIBO exon_start: ', $read->query->name, ' : ', $read->start, ' < ', $exon_start, ' @ ', $read->get_tag_values('FLAGS'), ' ? ', $read->get_tag_values('MD'), ' ! ', $read->query->start, ' + ', $read->target->start, ' $ ', $read->mate_start, "\n";
#        return;
#    }
    $self->_increase_read_count;
    # get rid of any reads overlapping out start point
    return if ($self->offset and $self->read_count <= $self->offset);
#    if ($self->offset and $self->read_count <= $self->offset) {
#        print STDERR 'TIBO offset: ', $read->query->name, ' : ', $self->offset, ' @ ', $read->get_tag_values('FLAGS'), "\n";
#        return;
#    }
    # calculate the offset 
#    print STDERR 'TIBO GOOD: ', $read->query->name,
#        ' start ', $read->start, ";\n",
#        ' end ', $read->end, ";\n",
#        ' length ', $read->length, ";\n",
#        ' mate_start ', $read->mate_start, ";\n",
#        ' mate_end ', $read->mate_end, ";\n",
#        ' mate_len ', $read->mate_len, ";\n",
#        ' strand ', $read->strand, ";\n",
#        ' mstrand ', $read->mstrand, ";\n",
#        ' query->name ', $read->query->name, ";\n",
#        ' query->start ', $read->query->start, ";\n",
#        ' query->end ', $read->query->end, ";\n",
#        ' query->length ', $read->query->length, ";\n",
#        ' target->name ', $read->target->name, ";\n",
#        ' target->start ', $read->target->start, ";\n",
#        ' target->end ', $read->target->end, ";\n",
#        ' primary_id ', $read->primary_id, ";\n",
#        ' get_all_tags ', $read->get_all_tags, ";\n",
#        ' cigar_str ', $read->cigar_str, ";\n",
#        ' aux ', $read->aux, ";\n",
#        ' tid ', $read->tid, ";\n",
#        ' qname ', $read->qname, ";\n",
#        ' pos ', $read->pos, ";\n",
#        ' flag ', $read->flag, ";\n",
#        ' mtid ', $read->mtid, ";\n",
#        ' mpos ', $read->mpos, ";\n",
#        ' n_cigar ', $read->n_cigar, ";\n",
#        ' aux_keys ', $read->aux_keys, ";\n",
#        ' paired ', $read->paired, ";\n",
#        ' proper_pair ', $read->proper_pair, ";\n",
#        ' unmapped ', $read->unmapped, ";\n",
#        ' munmapped ', $read->munmapped, ";\n",
#        ' reversed ', $read->reversed, ";\n",
#        ' mreversed ', $read->mreversed, ";\n",
#        ' isize ', $read->isize, ";\n",
#        ' < ', $exon_start, ";\n";
#        foreach my $tag ($read->get_all_tags) {
#            print STDERR "\t", $tag, ': ', $read->get_tag_values($tag);
#        }
#        print STDERR "\n";
    if ( $read->start == $self->_last_start ) {
        $self->_increase_read_offset;
    } else {
        $self->_reset_read_offset;
        $self->_last_start($read->start);
    }

    my $rg = $read->get_tag_values('RG') || '*';

    $self->_increase_batch;
    $self->_increase_count;

    if (  $self->count <= $self->BATCH_SIZE ) {
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
        $self->seq_hash($name,$bioseq);
    } else {
        # if running a batch finish here
#        print STDERR 'TIBO: entering iids generation ', $self->start, ' : ', $self->batch, "\n";
        if ($self->start) {
            $self->_kill_loop;
            return;
        }
        if ($self->BATCH_SIZE+1 == $self->batch ) {
#            print STDERR 'TIBO : creating iids', $self->batch, "\n";
            $self->_reset_batch;
            $self->_increase_batch;
            push(@$iids, $stable_id .":$i:" . $read->start .':'.$self->read_offset);
        }
        # otherwise figure out the ids for the rest of the batches
    }
}

sub _kill_loop {
    my $self = shift;

    $self->{'_kill_loop'} = 1;
}

sub kill_loop {
 my $self = shift;

 return $self->{'_kill_loop'};
}

sub _loop {
    my $self = shift;

    $self->{'_kill_loop'} = 0;
}

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
  my $bam = Bio::DB::Bam->open($self->BAM_FILE);
  my $header = $bam->header();
  my ($tid, $start, $end) = $header->parse_region($rough->seq_region_name);
  my $bam_index = Bio::DB::Bam->index_open($self->BAM_FILE);
  $self->throw("Bam file " . $self->BAM_FILE . "  not found \n") unless $bam_index; 
  $self->_reset_count;
  $self->_reset_batch;
  my $exon = 0;
  my @exons = sort { $a->start <=> $b->start } @{$rough->get_all_Exons};
  my @iids;
  $self->_loop;

 EXON: for ( my $i = $self->start_exon ; $i <= $#exons ; $i++ ){
    my $exon = $exons[$i];
    my $exon_start = $exon->start;
    $self->_reset_read_offset;
    $self->_last_start(0);
    my %split_points;
    $exon_start = $self->start if $self->start && $self->start_exon == $i;
    my @callback_data = ($self, $i, $exon_start, $rough->stable_id, \@reads, \@iids);
    $bam_index->fetch($bam, $tid, $exon_start-1, $exon->end-1, \&_process_reads, \@callback_data);
  }
  if (  scalar(@reads) == 0 ) {
    $self->delete_files();
    return 0;
  }
  if ( scalar(@iids) > 0  && !$self->start ) {
    print 'Making ',  scalar(@iids), ' New input ids from ', $self->count, " reads\n";
    $self->iids(\@iids);
  }
  $self->write_seq_file(\@reads,$query_seq);
  foreach my $val (@iids) {
      print STDERR $val, "\n";
  }
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

sub start_exon {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_start_exon'} = $value;
  }
  
  if (exists($self->{'_start_exon'})) {
    return $self->{'_start_exon'};
  } else {
    return 0;
  }
}

sub offset {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_offset'} = $value;
  }
  
  if (exists($self->{'_offset'})) {
    return $self->{'_offset'};
  } else {
    return 0;
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
