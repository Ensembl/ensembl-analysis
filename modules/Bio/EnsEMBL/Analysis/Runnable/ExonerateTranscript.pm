
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
								 -query_seqs     => \@q_seqs,
                                                             [or -query_file     => $q_file]   
								 -query_type     => 'dna',
								 -target_seqs    => \@t_seqs,
                                                             [or -target_file    => $t_file]   
                                                                 -exonerate      => $exonerate,
								 -options        => $options,
								);

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, namely 
the prediction of the transcript structure in a piece of genomic DNA by the alignment 
of a 'transcribed' sequence (ESt, cDNA or protein). The results is a set of 
Bio::EnsEMBL::Transcript objects


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($query_type, $query_seqs, $query_file, $q_chunk_num, $q_chunk_total,
      $target_file, 
      $verbose) = 
          rearrange([qw(
                        QUERY_TYPE
                        QUERY_SEQS                        
                        QUERY_FILE
                        QUERY_CHUNK_NUMBER
                        QUERY_CHUNK_TOTAL
                        TARGET_FILE
                        VERBOSE
                        )
                     ], @args);
  $self->_verbose($verbose) if $verbose;

  if (defined($query_seqs)) {     
    throw("You must supply an array reference with -query_seqs") 
        if ref($query_seqs) ne "ARRAY";
    $self->query_seqs($query_seqs);
  }
  elsif (defined $query_file) {
    throw("The given query file does not exist") if ! -e $query_file;
    $self->query_file($query_file);
  }

  if ($query_type){
    $self->query_type($query_type);
  } else{
    # default to DNA for backwards compatibilty
    $self->query_type('dna');
  }

  if (defined $target_file) {
    throw("The given database does not exist") if ! -e $target_file;
    $self->target_file($target_file);
  }

  ############################################################
  # We default exonerate-0.6.7
  if (not $self->program) {
    $self->program('/usr/local/ensembl/bin/exonerate-0.8.3');
  }

  ############################################################
  # options
  my $basic_options = "--ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" "; 
  if (defined $q_chunk_num and defined $q_chunk_total) {
    $basic_options .= "--querychunkid $q_chunk_num --querychunktotal $q_chunk_total ";
  }

  if ($self->options){
    $basic_options .= $self->options;
  }
  $self->options($basic_options);

  return $self;
}



############################################################
#
# Analysis methods
#
############################################################

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs exonerate script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output
=cut

sub run {
  my ($self) = @_;

  if ($self->query_seqs) {
    # Write query sequences to file if necessary
    my $query_file = $self->workdir . "/exonerate_q.$$";
    my $seqout = Bio::SeqIO->new('-format' => 'fasta',
                                 '-file'     => ">$query_file");   
    foreach my $seq ( @{$self->query_seqs} ) {
      $seqout->write_seq($seq);
    }
    # register the file for deletion
    $self->files_to_delete($query_file);
    $self->query_file($query_file);
  }

  # Build exonerate command

  my $command =$self->program . " " .$self->options .
      " --querytype "  . $self->query_type .
      " --targettype " . $self->target_type .
      " --query "  . $self->query_file .
      " --target " . $self->target_file;
  
  # Execute command and parse results

  print STDERR "Exonerate command : $command\n" if $self->_verbose;

  my $exo_fh;
  open( $exo_fh, "$command |" ) or throw("Error opening exonerate command: $? $!");
  $self->parse_results( $exo_fh );
  close( $exo_fh ) or throw ("Error closing exonerate command: $? $!");
  $self->delete_files;

  return 1;
}

sub parse_results {
  my ($self, $fh) = @_;

  my %strand_lookup = ('+' => 1, '-' => -1, '.' => 1);

  # Each alignment will be stored as a transcript with 
  # exons and supporting features.  Initialise our
  # transcript.

  my @transcripts;

  # Parse output looking for lines beginning with 'RESULT:'.
  # Each line represents a distinct match to one sequence
  # containing multiple 'exons'.

 TRANSCRIPT:
  while (<$fh>){
    print STDERR $_ if $self->_verbose;

    next unless /^RESULT:/;

    chomp;

    my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end,
	$t_strand, $score, $perc_id, $q_length, $t_length, $gene_orientation,
	@align_components) = split;
   
    $t_strand = $strand_lookup{$t_strand};
    $q_strand = $strand_lookup{$q_strand};
    $gene_orientation = $strand_lookup{$gene_orientation};
        
    my $coverage = sprintf("%.2f", 100 * (abs($q_end - $q_start) / $q_length));

    # Read vulgar information and extract exon regions.
    my $exons = $self->_parse_vulgar_block($t_start,
                                           $t_end,
                                           $t_strand,
                                           $t_length,
                                           $q_start, 
                                           $q_end,
                                           $q_strand,
                                           $q_length,
                                           \@align_components);

    # now we have extracted the exons and the coordinates are with 
    # reference to the forward strand of the query and target, we can 
    # use the gene_orienation to flip the strands if necessary
    if ($gene_orientation == -1 and $t_strand == 1) {
      $t_strand *= -1;
      $q_strand *= -1;
    }

    my $transcript = Bio::EnsEMBL::Transcript->new();

    # Build FeaturePairs for each region of query aligned to a single
    # Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
    # attach this to our Exon.

    foreach my $proto_exon (@$exons){
      
      # Build our exon and set its key values.
      my $exon = Bio::EnsEMBL::Exon->new();
      
      $exon->seqname($t_id);
      $exon->start($proto_exon->{exon_start});
      $exon->end($proto_exon->{exon_end});
      $exon->phase($proto_exon->{phase});
      $exon->end_phase($proto_exon->{end_phase});
      $exon->strand($t_strand);
            
      my @feature_pairs;
      foreach my $sf (@{$proto_exon->{sf}}){
        my $feature_pair = Bio::EnsEMBL::FeaturePair->new(-seqname    => $t_id,
                                                          -start      => $sf->{target_start},
                                                          -end        => $sf->{target_end},
                                                          -strand     => $t_strand,
                                                          -hseqname   => $q_id,
                                                          -hstart     => $sf->{query_start},
                                                          -hend       => $sf->{query_end},
                                                          -hstrand    => $q_strand,
                                                          -score      => $coverage,
                                                          -percent_id => $perc_id);

	push @feature_pairs, $feature_pair;

      }

      # Use our feature pairs for this exon to create a single 
      # supporting feature (with cigar line).
      my $supp_feature;

      eval{
        if ($self->query_type eq 'protein') {
          $supp_feature =
              Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@feature_pairs);
        } else {
          $supp_feature = 
              Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@feature_pairs);
        }
      };
      if ($@){
        warning($@);
        next TRANSCRIPT;
      }
      
      $exon->add_supporting_features($supp_feature);
      
      $transcript->add_Exon($exon);
    }

    my @exons = @{$transcript->get_all_Exons};
    if (scalar(@exons)) {

      if ($self->query_type eq 'protein') {
        # add a translation if this is a protein alignment

        my $translation = Bio::EnsEMBL::Translation->new();

        $translation->start_Exon($exons[0]);
        $translation->end_Exon  ($exons[-1]);
        
        # phase is relative to the 5' end of the transcript (start translation)
        if ($exons[0]->phase == 0) {
          $translation->start(1);
        } elsif ($exons[0]->phase == 1) {
          $translation->start(3);
        } elsif ($exons[0]->phase == 2) {
          $translation->start(2);
        }
        $translation->end($exons[-1]->end - $exons[-1]->start + 1);

        $transcript->translation($translation);
      }

      push @transcripts, $transcript;
    }

  }

  $self->output(\@transcripts);

  return 1;
}

sub _parse_vulgar_block {
  my ($self, 
      $target_start, $target_end, $target_strand, $target_length,
      $query_start, $query_end,  $query_strand, $query_length,
      $vulgar_components) = @_;

  # This method works along the length of a vulgar line 
  # exon-by-exon.  Matches that comprise an exon are 
  # grouped and an array of 'proto-exons' is returned.
  # Coordinates from the vulgar line are extrapolated 
  # to actual genomic/query coordinates.

  my @exons;
  my $exon_number = 0;


  # We sometimes need to increment all our start coordinates. Exonerate 
  # has a coordinate scheme that counts _between_ nucleotides at the start.
  # However, for reverse strand matches 
  
  my ($query_in_forward_coords, $target_in_forward_coords);
  my ($cumulative_query_coord, $cumulative_target_coord);

  if ($target_start > $target_end) {
    warn("For target, start and end are in thew wrong order for a reverse strand match")
        if $target_strand != -1;
    $cumulative_target_coord = $target_start;
    $target_in_forward_coords = 1;
  } else {
    $cumulative_target_coord = $target_start + 1;
    $target_in_forward_coords = 0;
  }
  if ($query_start > $query_end) {
    warn("For query, start and end are in thew wrong order for a reverse strand match")
        if $query_strand != -1;
    $cumulative_query_coord = $query_start;
    $query_in_forward_coords = 1;
  } else {
    $cumulative_query_coord = $query_start + 1;
    $query_in_forward_coords = 0;
  }


  while (@$vulgar_components){
    throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @$vulgar_components . "] items left to process.")
      unless scalar @$vulgar_components >= 3;

    my $type                = shift @$vulgar_components;
    my $query_match_length  = shift @$vulgar_components;
    my $target_match_length = shift @$vulgar_components;

    throw("Vulgar string does not start with a match.  Was not " . 
		 "expecting this.")
      if ((scalar @exons == 0) && ($type ne 'M'));

    if ($type eq 'M'){
      my %hash;

      if ($target_strand == -1) {
        if ($target_in_forward_coords) {
          $hash{target_start} = $cumulative_target_coord - ($target_match_length - 1);
          $hash{target_end}   = $cumulative_target_coord;
        } else {
          $hash{target_end}   = $target_length - ($cumulative_target_coord - 1);
          $hash{target_start} = $hash{target_end} - ($target_match_length - 1);
        }
      } else {
        $hash{target_start} = $cumulative_target_coord;
        $hash{target_end}   = $cumulative_target_coord + ($target_match_length - 1);
      }

      if ($query_strand == -1) {
        if ($query_in_forward_coords) {
          $hash{query_start} = $cumulative_query_coord - ($query_match_length - 1);
          $hash{query_end}   = $cumulative_query_coord;
        } else {
          $hash{query_end}   = $query_length - ($cumulative_query_coord - 1);
          $hash{query_start} = $hash{query_end} - ($query_match_length - 1);
        }
      } else {
        $hash{query_start} = $cumulative_query_coord;
        $hash{query_end}   = $cumulative_query_coord + ($query_match_length - 1);
      }

      # there is nothing to add if this is the last state of the exon
      $exons[$exon_number]->{gap_end}   = 0;
      push @{$exons[$exon_number]->{sf}}, \%hash;
    }
    elsif ($type eq "S") {
      if ($exons[$exon_number]) {
        # this is a split codon at the end of an exon
        $exons[$exon_number]->{split_end} = $target_match_length;
      } else {
        $exons[$exon_number]->{split_start} = $target_match_length;
      }
    }
    elsif ($type eq "G") {
      if (exists($exons[$exon_number]->{sf})) {
        # this is the gap in the middle of an exon, or at the end. Assume it is 
        # at the end, and then reset if we see another match state in this exon
        $exons[$exon_number]->{gap_end}   = $target_match_length;
      } else {
        # this is a gap at the start of an exon; 
        $exons[$exon_number]->{gap_start} = $target_match_length;
      }
    }
    elsif ($type eq "I" or
           $type eq "F") {

      # in protein mode, any insertion on the genomic side should be treated as 
      # an intron to ensure that the result translates. However, we allow for
      # codon insertions in the genomic sequence with respect to the protein. 
      # This introduces the possibility of in-frame stops, but I don't
      # think "introning over" these insertions is appropriate here. 

      # if we see a gap/intron immediately after an intron, the current exon is "empty"
      if ($exons[$exon_number]) {
        $exon_number++;
      }
    }

    if ($target_in_forward_coords and $target_strand == -1) {
      $cumulative_target_coord -= $target_match_length;
    } else {
      $cumulative_target_coord += $target_match_length;
    }
    if ($query_in_forward_coords and $query_strand == -1) {
      $cumulative_query_coord  -= $query_match_length;
    }
    else {
      $cumulative_query_coord  += $query_match_length;
    }

  }

  for(my $i = 0; $i < @exons; $i++) {
    my $ex = $exons[$i];
    my $ex_sf = $ex->{sf};

    $ex->{phase} = 0;
    $ex->{end_phase} = 0;
    
    if ($target_strand == -1) {
      $ex->{exon_start} = $ex_sf->[-1]->{target_start};
      $ex->{exon_end}   = $ex_sf->[0]->{target_end};

      if (exists $ex->{split_start}) {
        $ex->{exon_end} += $ex->{split_start};
        $ex->{phase} = 3 - $ex->{split_start};
      }
      if (exists $ex->{split_end}) {
        $ex->{exon_start} -= $ex->{split_end};
        $ex->{end_phase} = $ex->{split_end};
      }
      if (exists $ex->{gap_start}) {
        $ex->{exon_end} += $ex->{gap_start};
      }
      if (exists $ex->{gap_end}) {
        $ex->{exon_start} -= $ex->{gap_end};
      }

    } else {
      $ex->{exon_start} = $ex_sf->[0]->{target_start};
      $ex->{exon_end}   = $ex_sf->[-1]->{target_end};

      if (exists $ex->{split_start}) {
        $ex->{exon_start} -= $ex->{split_start};
        $ex->{phase} = 3 - $ex->{split_start};
      }
      if (exists $ex->{split_end}) {
        $ex->{exon_end} += $ex->{split_end};
        $ex->{end_phase} = $ex->{split_end};
      }
      if (exists $ex->{gap_start}) {
        $ex->{exon_start} -= $ex->{gap_start};
      }
      if (exists $ex->{gap_end}) {
        $ex->{exon_end} += $ex->{gap_end};
      }
    }
  }

  return \@exons;
}


############################################################
#
# get/set methods
#
############################################################


############################################################

sub query_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'protein' ){
      throw("not the right query type: $type");
    }
    $self->{_query_type} = $type;
  }
  return $self->{_query_type};
}

############################################################

sub query_seqs {
  my ($self, $seqs) = @_;
  if ($seqs){
    unless ($seqs->[0]->isa("Bio::PrimarySeqI") || $seqs->[0]->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_query_seqs} = $seqs;
  }
  return $self->{_query_seqs};
}


############################################################

sub query_file {
  my ($self, $file) = @_;
  
  if ($file) {
    $self->{_query_file} = $file;
  }
  return $self->{_query_file};
}

############################################################

sub target_type {
  my ($self) = @_;

  # the target type has to be DNA, because we are making transcripts

  return 'dna';
}

############################################################

sub target_file {
  my ($self, $file) = @_;
  
  if ($file) {
    $self->{_target_file} = $file;
  }
  return $self->{_target_file};
}

############################################################

sub _verbose {
  my ($self, $val) = @_;
  
  if ($val){
    $self->{_verbose} = $val;
  }
  
  return $self->{_verbose};
}


1;

