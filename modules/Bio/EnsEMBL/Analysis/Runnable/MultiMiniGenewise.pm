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

Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise - 

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise->new(-genomic  => $genseq,
								  -features => $features)

    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION


=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable::MiniGenewise;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
  qw(attach_Slice_to_Transcript);
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args); 
  my( $features,$seqfetcher, $cluster_start, $cluster_end, $full_seq,
      $genewise_options, $minigenewise_options, $min_length) =
        rearrange([qw( FEATURES SEQFETCHER CLUSTER_START CLUSTER_END 
                       FULLSEQ GENEWISE_OPTIONS MINIGENEWISE_OPTIONS
                       MINIMUM_FEATURE_LENGTH)], @args);

  #SETTING DEFAULTS
  $self->minimum_feature_length(50);
  #################
  $self->seqfetcher($seqfetcher);
  $self->features($features);
  $self->cluster_start($cluster_start);
  $self->cluster_end($cluster_end);
  $self->full_seq($full_seq);
  $self->genewise_options($genewise_options);
  $self->minigenewise_options($minigenewise_options);
  $self->minimum_feature_length($min_length);
  throw("MultiMiniGenewise needs a query sequence ") unless($self->query);
  throw("MultiMiniGenewise needs features to run across")
    unless($self->features && scalar(@{$self->features}));
  throw("MultiMiniGenewise needs a sequence fetcher ") unless($self->seqfetcher);


  
  return $self;
}



#ACCESSOR METHODS

=head2 genewise_options

    Title   :   genewise_options
    Usage   :   $self->genewise_options($genewise_options)
    Function:   Get/set method for genewise genewise_options
    Returns :   
    Args    :   

=cut

sub genewise_options {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      throw("Must pass genewise options a hash ref not a $arg") 
        unless(ref($arg) eq "HASH");
        $self->{'_genewise_options'} = $arg;
    }

    return $self->{'_genewise_options'};
}

=head2 minigenewise_options

    Title   :   minigenewise_options
    Usage   :   $self->minigenewise_options($minigenewise_options)
    Function:   Get/set method for minigenewise minigenewise_options
    Returns :   
    Args    :   

=cut

sub minigenewise_options {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      throw("Must pass minigenewise options a hash ref not a $arg") 
        unless(ref($arg) eq "HASH");
        $self->{'_minigenewise_options'} = $arg;
    }

    return $self->{'_minigenewise_options'};
}

=head2 features

  Arg [1]   : arrayref 
  Function  : sets varible to arrayref
  Returntype: arrayref
  Exceptions: throws if not given an arrayref or if elements of array aren't featurepairs'
  Caller    : $self
  Example   : my $features = $self->features();

=cut




sub features {
  my ($self,$features) = @_;
  
  if (!defined($self->{_features})) {
    $self->{_features} = [];
  }
  if (defined($features)) {
    if (ref($features) eq "ARRAY") {
      foreach my $f (@$features) {
	if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
	  push(@{$self->{_features}},$f);
	} else {
	  throw("Object [$f] is not a Bio::EnsEMBL::FeaturePair");
	}
      }
    } else {
      throw("[$features] is not an array ref.");
    }
  }
  return $self->{_features};
}



=head2 genomic_sequence

  Arg [1]   : Object of correct type
  Function  : get/set object
  Returntype: object
  Exceptions: throws if object not of correct type
  Caller    : $self
  Example   : my $genomic = $self->query;

=cut

 
sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
      #need to check if passed sequence is Bio::Seq object
      $value->isa("Bio::PrimarySeqI") || throw("Input isn't a Bio::PrimarySeqI");
      $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
  }


sub seqfetcher {

  my( $self, $value ) = @_;    
  
  if (defined($value)) {
    #need to check if we are being passed a Bio::DB::RandomAccessI object
    throw("[$value] is not a Bio::DB::RandomAccessI") unless $value->isa("Bio::DB::RandomAccessI");
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}



sub minimum_intron {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_minimum_intron'} = $arg;
  }

  return $self->{'_minimum_intron'};
}


sub exon_padding {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_padding'} = $arg;
  }

  return $self->{'_padding'};
  
}

sub terminal_padding {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_terminal_padding'} = $arg;
  }

  return $self->{'_terminal_padding'};
  
}

sub cluster_start {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{_cluster_start} = $arg;
  }
  return $self->{_cluster_start};
}

sub cluster_end {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{_cluster_end} = $arg;
  }
  return $self->{_cluster_end};
}

sub full_seq {
  my( $self, $value ) = @_;    

  if ($value) {
    $self->{'_full_seq'} = $value;
  }

  return $self->{'_full_seq'};
}

sub minimum_feature_length {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{_minimum_feature_length} = $arg;
  }
  return $self->{_minimum_feature_length};
}


#Functionality


=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise
  Function  : sort through the features passed in and create MiniGenewise
  objects for them and run them
  Returntype: 1
  Exceptions: warns if a single id failed to find a sequence, throws if
  all the ids fail to find a sequence
  Example   : 

=cut



sub run{
  my ($self) = @_;

  my ($fhash,$ids) = $self->get_all_features_by_id;
  my $failed_count = 0;
  my @output;
  foreach my $id(@$ids){
    my $features = $fhash->{$id};
    #printf STDERR "MMG doing $id (%d feats)\n", scalar(@$features);
    my $peptide_sequence = $self->get_Sequence($features->[0]->hseqname);
    if (!$peptide_sequence) {
      warning($id.' has produced no peptide sequence from '.$self->seqfetcher);
      $failed_count++;
    }
    else {
      my @forward;
      my @reverse;

      foreach my $feature(@$features) {
        if ($feature->strand == -1) {
          push(@reverse, $feature);
        }
        elsif ($feature->strand == 1) {
          push(@forward, $feature);
        }
        else {
          throw("MultiMiniGenewise $feature from id $id seems to have no strand defined");
        }
      }

      logger_info('MMG have '.@forward.' forward strand features and '.@reverse." reverse strand features\n");
      my $slice = $self->query;
      if ($self->cluster_end) {
        print "MMG making subseq based on ".$self->cluster_start." ".$self->cluster_end." from ".$self->query->name."\n";
        my $string_seq = ('N' x ($self->cluster_start - 1)).
          $self->query->subseq($self->cluster_start, $self->cluster_end).
            ('N' x ($self->query->length - ($self->cluster_end)));
        $slice = Bio::EnsEMBL::Slice->new(
           -seq => $string_seq,
           -seq_region_name  => $self->query->seq_region_name,
           -start => $self->query->start,
           -end => $self->query->end,
           -coord_system => $self->query->coord_system,
          );
      }else{
        print "MMG using whole query ".$self->query->name."\n";
      }
      if (@forward >= 1) {
        my $forward_output = $self->run_MiniGenewise(\@forward, $slice, $peptide_sequence);
        $self->output($forward_output);
      }
      if (@reverse >= 1) {
        my $reverse_output = $self->run_MiniGenewise(\@reverse, $slice, $peptide_sequence);
        $self->output($reverse_output);
      }
    }
  }
  if ($failed_count == @$ids) {
    throw("Can't find any sequences for the ".@$ids.' ids which match '.
          $self->query->name);
  }
  if($self->cluster_end){
    foreach my $output(@{$self->output}){
      attach_Slice_to_Transcript($output, $self->query);
    }
  }
  return 1;
}


=head2 run_MiniGenewise

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Tools::MultiMiniGenewise
  Arg [2]   : arrayref of Features
  Arg [3]   : Dna sequence (Should be a primaryseq)
  Arg [4]   : Protein sequence (should be a primary seq)
  Function  : Create and run a MiniGenewise object on the given features
  Returntype: arrayref of transcripts from MiniGenewise
  Exceptions: 
  Example   : 

=cut


sub run_MiniGenewise{
  my ($self, $features, $slice, $peptide) = @_;
  my @to_use = @{$self->find_extras($features)};
 
  return () if(@to_use == 0);
  if($self->full_seq){
    my $rangefeat = new Bio::EnsEMBL::FeaturePair
      (
       -start => $features->[0]->start,
       -end   => $features->[-1]->end,
       -strand=> 1,
       -slice => $slice
      );
    $rangefeat->strand(-1) if($features->[0]->strand == -1);
    $features = [$rangefeat];
  }
  my %params = %{$self->minigenewise_options} if($self->minigenewise_options); 
  my $runnable  = Bio::EnsEMBL::Analysis::Runnable::MiniGenewise->new
    (
     -query            => $slice,
     -protein          => $peptide,
     -features         => $features,
     -genewise_options => $self->genewise_options,
     -analysis         => $self->analysis,
     %params,
    );
  %params = ();
  $runnable->run;
  my @output = @{$runnable->output};
  return(\@output);
}


=head2 get_all_features_by_id

  Arg [1]   : none
  Function  : arranges all feature into hash keyed by hseqname, each element 
  containing an anonymous array of features with that name also produces an hash of 
  key hseqname and value of score which is used to sort an array of hseqname 
  Returntype: hasfref and array ref
  Exceptions: warns and skips if a feature doesn't have a hseqname'
  Caller    : $self
  Example   : my ($idhash, $idarray) = $self->get_all_features_by_id;

=cut

sub get_all_features_by_id {
  my ($self) = @_;

  my  %idhash;
  my  %scorehash;

 FEAT: foreach my $f (@{$self->features}) {
    if (!$f->hseqname) {
      warning("No hit name for " . $f->seqname . "\n");
      next FEAT;
    }

    $idhash{$f->hseqname} = [] if(!$idhash{$f->hseqname});
    push(@{$idhash{$f->hseqname}},$f);

    $scorehash{$f->hseqname} = 0 if(!$scorehash{$f->hseqname});
    $scorehash{$f->hseqname} = $f->score if($f->score > 
                                            $scorehash{$f->hseqname});;
  }

  my @ids = keys %idhash;
  @ids = sort {$scorehash{$b} <=> $scorehash{$a}} @ids;
  return (\%idhash,\@ids);
}



=head2 get_Sequence

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise
  Arg [2]   : string, peptide id
  Function  : fetch given peptide from sequence index
  Returntype: Bio::PrimarySeq
  Exceptions: warns if not passed an id or if the sequence cannot be found
  Example   :

=cut



sub get_Sequence {
  my ($self,$id) = @_;
  my $seqfetcher = $self->seqfetcher;
  my $seq;
  my $name;
  if($seqfetcher->can('db')){
    my @dbs = $seqfetcher->db;
    $name = $dbs[0];
  }
  if (!$id) {
    warning("No id input to get_Sequence");
    return undef;
  }

  eval {
    $seq = $seqfetcher->get_Seq_by_acc($id);
  };

  if($@) {
    warning("Problem fetching sequence for id [$id] with $seqfetcher $name  $@\n");
    return undef;
  }

  if(!$seq){
    warning("Could not find sequence for [$id] with $seqfetcher $name");
  }

  return $seq;
}


sub find_extras{
  my ($self, $features) = @_;
  my $output = $self->output;
  my @out;
  foreach my $f(@$features){
    my $found = 0;
    $f->slice($self->query);
    $f->seqname($f->slice->name);
    if($f->length <  $self->minimum_feature_length){
      next;
    }
    foreach my $transcript(@$output){
      foreach my $exon(@{$transcript->get_all_Exons}){
        $exon->slice($self->query);
        if($f->overlaps($exon)){
          $found = 1;
        }
      }
    }
    if($found == 0){
      push(@out, $f);
    }
  }
  return \@out;
}

1;

