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

Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise - 

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise->new
    ('-genomic'        => $genseq,
     '-features'       => $features,
     '-protein'        => $protein,
     '-seqfetcher'     => $seqfetcher,
     '-check_repeated' => 1);

    
    $obj->run

    my @newfeatures = $obj->output;

(where $protein and $genseq are Bio::Seq objects, 
 $features are X objects and $seqfetcher is a 
 SeqFetcher object.)


=head1 DESCRIPTION

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable::BlastMiniBuilder;
use Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise;
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::BlastDB;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info logger_verbosity);
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use Bio::EnsEMBL::Analysis::Tools::FeatureFilter;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BlastMiniBuilder);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my( $ids, $seqfetcher, $minigenewise_options, $check_repeated, 
      $full_seq,$exonerate,$exonerate_path, $exonerate_options, $genewise_options, 
      $blast_parser, $blast_filter, $post_filter, $score_cutoff, $blastdb, $mmg_options, $max_intron) =
        rearrange([qw(IDS SEQFETCHER MINIGENEWISE_OPTIONS
                      CHECK_REPEATED FULLSEQ EXONERATE EXONERATE_PATH 
                      EXONERATE_OPTIONS GENEWISE_OPTIONS BLAST_PARSER
                      BLAST_FILTER POST_BLAST_FILTER SCORE_CUTOFF BLASTDB 
                      MULTIMINIGENEWISE_OPTIONS MAX_INTRON)], @args);

  ####SETTING DEFAULTS#####
  $self->check_repeated(1);
  #########################
  
  logger_verbosity('INFO');
  $self->seqfetcher($seqfetcher);
  $self->ids($ids);
  $self->minigenewise_options($minigenewise_options);
  $self->full_seq($full_seq);
  $self->exonerate($exonerate);
  $self->exonerate_path($exonerate_path);
  $self->exonerate_options($exonerate_options);
  $self->check_repeated($check_repeated) if(defined($check_repeated));
  $self->genewise_options($genewise_options);
  $self->multiminigenewise_options($mmg_options);
  $self->blast_parser($blast_parser);
  $self->blast_filter($blast_filter);
  $self->post_blast_filter($post_filter);
  $self->score_cutoff($score_cutoff);
  $self->blastdb($blastdb);
  $self->max_intron($max_intron || 5000000);
  #note this is a default but due to the nature of it defaulting to on this
  #way is needed
  throw("No query sequence input")    unless ($self->query);
  throw("No seqfetcher provided")     unless ($self->seqfetcher);
  throw("No ids arrary ref provided") unless ($self->ids);
  throw("No analysis defined") unless($self->analysis);
  
  return $self;
}



sub ids{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw("BlastMiniGenewise::ids must be given an arrayref not a $arg")
      unless(ref($arg) eq "ARRAY");
    $self->{'ids'} = $arg;
  }
  return $self->{'ids'};
}


#Acessor methods



=head2 blast_filter

    Title   :   blast_filter
    Usage   :   $self->blast_filter($blast_filter)
    Function:   Get/set method for genewise blast_filter
    Returns :   
    Args    :   

=cut


sub blast_filter{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw($arg." must have method filter_results") 
      unless($arg->can("filter_results"));
    $self->{blast_filter} = $arg;
  }
  if(!$self->{blast_filter}){
    $self->{blast_filter} = Bio::EnsEMBL::Analysis::Tools::FeatureFilter
      ->new(
            -max_pvalue => 100,
            -prune => 1,
            -coverage => 100,
           );
  }
  return $self->{blast_filter};
}


=head2 minigenewise_options

    Title   :   minigenewise_options
    Usage   :   $self->minigenewise_options($minigenewise_options)
    Function:   Get/set method for minigenewise minigenewise_options
    Returns :   
    Args    :   

=cut


sub minigenewise_options{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw("BlastMiniGenewise::minigenewise_options must be given an hasref".
          "not a $arg") unless(ref($arg) eq "HASH");
    $self->{'minigenewise_options'} = $arg;
  }
  return $self->{'minigenewise_options'};
}

=head2 multiminigenewise_options

    Title   :   multiminigenewise_options
    Usage   :   $self->multiminigenewise_options($multiminigenewise_options)
    Function:   Get/set method for multiminigenewise multiminigenewise_options
    Returns :   
    Args    :   

=cut


sub multiminigenewise_options{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw("BlastMiniGenewise::multiminigenewise_options must be given an hasref".
          "not a $arg") unless(ref($arg) eq "HASH");
    $self->{'multiminigenewise_options'} = $arg;
  }
  return $self->{'multiminigenewise_options'};
}

=head2 genewise_options

    Title   :   genewise_options
    Usage   :   $self->genewise_options($genewise_options)
    Function:   Get/set method for genewise genewise_options
    Returns :   
    Args    :   

=cut


sub genewise_options{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw("BlastMiniGenewise::genewise_options must be given an hasref".
          "not a $arg") unless(ref($arg) eq "HASH");
    $self->{'genewise_options'} = $arg;
  }
  return $self->{'genewise_options'};
}


=head2 score_cutoff

    Title   :   score_cutoff
    Usage   :   $self->score_cutoff($score_cutoff)
    Function:   Get/set method for genewise score_cutoff
    Returns :   
    Args    :   

=cut

sub score_cutoff{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{score_cutoff} = $arg;
  }
  return $self->{score_cutoff};
}

sub post_blast_filter{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{post_blast_filter} = $arg;
  }
  return $self->{post_blast_filter};
}


=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::EnsEMBL::Analysis::Tools::SeqFetcher object
    Args    :   Bio::EnsEMBL::Analysis::Tools::SeqFetcher object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    
  if ($value) {
    throw("BlastMiniGenewise, seqfetcher must be a ".
          "Bio::EnsEMBL::Analysis::Tools::SeqFetcher")
      unless ($value->isa("Bio::DB::RandomAccessI"));
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}

=head2 exonerate

    Title   :   exonerate
    Usage   :   $self->exonerate
    Function:   Get/Set method for using exonerate rather than Blast
    Returns :   0 (False) or 1 (True)
    Args    :   0 (False) or 1 (True)

=cut

sub exonerate {
  my( $self, $value ) = @_;    

  if ($value) {
    $self->{'_exonerate'} = $value;
  }

  return $self->{'_exonerate'};
}

=head2 exonerate_path

    Title   :   exonerate_path
    Usage   :   $path = $self->exonerate_path
    Function:   Get/Set method for the path to the Exonerate executeable
    Returns :   String
    Args    :   String

=cut

sub exonerate_path {
  my( $self, $value ) = @_;    
  if($value){
    my $path = $self->locate_executable($value);
    $self->{'exonerate_path'} = $path;
  }
  throw($self->{'exonerate_path'}." is not executable") 
    if($self->{'exonerate_path'} && !(-x $self->{'exonerate_path'}));
  return $self->{'exonerate_path'};
}

=head2 exonerate_options

    Title   :   exonerate_options
    Usage   :   $options = $self->exonerate_options
    Function:   Get/Set method for the options for running Exonerate
    Returns :   String
    Args    :   String

=cut

sub exonerate_options {
  my( $self, $value ) = @_;    

  if ($value) {
    $self->{'_exonerate_options'} = $value;
  }

  return $self->{'_exonerate_options'};
}

=head2 full_seq

    Title   :   full_seq
    Usage   :   $self->full_seq
    Function:   Get/Set method for using the full genomic
            :   sequence rather than the mini seq
    Returns :   0 (False) or 1 (True)
    Args    :   0 (False) or 1 (True)

=cut

sub full_seq {
  my( $self, $value ) = @_;    

  if ($value) {
    $self->{'_full_seq'} = $value;
  }

  return $self->{'_full_seq'};
}


sub blast_parser{
  my ($self, $value) = @_;
  if($value){
    throw($value." must have a method called parse_file") 
      unless($value->can("parse_file"));
    $self->{blast_parser} = $value;
  }
  if(!$self->{blast_parser}){
    my $regex;
    if($self->query->name =~ /^\S+\:\S*\:(\S+)\:\S*:\S*\:\S*/){
      $regex = '^\S+\:\S*\:(\S+)\:\S*:\S*\:\S*';
    }elsif ($self->query->name =~ /^(.*)\|(.*)\|(.*)/) {
      $regex = '^.*\|(.*)\|.*';
    } elsif ($self->query->name =~ /^..\:(.*)/) {
      $regex = '^..\:(.*)';
    }else {
      $regex = '^(\w+)';
    }
    $self->{blast_parser}  = Bio::EnsEMBL::Analysis::Tools::FilterBPlite
      ->new(
            -regex => $regex,
            -query_type => "pep",
            -database_type => "dna",
           );
  }
  return $self->{blast_parser};
}


sub blastdb{
  my ($self, $arg) = @_;
  if($arg){
    throw($arg." must be a blastdb object") 
      unless($arg->isa("Bio::EnsEMBL::Analysis::Tools::BlastDB"));
    $self->{blastdb} = $arg;
  }
  if(!$self->{blastdb}){
    my $blastdb = Bio::EnsEMBL::Analysis::Tools::BlastDB
      ->new(
            -sequences => [$self->query],
            -output_dir => $self->workdir,
            -mol_type => "DNA",
           );
    $self->{blastdb} = $blastdb;
  }
  return $self->{blastdb};
}

sub blastdb_file{
  my ($self, $arg) = @_;
  if($arg){
    throw($arg." must exist as a file") unless(-e $arg);
    $self->{blastdb_file} = $arg;
  }
  if(!$self->{blastdb_file}){
    my $blastdb_file = $self->blastdb->create_blastdb;
    foreach my $file(@{$self->blastdb->list_dbfiles}){
      $self->files_to_delete($file);
    }
    $self->{blastdb_file} = $blastdb_file;
  }
  return $self->{blastdb_file};
}

#functional methods


=head2 run
  Title   : run.
  Usage   : $self->run()
  Function: Performs a Blast search or an Exonerate search 
          : with all the protiens, passes the resulting proteins
          : into a Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise
          : runnable
  Returns : none
  Args    : none

=cut

sub run {
  my ($self) = @_;

  my $features = $self->run_alignment;

  unless (@$features) {
    print $self->query->name." has no associated features. ".
           "Finishing run.\n";
    return;
  }

  my $mg_runnables;

  my @sorted_features = sort{$a->start <=> $b->start  
                     || $a->end <=> $b->end} @$features; 
  my %hash;
  foreach my $f(@sorted_features){
    if(!$hash{$f->hseqname}){
      $hash{$f->hseqname} = 1;
    }else{
      $hash{$f->hseqname}++;
    }
  }
  
  if ($self->check_repeated){ 
    $mg_runnables = $self->build_runnables(@sorted_features);
  } else {
    my $runnable = $self->make_object($self->query, \@sorted_features);
    push (@$mg_runnables, $runnable); 
  }
  my $output_count;
  foreach my $mg (@$mg_runnables){
    $mg->run;
    my @f = @{$mg->output};
    $self->output(\@f);
    $output_count += @f;
  }

  $self->delete_files;
  return 1;
  
}

sub run_alignment{
  my ($self) = @_;
  my $seqs = $self->get_Sequences;
  if(@$seqs != @{$self->ids}){
    warning("Managed to get only ".scalar(@$seqs)." sequences from ".scalar(@{$self->ids}).
            " for alignment run, check indices\n");
  }

  my $valid_seqs = $self->validate_Sequences($seqs);
  if(@$valid_seqs != @$seqs){
    warning("Only ".scalar(@$valid_seqs)." sequences were validated from ".@$seqs.
            " sequences");
  }
  logger_info("There are ".@$valid_seqs." sequences to align using blast/exonerate");
  my @features; 
  my @sorted_seqs = sort {$a->id cmp $b->id} @$valid_seqs;
  foreach my $seq(@sorted_seqs){
    if($self->exonerate){
      push(@features, @{$self->run_exonerate($seq)});
    }else{
      #push(@features, @{$self->run_blast($seq)});
      foreach my $f(@{$self->run_blast($seq)}){
        $f->invert;
        $f->slice($self->query);
        push(@features, $f);
      }
    }
  }
  return \@features;
}

sub run_blast {
  my ($self, $seq) = @_;
  my $blastdb_file = $self->blastdb_file;
  my $run = Bio::EnsEMBL::Analysis::Runnable::Blast
    ->new(
          -query => $seq,
          -program => 'tblastn',
          -database => $blastdb_file,
          -parser => $self->blast_parser,
          -filter => $self->blast_filter,
          -analysis => $self->analysis,
         );
  $run->run;
  my @blast_features;
  foreach my $f (@{$run->output}){
    my $id = $seq->id;
    if($f->seqname =~ /$id/){
      $f->seqname($id);
    }else{
      throw($f->seqname." doesn't contain ".$id);
    }
    $f->slice($self->query);
    push(@blast_features,$f);
  }
  my @filter_features;
  if($self->post_blast_filter){
    #print STDERR "Filtering blast results\n";
    #this code will through out any sets of features where
    #none of the blast scores are higher than score set in config 
    #on the grounds its unlikely in that case that it will produce 
    #a good gene if a gene at all
    my @fs;
    my %feature_hash;
    while(my $f = shift(@blast_features)){
      if(!$feature_hash{$f->seqname}){
	$feature_hash{$f->seqname} = [];
	push(@{$feature_hash{$f->seqname}}, $f);
      }else{
	push(@{$feature_hash{$f->seqname}}, $f);
      }
    }
  HIT: foreach my $hid(keys(%feature_hash)){
      my @hit_features = @{$feature_hash{$hid}};
      foreach my $f (@hit_features){
        if($f->score > $self->score_cutoff){
          push(@filter_features, @hit_features);
          next HIT;
        }
      }
    }
  }else{
    @filter_features = @blast_features;
  }
  return \@filter_features;
}


=head2 run_exonerate

  Title   : run_exonerate
  Usage   : @features = $self->run_exonerate()
  Function: Uses Bio::Ensembl::Analysis::ExonerateTranscript to align
          : the proteins to the slice. The supoporting features of 
          : the transcripts are returned
  Returns : none
  Args    : list of Bio::EnsEMBL::BaseAlignFeatue objects

=cut

sub run_exonerate {
  my ($self, $seq)= @_;

  my @features;
#  my @seqs = @{$self->get_Sequences};
#  if (@seqs != @{$self->ids}) {
#    warning("Managed to get only " . scalar(@seqs) . "  of ".
#            scalar(@{$self->ids}) ."for Exonerate run; check indices\n");
#  }
#  my @valid_seqs   = $self->validate_sequence(@seqs);
#  my @sorted_seqs = sort {$a->id cmp $b->id} @valid_seqs;
#  foreach my $seq (@sorted_seqs) {
    my $exonerate = new  Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
      (
       -program     => $self->exonerate_path,
       -analysis    => $self->analysis,
       -query_seqs  => ([$seq]),
       -query_type  => 'protein',
       -target_seqs => ([$self->query]),
       -target_type => 'dna',
       -options     => "-W 7 ". $self->exonerate_options,
      );
    eval {
      $exonerate->run;
    };
    if ($@){
      throw("Exonerate died on me$@\n");
    }
     my $transcripts = $exonerate->output;
     unless (scalar(@$transcripts > 0)){
      print "Didn't get a good exonerate alignment, trying again with a shorter word length\n";
      $exonerate = new  Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
	(
       -program     => $self->exonerate_path,
       -analysis    => $self->analysis,
       -query_seqs  => ([$seq]),
       -query_type  => 'protein',
       -target_seqs => ([$self->query]),
       -target_type => 'dna',
       -options     => "-W 5 ". $self->exonerate_options,
	);
      eval {
	$exonerate->run;
      };
      if ($@){
	throw("Exonerate died on me$@\n");
      }
      $transcripts = $exonerate->output;
    }
     unless (scalar(@$transcripts > 0)){
       print STDERR "Exonerate cannot align  ".$seq->display_id." \n";
      next;
    }
  
    foreach my $trans (@{$exonerate->output}){
      foreach my $exon (@{$trans->get_all_Exons}){
	push @features, @{$exon->get_all_supporting_features};
      }
    }
#  }
  return \@features;
}

#=head2 make_mmgw
=head2 make_object

  Args [1]   : $miniseq - a Bio::Seq object representing the
               target sequence for the genewise run.
  Args [2]   : $features - reference to a list of 
               FeaturePairs generated by a blast run.
  Example    : $self->make_object($miniseq, $features);
  Description: Takes a genomic sequence and builds a
               MultiMiniGenewise runnable object using the 
               list of FeaturePairs.
  Returntype : A list of 
               Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise
  Exceptions : none
  Caller     : $self->build_runnables

=cut


#sub make_mmgw {
sub make_object {
  my ($self, $miniseq, $features, $cluster_start, $cluster_end) = @_;
  
#  # Before we pass our blast generated features to 
#  # MultiMiniGenewise we must first convert them from 
#  # PepDnaAlignFeatures to FeaturePairs.
#
#  my @newf;
#  foreach my $f (@$features){
#    my $newf = new Bio::EnsEMBL::FeaturePair(-feature1 => $f->feature2,
#					     -feature2 => $f->feature1);
#    push(@newf,$newf);
#  }

  # Create a MultiMiniGenewise object with the features we've
  # just converted.
  
  my %pars = %{$self->multiminigenewise_options} if($self->multiminigenewise_options);
  if (defined($cluster_end)) {
    $pars{-cluster_start} = $cluster_start;
    $pars{-cluster_end}   = $cluster_end;
  }
  my $mg = Bio::EnsEMBL::Analysis::Runnable::MultiMiniGenewise->new
    (     
     -query            => $miniseq,
     -features         => $features,
     -seqfetcher       => $self->seqfetcher,
     -minigenewise_options => $self->minigenewise_options,
     -fullseq          => $self->full_seq,
     -genewise_options => $self->genewise_options,
     -analysis         => $self->analysis,
     %pars,
    );
  %pars = ();
  return $mg;
}

sub get_Sequences {
    my ($self) = @_;

    my @seq;

    foreach my $id (@{$self->ids}) {
      logger_info("Fetching ".$id." sequence\n");
        my $seq = $self->get_Sequence($id);

        if ($seq && $seq->length > 0) {
            push(@seq,$seq);
        } else {
            print "Invalid sequence for $id - skipping\n";
        }
    }

    return \@seq;

}

sub validate_Sequences {
  my ($self, $seqs) = @_;
  my @validated;
  
  foreach my $seq (@$seqs) {
    
    my $sequence = $seq->seq;
    if ($sequence !~ /[^acgtn]/i) {
      push (@validated, $seq);
    } else {
      $_ = $sequence;
      my $len = length ($_);
      my $invalidCharCount = tr/bB/xX/;
      
      if ($invalidCharCount / $len > 0.05) {
        warning("Ignoring ".$seq->display_id()
                ." contains more than 5% ($invalidCharCount) "
                ."odd nucleotide codes ($sequence)\n Type returns "
                .$seq->moltype().")\n");
      } else {
        $seq->seq($_);
        push (@validated, $seq);
      }
    }
  } 
  return \@validated;  
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

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
  my ($p, $f, $l) = caller;
  if (!defined($id)) {
    warning("No id input to get_Sequence");
  }  
  
  eval {
    $seq = $seqfetcher->get_Seq_by_acc($id);
  };
  
  if($@) {
    warning("Problem fetching sequence for id [$id] with $seqfetcher $name  $@\n");
    return undef;
  }
  
  if(!defined($seq)){
    warning("Could not find sequence for [$id] with $seqfetcher $name called by $f:$l");
  }
  
  return $seq;
}


1;
