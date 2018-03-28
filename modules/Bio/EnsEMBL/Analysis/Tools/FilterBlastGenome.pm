# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome->
  new(
      -regex => '^\w+\s+(\w+)'
      -query_type => 'pep',
      -database_type => 'dna',
      -threshold_type => 'PVALUE',
      -threshold => 0.01,
     );
 my @results = @{$parser->parse_results('blast.out')};

=head1 DESCRIPTION

This is an attempt at making a blast parser/filter that uses the protein sequence as the
query and the genome as the db. Typically the blast modules expect the db to be a protein
db. There are some changes in this module to account for the switch

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::FeatureFilter;
use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Tools::BPliteWrapper);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  Arg [THRESHOLD_TYPE] : string, threshold type
  Arg [THRESHOLD] : int, threshold
  Arg [COVERAGE] : int, coverage value
  Arg [FILTER] : int, boolean toggle as whether to filter
  Function  : create a Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  object
  Returntype: Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  Exceptions: 
  Example   : 

=cut


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($threshold_type, $threshold, $coverage, $filter) = rearrange
    (['THRESHOLD_TYPE', 'THRESHOLD', 'COVERAGE', 'FILTER'], @args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->coverage(10);
  $self->filter(1);
  ######################

  $self->threshold_type($threshold_type);
  $self->threshold($threshold);
  $self->coverage($coverage) if(defined($coverage));
  $self->filter($filter) if(defined($filter));
#  my $slice;
#  my $analysis;
#  $self->slice($slice) if(defined($slice));
#  $self->analysis($analysis) if(defined($analysis));
  return $self;
}


=head2 threshold_type 

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  Arg [2]   : string/int
  Function  : container methods, this documents the 4 methods
  below threshold_type, threshold, coverage, filter
  Returntype: string/int
  Exceptions: 
  Example   : 

=cut


sub threshold_type{
  my $self = shift;
  $self->{'threshold_type'} = shift if(@_);
  return $self->{'threshold_type'};
}

sub threshold{
  my $self = shift;
  $self->{'threshold'} = shift if(@_);
  return $self->{'threshold'};
}

sub coverage{
  my $self = shift;
  $self->{'coverage'} = shift if(@_);
  return $self->{'coverage'};
}

sub filter{
  my $self = shift;
  $self->{'filter'} = shift if(@_);
  return $self->{'filter'};
}



=head2 get_hsps

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::BPlite
  Function  : prefilter the hsps then parser then and turn them into
  features
  Returntype: none 
  Exceptions: throw if no name can be parser from the subject
  Example   : 

=cut



sub get_hsps {
  my ($self, $parsers) = @_;
  my $regex = $self->regex;
  my @output;
  my $ids;
  if($self->filter){
    $ids = $self->filter_hits($parsers);
  }
  my $seconds = $self->get_parsers($self->filenames);
  PARSER:foreach my $second (@$seconds) {
  NAME: while(my $sbjct = $second->nextSbjct) {

      # First parse the slice name out of the hit name, since the subject is the genome itself
      my $slice_regex = '\S+\:\S+\:\S+\:\d+\:\d+\:1/';
      my $slice_name = $sbjct->name;
      unless($slice_name =~ /\S+\:\S+\:\S+\:\d+\:\d+\:1/) {
        $self->throw("Could not find a slice name in the hit name.\nRegex used:\n".$slice_regex."\nLine searched:\n".$slice_name);
      }

      $slice_name = $&;
      $sbjct->slice_name($slice_name);

      # Retrieve the query name of the protein and then swap it into the the name field
      my $query_name = $sbjct->query_name;
      $sbjct->name($query_name);

    HSP: while (my $hsp = $sbjct->nextHSP) {
        if($self->is_hsp_valid($hsp)) {
          say "FM2 SPLIT: ".Dumper($self->split_hsp($hsp, $query_name));
          my $fp = $self->split_hsp($hsp, $query_name);
          $fp->{'SLICE_NAME'} = $slice_name;
          push(@output, $fp);
        }
      }
    }
  }
  $parsers = [];
  $self->output(\@output);
}



=head2 filter_hits

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::BPlite
  Function  : prefilter the blast results using specified thresholds
  and FeatureFilter
  Returntype: hashref
  Exceptions:
  Example   :

=cut



sub filter_hits {
  my ($self, $parsers) = @_;
  my %ids;
  my @features;
 PARSER:foreach my $parser(@$parsers){
  SUB:while(my $sbjct = $parser->nextSbjct) {
    my $name = $sbjct->name;
    HSP:while (my $hsp = $sbjct->nextHSP) {
        if($self->is_hsp_valid($hsp)) {
#          my $qstart = $hsp->query->start();
#          my $hstart = $hsp->subject->start();
#          my $qend   = $hsp->query->end();
#          my $hend   = $hsp->subject->end();
#          my $qstrand = $hsp->query->strand();
#          my $hstrand = $hsp->subject->strand();
#          my $score  = $hsp->score;
#          my $p_value = $hsp->P;
#          my $percent = $hsp->percent;

          my $qstart = $hsp->subject->start();
          my $hstart = $hsp->query->start();
          my $qend   = $hsp->subject->end();
          my $hend   = $hsp->query->end();
          my $qstrand = $hsp->subject->strand();
          my $hstrand = $hsp->query->strand();
          my $score  = $hsp->score;
          my $p_value = $hsp->P;
          my $percent = $hsp->percent;

          my $fp = $self->feature_factory->create_feature_pair
            ($qstart, $qend, $qstrand, $score, $hstart,
             $hend, $hstrand, $name, $percent, $p_value);

#          say "FM2 Dumper FP:".Dumper($fp);
          push(@features,$fp);
        }
      }
    }
  }

  my $search = Bio::EnsEMBL::Analysis::Tools::FeatureFilter->new(-coverage => $self->coverage);

  my @newfeatures = @{$search->filter_results(\@features)};

  foreach my $f (@newfeatures) {
    my $id = $f->hseqname;
    $ids{$id} = 1;
  }
  return \%ids;
}



=head2 is_hsp_valid

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FilterBPlite
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::BPlite::HSP
  Function  : checks hsp against specified threshold returns hsp
  if above value 0 if not
  Returntype: Bio::EnsEMBL::Analysis::Tools::BPlite::HSP/0
  Exceptions:
  Example   :

=cut



sub is_hsp_valid{
  my ($self, $hsp) = @_;
  if($self->threshold_type){
    if ($self->threshold_type eq "PID") {
      return 0 if ($hsp->percent < $self->threshold);
    } elsif ($self->threshold_type eq "SCORE") {
      return 0 if ($hsp->score < $self->threshold);
    } elsif ($self->threshold_type eq "PVALUE") {
      return 0 if($hsp->P > $self->threshold);
    } 
  }
  return $hsp;
}



sub split_hsp {
    my ($self,$hsp,$name) = @_;
    my $qstrand = $hsp->subject->strand;
    my $hstrand = $hsp->query->strand;
    my ($qinc,   $hinc)    = $self->find_increments($qstrand,$hstrand);
    my @qchars = split(//,$hsp->querySeq);  # split alignment into array of
                                            # chars
    my @hchars = split(//,$hsp->sbjctSeq);  # ditto for hit sequence
    my $qstart = $hsp->subject->start(); # Start off the feature pair start
    my $hstart = $hsp->query->start(); # ditto
    my $qend   = $hsp->subject->start(); # Set the feature pair end also
    my $hend   = $hsp->query->start(); # ditto
    if ($qstrand == -1) {
      $qstart = $hsp->subject->end;
      $qend   = $hsp->subject->end;
    }
    if ($hstrand == -1) {
      $hstart = $hsp->query->end;
      $hend   = $hsp->query->end;
    }

    my $count = 0; # counter for the bases in the alignment
    my $found = 0; # flag saying whether we have a feature pair


    my @tmpf;

    while ($count <= $#qchars) {
      # We have hit an ungapped region.  Increase the query and hit
      #counters and flag that we have a feature pair.

      if ($qchars[$count] ne '-' &&
          $hchars[$count] ne '-') {

        $qend += $qinc;
        $hend += $hinc;

        $found = 1;
      } else {

        # We have hit a gapped region.  If the feature pair flag is set
        # ($found) then make a feature pair, store it and reset the start
        # and end variables.

        my $query_seqname;
        if ($hsp->query->can('seq_id')) {
          $query_seqname = $hsp->query->seq_id;
        } else {
          $query_seqname = $hsp->query->seqname;
        }

        if ($found == 1) {
          my $fp = $self->convert_to_featurepair($qstart, $qend, $qstrand,
                                                 $qinc, $hstart, $hend,
                                                 $hstrand, $hinc, $name,
						  $query_seqname,
                                                 $hsp->score,
                                                 $hsp->percent, $hsp->P,
                                                 $hsp->positive,
                                                 $hsp->match);
          push(@tmpf,$fp);
        }

        # We're in a gapped region.  We need to increment the sequence that
        # doesn't have the gap in it to keep the coordinates correct.
        # We also need to reset the current end coordinates.

        if ($qchars[$count] ne '-') {
          $qstart = $qend   + $qinc;
        } else {
          $qstart = $qend;
        }
        if ($hchars[$count] ne '-') {
          $hstart = $hend   + $hinc;
        } else {
          $hstart = $hend;
        }

        $qend = $qstart;
        $hend = $hstart;

        $found = 0;
      }
      $count++;
    }
    # Remember the last feature
    if ($found == 1) {
      my $query_seqname;
      if ($hsp->query->can('seq_id')) {
        $query_seqname = $hsp->query->seq_id;
      } else {
        $query_seqname = $hsp->query->seqname;
      }

      my $fp = $self->convert_to_featurepair($qstart, $qend, $qstrand,
                                             $qinc, $hstart, $hend,
                                             $hstrand, $hinc, $name,
					          $query_seqname,
                                             $hsp->score,
                                             $hsp->percent, $hsp->P,
                                             $hsp->positive,
                                             $hsp->match);
      push(@tmpf,$fp);
    }
    my $fp;


    $qinc = abs( $qinc );
    $hinc = abs( $hinc );

    if( $qinc == 3 && $hinc == 1 ) {
      $fp = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tmpf, -align_type => 'ensembl');
    } elsif( $qinc == 1 && $hinc == 3 ) {
      $fp = Bio::EnsEMBL::PepDnaAlignFeature->new(-features => \@tmpf, -align_type => 'ensembl');
    } elsif( $qinc == 1 && $hinc == 1 ) {
      $fp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tmpf, -align_type => 'ensembl');
    } else {
      throw( "Hardcoded values wrong?? " );
    }

    # helps debugging subsequent steps
    $fp->{'qseq'} = $hsp->querySeq();
    $fp->{'sseq'} = $hsp->sbjctSeq();

    # for compara
    $fp->positive_matches($hsp->positive);
    $fp->identical_matches($hsp->match);
    if($fp->hstart > $fp->hend){
      throw("Failed start ".$fp->hstart." is greater than end ".$fp->hend." ".
            "for ".$fp->hseqname."\n");
    }
    return $fp;
  }

