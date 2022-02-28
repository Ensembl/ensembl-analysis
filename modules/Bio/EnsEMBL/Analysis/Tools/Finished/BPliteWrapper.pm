# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

  Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper

=head1 SYNOPSIS

no synopsis

=head1 DESCRIPTION

This module is a wrapper for Tools::BPliteWrapper to provide an interface between
Bio::EnsEMBL::Analysis::Runnable::Finished::Blast and BPlite. This method fits model
for standard blast parsers as it provides the parse_file method which
returns an array of results. This method uses BPliteWrapper to parse the
file and contains methods to do filtering which is called by the Runnable::Finished::Blast
It contains the actual methods to create FeaturePairs from HSPs after
doing any depth filtering to save time and memory when searching genomic sequences that generate
large numbers of blast matches.

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::PepDnaAlignFeature;
use base 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper';

=head2 new

  Arg [1]               : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [FILENAME]        : string, filename
  Arg [REGEX]           : string, regex
  Arg [QUERY_TYPE]      : string, query sequence type, should be pep or dna
  Arg [DATABASE_TYPE]   : string, database sequence type as above
  Arg [ANALYSIS]        : Bio::EnsEMBL::Analysis object
  Arg [COVERAGE]        : integer, coverage , defined in config file
  Arg [THRESHOLD_TYPE]  : string, threshold_type defined in config file
  Arg [THRESHOLD]       : integer, threshold defined in config file
  Arg [DISCARD_OVERLAPS]: flag, 1 or 0 defined in config file
  Function  :
  Returntype:
  Exceptions:
  Example   :

=cut


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($regex, $skip_name, $query, $target, $analysis,$coverage,$threshold_type,$threshold,$discard_overlaps) =
    rearrange(['REGEX', 'SKIP_NAME', 'QUERY_TYPE', 'DATABASE_TYPE',
               'ANALYSIS','COVERAGE','THRESHOLD_TYPE','THRESHOLD','DISCARD_OVERLAPS'], @args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->regex('(\w+)\s+');
  ######################
  $self->regex($regex) if(defined $regex);
  $self->skip_name($skip_name) if(defined $skip_name);
  $self->query_type($query) if($query);
  $self->database_type($target) if($target);
  $self->analysis($analysis);
  $self->threshold_type($threshold_type);
  $self->threshold($threshold);
  $self->coverage($coverage);
  $self->discard_overlaps($discard_overlaps);
  return $self;
}

=head2 skip_name

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : string, regex
  Function  : skip entry name that matches regex
  Returntype: string, regex
  Exceptions:
  Example   :

=cut


sub skip_name{
  my $self = shift;
  $self->{'skip_name'} = shift if(@_);
  return $self->{'skip_name'};
}


=head2 threshold_type

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : string, threshold type
  Function  : container
  Returntype: string
  Exceptions: throws if string is not either PVALUE or PID
  Example   :

=cut

sub threshold_type{

  my ($self, $ttype) = @_;
  if ( $ttype ) {
      throw ("Threshold_type should be either PVALUE or PID") unless ($ttype eq 'PVALUE' || $ttype eq 'PID');
      $self->{'threshold_type'} = $ttype;
  }
  return $self->{'threshold_type'};
}

=head2 threshold

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : integer, threshold
  Function  : container
  Returntype: integer
  Example   :

=cut

sub threshold{
  my ($self, $tvalue) = @_;
  if ( $tvalue ) {
      $self->{'threshold'} = $tvalue;
  }
  return $self->{'threshold'};
}

=head2 coverage

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : integer, coverage
  Function  : container
  Returntype: integer
  Exceptions: throws if string is not >= 0 or <= 255
  Example   :

=cut

sub coverage{
  my ($self, $coverage) = @_;
  if ( $coverage ) {
      throw ("Coverage set should be between 0 and 255") unless ($coverage >= 0 || $coverage <= 255 );
      $self->{'coverage'} = $coverage;
  }
  return $self->{'coverage'};
}

=head2 discard_overlaps

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : string, discard_overlap
  Function  : container
  Returntype: string
  Example   :

=cut

sub discard_overlaps{
  my ($self, $dtype) = @_;
  if ( $dtype ) {
      $self->{'discard_overlaps'} = $dtype;
  }
  return $self->{'discard_overlaps'};
}


=head2 parse_files

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : string filename
  Function  : using BPlite to parse the blast output
  Returntype: arrayref
  Exceptions: throws if file does not exist
  Example   :

=cut



sub parse_files{

  my ($self, $files) = @_;
  $self->clean_output;
  $self->clean_filenames;
  $self->filenames($files);
  my $bplites = $self->get_parsers($files);
  return $bplites;

}



=head2 get_parsers

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : string, filename
  Function  : opens file using Filehandle and passes filehandle to BPlite
  Returntype: Bio::EnsEMBL::Analysis::Tools::BPlite
  Exceptions: none
  Example   :

=cut



sub get_parsers {
  my ($self, $files)  = @_;
  if(!$files){
    $files = $self->filenames;
  }
  my @parsers;
  foreach my $file (@$files) {
    my $fh = new FileHandle;
    $fh->open("<".$file);
    my $parser = Bio::EnsEMBL::Analysis::Tools::BPlite->new('-fh' => $fh);
    push(@parsers,$parser);
  }

  return \@parsers;
}


=head2 get_hsps

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::BPlite
  Function  : get the hsps from bplite and return the best hits
  Returntype: none
  Exceptions: throws if regex does not produce a result
  Example   :

=cut


sub get_best_hits{

  my ($self, $parsers,$thresh_type,$threshold) = @_;
  my $regex = $self->regex;
  my $skip_name = $self->skip_name;
  my $best_hits = {};

  PARSER:
  foreach my $parser(@$parsers){

      HIT:
      while(my $sbjct = $parser->nextSbjct){

	  my ($name) = $sbjct->name =~ /$regex/;
	  throw("Error parsing name from ".$sbjct->name." check your ".
            "blast setup and blast headers") unless($name);
	  # ignore blast db entries like...
 	  next HIT if $skip_name && $name =~ /$skip_name/;

	  my $hsps=[];
	  my $above_thresh=0;
	  my $best_value=1e6;
	  my $top_score = 0;

	  HSP:
	  while (my $hsp = $sbjct->nextHSP) {

	      push(@$hsps, $hsp);
	      my $score=$hsp->score;
	      my $pc=$hsp->percent;
	      my $p =$hsp->P;

	      if ( $thresh_type eq 'PVALUE' ) {
		if ( $p <= $threshold ) {
		  $above_thresh = 1;
		}
		if ( $p < $best_value ) {
		  $best_value = $p;
		}
	      }
	      elsif ( $thresh_type eq 'PID' ) {

		$best_value = 0;
		if ( $pc >= $threshold ) {
		  $above_thresh = 1;
		}
		if ( $pc > $best_value ) {
		  $best_value = $pc;
		}
	      }
	      else {
		  throw("Unknown threshold type '$thresh_type'");
	      }

	      #get top_score

	      if ( $score > $top_score ) {
		  $top_score = $score;
	      }

	  }

	  next unless $above_thresh;
          my $best = $best_hits->{$best_value}{$top_score} ||= [];

          #Put as first element of hsps array
          unshift ( @$hsps, $name );

          push ( @$best, $hsps );

	}#while next->sbjt

    }#foreach parser

  $parsers = [];
  return $best_hits;

}

sub _apply_coverage_filter {

    my ( $self, $query_length, $best_hits,$threshold_type,$threshold,$coverage,$discard_overlaps ) = @_;
    my %id_list;
    my @output;
    my $max_coverage = defined($coverage) ? $coverage : 0;
    $discard_overlaps = $discard_overlaps || 0;
    throw("Max coverage '$max_coverage' is beyond limit of method '255'") if $max_coverage > 255;

    # Make a string of nulls (zeroes) the length of the query
    my $coverage_map = "\0" x ( $query_length + 1 );
    my $ana = $self->analysis;

    # Loop through from best to worst according
    # to our threshold type.

    my (@bin_numbers);

    if ( $threshold_type eq 'PID' ) {
	@bin_numbers = sort { $b <=> $a } keys %$best_hits;
    }
    elsif ( $threshold_type eq 'PVALUE' ) {
	@bin_numbers = sort { $a <=> $b } keys %$best_hits;
    }
    else {
	throw("Unknown threshold type '$threshold_type'");
    }

    foreach my $bin_n (@bin_numbers) {
	my $score_hits = $best_hits->{$bin_n};

	# Loop through from best to worst according
	# to score within threshold bin
	foreach my $top_score ( sort { $b <=> $a } keys %$score_hits ) {
	    my $hits = $score_hits->{$top_score};

	    # Loop through each hit with this score
	    foreach my $hit (@$hits) {
      		my $keep_hit = 0;
      		my ( $name, @hsps ) = @$hit;
		my $ct=scalar(@hsps);

		# Don't keep multiple matches to the same sequence
		# at the same genomic location.
		@hsps = $self->_discard_worst_overlapping(@hsps) if $discard_overlaps;
		unless ($max_coverage == 0){
		    foreach my $hsp (@hsps) {
			my $q = $hsp->query;
			foreach my $i ( $q->start .. $q->end ) {
			    my $depth = unpack( 'C', substr( $coverage_map, $i, 1 ) );
			    if ( $depth < $max_coverage ) {
				$keep_hit = 1;
				$depth++;

				# Increment depth at this position in the map
				substr( $coverage_map, $i, 1 ) = pack( 'C', $depth );
			    }
			}
		    }
		}

		# Make FeaturePairs if we want to keep this hit
		if ($keep_hit || $max_coverage == 0) {
		    foreach my $hsp (@hsps) {
			push (@output,$self->split_hsp( $hsp, $name, $ana ));
		    }
		}

	    }
	}
    }
    print "*** FINISHED PARSING ***\n";

    $self->output(\@output);
}




sub _discard_worst_overlapping {

    my $self = shift;
    my @hsps = sort { $a->query->start <=> $b->query->start } @_;

    # Put all the hsps hits into overlapping bins
    my $first   = shift @hsps;
    my $start   = $first->start;
    my $end     = $first->end;
    my $current = [$first];        # The current bin
    my @bins    = ($current);
    while ( my $hsp = shift @hsps ) {
        my $q = $hsp->query;

        # Does this hsp overlap the current bin?
	# compare the first hsp with the second hsp, if both overlap then the current array will have both the first hsp and the second hsp
        if ( $q->start <= $end and $q->end >= $start ) {
            push ( @$current, $hsp );
            if ( $q->end > $end ) {
                $end = $q->end;
            }
        }
        else {
            $current = [$hsp];
            push ( @bins, $current );
            $start = $hsp->start;
            $end   = $hsp->end;
        }
    }

    foreach my $bin (@bins) {
        if ( @$bin == 1 ) {
            push ( @hsps, $bin->[0] );
        }
        else {

            # Remove the hsp with the lowest percent identity
            # (or lowest score or highest P value)
            my @bin_hsps = sort { $a->percent <=> $b->percent || $a->score <=> $b->score || $b->P <=> $a->P } @$bin;
            shift (@bin_hsps);

            if ( @bin_hsps == 1 ) {

                # We are left with 1 hsp, so add it.
                push ( @hsps, $bin_hsps[0] );
            }
            else {

                # Remaining hsps may not be overlapping
                push ( @hsps, $self->_discard_worst_overlapping(@bin_hsps) );
            }
        }
    }
    return @hsps;
}



1;
