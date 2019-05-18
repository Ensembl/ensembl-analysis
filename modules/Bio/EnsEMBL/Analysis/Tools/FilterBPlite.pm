# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

  Bio::EnsEMBL::Analysis::Tools::FilterBPlite

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Analysis::Tools::FilterBPlite->
  new(
      -regex => '^\w+\s+(\w+)'
      -query_type => 'dna',
      -database_type => 'pep',
      -threshold_type => 'PVALUE',
      -threshold => 0.01,
     );
 my @results = @{$parser->parse_results('blast.out')};

=head1 DESCRIPTION

This module inherits from BPliteWrapper so follows the same basic 
methodology but it implements some prefiltering of the HSPs to mimic how
the old pipeline blast runnable was used in the raw computes

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::FilterBPlite;

use strict;
use warnings;

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



sub get_hsps{
  my ($self, $parsers) = @_;
  my $regex = $self->regex;
  my @output;
  my $ids;
  if($self->filter){
    $ids = $self->filter_hits($parsers);
  }
  my $seconds = $self->get_parsers($self->filenames);
 PARSER:foreach my $second(@$seconds){
  NAME: while(my $sbjct = $second->nextSbjct){
      if($self->filter && !($ids->{$sbjct->name})){
        next NAME;
      }

      # Input IDs which are longer than 78 characters would be printed
      # over two rows (80 chars per row) in the BLAST temporary output 
      # (to be parsed into Ensembl DB). (78 chars is the threshold, 
      # allowing for space for "> ".) Therefore, a single whitespace 
      # would be introduced after the 78th character of long input_IDs
      # ($sbjct->name). This is bad because with an extra whitespace, 
      # the name won't be parsed.  Therefore, a new variable 
      # $name_may_need_fix was introduced to fix the name by subsituting 
      # the whitespace away. For input IDs which are shorter than 78 
      # characters, the substitution has no effect. 
      # 
      # This heaeder-reformatting ( insertion of white spaces + linebreaks ) 
      # done by wublast can't be undone and the original heder can't be restored. 
      #
      my $name_may_need_fix = $sbjct->name;  
      #print "BEFORE FIXING NAME : $name_may_need_fix \n" ; 
      # $name_may_need_fix =~ s/\s//; 
      #print "AFTER  FIXING NAME : $name_may_need_fix \n" ; 
      
      my ($name) = $name_may_need_fix =~ /$regex/;   

      unless($name) { 
         throw("Error parsing name from ".$sbjct->name."\nCheck your ".
            "blast setup and blast headers ( your regex in your Blast config file)\n\t".
            "If your identifiers in the blast database are > 78 chars, you have to \n\t".
            "edit FilterBPLite.pm manually \n"); 
      }
    HSP: while (my $hsp = $sbjct->nextHSP) {
        if($self->is_hsp_valid($hsp)){     
          push(@output, $self->split_hsp($hsp, $name));
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



sub filter_hits{
  my ($self, $parsers) = @_;
  my %ids;
  my @features;
 PARSER:foreach my $parser(@$parsers){
  SUB:while(my $sbjct = $parser->nextSbjct){
      my $name = $sbjct->name;
    HSP:while (my $hsp = $sbjct->nextHSP) {
        if($self->is_hsp_valid($hsp)){
          my $qstart = $hsp->query->start();
          my $hstart = $hsp->subject->start();
          my $qend   = $hsp->query->end();
          my $hend   = $hsp->subject->end();
          my $qstrand = $hsp->query->strand();
          my $hstrand = $hsp->subject->strand();
          my $score  = $hsp->score;
          my $p_value = $hsp->P;
          my $percent = $hsp->percent;
          
          my $fp = $self->feature_factory->create_feature_pair
            ($qstart, $qend, $qstrand, $score, $hstart,
             $hend, $hstrand, $name, $percent, $p_value);
          
          push(@features,$fp);
        }
      }
    }
  }
 
  my $search = Bio::EnsEMBL::Analysis::Tools::FeatureFilter->new
    (
     -coverage => $self->coverage,
    );

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
