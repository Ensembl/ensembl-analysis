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

  Bio::EnsEMBL::Analysis::Tools::BPliteWrapper

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper->
  new(
      -regex => '^\w+\s+(\w+)'
      -query_type => 'dna',
      -database_type => 'pep',
     );
 my @results = @{$parser->parse_results('blast.out')};

=head1 DESCRIPTION

This module is a wrapper for BPlite to provide an interface between
Bio::EnsEMBL::Analysis::Runnable::Blast and BPlite. This method fits model
for standard blast parsers as it provides the parse_file method which
returns an array of results. This method just uses BPlite to parse the
file it does no pre or post filtering and as such will not mimic the 
behaviour or the current blast runnable in the pipeline code

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;

use strict;
use warnings;
use FileHandle;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Analysis::Tools::BPlite;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::PepDnaAlignFeature;

use vars qw (@ISA);

@ISA = qw();


=head2 new

  Arg [1]             : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [FILENAME]      : string, filename
  Arg [REGEX]         : string, regex
  Arg [QUERY_TYPE]    : string, query sequence type, should be pep or dna
  Arg [DATABASE_TYPE] : string, database sequence type as above
  Arg [ANALYSIS]      : Bio::EnsEMBL::Analysis object
  Function  : 
  Returntype: Bio::EnsEMBL::Analysis::Tools::BPliteWrapper 
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($caller,@args) = @_;
  
  my $class = ref($caller) || $caller;
  my $self = bless({}, $class);

  &verbose('WARNING');
  my ($regex, $query, $target, $analysis) = 
    rearrange(['REGEX', 'QUERY_TYPE', 'DATABASE_TYPE',
               'ANALYSIS'], @args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->regex('(\w+)\s+');
  ######################
  $self->regex($regex) if(defined $regex);
  $self->query_type($query) if($query);
  $self->database_type($target) if($target);
  $self->analysis($analysis);
  return $self;
}




=head2 regex

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : string, regex
  Function  : container
  Returntype: string, regex
  Exceptions: 
  Example   : 

=cut


sub regex{
  my $self = shift;
  $self->{'regex'} = shift if(@_);
  return $self->{'regex'};
}




sub filenames{
  my ($self, $files) = @_;
  if(!$self->{'filenames'}){
    $self->{'filenames'} = [];
  }
  if($files){
    throw($files." must be an arrayref BPliteWrapper:filenames ") 
      unless(ref($files) eq 'ARRAY');
    push(@{$self->{'filenames'}}, @{$files});
  }
  return $self->{'filenames'};
}


=head2 query_type

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : string, query sequence type
  Function  : container
  Returntype: string
  Exceptions: throws if string is not either dna or pep
  Example   : 

=cut



sub query_type{
  my ($self, $dtype) = @_; 
  if($dtype){
    $dtype = lc($dtype);
    throw("Query type must be either dna or pep not ".$dtype)
      unless($dtype eq 'pep' || $dtype eq 'dna');
    $self->{'query_type'} = $dtype;
  }
  return $self->{'query_type'};
}

=head2 database_type

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : string, database sequence type
  Function  : container
  Returntype: string
  Exceptions: throws if string is not either dna or pep
  Example   : 

=cut

sub database_type {
  my ($self, $dtype) = @_; 
  if($dtype){
    $dtype = lc($dtype);
    throw("Database type must be either dna or pep not ".$dtype)
      unless($dtype eq 'pep' || $dtype eq 'dna');
    $self->{'database_type'} = $dtype;
  }
  return $self->{'database_type'};
}


=head2 output

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : arrayref of output features
  Function  : store the output features in an array
  Returntype: arrayref
  Exceptions: throws if not passed an arrayref
  Example   : 

=cut



sub output {
  my ($self, $output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($output){
    throw("Must pass and arrayref not a ".$output." BPliteWrapper:output")
      unless(ref($output) eq 'ARRAY');
    push(@{$self->{'output'}}, @$output);
  }
  return $self->{'output'};
}


=head2 clean_output

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Function  : empties the internal output array
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub clean_output {
  my ($self) = @_;
  $self->{'output'} = [];
}
sub clean_filenames{
  my ($self) = @_;
  $self->{'filenames'} = [];
}

=head2 feature_factory

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Function  : container for feature factory creates one if one is requested
  but one does not currently exist
  Returntype: Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Exceptions: 
  Example   : 

=cut



sub feature_factory{
  my ($self, $feature_factory) = @_;
  if($feature_factory){
    $self->{'feature_factory'} = $feature_factory;
  }
  if(!$self->{'feature_factory'}){
    $self->{'feature_factory'} = 
      Bio::EnsEMBL::Analysis::Tools::FeatureFactory->new();
  }
  return $self->{'feature_factory'};
}

=head2 analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Analysis
  Function  : container for analysis object
  Returntype: Bio::EnsEMBL::Analysis
  Exceptions: throws passed incorrect object type
  Example   : 

=cut



sub analysis{
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis".
          "not a ".$analysis) unless($analysis->isa
                                     ('Bio::EnsEMBL::Analysis'));
    $self->{'analysis'} = $analysis;
  }
  return $self->{'analysis'};
}



=head2 parse_files

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
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
  $self->get_hsps($bplites);
  return $self->output;
}



=head2 get_parsers

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
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

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::BPlite
  Function  : get the hsps from bplite and turn then into features
  Returntype: none
  Exceptions: throws if regex does not produce a result
  Example   : 

=cut



sub get_hsps{
  my ($self, $parsers) = @_;
  my $regex = $self->regex;
  my @output;
 PARSER:foreach my $parser(@$parsers){
  NAME: while(my $sbjct = $parser->nextSbjct){
      my ($name) = $sbjct->name =~ /$regex/;
      throw("Error parsing name from ".$sbjct->name." check your ".
            "blast setup and blast headers") unless($name);
    HSP: while (my $hsp = $sbjct->nextHSP) {
        push(@output, $self->split_hsp($hsp, $name));
      }
    }
  }
  $parsers = [];
  $self->output(\@output);
}



=head2 split_hsp

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::BPlite::HSP
  Function  : turn hsp into alignfeature
  Returntype: Bio::EnsEMBL::BaseAlignFeature
  Exceptions: 
  Example   : 

=cut


sub split_hsp {
    my ($self,$hsp,$name) = @_;
    my $qstrand = $hsp->query->strand;
    my $hstrand = $hsp->subject->strand;
    my ($qinc,   $hinc)    = $self->find_increments($qstrand,$hstrand);
    my @qchars = split(//,$hsp->querySeq);  # split alignment into array of
                                            # chars
    my @hchars = split(//,$hsp->sbjctSeq);  # ditto for hit sequence
    my $qstart = $hsp->query->start(); # Start off the feature pair start
    my $hstart = $hsp->subject->start(); # ditto
    my $qend   = $hsp->query->start(); # Set the feature pair end also
    my $hend   = $hsp->subject->start(); # ditto
    if ($qstrand == -1) {
      $qstart = $hsp->query->end;
      $qend   = $hsp->query->end;
    }
    if ($hstrand == -1) {
      $hstart = $hsp->subject->end;
      $hend   = $hsp->subject->end;
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



=head2 find_increments

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : int, query strand
  Arg [3]   : int, hit strand
  Function  : work out the query and hit increments 
  Returntype: int, int query inc, hit inc
  Exceptions: none
  Example   : 

=cut



sub find_increments{
  my ($self, $qstrand, $hstrand) = @_;
  my $qinc   = 1 * $qstrand;
  my $hinc   = 1 * $hstrand;

  my $qtype = lc($self->query_type);
  my $htype = lc($self->database_type);

  if ($qtype eq 'dna' && $htype eq 'pep') {
    $qinc = 3 * $qinc;
  } 
  if ($qtype eq 'pep' && $htype eq 'dna') {
    $hinc = 3 * $hinc;
  }
  
  return ($qinc,$hinc);
}



=head2 convert_to_featurepair

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
  Arg [2]   : int, query start
  Arg [3]   : int, query end
  Arg [4]   : int, query strand
  Arg [5]   : int, query increment
  Arg [6]   : int, hit start
  Arg [7]   : int, hit end
  Arg [8]   : int, hit strand
  Arg [9]   : int, hit increment
  Arg [10]  : string, name
  Arg [11]  : string, query name
  Arg [12]  : int, bit score
  Arg [13]  : int, percent identity
  Arg [14]  : int, p value
  Arg [15]  : int, positive matches
  Arg [16]  : int, matches
  Function  : take values and taking account for inc values make feature
  pair 
  Returntype:  Bio::EnsEMBL::FeaturePair
  Exceptions: 
  Example   : 

=cut


sub convert_to_featurepair{
  my ($self, $qstart, $qend, $qstrand, $qinc, $hstart, $hend, $hstrand, 
      $hinc, $name, $seqname,$score, $percent, $pvalue, $positive, $matches) = @_;
  my $tmpqend = $qend; $tmpqend -= $qinc;
  my $tmphend = $hend; $tmphend -= $hinc;
    
  my $tmpqstart = $qstart;
  my $tmphstart = $hstart;
  
  # This is for dna-pep alignments.  The actual end base
  # will be +- 2 bases further on.
  if (abs($qinc) > 1) {
    $tmpqend += $qstrand * 2;
  }
  if (abs($hinc) > 1) {
    $tmphend += $hstrand * 2;
  }
  # Make sure start is always < end
  if ($tmpqstart > $tmpqend) {
    my $tmp    = $tmpqstart;
    $tmpqstart = $tmpqend;
    $tmpqend   = $tmp;
  }
  if ($tmphstart > $tmphend) {
    my $tmp    = $tmphstart;
    $tmphstart = $tmphend;
    $tmphend   = $tmp;
  }
  my $fp = $self->feature_factory->create_feature_pair($tmpqstart, 
                                                       $tmpqend, $qstrand,
                                                       $score, $tmphstart,
                                                       $tmphend, $hstrand,
                                                       $name, $percent, 
                                                       $pvalue, $seqname, 
                                                       undef, 
                                                       $self->analysis, 
                                                       $positive, 
                                                       $matches);

  return $fp;

}


1;
