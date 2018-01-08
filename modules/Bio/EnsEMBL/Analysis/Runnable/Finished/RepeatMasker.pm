
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

Bio::EnsEMBL::Analysis::Runnable::Finished::RepeatMasker

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::Runnable::Finished::RepeatMasker->
  new(
      -query => $slice,
      -program => 'repeatmasker',
      -options => '-low'
     );
  $repeat_masker->run;
  my @repeats = @{$repeat_masker->output};

=head1 DESCRIPTION

RepeatMasker expects to run the program RepeatMasker
and produce RepeatFeatures which can be stored in the repeat_feature
and repeat_consensus tables in the core database


=head1 CONTACT

Post questions to : anacode@sanger.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::Finished::RepeatMasker;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::RepeatMasker);



=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Finished::RepeatMasker
  Arg [2]   : string, filename
  Function  : open and parse the results file into repeat
  features
  Returntype: none
  Exceptions: throws on failure to open or close output file
  Example   :

=cut



sub parse_results{
  my ($self, $results) =  @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $feature_factory = $self->feature_factory;
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "RepeatMasker:parse_results");
  REPEAT:while(<OUT>){
      if (/no repetitive sequences detected/) {
        print "RepeatMasker didn't find any repetitive sequences";
        last REPEAT; #no repeats found no point carrying on
      }
      if(/only contains ambiguous bases/){
        print "Sequence contains too many N's \n";
        last REPEAT;
      }
      my @columns;
      if(/\d+/){ #ignoring introductory lines
        chomp;
        @columns = split;
        pop @columns if $columns[-1] eq '*';
        #if (@columns != 15 || @columns != 14) {
        #  throw("Can't parse repeatmasker output unexpected number ".
        #        "of columns in the output ".@columns." in ".$_." ".
        #        "RepeatMasker:parse_results");
        #}

        my ($score, $name, $start, $end, $strand,
            $repeat_name, $repeat_class, $repeatmunge, $repeatmunge2);
        if(@columns == 15){
          ($score, $name, $start, $end, $strand,
           $repeat_name, $repeat_class) =  @columns[0, 4, 5, 6, 8, 9, 10];
        }elsif(@columns == 14){
          ($score, $name, $start, $end, $strand,
           $repeatmunge,$repeatmunge2) =  @columns[0, 4, 5, 6, 8, 9, 10];
          if ($repeatmunge =~ /(\S+)(LINE\S+)/) {
            $repeatmunge =~ /(\S+)(LINE\S+)/;
            $repeat_name = $1;
            $repeat_class = $2;
          } elsif ($repeatmunge2 eq 'Unknown') {
            print "Unknown repeat name\n";
            $repeat_name = 'Unknown';
          } else {
            $repeat_name = $repeatmunge;
            $repeat_class = $repeatmunge2;
            #throw("Can't parse repeatmasker output for line = $_\n");
          }
          if(!$repeat_class){
            $repeat_class = 'UNK';
          }
        }else{
          throw("Can't parse repeatmasker output unexpected number ".
                "of columns in the output ".@columns." in ".$_." ".
                "RepeatMasker:parse_results");
        }
        my $start_column;
        if($strand eq '+'){
          $start_column = 11;
          $strand = 1;
        }
        if($strand eq 'C'){
          $start_column = 13;
          $strand = -1;
        }#the location of the start and end inside the repeat
        #is different depending on the strand
        my $repeat_start = $columns[$start_column];
        my $repeat_end = $columns[12];
        if($repeat_end < $repeat_start){
          my $temp_end = $repeat_start;
          $repeat_start = $repeat_end;
          $repeat_end = $temp_end;
        }
        my $rc = $feature_factory->create_repeat_consensus($repeat_name,
                                                           $repeat_class);
        my $rf = $feature_factory->create_repeat_feature
          ($start, $end, $strand, $score,
           $repeat_start, $repeat_end,
           $rc, $self->query);

        $self->output([$rf]);
      }
    }
  close(OUT) or throw("FAILED to close ".$results.
                      "RepeatMasker:parse_results");
}

sub write_seq_file{
  my ($self, $seq, $filename) = @_;

  if(!$seq){
    my $slice = $self->query;
    $seq = Bio::PrimarySeq->new(
		-display_id => $slice->seq_region_name(),
		-seq        => $slice->seq()
	);

  }
  if(!$filename){
    $filename = $self->queryfile;
  }
  $filename = write_seqfile($seq, $filename);
  return $filename;
}


1;
