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

package Bio::EnsEMBL::Analysis::Tools::Logger;

use strict;
use warnings;
use Exporter;
use vars qw(@ISA @EXPORT);
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);


@ISA = qw(Exporter);
@EXPORT = qw(logger_verbosity logger_info);


my $DEFAULT_OFF             = 0;
my $DEFAULT_INFO            = 4000;
my $DEFAULT_INFO_WITH_TRACE = 5000;

my $VERBOSITY = $DEFAULT_OFF;

sub logger_verbosity{
  if(@_) {
    my $verbosity = shift;
    $verbosity = shift 
      if($verbosity && 
         $verbosity eq "Bio::EnsEMBL::Utils::Exception");
    $verbosity = $VERBOSITY if(!$verbosity);
    if($verbosity =~ /\d+/) { #check if verbosity is an integer
      $VERBOSITY = $verbosity;
    } else {
      if($verbosity eq 'OFF' || $verbosity eq 'NOTHING' ||
         $verbosity eq 'NONE'){
        $VERBOSITY = $DEFAULT_OFF;
      }elsif($verbosity eq 'INFO' || 
             $verbosity eq 'LOGGER_INFO'){
        $VERBOSITY = $DEFAULT_INFO;
      } elsif ($verbosity eq 'INFO_STACK_TRACE' ||
               $verbosity eq 'LOGGER_INFO_STACK_TRACE') {
        $VERBOSITY = $DEFAULT_INFO_WITH_TRACE;
      }
      elsif($verbosity eq 'ALL' || 'ON'){
        $VERBOSITY = 1e6;
      } else {
        $VERBOSITY = $DEFAULT_OFF;
        warning("Unknown level or verbosity :".$verbosity);
      }
    }
  }
  return $VERBOSITY;
}

sub logger_info {
  my $string = shift;

  if ($VERBOSITY < $DEFAULT_INFO) {
    return;
  } elsif ($VERBOSITY < $DEFAULT_INFO_WITH_TRACE) {
    print STDERR "INFO: $string\n";
    return;
  }

  my @caller = caller;
  my $line = $caller[2] || '';

  #use only 2 subdirs for brevity when reporting the filename
  my $file;
  my @path = split(/\//, $caller[1]);
  $file = pop(@path);
  my $i = 0;
  while(@path && $i < 2) {
    $i++;
    $file = pop(@path) ."/$file";
  }

  @caller = caller(1);
  my $caller_line;
  my $caller_file;
  $i=0;
  if(@caller) {
     @path = split(/\//, $caller[1]);
     $caller_line = $caller[2];
     $caller_file = pop(@path);
     while(@path && $i < 2) {
       $i++;
       $caller_file = pop(@path) ."/$caller_file";
     }
  }

  
  my $out = "\n-------------------- LOG INFO ---------------------\n".
              "MSG: $string\n".
              "FILE: $file LINE: $line\n";
  $out .=     "CALLED BY: $caller_file  LINE: $caller_line\n" if($caller_file);
  $out .=     "---------------------------------------------------\n";
  print STDERR $out;
}

1;
