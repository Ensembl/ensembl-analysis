package Bio::EnsEMBL::Analysis::Tools::Logger;

use strict;
use warnings;
use Exporter;
use vars qw(@ISA @EXPORT);
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);


@ISA = qw(Exporter);
@EXPORT = qw(logger_verbosity logger_info 
             logger_warning);


my $VERBOSITY = 3000;
my $DEFAULT_INFO = 4000;
my $DEFAULT_WARNING = 2000;


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
        $VERBOSITY = 0;
      }elsif($verbosity eq 'WARN' || $verbosity eq 'WARNING' ||
             $verbosity eq 'LOGGER_WARNING'){
        $VERBOSITY = $DEFAULT_WARNING;
      }elsif($verbosity eq 'INFO' || 
             $verbosity eq 'LOGGER_INFO'){
        $VERBOSITY = $DEFAULT_INFO;
      }elsif($verbosity eq 'ALL' || 'ON'){
        $VERBOSITY = 1e6;
      }else{
        warning("Unknown level or verbosity :".$verbosity);
      }
    }
  }
}

sub logger_info{
  my $string = shift;
  my $level = shift;
  $level = $DEFAULT_INFO if(!defined($level));
  return if($VERBOSITY < $level);
  print STDERR "INFO: $string\n";
}

sub logger_warning{
  my $string = shift;
  my $level = shift;
  $level = $DEFAULT_WARNING if(!defined($level));
  return if($VERBOSITY < $level);
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

  
  my $out = "\n-------------------- WARNING ----------------------\n".
              "MSG: $string\n".
              "FILE: $file LINE: $line\n";
  $out .=     "CALLED BY: $caller_file  LINE: $caller_line\n" if($caller_file);
  $out .=     "---------------------------------------------------\n";
  print STDERR $out;
}
