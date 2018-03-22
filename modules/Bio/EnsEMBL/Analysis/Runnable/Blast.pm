=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Blast - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::Blast->
  new(
      -query => $slice,
      -program => 'wublastn',
      -database => 'embl_vertrna',
      -options => 'hitdist=40 -cpus=1',
      -parser => $bplitewrapper,
      -filter => $featurefilter,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

  This module is a wrapper for running blast. It knows how to construct
  the commandline and can call to other modules to run the parsing and 
  filtering. By default is constructs wublast commandlines but it can be
  told to construct ncbi command lines. It needs to be passed a Bio::Seq
  and a database name (this database should either have its full path 
  given or it should live in the location specified by the $BLASTDB 
  environment variable). It should also be given a parser object which has
  the method parse_file which takes a filename and returns an arrayref of
  results and optionally it can be given a filter object which has the 
  method filter_results which takes an arrayref of results and returns the
  filtered set of results as an arrayref. For examples of both parser
  objects and a filter object look in Bio::EnsEMBL::Analysis::Tools for
  BPliteWrapper, FilterBPlite and FeatureFilter


=cut


package Bio::EnsEMBL::Analysis::Runnable::Blast;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_timer);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [Parser]               : Object A blast parser object must meet specified interface
 Arg [Filter]               : Object A Filter object must meet specified interface
 Arg [Database]             : String, database name/path
 Arg [Type]                 : String, wu or ncbi to specify which type of input
 Arg [Unknown_error_string] : String to throw if the blast runs fails
                              with an unexpected error 4
 Description                : create a Blast runnable 
 Returntype                 : Bio::EnsEMBL::Analysis::Runnable::Blast
 Exceptions                 : throws if not given a database name or if not given a parser object

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($parser, $filter, $database, $type,
      $unknown_error ) = rearrange(['PARSER', 'FILTER', 'DATABASE',
                                    'TYPE', 'UNKNOWN_ERROR_STRING',
                                   ], @args);
  $type = undef unless($type);
  $unknown_error = undef unless($unknown_error);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->type('ncbi');
  $self->unknown_error_string('FAILED');
  $self->options('-num_threads=1') if(!$self->options);
  ######################
  $self->databases($database);
  $self->parser($parser);
  $self->filter($filter);
  $self->type($type) if($type);
  $self->unknown_error_string($unknown_error) if($unknown_error);

  throw("No valid databases to search")
      unless(@{$self->databases});

  throw("Must pass Bio::EnsEMBL::Analysis::Runnable::Blast ".
        "a parser object ") 
      unless($self->parser);

  return $self;
}



=head2 databases

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub databases{
  my ($self, @vals) = @_;

  if (not exists $self->{databases}) {
    $self->{databases} = [];
  }

  foreach my $val (@vals) {
    my $dbname = $val;

    my @dbs;

    $dbname =~ s/\s//g;

    # prepend the variable $BLASTDB from Config/Blast.pm
    # if database name is not an absolute path
  
    unless ($dbname =~ m!^/!) {
        warning('Specify an full path to the database, we are NOT using Bio::EnsEMBL::Analysis::Config::Blast anymore');
    }
  
    # If the expanded database name exists put this in
    # the database array.
    #
    # If it doesn't exist then see if $database-1,$database-2 exist
    # and put them in the database array
    
    if (-f "$dbname" || -f "$dbname.fa" || -f "$dbname.xpd" || -f "$dbname.phr" || -f "$dbname.nhr" || -f "$dbname.nal") {
      push(@dbs,$dbname);
    } else {
      my $count = 1;
      while (-f "$dbname-$count") {
        push(@dbs,"$dbname-$count");
        $count++; 	 
      }
    }

    if (not @dbs) {
      warning("Valid BLAST database could not be inferred from '$val'");
    } else {
      push @{$self->{databases}}, @dbs;
    }
  }

  return $self->{databases};
}


=head2 parser 

  Arg [1]    : Object
  Description: Getter/setter
  Returntype : Object
  Exceptions : None

=cut
sub parser{
  my $self = shift;
  $self->{'parser'} = shift if(@_);
  return $self->{'parser'};
}


=head2 filter

 Arg [1]    : Object
 Description: Getter/setter
 Returntype : Object
 Exceptions : None

=cut
sub filter{
  my $self = shift;
  $self->{'filter'} = shift if(@_);
  return $self->{'filter'};
}


=head2 type

 Arg [1]    : String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub type{
  my $self = shift;
  $self->{'type'} = shift if(@_);
  return $self->{'type'};
}


=head2 unknown_error_string

 Arg [1]    : String
 Description: container for given value, this describes the 5 methods
 Returntype : String
 Exceptions : None

=cut
sub unknown_error_string{
  my $self = shift;
  $self->{'unknown_error_string'} = shift if(@_);
  return $self->{'unknown_error_string'};
}


=head2 results_files

 Arg [1]    : String, filename
 Description: holds a list of all the output files from the blast runs
 Returntype : Arrayref
 Exceptions : None

=cut

sub results_files{
  my ($self, $file) = @_;
  if(!$self->{'results_files'}){
    $self->{'results_files'} = [];
  }
  if($file){
    push(@{$self->{'results_files'}}, $file);
  }
  return $self->{'results_files'};
}


=head2 run_analysis

 Arg [1]    : None
 Description: Gets a list of databases to run against and constructs
              commandlines against each one and runs them
 Returntype : None
 Exceptions : Throws if there is a problem opening the commandline or
              if blast produces an error

=cut

sub run_analysis {
  my ($self) = @_;

  foreach my $database (@{$self->databases}) {

    my $db = $database;
    $db =~ s/.*\///;
    #allow system call to adapt to using ncbi blastall. 
    #defaults to WU blast
    my $command  = $self->program;
    my $filename = $self->queryfile;
    my $results_file = $self->create_filename($db, 'blast.out');
    $self->files_to_delete($results_file);
    $self->results_files($results_file);
    if ($self->type eq 'legacy_ncbi') {
      if (!exists $ENV{BLASTMAT} or ! -e $ENV{BLASTMAT}) {
        throw(' Your environment variable $BLASTMAT is not set !!! '.
              " Point it to /usr/local/ensembl/data/blastmat/ or where your BLOSUM62 matrices live\n");
      }
      $command .= " -d $database -i $filename ";
    }
    elsif ($self->type eq 'ncbi') {
      $command .= " -db $database -query $filename ";
    }
    else {
      if (!exists $ENV{WUBLASTMAT} or ! -e $ENV{WUBLASTMAT}) {
        throw(' Your environment variable $BLASTMAT is not set !!! '.
              " Point it to /usr/local/ensembl/data/blastmat/ or where your BLOSUM62 matrices live\n");
      }
      $command .= " $database $filename -gi ";
    }
    $command .= $self->options. ' 2>&1 > '.$results_file;

    print "Running blast ".$command."\n";
    info("Running blast ".$command); 


    # I don't thing the vast majority of error parsing is needed, so I will not integrate it for the timer
    # assuming we don't run into problems in the future with BLAST the else can probably be deleted
    if($self->timer) {
      execute_with_timer($command, $self->timer);
    } else {
      open(my $fh, "$command |") || 
        throw("Error opening Blast cmd <$command>." .
              " Returned error $? BLAST EXIT: '" . 
              ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) . 
              "', There was " . ($? & 128 ? 'a' : 'no') . 
              " core dump");
      # this loop reads the STDERR from the blast command
      # checking for FATAL: messages (wublast) [what does ncbi blast say?]
      # N.B. using simple die() to make it easier for RunnableDB to parse.
      while(<$fh>){
        if(/FATAL:(.+)/){
          my $match = $1;
          print $match;
          # clean up before dying
          $self->delete_files;
          if($match =~ /no valid contexts/){
            die qq{"VOID"\n}; # hack instead
          }elsif($match =~ /Bus Error signal received/){
            die qq{"BUS_ERROR"\n}; # can we work out which host?
          }elsif($match =~ /Segmentation Violation signal received./){
            die qq{"SEGMENTATION_FAULT"\n}; # can we work out which host?
          }elsif($match =~ /Out of memory;(.+)/){
            # (.+) will be something like "1050704 bytes were last 
            #requested."
            die qq{"OUT_OF_MEMORY"\n}; 
            # resenD to big mem machine by rulemanager
          }elsif($match =~ /the query sequence is shorter than the word length/){
            #no valid context 
            die qq{"VOID"\n}; # hack instead
          }elsif($match =~ /External filter/){
            # Error while using an external filter
            die qq{"EXTERNAL_FITLER_ERROR"\n}; 
          }else{
            warning("Something FATAL happened to BLAST we've not ".
                    "seen before, please add it to Package: " 
                    . __PACKAGE__ . ", File: " . __FILE__);
            die ($self->unknown_error_string."\n"); 
            # send appropriate string 
            #as standard this will be failed so job can be retried
            #when in pipeline
          }
        }elsif(/WARNING:(.+)/){
          # ONLY a warning usually something like hspmax=xxxx was exceeded
          # skip ...
        }elsif(/^\s{10}(.+)/){ # ten spaces
          # Continuation of a WARNING: message
          # Hope this doesn't catch more than these.
          # skip ...
        }
      }
      unless(close $fh){
        # checking for failures when closing.
        # we should't get here but if we do then $? is translated 
        #below see man perlvar
        warning("Error running Blast cmd <$command>. Returned ".
                "error $? BLAST EXIT: '" . ($? >> 8) . 
                "', SIGNAL '" . ($? & 127) . "', There was " . 
                ($? & 128 ? 'a' : 'no') . " core dump");
        die ($self->unknown_error_string."\n"); 
      }
    }
  }
}


=head2 parse_results

 Arg [1]    : None
 Description: call to parser to get results from output file
              and filter those results if there is a filter object
 Returntype : None
 Exceptions : None

=cut

sub parse_results{
  my ($self) = @_;
  my $results = $self->results_files;
  my $output = $self->parser->parse_files($results);
  my $filtered_output;
  #print "Have ".@$output." features to filter\n";
  if($self->filter){
    $filtered_output = $self->filter->filter_results($output);
  }else{
    $filtered_output = $output;
  }
  $self->output($filtered_output);
}


=head2 store_slice

 Arg [1]    : Bio::EnsEMBL::Slice
 Description: used to explictly get/set the slice since query is used inconsistently
 Returntype : Bio::EnsEMBL::Slice
 Exceptions : throw if not a slice ref

=cut

sub store_slice {
  my ($self,$slice) = @_;
  if($slice) {
    unless(ref($slice) eq 'Bio::EnsEMBL::Slice') {
      throw("Must pass in a slice ref. Ref type found: ".ref($slice));
    }
    $self->{'_store_slice'} = $slice;
  }

  return($self->{'_store_slice'});
}

1;
