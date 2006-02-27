# Cared for by Abel Ureta-Vidal <abel@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Blastz

=head1 SYNOPSIS

  # To run a blastz job from scratch do the following.

  my $query = new Bio::SeqIO(-file   => 'somefile.fa',
                           -format => 'fasta')->next_seq;

  my $database = 'multifastafile.fa';

  my $blastz =  Bio::EnsEMBL::Analysis::Runnable::Blastz->new 
    ('-query'     => $query,
     '-database'  => $database,
     '-options'   => 'T=2');

  $blastz->run();

  @featurepairs = $blast->output();

  foreach my $fp (@featurepairs) {
      print $fp->gffstring . "\n";
  }

  # Additionally if you have blast runs lying around that need parsing
  # you can use the EnsEMBL blastz parser module 
  # perldoc Bio::EnsEMBL::Analysis::Runnable::Parser::Blastz


=head1 DESCRIPTION

Blastz takes a Bio::Seq (or Bio::PrimarySeq) object and runs blastz with against 
the specified multi-FASTA file database. Tthe output is parsed by 
Bio::EnsEMBL::Analysis::Runnable::Parser::Blastz and stored as Bio::EnsEMBL::DnaDnaAlignFeature 

Other options can be passed to the blast program using the -options method

=head1 CONTACT

Describe contact details here

=head1 APPENDIX


=cut

package Bio::EnsEMBL::Analysis::Runnable::Blastz;


use vars qw(@ISA);
use strict;

# Object preamble

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::Blastz;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Analysis::Runnable::Blast->new 
    (-query    => $seq,
     -database => $database,
     -options   => 'C=2 K=3000 H=2200');

    Function:   Initialises Blastz object
    Returns :   a Blastz Object
    Args    :   A Bio::Seq object (-query)
                A database file (-database)
                The blastz executable (-program)
                Options (-options)

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'}     = undef;     # file location of query sequence
    $self->{'_program'}   = "blastz";  # location of Blast
    $self->{'_database'}  = undef;     # name of database filename
    $self->{'_options'}   = "";        # options for blastz
    $self->{'_fplist'}    = [];        # an array of feature pairs (the output)
    $self->{'_workdir'}   = "/tmp";    # location of temp directory
    $self->{'_results'}   = $self->{'_workdir'}."/results.".$$; # location of result file
    $self->{'_results_to_tmp_file'} = 0;  # switch on whether to use pipe or /tmp file
    $self->{'_delete_results'} = 1;       # switch on whether to delete /tmp/results file or not
    $self->{'_verbose_debug'} = 0;
    
    # Now parse the input options and store them in the object
    my($program,$query,$database,$options) = rearrange([qw(PROGRAM
                                                           QUERY 
                                                           DATABASE 
                                                           OPTIONS)], 
                                                       @args);

    if ($query) {
      $self->query($query);
    } else {
      throw("No query sequence input.");
    }

    if ($database) {
      $self->database($database);
    } else {
      throw("No database input");
    }
    
    if ($options) {
       $self->options($options);
    } 
    # $options ||= ' C=2 ';
    
    return $self; # success - we hope!
}

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
  my ($self) = @_;
  
  $self->checkdir();
  
  $self->run_analysis();
  
  #parse output and create features
   #kfb 11.01.2006 changed to Analysis/Runnable method 
#  $self->deletefiles();
  $self->delete_files();
}


sub run_analysis {
  my $self = shift;

  my $cmd = $self->program  ." ".
            $self->query ." ".
            $self->database ." ".
            $self->options;

  my $BlastzParser;
  my $blastz_output_pipe = undef;
              
  if($self->{'_results_to_tmp_file'}) {
    $cmd .=  " > ". $self->results;
    print STDERR "Running blastz...\n$cmd\n" if($self->{'_verbose_debug'});
    throw("Error runing blastz cmd\n$cmd\n." .
                 " Returned error $? BLAST EXIT: '" .
                 ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
                 "', There was " . ($? & 128 ? 'a' : 'no') .
                 " core dump") unless(system($cmd) == 0);
    #throw("Failed during blastz run, $!\n") unless (system ($cmd));
    if($self->{'_delete_results'}) {
      $self->file($self->results);
    }
    $BlastzParser = Bio::EnsEMBL::Analysis::Tools::Blastz->new('-file' => $self->results);
  } else {
    print STDERR "Running blastz to pipe...\n$cmd\n" if($self->{'_verbose_debug'});
    open($blastz_output_pipe, "$cmd |") ||
      throw("Error opening Blasts cmd <$cmd>." .
                   " Returned error $? BLAST EXIT: '" .
                   ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
                   "', There was " . ($? & 128 ? 'a' : 'no') .
                   " core dump");
    $BlastzParser = Bio::EnsEMBL::Analysis::Tools::Blastz->new('-fh' => $blastz_output_pipe);
  }

  my $count=0;
  while (defined (my $alignment = $BlastzParser->nextAlignment)) { # nextHSP-like
    $self->_add_fp($alignment);
    $count++;
  }

  close($blastz_output_pipe) if(defined($blastz_output_pipe));
}


#################
# get/set methods 
#################

=head2 query

    Title   :   query
    Usage   :   $self->query($seq)
    Function:   Get/set method for query.  If set with a Bio::Seq object it
                will get written to the local tmp directory
    Returns :   filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub query {
  my ($self, $value) = @_;

  if (defined $value) {

     if (! ref($value)) {
      # assume it's a filename - check the file exists
      throw("[$value] : file does not exist\n") unless -e $value;
      $self->{'_query'} = $value;
    }
    elsif ($value->isa("Bio::PrimarySeqI") || $value->isa("Bio::Seq")) {
      my $filename = "/tmp/genfile_$$.".rand(time()).".fn";
      my $genOutput = Bio::SeqIO->new(-file => ">$filename" , '-format' => "Fasta")
	or throw("Can't create new Bio::SeqIO from $filename '$' : $!");

      throw("problem writing genomic sequence to $filename\n" ) unless $genOutput->write_seq($value);
      $self->{'_query'} = $filename;
      $self->file($filename);
    }
    else {
      throw("$value is neither a Bio::Seq  nor a filename\n");
    }
  }

  return $self->{'_query'};
}

=head2 database
  
    Title   :   database
    Usage   :   $self->database($seq)
    Function:   Get/set method for database.  If set with a Bio::Seq object it
                will get written to the local tmp directory
    Returns :   filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub database {
  my($self, $value) = @_;    
  
  if (defined $value) {

     if (! ref($value)) {
      # assume it's a filename - check the file exists
      throw("[$value] : file does not exist\n") unless -e $value;
      $self->{'_database'} = $value;
    }
    elsif ($value->isa("Bio::PrimarySeqI")) {
      my $filename = "/tmp/genfile_$$.".rand(time()).".fn";
      my $genOutput = Bio::SeqIO->new(-file => ">$filename" , '-format' => "Fasta")
	or throw("Can't create new Bio::SeqIO from $filename '$' : $!");
    
      throw("problem writing genomic sequence to $filename\n" ) unless $genOutput->write_seq($value);
      $self->{'_database'} = $filename;
      $self->file($filename);
    }
    else {
      throw("$value is neither a Bio::Seq  nor a filename\n");
    }
  }
  
  return $self->{'_database'};
}

=head2 options

    Title   :   options
    Usage   :   $obj->options(' -I ');
    Function:   Get/set method for blast arguments
    Args    :   File path (optional)

=cut

sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_options'} = $args ;
  }
  return $self->{'_options'};
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   None

=cut

sub output {
  my ($self) = @_;
   #kfb 11.01.2006 changed to Analysis/Runnable method 
#  return @{$self->{'_fplist'}};
  return $self->{'_fplist'};
}

sub _add_fp {
  my ($self,@args) = @_;
  if (@args) {
    push(@{$self->{'_fplist'}},@args);
  } else {
    warning("Bio::EnsEMBL::Analysis::Runnable::Blastz->_add_fp should have an argument\n");
  }
}

=head2 workdir

 Title   : workdir
 Usage   : $obj->workdir($newval)
 Function: 
 Example : 
 Returns : value of workdir
 Args    : newvalue (optional)


=cut

sub workdir{
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'_workdir'} = $value;
   }
   return $self->{'_workdir'};
}

