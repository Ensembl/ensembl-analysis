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


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($database) = rearrange(['DATABASE'], @args);
  $self->database($database) if defined $database;

  throw("You must supply a database") if not $self->database; 
  throw("You must supply a query") if not $self->query;

  $self->program("blastz") if not $self->program;

  return $self;
}

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run{
  my ($self, $dir) = @_;

  $self->workdir($dir) if($dir);

  throw("Can't run ".$self." without a query sequence")
    unless($self->query);

  $self->write_seq_files();
  $self->run_analysis();

  $self->delete_files;
  return 1;
}



sub run_analysis {
  my $self = shift;

  my $cmd = $self->program  ." ".
            $self->query ." ".
            $self->database ." ".
            $self->options;

  my $BlastzParser;
  my $blastz_output_pipe = undef;
              
  if($self->results_to_file) {
    if (not $self->resultsfile) {
      my $resfile = $self->create_filename("blastz", "results");
      $self->resultsfile($resfile);
      $self->files_to_delete($resfile);
    }

    $cmd .=  " > ". $self->resultsfile;
    info("Running blastz...\n$cmd\n");

    throw("Error runing blastz cmd\n$cmd\n." .
                 " Returned error $? BLAST EXIT: '" .
                 ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
                 "', There was " . ($? & 128 ? 'a' : 'no') .
                 " core dump") unless(system($cmd) == 0);

    $BlastzParser = Bio::EnsEMBL::Analysis::Tools::Blastz->
        new('-file' => $self->resultsfile);
  } else {
    info("Running blastz to pipe...\n$cmd\n");

    open($blastz_output_pipe, "$cmd |") ||
      throw("Error opening Blasts cmd <$cmd>." .
                   " Returned error $? BLAST EXIT: '" .
                   ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
                   "', There was " . ($? & 128 ? 'a' : 'no') .
                   " core dump");
    $BlastzParser = Bio::EnsEMBL::Analysis::Tools::Blastz->
        new('-fh' => $blastz_output_pipe);
  }

  my @results;

  while (defined (my $alignment = $BlastzParser->nextAlignment)) { # nextHSP-like
    push @results, $alignment;
  }
  close($blastz_output_pipe) if(defined($blastz_output_pipe));

  $self->output(\@results);
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
  my ($self, $val) = @_;

  if (defined $val) {
    if (not ref($val)) {   
      throw("[$val] : file does not exist\n") unless -e $val;
    } elsif (not $val->isa("Bio::PrimarySeqI")) {
      throw("[$val] is neither a Bio::Seq not a file");
    }
    $self->{_query} = $val;
  }

  return $self->{_query}
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
  my ($self, $val) = @_;

  if (defined $val) {
    if (not ref($val)) {   
      throw("[$val] : file does not exist\n") unless -e $val;
    } else {
      if (ref($val) eq 'ARRAY') {
        foreach my $el (@$val) {
          throw("All elements of given database array should be Bio::PrimarySeqs")
              if not ref($el) or not $el->isa("Bio::PrimarySeq");
        }
      } elsif (not $val->isa("Bio::PrimarySeq")) {
        throw("[$val] is neither a file nor array of Bio::Seq");
      } else {
        $val = [$val];
      }
    }
    $self->{_database} = $val;
  }

  return $self->{_database};
}
  

sub write_seq_files {
  my ($self) = @_;

  if (ref($self->query)) {
    # write the query
    my $query_file = $self->create_filename("blastz", "query");
    my $seqio = Bio::SeqIO->new(-format => "fasta",
                                -file   => ">$query_file");
    $seqio->write_seq($self->query);
    $seqio->close;

    $self->query($query_file);
    $self->files_to_delete($query_file);
  }
  if (ref($self->database)) {
    my $db_file = $self->create_filename("blastz", "database");    
    my $seqio = Bio::SeqIO->new(-format => "fasta",
                                -file   => ">$db_file");
    foreach my $seq (@{$self->database}) {
      $seqio->write_seq($seq);
    }
    $seqio->close;

    $self->database($db_file);
    $self->files_to_delete($db_file);
  }
}


sub results_to_file {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_results_to_file} = $val;
  }

  return $self->{_results_to_file};
}

1;
