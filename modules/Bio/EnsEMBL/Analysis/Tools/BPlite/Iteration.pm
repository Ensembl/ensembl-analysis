=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
BPlite.pm is copyright (C) 1999 by Ian Korf. 
# Copyright Peter Schattner
# You may distribute this module under the same terms as perl itself
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 AUTHORS

  Peter Schattner - schattner@alum.mit.edu

=head1 NAME

Bio::Tools:: BPlite::Iteration - object for parsing single iteration
of a PSIBLAST report

=head1 SYNOPSIS

   use Bio::Tools:: BPpsilite;

   open FH, "t/psiblastreport.out";
   $report = Bio::Tools::BPpsilite->new(-fh=>\*FH);

   # determine number of iterations executed by psiblast
   $total_iterations = $report->number_of_iterations;
   $last_iteration = $report->round($total_iterations);

   # Process only hits found in last iteration ...
   $oldhitarray_ref = $last_iteration->oldhits;
   HIT: while($sbjct = $last_iteration->nextSbjct) {
       $id = $sbjct->name;
       $is_old =  grep  /\Q$id\E/, @$oldhitarray_ref;
       if ($is_old ){next HIT;}
   #  do something with new hit...
   }

=head1 DESCRIPTION

See the documentation for BPpsilite.pm for a description of the
Iteration.pm module.

=head1 METHODS

=cut

=head1 ACKNOWLEDGEMENTS

Based on work of:
Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
Lorenz Pollak (lorenz@ist.org, bioperl port)

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

package Bio::EnsEMBL::Analysis::Tools::BPlite::Iteration;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::BPlite; #
use Bio::EnsEMBL::Analysis::Tools::BPlite::Sbjct;

@ISA = qw();

sub new {
    my ($caller, @args) = @_;
    # my $self = $class->SUPER::new(@args);

    my $class = ref($caller) || $caller;
    my $self = bless({}, $class);


    ($self->{'FH'},$self->{'PARENT'},$self->{'ROUND'}) =
           rearrange([qw(FH
			      PARENT
			      ROUND
			      )],@args);
    
    if((! ref($self->{'FH'})) ||
       ((ref($self->{'FH'}) ne 'GLOB') &&
	(! $self->{'FH'}->isa('IO::Handle')))) {
	throw("Expecting a GLOB reference, not $self->{'FH'} !");
    }
    
    $self->{'LASTLINE'} = "";
    $self->{'QUERY'} = $self->{'PARENT'}->{'QUERY'};
    $self->{'LENGTH'} = $self->{'PARENT'}->{'LENGTH'};

    if ($self->_parseHeader) {$self->{'REPORT_DONE'} = 0} # there are alignments
    else                     {$self->{'REPORT_DONE'} = 1} # empty report
  
    return $self; # success - we hope!
}

=head2 query

 Title    : query
 Usage    : $query = $obj->query();
 Function : returns the query object
 Example  :
 Returns  : query object
 Args     :

=cut

sub query    {shift->{'QUERY'}}

=head2 qlength

 Title    : qlength
 Usage    : $len = $obj->qlength();
 Returns  : length of query
 Args     : none

=cut

sub qlength  {shift->{'LENGTH'}}

=head2 newhits

 Title    :  newhits
 Usage    : $newhits = $obj->newhits();
 Returns  : reference to an array listing all the hits 
            from the current iteration which were not identified 
            in the previous iteration
 Args     : none

=cut

sub newhits  {shift->{'NEWHITS'}}

=head2 oldhits

 Title    :  oldhits
 Usage    : $oldhits = $obj->oldhits();
 Returns  : reference to an array listing all the hits from 
            the current iteration which were identified and 
            above threshold in the previous iteration
 Args     : none

=cut

sub oldhits  {shift->{'OLDHITS'}}


=head2 nextSbjct

 Title    : nextSbjct
 Usage    : $sbjct = $obj->nextSbjct();
 Function : Method of iterating through all the Sbjct retrieved
            from parsing the report 
 Example  : while ( my $sbjct = $obj->nextSbjct ) {}
 Returns  : next Sbjct object or undef if finished
 Args     :

=cut

sub nextSbjct {
  my ($self) = @_;
  $self->_fastForward or return undef;

  #######################
  # get all sbjct lines #
  #######################
  my $def = $self->{'LASTLINE'};
  my $FH = $self->{'FH'};
  while(<$FH>) {
    if    ($_ !~ /\w/)            {next}
    elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
    elsif ($_ =~ /^\s{0,2}Score/) {$self->{'LASTLINE'} = $_; last}
    else                          {$def .= $_}
  }
  $def =~ s/\s+/ /g;
  $def =~ s/\s+$//g;
  $def =~ s/Length = ([\d,]+)$//g;
  my $length = $1;
  return 0 unless $def =~ /^>/;
  $def =~ s/^>//;

  ####################
  # the Sbjct object #
  ####################
  my $sbjct = new Bio::EnsEMBL::Analysis::Tools::BPlite::Sbjct('-name'=>$def,
					    '-length'=>$length,
                                            '-fh'=>$self->{'FH'}, 
					    '-lastline'=>$self->{'LASTLINE'}, 
					    '-parent'=>$self);
  return $sbjct;
}

sub _parseHeader {
  my ($self) = @_;
  my (@old_hits, @new_hits);
  my $FH = $self->{'FH'};
  my $newhits_true = ($self->{'ROUND'} < 2) ? 1  : 0 ;
  while(<$FH>) {
    if ($_ =~ /(\w\w|.*|\w+.*)\s\s+(\d+)\s+([-\.e\d]+)$/)    {
	my $id= $1;
	my $score= $2;	#not used currently
	my $evalue= $3; 	#not used currently
    	if ($newhits_true) { push ( @new_hits, $id);}
    	else { push (@old_hits, $id);}
    }
    elsif ($_ =~ /^Sequences not found previously/)  {$newhits_true = 1 ;}
    elsif ($_ =~ /^>/)
        {$self->{'LASTLINE'} = $_;
	 $self->{'OLDHITS'} = \@old_hits;
	 $self->{'NEWHITS'} = \@new_hits;
	 $self->{'LASTLINE'} = $_;	
	 return 1;
    }
    elsif ($_ =~ /^Parameters|^\s+Database:|^\s*Results from round\s+(d+)/) {
      	$self->{'LASTLINE'} = $_;
      	return 0; #  no sequences found in this iteration
    }
  }
  return 0; # no sequences found in this iteration
}

sub _fastForward {
  my ($self) = @_;
  return 0 if $self->{'REPORT_DONE'}; # empty report
  return 1 if $self->{'LASTLINE'} =~ /^>/;

  my $FH = $self->{'FH'};
  while(<$FH>) {
    if ($_ =~ /^>|^Parameters|^\s+Database:/) {
      $self->{'LASTLINE'} = $_;
      return 1;
    }
  }
  warning("Possible error while parsing BLAST report!");
}

1;
__END__
