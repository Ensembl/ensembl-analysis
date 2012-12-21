# Author: Gary Williams (gw3@sanger.ac.uk)
# Copyright (c) Gary Williams (2005)
# based on Michelle Clamp's Blast.pm
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Panther

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Panther->new 
    ('-query'          => $query,
     '-program'        => 'Panther' or '/usr/local/pubseq/bin/Panther',
     '-database'       => 'sf_hmm');

  $hmm->workdir ($workdir);
  $hmm->run;
  my @results = $hmm->output;

=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs Panther.
  The resulting output file is parsed to produce a set of Bio::EnsEMBL::FeaturePairs.

=head1 CONTACT

   Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

=cut

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Runnable/ProteinAnnotation/Panther.pm,v $
# $Revision: 1.4 $
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Panther;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub run_analysis {
  my ($self) = @_;

  # run program
  print STDERR "running ".$self->program." against ".$self->database."\n";

  # some of these options require HMMER 2.2g (August 2001)
    
  my @dbfiles = split(/;/,$self->analysis->db_file);

  print STDERR "FILENAME: ".$self->queryfile."\n";
  
                                                                                                                   
  my $cmd = $self->program .' '.
	$self->options .' '.
  	'-n '.
	'-l ' . $dbfiles[0] . ' '.
	'-i ' . $self->queryfile . ' '.
	'-H /software/worm/iprscan/bin/binaries/hmmsearch ' . 
	'-B /software/worm/iprscan/bin/binaries/blast/blastall ' .
	'-T /tmp ' .
	'-D I -z -E 1e-11 ' .
	'-o ' . $self->resultsfile;

  print STDERR "$cmd\n";   
    $self->throw ("Error running ".$self->program." on ".$self->queryfile." against ".$dbfiles[0]) 
    unless ((system ($cmd)) == 0);
}


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses program output to give a set of features
 Example  :
 Returns  : 
 Args     : filename (optional, can be filename, filehandle or pipe, not implemented)
 Throws   :

=cut

sub parse_results {
  my ($self) = @_;

  my $filehandle;
  my $resfile = $self->resultsfile();
  my @fps;
 

    if (-e $resfile) {
	if (-z $self->resultsfile) {  
	    print STDERR "Panther didn't find any hits\n";
	    return; }       
	else {
	    open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
	    }
    }

   
# example output line:
#R186.4  PTHR10192       MOLYBDOPTERIN BIOSYNTHESIS PROTEIN      3.7e-116        396.7   12-391
#5       PTHR11865       NUCLEAR HORMONE RECEPTOR        3.2e-13     54.8         1-40
#^(\S+)\s+(\S+)      \s+          ([A-Z]+\s)+              (\S+)  \s+(\S+)\s+ (\S+)-(\S+)/


  while (<CPGOUT>) {
    chomp;
    next if (/^\/\//);

    print "$_\n";
    if (my ($id, $hid,$description,$evalue, $score, $start, $end) = /^(\S+)\s+(\S+)\s+(.+)?\s+(\S+)\s+(\S+)\s+(\S+)-(\S+)/) {
      print "matched ($id, $hid,$description,$evalue, $score, $start, $end)\n";
      $hid =~ s/:SF\d+//;

      my $fp= $self->create_protein_feature($start,$end,$score,$id,0,0,$hid,$self->analysis, sprintf("%.3e", $evalue),0);
      push @fps,$fp;
    }
  }
  close (CPGOUT);
  $self->output(\@fps);
}

1;
