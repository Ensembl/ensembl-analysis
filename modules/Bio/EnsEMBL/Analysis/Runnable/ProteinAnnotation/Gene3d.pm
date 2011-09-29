=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Gene3d - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# based on Michelle Clamp's Blast.pm
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Genes3d

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Genes3d->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => 'Pfam');

  $hmm->workdir ($workdir);
  $hmm->run;
  my @results = $hmm->output;

=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs hmmpfam.
  The resulting output file is parsed to produce a set of Bio::EnsEMBL::FeaturePairs.

=head1 CONTACT

   Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Gene3d;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

###################
# analysis methods
###################

=head2 run_analysis

 Title    : run_program
 Usage    : $self->program
 Function : makes the system call to program
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut

sub run_analysis {
    my ($self) = @_;

    # run program
    print STDERR "running ".$self->program." against ".$self->database."\n";

    # some of these options require HMMER 2.2g (August 2001)
    
    my @dbfiles = split(/;/,$self->analysis->db_file);

    print STDERR "FILENAME: ".$self->queryfile."\n";
 
# /software/worm/iprscan/bin/binaries/hmmer3/hmmscan --domE 1 -Z 100000 /data/blastdb/Worms/interpro_scan/iprscan/data/gene3d.lib inputfile > resultsfile

    my $cmd = $self->program .' '. 
	$self->analysis->parameters .' '.
	  $self->database .' '.
	    $self->queryfile.' > '.
	      $self->resultsfile;
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
    if (-z $resfile) {  
      print STDERR "Genes3D didn't find any hits\n";
      return; }       
    else {
      open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
    }
  }

   
  my %seen_before;




#Domain annotation for each model (and alignments):
#>> 1ukxA00
#   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
# ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#   1 ?   10.5   0.0   2.1e-05       2.1      46      81 ..      23      58 ..      20      88 .. 0.90
#   2 ?   -2.7   0.0      0.24   2.4e+04      67      83 ..     136     152 ..     108     168 .. 0.77
#
#  Alignments for each domain:
#  == domain 1    score: 10.5 bits;  conditional E-value: 2.1e-05
#                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX RF
#      1ukxA00 46 sfkikikpesdeeeeqkvsltLqvelpetYPdeaPe 81
#                  fk+++++  de++ + + + L+++l ++ Pde+P
#  temp_gene_1 23 AFKFHLSTLVDEDDPTPFDFALSFQLKPELPDELPA 58
#                 5999999999************************96 PP
#
#  == domain 2    score: -2.7 bits;  conditional E-value: 0.24
#                  XXXXXXXXXXXXXXXXX RF
#      1ukxA00  67 LqvelpetYPdeaPeie 83
#                  + +  +  YPd+   ++
#  temp_gene_1 136 MLIYKTRYYPDTKCLVK 152
#                  55556677888887765 PP






#First parse what comes from the ls mode matches. Every match in that case is taken
  my $id;
  my $hid;
  while (my $line = <CPGOUT>) {
    chomp;
    if ($line =~ /Query:\s+(\S+)/) {$id = $1}
    if ($line =~ />>\s+(\S+)/) {$hid = $1} 
#   1 ?   10.5   0.0   2.1e-05       2.1      46      81 ..      23      58 ..      20      88 .. 0.90
    if (my ($evalue, $hstart, $hend, $start, $end) = ($line =~ /^\s+\d+\s+\?\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\.\.\s+(\d+)\s+(\d+)/)) {
      $evalue *= 1;           # force evalue to be a float not a string
      my $score = 0;
      my $percentIdentity = 0;
      my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);

      push @fps,$fp;      
    }
	

  }
  close (CPGOUT); 
  $self->output(\@fps);

}

1;
