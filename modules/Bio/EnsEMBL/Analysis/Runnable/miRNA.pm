=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::miRNA - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::miRNA->new
    (
     -queries => \%families,
     -analysis => $self->analysis,
    );
  $self->runnable($runnable);

=head1 DESCRIPTION

Takes blast alignments in the form of a  hash ref of Bio::EnsEMBL::DnaDnaAlignFeatures
 grouped by miRNA family.
Opens an EMBL format file of miRNA sequences used to create the blast database and
parses out the positions of the mature RNA for each dna align feature.
Determines if the mature sequence is wholly contained within the dna align feature.
Applies a percent identity filter to the blast alignments.
If the alignment meet these criteria it then runs RNAfold and parses the results.
If the results of RNAfold indicate that the sequence can form a hairpin structure with 
a score < -20 it is classed as a miRNA.
Single exon gene objects are made and the  DnaDnaAlignFeature is added as supporting evidence.
The predicted structure and coordinates of the mature sequence are added to the transcript as
transcript attributes
=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::miRNA;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis::Runnable::RNAFold;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

my $verbose = "yes";

=head2 new

  Title      : new
  Usage      :   my $runnable = Bio::EnsEMBL::Analysis::Runnable::miRNA->new
             :    (
             :     -queries => \%families,
             :     -analysis => $self->analysis,
             :    );
  Function   : Instantiates new miRNA runnable
  Returns    : Bio::EnsEMBL::Analysis::Runnable::miRNA
  Exceptions : none
  Args       : Hash ref of Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
#  my ($queries,$thresholds) = rearrange(['QUERIES'], @args);
#  $self->queries($queries); 
#  $self->throw("miRNA: dying because cannot find database".$self->analysis->db_file."\n")
#    unless (-e $self->analysis->db_file);
  
  my ($outdir, $queries,$thresholds) = rearrange(['OUTPUT_DIR', 'QUERIES'], @args);
  $self->queries($queries); 
  $self->throw("miRNA: dying because cannot find database ".$self->analysis->db_file."\n")
    unless (-e $self->analysis->db_file);
  $self->outdir($outdir);

  return $self;
}

=head2 new

  Title      : run
  Usage      : my $runnable->run;
  Function   : Run method.
  Returns    : None
  Exceptions : Throws if no query sequences found
  Args       : None

=cut

sub run{
  my ($self) = @_;
  print STDERR "get miRNAs\n"  if $verbose;
  # fetch the coordinates of the mature miRNAs
  $self->get_miRNAs;
  my %queries = %{$self->queries};
  $self->throw("Cannot find query sequences $@\n") unless %queries;
  print STDERR "Run analysis\n" if $verbose;
 FAM:  foreach my $family (keys %queries){
  DAF:foreach my $daf (@{$queries{$family}}){    
      # Does the alignment contain the mature sequence?
      my ($align,$status) = $self->run_analysis($daf);
      next DAF unless ($align && scalar @$align > 0);
      # does the sequence fold into a hairpin ?
      my $seq = Bio::PrimarySeq->new
	(
	 -display_id => $daf->hseqname,
	 -seq => $daf->seq,
	);
      my $RNAfold = Bio::EnsEMBL::Analysis::Runnable::RNAFold->new
	(
	 -analysis  => $self->analysis,
	 -sequence  => $seq
	);
      $RNAfold->run;
      # get the final structure encoded by run length
      my $structure = $RNAfold->encoded_str;
      # next DAF unless ($RNAfold->score < -20); # uncommented to report all MFEs from RNAFold, hard threshold no longer needed - osagie 10/2017
      next DAF unless ($RNAfold->structure);
      $self->display_stuff($daf,$RNAfold->structure,$align,$RNAfold->score) if $verbose;
      # create the gene
      $self->make_gene($daf,$structure,$align);
    }
  }
  print STDERR "delete temp files\n"  if $verbose;
  $self->delete_files;
}

=head2 get_miRNAs

  Title      : get_miRNAs
  Usage      : $self->get_miRNAs;
  Function   : Opens and reads EMBL formatted file of miRNA sequences used to make 
             : the miRNA blast database
  Returns    : None
  Exceptions : Throws if sequence file is not found / cannot be opened
  Args       : None

=cut

sub get_miRNAs{
  my ($self)=@_;
  my %miRNA;
  my $seqIO = Bio::SeqIO->new(
			      -file   => $self->analysis->db_file,
			      -format => "embl"
			     );
  $self->throw("unable to open file ".$self->analysis->db_file."\n")
    unless $seqIO;
  while (my $seq = $seqIO->next_seq){
    $miRNA{$seq->accession} = $seq;
  }
  $self->miRNAs(\%miRNA);
}

=head2 run_analysis

  Title      : run_analysis
  Usage      : my $align = $self->run_analysis($dna_align_feature);
  Function   : makes an alignment that represents the mature sequence of the miRNA
             : tests to check that it contains the entire mature sequence and  is of 
             : high enough percent identity
  Returns    : Bio::SimpleAlign object
  Exceptions : Throws if sequence file is not found / cannot be opened
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub run_analysis{
  my ($self,$daf)=@_;
  my %miRNAs = %{$self->miRNAs};
  my ($all_mature, $status);
  my @mature_aligns;
  if ($self->get_mature($daf)){
    ($all_mature,$status) = $self->get_mature($daf)
  } else {
    # ignore if you dont have a mature form identified.
    $self->warning("No mature sequence identified for sequence ".$daf->hseqname." $@");
    return 0;
  }
  # can have more than 1 mature sequence per hairpin
  foreach my $mature(@$all_mature){
    my $miRNA = $miRNAs{$daf->hseqname};
    # length is 1 longer than end - start because it includes both the start and end base
    my $miRNA_length = $mature->{'end'} - $mature->{'start'} +1;
    # change the U to T so as not to confuse the BioPerl
    my $seq = $miRNA->seq;
    $seq =~ s/[uU]/T/g;
    # make a fake slice to use as the hit sequence
    my $slice = Bio::EnsEMBL::Slice->new
      (
       -seq_region_name => $daf->hseqname,
       -start           => 1,
       -end             => $miRNA->length,
       -seq             => $seq,
       -strand          => 1,
       -coord_system    => $daf->slice->coord_system,
      );
    $daf->hslice($slice);
    # make an alignment of just the mature sequence
    my $temp_align = $daf->restrict_between_positions($mature->{'start'},$mature->{'end'},"HSEQ");
    unless ($temp_align)
      {
#       print "rejecting ".$daf->p_value." cos of not getting mature alignment\n";
       next;
      }
    my $align = $temp_align->transfer($daf->feature_Slice);
    # is the enitre mature sequence present in the alignment, no indels?
    unless ( $align->alignment_length == $miRNA_length && $align->alignment_length == $align->length)
      {
#	print "rejecting ".$daf->dbID." coz of not getting aligned length\n";
	next;
      }
    # is the mature alignment of high enough percent identity?
    # also no gaps in the mature alignment
    unless ( $align->get_SimpleAlign->overall_percentage_identity >= 80 ) 
      {
#	print "rejecting ".$daf->p_value." cos of aligned ID\n";
	next;
      }
    push @mature_aligns,$align;
  }
  return \@mature_aligns,$status;
}


=head2 get_mature

  Title      : get_mature
  Usage      : my @mature = $self->get_mature($dna_align_feature);
  Function   : Returns the coordinates of the mature miRNA on the target sequence
  Returns    : Array containing 2 elements, the start and stop position (inclusive)
  Exceptions : Throws if coordinates cannot be parsed or the miRNA cannot be found
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub get_mature{
  my ($self,$query)=@_;
  my %miRNAs = %{$self->miRNAs};
  my $miRNA = $miRNAs{$query->hseqname};
  my @mature;
  my $status = "PUTATIVE";
  unless ($miRNA){
  print STDERR "Unable to locate miRNA fom embl file that corresponds to ".
	       $query->hseqname."$@ \n";
    return 0;
    }       
  my @feats = $miRNA->can('top_SeqFeatures') ? $miRNA->top_SeqFeatures : ();
  foreach my $sf ( @feats ) {
    my @fth = Bio::SeqIO::FTHelper::from_SeqFeature($sf);
    foreach my $fth ( @fth ) {
      if( ! $fth->isa('Bio::SeqIO::FTHelper') ) {
	$self->throw("Cannot process FTHelper... $fth $@\n");
      }
      my $location =  $fth->loc;
      $status = "KNOWN" 
	if ($fth->field->{"evidence"}[0] && $fth->field->{"evidence"}[0] eq "experimental");
      if ($location =~ /(\d+)\.+(\d+)/){
	push @mature,{ 'start' => $1,
		       'end'   => $2
		     };
      } else {
	return 0;
      }
    }
  }
  return \@mature,$status;
}


=head2 make_gene

  Title      : make_gene
  Usage      : my $gene = $runnable->make_gene($dna_align_feature,$structure,$simple_alignment)
  Function   : Creates the non coding gene object from the parsed result file.
  Returns    : Hashref of Bio::EnsEMBL::Gene, Bio::EnsEMBL::Attribute and Bio::EnsEMBL::DBEntry
  Exceptions : None
  Args       : dna_align_feature (Bio::EnsEMBL::DnaDnaAlignFeature)
             : structure (String)
             : alignment (Bio::SimpleAlign)

=cut

sub make_gene{
  my ($self,$daf,$structure,$aligns) = @_;
  my %miRNAs = %{$self->miRNAs};
  my @mature;
  my %gene_hash;
  my $description = $miRNAs{$daf->hseqname}->display_id;
  my @attributes;
  # exons
  my $slice = $daf->slice;
  my $exon = Bio::EnsEMBL::Exon->new
    (
     -start => $daf->start,
     -end   => $daf->end,
     -strand => $daf->strand,
     -slice => $slice,
     -phase => -1,
     -end_phase => -1
    );
  # dont let it fall of the slice 
  return if ($exon->start < 1 or $exon->end > $slice->length);
  $exon->add_supporting_features($daf);
  # transcripts
  my $transcript = Bio::EnsEMBL::Transcript->new;
  $transcript->add_Exon($exon);
  $transcript->start_Exon($exon);
  $transcript->end_Exon($exon);
  $transcript->biotype("miRNA");
  $transcript->source("ensembl");
  # add the transcript attributes for the position of the mature miRNA as well as
  # the predicted structure
  foreach my $code(@$structure){
    my $str_attrib = Bio::EnsEMBL::Attribute->new
      (-CODE => 'ncRNA',
       -NAME => 'Structure',
       -DESCRIPTION => 'RNA secondary structure line',
       -VALUE => $code
      );
    push @attributes,$str_attrib;
  }
  foreach my $align(@$aligns){
    my $miRNA_attrib = Bio::EnsEMBL::Attribute->new
      (-CODE => 'miRNA',
       -NAME => 'Micro RNA',
       -DESCRIPTION => 'Coordinates of the mature miRNA',
       -VALUE => $align->start."-".$align->end,
      );
    push @attributes,$miRNA_attrib;
  }
  # gene
  my $gene = Bio::EnsEMBL::Gene->new;
  $gene->biotype('miRNA');
  $gene->source("ensembl");
  $gene->analysis($self->analysis);
  $gene->add_Transcript($transcript);
  $gene_hash{'gene'} = $gene;
  $gene_hash{'attrib'} = \@attributes;
  $self->output(\%gene_hash);
}

sub display_stuff{
  my ($self,$daf,$structure,$aligns,$score)=@_;
  my $all_mature;
  my $status;
  my %miRNAs = %{$self->miRNAs};
  my $description = $miRNAs{$daf->hseqname}->display_id;
  ($all_mature,$status) = $self->get_mature($daf);
  print STDERR "DAF ".$daf->dbID." ".$daf->seq_region_name." ".$daf->seq_region_start." ".$daf->seq_region_end."\n";
  foreach my $mature (@$all_mature){
    print STDERR $daf->hseqname." $description miRNA at ".$mature->{'start'}." ".$mature->{'end'}."\n ";
  }
  print STDERR "target strand ".$daf->strand."\n";
  print STDERR "Query start end ".$daf->hstart.":".$daf->hend.
    " strand ".$daf->hstrand."\n";
  print STDERR "Start ".$daf->start." end ".$daf->end.
    " Hit end ".$daf->hstart." hit end ".$daf->hend."\n";
  print STDERR "Score $score\n";
  foreach my $align(@{$aligns}){
    my $simple = $align->get_SimpleAlign;
    foreach my $seq ( $simple->each_seq() ) {
      print STDERR $seq->seq."\n";
    }
    print STDERR $simple->match_line."\n";
    print STDERR "\n";
    $simple = $daf->get_SimpleAlign;
    foreach my $seq ( $simple->each_seq() ) {
      print STDERR $seq->seq."\n";
    }
    print STDERR $simple->match_line."\n";
    print STDERR "\n";
    print STDERR "$structure\n";

    for (my $i=1 ; $i< $align->start;$i++) {
      print STDERR ".";
    }
    print STDERR substr($daf->seq,$align->start-1,$align->end-$align->start+1);
    print STDERR "\nMirna position = ".$align->start." ". $align->end."\n";
  }

  my $fn = $self->outdir . "/rna_fold_results.txt";
  open(FH, '>>', $fn) or die "Could not write to $fn ; please check basedir exists";
  print FH $daf->seq_region_name . "\t" . $daf->seq_region_start . "\t" . $daf->seq_region_end . "\t" .
    $daf->seq_region_name . ":" . $daf->seq_region_start . "-" . $daf->seq_region_end . "\t" . $score . "\t" .
    ($daf->strand > 0 ? "+" : "-") . "\t" . $daf->dbID . "\n";
}

##################################################################################
# Containers

=head2 queries

  Title      : queries
  Usage      : my %queries = %$runnable->queries
  Function   : 
  Returns    : Hash reference
  Exceptions : None
  Args       : 

=cut

sub queries {
  my ($self, $queries) = @_;
  if ($queries){
    $self->{'_queries'} = $queries;
  }
  return $self->{'_queries'};
}

=head2 miRNAs

  Title      : miRNAs
  Usage      : my %miRNAs = %{$self->miRNAs};
  Function   : Container for storing miRNA BioSeq objects
  Returns    : Hash reference to Bio::Seq object
  Exceptions : None
  Args       : Hash reference to Bio::Seq object

=cut


sub miRNAs {
  my ($self, $miRNAs) = @_;
  if ($miRNAs){
    $self->{'_miRNAs'} = $miRNAs;
  }
  return $self->{'_miRNAs'};
}


=head2 output

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : hashref of output
  Function  : overrides runnable output method to allow storing of hash refs
  Exceptions: throws if not passed an hashref
  Example   : 

=cut


sub output{
  my ($self, $output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($output){
    throw("Must pass Runnable:output an hashref not a ".$output)
      unless(ref($output) eq 'HASH');
    push(@{$self->{'output'}}, $output);
  }
  return $self->{'output'};
}

sub outdir {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_outdir'} = $value;
  }
  if (exists($self->{'_outdir'})) {
    return $self->{'_outdir'};
  } else {
    return undef;
  }
}
1;
