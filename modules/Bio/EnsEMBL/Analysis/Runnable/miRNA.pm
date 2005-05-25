
# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::miRNA
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::miRNA

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::miRNA->new
    (
     -queries => \%families,
     -analysis => $self->analysis,
    );
  $self->runnable($runnable);

=head1 DESCRIPTION

Runnable to wrap RNAfold, part of the Vienna RNA package.
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

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

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
  my ($queries,$thresholds) = rearrange(['QUERIES'], @args);
  $self->queries($queries);
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
  my $loop = 0;
  print STDERR "get miRNAs\n"  if $verbose;
  # fetch the coordinates of the mature miRNAs
  $self->get_miRNAs;
  my %queries = %{$self->queries};
  $self->throw("Cannot find query sequences $@\n") unless %queries;
  print STDERR "Run analysis\n" if $verbose;
 FAM:  foreach my $family (keys %queries){
  DAF:foreach my $daf (@{$queries{$family}}){
      # Does the alignment contain the mature sequence?
      my $align = $self->run_analysis($daf);
      next DAF unless (scalar @$align);
      my $filename = $self->write_seq($daf);
      # does the sequence fold ino a hairpin ?
      next DAF unless $self->RNAfold($daf,$filename);
      # Parse the structure
      my $structure = $self->parse_miRNA;
      next DAF unless ($structure);
      $self->display_stuff($daf,$structure,$align) if $verbose;
      # create the gene
      $self->make_gene($daf,$structure,$align);
#      last FAM if $loop == 2;
      $loop++;
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
			      -file => "/ecs2/work2/sw4/miRNA/all_mirnas.embl",
			      -format => "embl"
			     );
  $self->throw("unable to open file /ecs2/work2/sw4/miRNA/all_mirnas.embl\n") 
    unless $seqIO;
  while (my $seq = $seqIO->next_seq){
    $miRNA{$seq->accession} = $seq;
  }
  $self->miRNAs(\%miRNA);
}

=head2 run_analysis

  Title      : run_analysis
  Usage      : my $align = $self->run_analysis($daf);
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
  my @all_mature;
  my @mature_aligns;
  if ($self->get_mature($daf)){
    @all_mature = @{$self->get_mature($daf)}
  } else {
    # ignore if you dont have a mature form identified.
    $self->warn("No mature sequence identified for sequence ".$daf->hseqname." $@\n");
    return 0;
  }
  # can have more than 1 mature sequence per hairpin
  foreach my $mature(@all_mature){
    my $miRNA = $miRNAs{$daf->hseqname};
    my $miRNA_length = $mature->{'end'} - $mature->{'start'};
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
    next unless $temp_align;
    $slice = $daf->feature_Slice;
    my $align = $temp_align->transfer($slice);
    my $align_length = $align->hend - $align->hstart;
    # is the enitre mature sequence present in the alignment?
    next unless ( $align_length >= $miRNA_length );
    # is the mature alignment of high enough percent identity?
    next unless ( $align->percent_id >= 90 );
    push @mature_aligns,$align;
  }
  return \@mature_aligns;
}

=head2 RNAfold

  Title      : RNAfold
  Usage      : $self->RNAfold($dna_align_feature,$filename);
  Function   : Wrapper for RNAfold
  Returns    : none
  Exceptions : Throws if RNAfold fails to run
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature, Scalar (filename)

=cut

sub RNAfold{
  my ($self,$daf,$filename)=@_;
  my $command  = "/usr/local/ensembl/bin/RNAfold < ";
  my $options ="";
  my $results_file = $self->create_filename("RNAfold","txt");
  $self->files_to_delete($results_file);
  # delete the postcript file that RNAfold generates
  $self->files_to_delete("/tmp/".substr($daf->hseqname,0,7)."_ss.ps");
  $self->resultsfile($results_file);
  $command .= " $filename ";
  $command .= "$options 2>&1 > ".$results_file;
  print STDERR "Running RNAfold ".$command."\n";
  open(my $fh, "$command |") || 
    throw("Error opening RNAfold cmd <$command>." .
	  " Returned error $? RNAfold EXIT: '" . 
	  ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) . 
	  "', There was " . ($? & 128 ? 'a' : 'no') . 
	  " core dump");
  # this loop reads the STDERR from the blast command
  # checking for FATAL: messages (wublast) [what does ncbi blast say?]
    # N.B. using simple die() to make it easier for RunnableDB to parse.
  while(<$fh>){
    if(/FATAL:(.+)/){
      my $match = $1;
      $self->throw("miRNA: RNAfold failed to run: $match$@\n");
    }
  }
  unless(close $fh){
    # checking for failures when closing.
    # we should't get here but if we do then $? is translated 
    #below see man perlvar
    warning("Error running RNAfold cmd <$command>. Returned ".
	    "error $? BLAST EXIT: '" . ($? >> 8) . 
	    "', SIGNAL '" . ($? & 127) . "', There was " . 
	    ($? & 128 ? 'a' : 'no') . " core dump");
    die ($self->unknown_error_string."\n"); 
  }
  return 1;
}

=head2 parse_miRNA

  Title      : parse_miRNA
  Usage      : my $structure = $self->parse_miRNA;
  Function   : Parses the results of RNAfold
  Returns    : Scalar (string)
  Exceptions : Throws if results file cannot be opened or closed
  Args       : None

=cut

sub parse_miRNA{
  my ($self)=@_;  
  my $results = $self->resultsfile;
  my $structure;
  my $score;
  open(RNAFOLD, $results) or $self->throw("FAILED to open ".$results.
					  " miRNA:parse_results\n$@\n");
 LINE: while(<RNAFOLD>){
    chomp;
    if ($_ =~ /([().]+)\s\((-\d+.\d+)\)$/){
      $structure = $1;
      $score = $2;
    }
  }
  close(RNAFOLD) or $self->throw("FAILED to close ".$results.
				 " miRNA:parse_results\n$@\n");
  if ($score && $score < -20){
    return $structure;
  } else {
    return 0;
  }
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
  $self->throw("Unable to locate miRNA fom embl file that corresponds to ".
	       $query->hseqname."$@ \n") unless $miRNA;
  my @feats = $miRNA->can('top_SeqFeatures') ? $miRNA->top_SeqFeatures : ();
  foreach my $sf ( @feats ) {
    my @fth = Bio::SeqIO::FTHelper::from_SeqFeature($sf,$miRNA);
    foreach my $fth ( @fth ) {
      if( ! $fth->isa('Bio::SeqIO::FTHelper') ) {
	$self->throw("Cannot process FTHelper... $fth $@\n");
      }
      my $location =  $fth->loc;
      if ($location =~ /(\d+)\.+(\d+)/){
	push @mature,{ 'start' => $1,
		       'end'   => $2
		     };
      } else {
	$self->throw("Cannot parse mature coordinates\n");
      }
    }
  }
  return \@mature;
}


=head2 make_gene

  Title      : make_gene
  Usage      : my $gene = $runnable->make_gene($start,$end,$strand,$dna_align_feature,$structure_line,$simple_alignment)
  Function   : Creates the non coding gene object from the parsed result file.
  Returns    : Hashref of Bio::EnsEMBL::Gene and Bio::EnsEMBL::Attribute
  Exceptions : None
  Args       : Scalar,Scalar,Scalar,Bio::EnsEMBL::DnaDnaAlignFeature,Scalar,Bio::SimpleAlign

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
     -phase => 0,
     -end_phase => (($daf->end - $daf->start + 1)%3)
    );

  $exon->add_supporting_features($daf);
  # transcripts
  my $transcript = Bio::EnsEMBL::Transcript->new;
  $transcript->add_Exon($exon);
  $transcript->start_Exon($exon);
  $transcript->end_Exon($exon);
  # add the transcript attributes for the position of the mature miRNA as well as
  # the predicted structure
  my @codes = @{$self->encode_str($structure)};
  foreach my $code(@codes){
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
  $gene->type('miRNA');
  $gene->description($description);
  $gene->analysis($self->analysis);
  $gene->add_Transcript($transcript);
  $gene_hash{'gene'} = $gene;
  $gene_hash{'attrib'} = \@attributes;
  $self->output(\%gene_hash);
}

=head2 write_seq

  Title      : write_seq
  Usage      : my $filename = $self->write_seq($daf);
  Function   : Writes the dna sequence file of the region covered by the DnaDnaAlignFeature
  Returns    : Bio::EnsEMBL::Gene
  Exceptions : Throws if it cannot write to the file
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature 

=cut

sub write_seq{
  my ($self,$daf)=@_;
  my $filename = $self->create_filename("miRNA","seq");
  # have to write file so the sequence is all on a single line 
  # cos thats the way RNAfold likes it
  my $seq = Bio::Seq->new(
			  -display_id => "/tmp/".$daf->hseqname,
			  -seq        => $daf->seq,
			 );
  $self->files_to_delete("/tmp/".$daf->hseqname);
  my $seqIO = Bio::SeqIO->new(
			     -file => ">$filename",
			     -format   => 'fasta',
			     -width    => length($daf->seq),
			     );
  $seqIO->write_seq($seq)or $self->throw
    ("FAILED to write to  $filename miRNA:run:write_seq $@\n");
  $self->files_to_delete($filename);
  return $filename;
}


=head2 encode_str

  Title      : encode_str
  Usage      : my $encoded_str = $runnable->encode_string($string)
  Function   : Does string length encoding to reduce size of structure string
             : splits strings if they are longer then 200 characters so they 
             : will fit in the transcript attribute table, gives a range value
             : at the start of the string indicating the start and stop positions of the 
             : structure on the transcript
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub encode_str{
  my ($self,$string)= @_;
  my @codes;
  my $start = 1;
  my $count=0;
  my $code;
  my @elements = split //,$string;
  my $last_chr = "";
  my @array =[];
  foreach my $chr (@elements){
    $count++;
    if ($chr eq $last_chr){
	push @array,$chr;
      }
    else {
      if ($code && length($code) > 200 && scalar(@array) == 1){
	push @codes,"$start:$count\t$code";
	$code = undef;
	$start = $count+1;
      }
      # Character has changed print STDERR it and the associated array length
      if (scalar(@array) > 1){
	$code .= scalar(@array);
	@array = [];
      }
      $code .= "$chr";
      $last_chr = $chr;
    }
  }
# last element
  if (scalar(@array) > 1){
    $code .= scalar(@array);
  }
  push @codes,"$start:$count\t$code";
  return \@codes;
}

sub display_stuff{
  my ($self,$daf,$structure,$aligns)=@_;
  my @all_mature;
  @all_mature = @{$self->get_mature($daf)};
  foreach my $mature (@all_mature){
    print STDERR $daf->hseqname." miRNA at ".$mature->{'start'}." ".$mature->{'end'}."\n ";
  }
  print STDERR "target strand ".$daf->strand."\n";
  print STDERR "Query start end ".$daf->hstart.":".$daf->hend.
    " strand ".$daf->hstrand."\n";
  print STDERR "Start ".$daf->start." end ".$daf->end.
    " Hit end ".$daf->hstart." hit end ".$daf->hend."\n";
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

sub  queries {
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


sub  miRNAs {
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
1;
