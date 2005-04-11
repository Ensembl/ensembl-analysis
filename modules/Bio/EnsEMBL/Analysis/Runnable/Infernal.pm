# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Blast
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::Infernal

=head1 SYNOPSIS

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Infernal->new
            (
            '-queries'  => \@array_ref_of_dna_align_features,
            '-analysis' => $self->analysis,
            );
    $runnable->run;
    $output = $runnable->output;

=head1 DESCRIPTION

Runnable for Infernal (Runs ncRNA analysis on blast hits).
Wraps cmsearch, part of the Infernal suite of programs by Sean Eddy.
Parses results to build non-coding gene objects and a representation
of secondary structure which is string length encoded and stored as a 
transcript attribute (currently)

=head1 CONTACT

Simon White

sw4@sanger.ac.uk

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::Infernal;

use strict;
use warnings;


use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Gene;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

  Title      : new
  Usage      : my $runnable = Bio::EnsEMBL::Analysis::Runnable::Infernal->new
  Function   : Instantiates new Infernal runnable
  Returns    : Infernal object
  Exceptions : none
  Args       : Array ref of dna align features

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($queries,$thresholds) = rearrange(['QUERIES'], @args);
  $self->queries($queries);
  $self->get_thresholds;
  $self->get_descriptions;
  return $self;
}

=head2 run

  Title      : run
  Usage      : my $runnable->run;
  Function   : Runs the runnable
  Returns    : None
  Exceptions : none
  Args       : None

=cut

sub run{
  my ($self) = @_;
  my $queries = $self->queries;
  foreach my $query (@$queries){
    $self->query($self->get_sequence($query));
    throw("Can't run ".$self." without a query sequence") 
      unless($self->query);
    $self->checkdir();
    my $filename = $self->write_seq_file();
    $self->files_to_delete($filename);
    $self->files_to_delete($self->resultsfile);
    $self->run_analysis();
    $self->parse_results($query);
  }
  $self->delete_files;
  return 1;
}

=head2 run_analysis

  Title      : run_analysis
  Usage      : $runnable->run_analysis;
  Function   : Runs the analysis (cmsearch)
  Returns    : None
  Exceptions : Throws in the event of cmsearch returning an error code
  Args       : None

=cut

sub run_analysis{
  my ($self) = @_;
  my $db =  $self->analysis->db_file.$self->query->display_id.".cm";
  my $command  =  '/nfs/acari/sw4/local/bin/cmsearch';
  my %thresholds = %{$self->thresholds};
  my $domain = $self->query->display_id;
  my $filename = $self->queryfile;
  my $results_file = $self->create_filename("Infernal.$domain.out");
  $self->files_to_delete($results_file);
  $self->resultsfile($results_file);
  $command .= " -W ".$thresholds{$domain}{'win'};
  if($thresholds{$domain}{'mode'} =~ /local/) {
	$command .= " --local";
    }
  $command .= " $db  $filename > $results_file";
  print "Running infernal ".$command."\n";
  open(my $fh, "$command |") || 
    $self->throw("Error opening Infernal cmd <$command>." .
	  " Returned error $? Infernal EXIT: '" .
	  ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
	  "', There was " . ($? & 128 ? 'a' : 'no') .
	  " core dump");
}

=head2 parse_results

  Title      : parse_results
  Usage      : $runnable->parse_results(\@dna_align_features)
  Function   : Parses cmsearch output: rejects alignments where score is below 
             : the threshold defined for the RFAM familly in question.
             : The thresholds are defined in the file /data/blastdb/Rfam/Rfam.thr.
             : A dna_align_feature representing the initial RFAM blast hit is used to determine 
             : the familly and position on the genome.
  Returns    : None
  Exceptions : Throws if it cannot find, open or close the cmsearch output
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

  # hit 1   :    201    266    40.08 bits
sub parse_results{
  my ($self, $daf) = @_;
  my @dafs;
  my $align;
  my $line =0;
  my $domain = substr($daf->hseqname,0,7);
  my %thresholds = %{$self->thresholds};
  my $threshold = $thresholds{$domain}{'thr'};
  print STDERR "Domain $domain has threshold $threshold\n";
  my $results = $self->resultsfile;
  my $ff = $self->feature_factory;
  if(!-e $results){
    $self->throw("Can't parse an no existance results file ".$results.
          " Infernal:parse_results");
  }
  my @output;
  my ($hit,$start,$end,$score,$strand,$str) =0;
  open(INFERNAL, $results) or $self->throw("FAILED to open ".$results.
				    " INFERNAL:parse_results");
 LINE: while(<INFERNAL>){
    chomp;
    if ($_ =~ /^hit/){
      if ($score > $threshold && $align->{'name'}){
	$str = $self->parse_structure($align);
	push @output, $self->make_gene($start,$end,$strand,$daf,$str);
      }	
      $align = {};
      $align->{'name'} = $domain ;
      $line = -1;
    }
    if ($_ =~ /^CPU/){
      $line = -1;
    }
    print STDERR "$line\t$_\n";
#    print "$_\n";
    if ($score && $score >= $threshold && $line >= 0){
      # parsing goes in 5 lines
      # split the lines and store each element in a anonymous hash array
      $line++;
      if ($_){
	push @{$align->{'str'}}   ,  split//,substr($_,11,60) if ($line == 1);
	push @{$align->{'target'}},  split//,substr($_,11,-9) if ($line == 2);
	push @{$align->{'match'}} ,  split//,substr($_,11,60) if ($line == 3);
	push @{$align->{'query'}} ,  split//,substr($_,11,-9) if ($line == 4);
      }
      else{
	$line =0;
      }
    }
    if ($_ =~ /^hit\s+(\d+)\s+:\s+(\d+)\s+(\d+)\s+(\d+).(\d+)\s+bits/){
      $hit = $1;
      $start = $2;
      $end = $3;
      $score = $4+$5/100;
      print "hit - $hit\nscore - $score\n";
      if ($score >= $threshold){ 
	if ($end < $start){
	  $strand = -1;
	  my $temp = $end;
	  $end = $start;
	  $start = $temp;
	}
	else {
	  $strand =1;
	}
	$line = 0;
      }
      else{
	$line = -1;
      }
    }
  }
  if ($align->{'name'} && $score > $threshold){
    $str = $self->parse_structure($align);
    print "$str\n";
    push @output, $self->make_gene($start,$end,$strand,$daf,$str);
  }
  $self->output(\@output);
  close(INFERNAL) or $self->throw("FAILED to close ".$results.
			   " INFERNAL:parse_results");
}

=head2  parse_structure

  Title      : parse_structure
  Usage      : my $structure = $runnable->parse_structure(\%hash_ref_containing_parsed_cmsearch_results)
  Function   : Parses cmsearch alignment to create a structure line which represents the
             : predicted secondary structure of the ungaped sequence
  Returns    : Bio::EnsEMBL::Attribute object representing a string length encoded form of the structure line
  Exceptions : Throws if it cannot find, open or close  the cmsearch output
  Args       : None

=cut

sub parse_structure{
  my ($self,$align)=@_;
  my @fwd_matches;
  my @stack;
  my @big_gaps;
  my $matchstring;
  my @matches;
  my $big_gap=0;
  print $align->{'name'}."\n";  
  print $align->{'str'}."\n";
  #Brace matching
  # sequence postions start at 1 not 0 so you dont get false booleans for real sequence locations, if you see what I mean
  for (my $i=0 ; $i< scalar(@{$align->{'str'}}); $i++){
    if ($align->{'str'}[$i] eq '(' or
	$align->{'str'}[$i] eq '<' or
	$align->{'str'}[$i] eq '[' or
	$align->{'str'}[$i] eq '{'){
      push @stack,$i+1;
    }
    if ($align->{'str'}[$i] eq ')' or
	$align->{'str'}[$i] eq '}' or
	$align->{'str'}[$i] eq ']' or
	$align->{'str'}[$i] eq '>'){
      $fwd_matches[$i+1] = pop @stack;
    }
  }
  for (my $i=0 ; $i< scalar(@{$align->{'str'}}); $i++){
    # Parse out large gaps by looking for ~ on the str line;
    if ($align->{'query'}[$i] eq '*'){
      my $string;
      for (my $j=$i+1 ; $j <= $i+10 ; $j++){	
	last if ($align->{'query'}[$j] eq '*');
	$string .= $align->{'query'}[$j];
      }
      if ($string =~ /\[(.+)\]/){
	my $gap = $1;
	$gap =~ s/\D//g;
	$big_gaps[$i] = $gap;
      }
    }
    #skip over if you have a gap - beware gaps at the other end of the alignment;
    if ($align->{'query'}[$i] eq '-'){
      next;
    }
    # skip over if you have a missmatch
    if ($align->{'match'}[$i] eq ' '){
      $matches[$i+1] = '.';
      next;
    }
    if ($fwd_matches[$i+1]){
      # check there isnt a gap at the other end
      if ($align->{'query'}[$fwd_matches[$i+1]-1] eq '-'){
	$matches[$i+1] = $align->{'str'}[$i];
	next;
      }
      $matches[$i+1] = $align->{'str'}[$i];
      $matches[$fwd_matches[$i+1]] = $align->{'str'}[$fwd_matches[$i+1]-1];
    }
    else {
      $matches[$i+1] = $align->{'str'}[$i];
    }
  }
  for (my $i=0 ; $i< scalar(@{$align->{'str'}}); $i++){
    if ($big_gaps[$i]){
      for (my $j = 0 ; $j < $big_gaps[$i] ; $j++){
	$matchstring .= ".";
      }
      $i=$i+5;
      next;
    }
    if ($matches[$i+1]){
      $matchstring.= $matches[$i+1];
    }
  }
  my $attribute = Bio::EnsEMBL::Attribute->new
    (-CODE => 'ncRNA',
     -NAME => 'Structure',
     -DESCRIPTION => 'RNA secondary structure line',
     -VALUE => $self->encode_string($matchstring)
    );
  print "$matchstring\n";

  return $attribute;
}

=head2 encode_string

  Title      : encode_string
  Usage      : my $encoded_string = $runnable->encode_string($string)
  Function   : Does simple string length encoding to reduce size of structure sting
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub encode_string{
  my ($self,$string)= @_;
  my $code;
  my @elements = split //,$string;
  my $last_chr;
  my @array =[];
  foreach my $chr (@elements){
    if ($chr eq $last_chr){
	push @array,$chr;
      }
    else {
      if (scalar(@array) > 1){
	$code .= scalar(@array);
	@array = [];
      }
      $code .= "$chr";
    }
    $last_chr = $chr;
  }
# last element
  if (scalar(@array) > 1){
    $code .= scalar(@array);
  }
  return $code;
}

=head2 decode_string

  Title      : decode_string
  Usage      : my $decoded_string = $runnable->decode_string($string)
  Function   : Does simple string length decoding.
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub decode_string{
  my ($self,$string)= @_;
  my $code;
  my $offset = 0;
  my $chr;
  my $num = 0;
  my $str = "";
  for (my $pos =0; $pos <= length($string) ; $pos++){
    $chr =  substr($string,$pos,1);
    if ($chr =~ /\D/){
      print $str unless($num);
      if ($num){
	for (my$i =1 ; $i <= $num ; $i++){
	  print $str;
	}
      }
      $str = $chr;
      $num = "";
      next;
    }
    if ($chr =~ /\d/){
      $num .= $chr;
    }
  }
  # empty array 
  if ($num){
    for (my$i =1 ; $i <= $num ; $i++){
      print $str;
    }
  }
  else{
    print $str;
  }
}

=head2 make_gene

  Title      : make_gene
  Usage      : my $gene = $runnable->make_gene($start,$end,$strand,$dna_align_feature,$structure_line)
  Function   : Creates the non coding gene object from the parsed result file.
  Returns    : Bio::EnsEMBL::Gene
  Exceptions : None
  Args       : scalar ($start,$end,$strand,$structure_line) and
             : Bio::EnsEMBL::DnaDnaAlignFeature ($dna_align_feature)

=cut

sub make_gene{
  my ($self,$start,$end,$strand,$daf,$str) = @_;
  my $domain = substr($daf->hseqname,0,7);
  my %descriptions = %{$self->descriptions};
  # exons
  my $slice = $daf->feature_Slice;
  my $exon = Bio::EnsEMBL::Exon->new
    (
     -start => $start,
     -end   => $end,
     -strand => $strand,
     -slice => $slice,
     -phase => 0,
     -end_phase => (($end - $start + 1)%3)
    );

  $exon->add_supporting_features($daf->transfer($slice));
  # transcripts
  my $transcript = Bio::EnsEMBL::Transcript->new;
  $transcript->add_Exon($exon);
  $transcript->start_Exon($exon);
  $transcript->end_Exon($exon);
  # dodgy
  $transcript->{'_structure'} = $str;
  # gene
  my $gene = Bio::EnsEMBL::Gene->new;
  $gene->type('ncRNA');
  $gene->description($descriptions{$domain});
  print "Rfam_id $domain ".$descriptions{$domain}."\n";
  $gene->analysis($self->analysis);
  $gene->add_Transcript($transcript);
  return $gene;
}

=head2 get_sequence

  Title      : get_sequence
  Usage      : my $seq = $runnable->get_sequence($dna_align_feature)
  Function   : Makes a BioPerl seq object containing the sequecnce for cmsearch to run on.
             : Uses the original RFAM blast hit represented by a dna alignment feature and adds 
             : 5' and 3' flanking sequence based on predefined values set for each RFAM familly
             : in the /data/blastdb/Rfam/Rfam.thr file
  Returns    : Bio::Seq object
  Exceptions : None
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub get_sequence{
  my ($self,$daf)=@_;
  my $domain = substr($daf->hseqname,0,7);
  my %thresholds = %{$self->thresholds};
  my $padding = $thresholds{$domain}{'win'};
  print STDERR "Using padding of $padding\t";
  $daf->start($daf->start - $padding);
  $daf->end($daf->end + $padding);
  my $seq = Bio::Seq->new(
			    -seq => $daf->seq,
			    -display_id => $domain,
			   );
  print STDERR " total seq length = ";
  print $seq->length."\n";
  return $seq;
}

=head2 get_thresholds

  Title      : get_thresholds
  Usage      : $runnable->get_thresholds
  Function   : Fetches and stores predefined thresholds set for each RFAM familly
             : in the /data/blastdb/Rfam/Rfam.thr file
  Returns    : None
  Exceptions : Throws if it cannot open the thresholds file
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub get_thresholds{
  my($self)= @_;
  # read threshold file
  my %thr;
  open( T, "/data/blastdb/Rfam/Rfam.thr" ) or $self->throw("can't file the Rfam.thr file");
  while(<T>) {
    if( /^(RF\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s*$/ ) {
      $thr{ $1 } = { 'id' => $2, 'thr' => $3, 'win' => $4, 'mode' => $5 };
    }
  }
  close T;
  $self->thresholds(\%thr);
}


# LOOK AT THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# I should get the descriptions file pushed out across the farm along with
# the other files, In fact maybe sam should do it as he periodically updates
# these files anyway

=head2 get_descriptions

  Title      : get_descriptions
  Usage      : $runnable->get_descriptions
  Function   : Fetches and stores predefined descriptions for each RFAM familly
             : in the /ecs2/work2/sw4/RFAM/Rfam.descriptions file
  Returns    : none
  Exceptions : Throws if it cannot open the thresholds file
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub get_descriptions{
  my($self)= @_;
  # read descriptions file
  my %desc;
  open( T, "/ecs2/work2/sw4/RFAM/Rfam.descriptions" ) or 
    $self->throw("can't file the Rfam.descriptions file");
  while(<T>) {
    if( /^(RF\d+)\s+\S+\s+(.+)$/ ) {
      $desc{ $1 } = $2;
    }
  }
  close T;
  $self->descriptions(\%desc);
}

#################################################################
# Containers

=head2 queries

  Title      : queries
  Usage      : my %queries = %$runnable->queries
  Function   : Get/ set for the dna alignfeatures defining as the query sequences for cmsearch
  Returns    : Hash reference
  Exceptions : None
  Args       : Array reference of Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub  queries {
  my ($self, $queries) = @_;
  if ($queries){
    $self->{'_queries'} = $queries;
  }
  return $self->{'_queries'};
}

=head2 thresholds

  Title      : thresholds
  Usage      : my %thresholds = %$runnable->thresholds
  Function   : Get/ set for the cmsearch thresholds for each RFAM familly
  Returns    : Hash reference
  Exceptions : None
  Args       : Hash Reference

=cut

sub  thresholds {
  my ($self, $thresholds) = @_;
  if ($thresholds){
    $self->{'_thresholds'} = $thresholds;
  }
  return $self->{'_thresholds'};
}

=head2 descriptions

  Title      : descriptions
  Usage      : my %descriptions = %$runnable->descriptions
  Function   : Get/ set descriptions for each RFAM familly
  Returns    : Hash reference
  Exceptions : None
  Args       : Hash Reference

=cut

sub  descriptions {
  my ($self, $descriptions) = @_;
  if ($descriptions){
    $self->{'_descriptions'} = $descriptions;
  }
  return $self->{'_descriptions'};
}

1;
