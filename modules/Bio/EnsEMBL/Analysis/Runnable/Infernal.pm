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

Bio::EnsEMBL::Analysis::Runnable::Infernal - 

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
transcript attribute

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Infernal;

use strict;
use warnings;


use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis::Runnable::RNAFold;
use IO::File;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

my $verbose = 1;

=head2 new

  Title      : new
  Usage      : my $runnable = Bio::EnsEMBL::Analysis::Runnable::Infernal->new
  Function   : Instantiates new Infernal runnable
  Returns    : Bio::EnsEMBL::Analysis::Runnable::Infernal object
  Exceptions : none
  Args       : Array ref of Bio::EnsEMBL::DnaDnaAlignFeature objects

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($queries,$thresholds) = rearrange(['QUERIES'], @args);
  $self->queries($queries);
#  $self->get_thresholds;
 return $self;
}

=head2 run

  Title      : run
  Usage      : my $runnable->run;
  Function   : Runs the runnable
  Returns    : None
  Exceptions : Throws if no query sequence is provided
  Args       : None

=cut

sub run{
  my ($self) = @_;
  my $queries = $self->queries;
  my @attributes;
  my $descriptions = $self->get_descriptions;
  my $cmfh = IO::File->new( $self->analysis->db_file."/Rfam.cm"  ) or $self->throw("can't file the Rfam.cm file"); 
  my $cm = $self->read_cm_library($cmfh);
  foreach my $query (@$queries){
    $self->query($self->get_sequence($query,$cm));
    $self->throw("Can't run ".$self." without a query sequence") 
      unless($self->query);
    $self->checkdir();
    # write the seq file
    my $filename = $self->write_seq_file();
    # write the cm file
    my $cmfile   = $self->write_cm_file($cm,$cmfh);
    $self->files_to_delete($cm);
    $self->files_to_delete($filename);
    $self->files_to_delete($self->resultsfile);
    # run cmsearch
    $self->run_analysis($cmfile,$cm);
    # parse the cmsearch results file - make a hash containing all 
    # the results ordered by score;
    my $results = $self->parse_results($query,$cm);
    next unless($results);
    # make the gene object out of the highest scoring result that
    # overlaps the blast hit
    my $gene = $self->make_gene($results,$query,$descriptions,$cm);
    $self->output($gene) if $gene;
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
  my ($self,$cmfile, $cm) = @_;
  my $command  =  $self->program;
  my %thresholds;# = %{$self->thresholds};
  my $domain = $self->query->display_id;
  my $filename = $self->queryfile;
  my $results_file = $self->create_filename("Infernal.$domain.out");
  my $options = " --ga " . $cm->{$domain}->{-options};
  print "OPTIONS $options\n";
  $options =~ s/\-\-toponly//;

  $self->files_to_delete($results_file);
  $self->resultsfile($results_file);
  $command .= "   $options $cmfile $filename > $results_file ";
  print STDOUT "Running infernal ".$command."\n" if $verbose;
  open(my $fh, "$command |") || 
    $self->throw("Error opening Infernal cmd <$command>." .
		 " Returned error $? Infernal EXIT: '" .
		 ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
		 "', There was " . ($? & 128 ? 'a' : 'no') .
		 " core dump");
  # this loop reads the STDERR from the blast command
  while (<$fh>) {
    if (/FATAL:(.+)/) {
      my $match = $1;
      $self->throw("Infernal failed to run - $match \n");
    }
  }
}

=head2 parse_results

  Title      : parse_results
  Usage      : $runnable->parse_results(\@dna_align_features)
  Function   : Parses all cmsearch output: rejects alignments where score is below 
             : the threshold defined for the RFAM family in question.
             : The thresholds are defined in the file /data/blastdb/Rfam/Rfam.thr.
             : A dna_align_feature representing the initial RFAM blast hit is used to determine 
             : the family and position on the genome.
             : Parses all the results (there can be more than one) and stores them in
             : an array of hashes which is then sorted by the cmsearch score
  Returns    : None
  Exceptions : Throws if it cannot find, open or close the cmsearch output
  Args       : ref to array of hashes

=cut

sub parse_results{
  my ($self, $daf,$cm) = @_;
  my @dafs;
  my $align;
  my $line =0;
  my $domain = substr($daf->hseqname,0,7);
  my $threshold = $cm->{$domain}->{-threshold};
  print STDERR "Domain $domain has threshold $threshold\n" if $verbose;
  my $results = $self->resultsfile;
  my $ff = $self->feature_factory;
  if(!-e $results){
    $self->throw("Can't parse a  results file that dosen't exist ".$results.
		 " Infernal:parse_results");
  }
  my @results;
  my ($hit,$start,$end,$hstart,$hend,$score,$evalue,$strand,$str) =0;

  open(INFERNAL, $results) or $self->throw("FAILED to open ".$results.
					   " INFERNAL:parse_results");
 LINE: while(<INFERNAL>){
   chomp;     
    if ($_ =~ /results:/){
      if ($score &&  $threshold && $score > $threshold && $align->{'name'}){
	$str = $self->parse_structure($align);
	print STDERR "positve result at ".$daf->seq_region_name.":".$daf->seq_region_start."-".$daf->seq_region_end." strand ".
           $daf->strand."\n" if $verbose;
	push @results, {str => $str,
			score => $score,
			start => $start,
			end   => $end,
			hstart => $hstart,
			hend => $hend,
			evalue => $evalue,
			strand=> $strand};		
      }	
      $align = {};
      $align->{'name'} = $domain ;
      $line = -1;
    }
    if ($_ =~ /^#/){
      $line = -1;
    }
    print STDERR "$_\t$line\n"if $verbose;
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
   if ($_ =~ />(RF\d+)/ ) {
     $hit = $1;
   }
   if ($_ =~ / Query =\s+(\d+)\s+-\s+(\d+),\s+Target =\s+(\d+)\s+-\s+(\d+)/) {
     $start = $3;
     $end = $4;
     $hstart = $1;
     $hend = $2;
   }
   if ($start && $end && $hstart && $hend && $_ =~ / Score =\s+(\S+),\s+E =\s+(\S+),/) {
     $score = $1;
     $evalue = $2;
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
    print STDERR "positve result at ".$daf->seq_region_name.":".$daf->seq_region_start."-".$daf->seq_region_end." strand ".
      $daf->strand."\n" if $verbose;
    $str = $self->parse_structure($align);
    push @results, {str => $str,
		    score => $score,
		    start => $start,
		    end   => $end,
		    strand=> $strand};
  }
  close(INFERNAL) or $self->throw("FAILED to close ".$results.
			   " INFERNAL:parse_results");
  # sort the results to get the highest scoring infernal alignment
  @results = sort {$b->{'score'} <=> $a->{'score'}} @results;
  if (@results){
    return \@results;
  }
  else {
    return;
  }
}

=head2  parse_structure

  Title      : parse_structure
  Usage      : my $structure = $runnable->parse_structure(\%hash_ref_containing_parsed_cmsearch_results)
  Function   : Parses cmsearch alignment to create a structure line which represents the
             : predicted secondary structure of the ungaped sequence
  Returns    : String representing a run length encoded form of the structure line
  Exceptions : Throws if it cannot find, open or close  the cmsearch output
  Args       : ref to an array of hashes wrapping the cmsearch output alignment

=cut

sub parse_structure{
  my ($self,$align)=@_;
  my @all_matches;
  my @stack;
  my @big_gaps;
  my $matchstring;
  my @matches;
  my $big_gap=0;
  my @attributes;

  # Brace matching
  # push open braces on to the stack
  for (my $i=0 ; $i< scalar(@{$align->{'str'}}); $i++){
    if ($align->{'str'}[$i] eq '(' or
	$align->{'str'}[$i] eq '<' or
	$align->{'str'}[$i] eq '[' or
	$align->{'str'}[$i] eq '{'){
      push @stack,$i;
    }
    # pop the positions of the open brace off the stack as you find close braces
    if ($align->{'str'}[$i] eq ')' or
	$align->{'str'}[$i] eq '}' or
	$align->{'str'}[$i] eq ']' or
	$align->{'str'}[$i] eq '>'){
      $all_matches[$i] = pop @stack;
    }
  }
  @stack = [];
# Need to do the reverse proces to get all matches
  for (my $i = scalar(@{$align->{'str'}}-1); $i >=0 ; $i--){
    if ($align->{'str'}[$i] eq ')' or
	$align->{'str'}[$i] eq '}' or
	$align->{'str'}[$i] eq ']' or
	$align->{'str'}[$i] eq '>'){
      push @stack,$i;
    }
    # pop the positions of the close brace off the stack as you find open braces
    if ($align->{'str'}[$i] eq '(' or
	$align->{'str'}[$i] eq '<' or
	$align->{'str'}[$i] eq '[' or
	$align->{'str'}[$i] eq '{'){
      $all_matches[$i] = pop @stack;
    }
  }
 for (my $i=0 ; $i< scalar(@{$align->{'str'}}); $i++){
    # Parse out large gaps by looking for ~ on the str line;
    if ($align->{'query'}[$i] eq '*'){
      my $string;
      for (my $j=$i+1 ; $j < scalar(@{$align->{'str'}}) ; $j++){	
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
    if ($align->{'match'}[$i] eq ' ') {
      $matches[$i] = '.';
      next;
    }
    # Found a match
    if (defined $all_matches[$i]){
      # check there isnt a gap at the other end
      if ($align->{'query'}[$all_matches[$i]] eq '-'){
	$matches[$i] = '.';
	next;
      }
      $matches[$i] = $align->{'str'}[$i];
      $matches[$all_matches[$i]] = $align->{'str'}[$all_matches[$i]];
    }
    else {
      $matches[$i] = $align->{'str'}[$i];
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
    if ($matches[$i]){
      $matchstring.= $matches[$i];
    }
  }
  # make all characters into either (,),or .
  $matchstring =~  s/[\<\[\{]/(/g;
  $matchstring =~  s/[\>\]\}]/)/g;
  $matchstring =~  s/[,:_-]/./g;
  return $matchstring;
}


=head2 make_gene

  Title      : make_gene
  Usage      : my $gene = $runnable->make_gene($result_hash,$dna_align_feature)
  Function   : Creates the non coding gene object from the parsed result file.
             : Takes all the results from the result file sorted by score and 
             : selects the highest scoring result that also overlaps the original BLAST hit
             : Uses Bio::EnsEMBL::Analysis::Runnable::RNAfold to create a structure
             : line for the gene. The structure obtained from the cmsearch alignment
             : is used to constrain the folding prediction made by RNAfold.
             : Adds the encoded structure as a transript attribute and xrefs are 
             : made based on the original BLAST hit
  Returns    : Bio::EnsEMBL::Gene
  Exceptions : None
  Args       : ref to array of hashes containing parsed cmsearch results
             : Bio::EnsEMBL::DnaDnaAlignFeature ($dna_align_feature)

=cut

sub make_gene{
  my ($self,$results,$daf,$descriptions,$cm) = @_;
  my $domain = substr($daf->hseqname,0,7);
  my $description = $descriptions->{$domain}->{'description'};
  my $name =     $descriptions->{$domain}->{'name'};
  my $padding =  $cm->{$domain}->{-length};
  my $type = $descriptions->{$domain}->{'type'};
  my %gene_hash;
  my $slice = $daf->slice;
  my @attributes;
  my ($start,$end,$str,$score,$exon);
  # step through all the results in order of score, highest first until
  # you find a model that overlapping the blast hit
  foreach my $result (@$results){
    if ($daf->strand == 1){
      $start = $daf->start - $padding + $result->{'start'} - 1 ;
      $end   = $daf->start - $padding + $result->{'end'}   - 1;
    } else {
      $start = $daf->end + $padding - $result->{'end'}   + 1 ;
      $end   = $daf->end + $padding - $result->{'start'} + 1;
    }
 
    $str = $result->{'str'};
    $score = $result->{'score'};
    # exons
    # need to remove padding from exon coordinates to put them in the correct place
    $exon = Bio::EnsEMBL::Exon->new
      (
       -start => $start,
       -end   => $end,
       -strand => $daf->strand,
       -slice => $slice,
       -phase => -1,
       -end_phase => -1
      );
    # reject if it falls of the start of the slice
    next if ($exon->start < 1);
    # reject if it falls of the end of the slice
    next if ($exon->end > $slice->length);
    # Only allow exons that overlap the origional dna_align_feature and 
    # have a secondary structure that is possible to parse
    last if ($exon->overlaps($daf));
  }
  # return undef if no suitable candidates are found
  return unless ($exon->start >= 1);
  return unless ($exon->overlaps($daf));
#  $daf->score($score);
  $exon->add_supporting_features($daf);
  # transcripts
  my $transcript = Bio::EnsEMBL::Transcript->new;
  $transcript->add_Exon($exon);
  $transcript->start_Exon($exon);
  $transcript->end_Exon($exon);
  $transcript->source("ensembl");
  my $gene = Bio::EnsEMBL::Gene->new;
  #Biotypes
  $gene->biotype("misc_RNA");
  $gene->biotype("snRNA")  if($type =~ /^snRNA;/ );
  $gene->biotype("snoRNA") if($type =~ /^snRNA; snoRNA;/);
  $gene->biotype("scaRNA") if($type =~ /^snRNA; snoRNA; scaRNA;/);
  $gene->biotype("rRNA")   if($type =~ /rRNA;/);
  $gene->biotype("tRNA")   if($type =~ /tRNA;/);
  $gene->biotype("sRNA")   if($type =~ /sRNA;/);
  $gene->biotype("miRNA")  if($type =~ /miRNA;/);
  $gene->biotype("CRISPR")      if($type =~ /CRISPR;/);
  #$gene->biotype("lncRNA")      if($type =~ /lncRNA;/);
  $gene->biotype("antisense")   if($type =~ /antisense;/);
  $gene->biotype("antitoxin")   if($type =~ /antitoxin;/);
  $gene->biotype("ribozyme")    if($type =~ /ribozyme;/);

  $gene->description($description." [Source: RFAM;Acc:$domain]");
  print STDERR "Rfam_id $domain ".$description."\n"if $verbose;;
  $gene->analysis($self->analysis);
  $gene->add_Transcript($transcript);
  $transcript->biotype($gene->biotype);
  $gene->source("ensembl");
  # XREFS
  my $xref = Bio::EnsEMBL::DBEntry->new
    (
     -primary_id => $domain,
     -display_id => $name,
     -dbname => 'RFAM',
     -version => 0,
     -description => $description." [Source: RFAM;Acc:$domain]",
    );
  # Use RNA fold to tidy up the structure parsed from cmsearch results
  my $seq = Bio::PrimarySeq->new(
				 -display_id => $domain,
				 -seq => $gene->seq,
				);
  my $RNAfold = Bio::EnsEMBL::Analysis::Runnable::RNAFold->new
    (
     -analysis  => $self->analysis,
     -sequence  => $seq,
     -structure => $str,
    );
  $RNAfold->run;
  # return if there is no structure to display
  return unless $RNAfold->structure;
  # get the final structure encoded by run length
  my @final_str = @{$RNAfold->encoded_str};
  foreach my $str (@final_str){
    # add the transcript attribute to the gene hash
    my $attribute = Bio::EnsEMBL::Attribute->new
      (-CODE => 'ncRNA',
       -NAME => 'Structure',
       -DESCRIPTION => 'RNA secondary structure line',
       -VALUE => $str
      );  
    push @attributes,$attribute;
  }    
  # add the final structure to the gene as a transcript attribute
  $gene_hash{'attrib'} = \@attributes;
  $gene_hash{'gene'} = $gene;
  $gene_hash{'xref'} = $xref;
  print "Chosen hit and structure constraint : $start $end " . $gene->strand ." $description $name\n$str\n";
  return \%gene_hash;
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
  my ($self,$daf,$cm)=@_;
  my $domain = substr($daf->hseqname,0,7);
  my $slice = $daf->slice;
  my $padding = $cm->{$domain}->{-length};
  my $old_start = $daf->start;
  my $old_end = $daf->end;
  print STDERR "Using padding of $padding\t"if $verbose;
  # add padding
  $daf->start($daf->start-$padding);
  $daf->end($daf->end+$padding);
  my $seq = Bio::Seq->new(
			    -seq => $daf->seq,
			    -display_id => $domain,
			   );
  # remove padding
  $daf->start($old_start);
  $daf->end($old_end);
  print STDERR " total seq length = "if $verbose;;
  print $seq->length."\n"if $verbose;;
  return $seq;
}


=head2 get_descriptions

  Title      : get_descriptions
  Usage      : $runnable->get_descriptions
  Function   : Fetches and stores predefined descriptions for each RFAM familly
             : in the Rfam.descriptions file
  Returns    : none
  Exceptions : Throws if it cannot open the descriptions file
  Args       : Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub get_descriptions{
  my($self)= @_;
  my %descriptions;
  my $domain;
  my $name;
  # read descriptions file
  return undef unless -e $self->analysis->db_file."/Rfam.seed";
  open( T,$self->analysis->db_file."/Rfam.seed") or
    $self->throw("can't file the ".$self->analysis->db_file."/Rfam.seed file");
  while(<T>) {
    chomp;
    if ($_ =~ /^\#=GF AC   (.+)/){
      $domain = $1;
     }
    if ($_ =~ /^\#=GF DE   (.+)/){
      $descriptions{$domain}{'description'} = $1;
    }
    if ($_ =~ /^\#=GF ID   (.+)/){
      $descriptions{$domain}{'name'} = $1;      
    }
    if ($_ =~ /^\#=GF TP   Gene; (.+)/){
      $descriptions{$domain}{'type'} = $1;      
    }
  }
  close T;
  return \%descriptions if scalar(keys %descriptions > 0);
  $self->throw("Unable to find descriptions");
  return undef;
}

sub read_cm_library {
  # rfam cm file is now one big file with lots of entries this code is from rfam scan.pl and is
  # used to parse the large file
  my ( $self,$fh) = @_;
 
  my $off = 0;
  my %cm;
  my $name;
  while(<$fh>) {
    if( /^ACC\s+(\S+)/ ) {
      $name = $1;
      $cm{$name}->{-accession} = $1; #RF number
      $cm{$name}->{-offset} = $off;
    }
    if( /^CLEN\s+(\d+)/ ) {
      $cm{$name}->{-length} = $1;
    }
    if( /^ACC\s+(\S+)/ ) {
      $cm{$name}->{-accession} = $1; #RF number
    }
    if( /^SCOM\s+cmsearch\s+\-Z.*\s+(\-\-toponly.*)\s+CM/ ) {
      $cm{$name}->{-options} = $1;
    }
    if( /^GA\s+(\d+)/ ) {
      $cm{$name}->{-threshold} = $1;
    }  
    if( /^\/\// ) {
      $off = tell($fh);
    }
  }
  return \%cm;
}

sub write_cm_file {
  my ($self,$cm,$cmfh) =@_;
  my $id = $self->query->display_id;
  my $outfile = $self->create_filename($id,'cm');
  seek($cmfh,$cm->{$id}->{-offset},0);
  open( F, ">$outfile" ) or die "FATAL: Failed to retrieve CM [$id] from handle [$cm]\n";
  while(<$cmfh>) {
    print F $_;
    last if( /^\/\// );
  }
  return $outfile;
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
