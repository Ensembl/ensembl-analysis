=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


package Bio::EnsEMBL::Analysis::Runnable::CMSearch;

use strict;
use warnings;
use Getopt::Long;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis::Runnable::RNAFold;
use IO::File;
use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name write_seqfile);
use Path::Tiny qw(path);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

my $verbose = 1;


=head2 new

  Title      : new
  Usage      : my $runnable = Bio::EnsEMBL::Analysis::Runnable::CMSearch->new
  Function   : Instantiates new CMSearch runnable
  Returns    : Bio::EnsEMBL::Analysis::Runnable::CMSearch object
  Exceptions : none
  Args       : Array ref of Bio::EnsEMBL::DnaDnaAlignFeature objects

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($queries,$thresholds) = rearrange(['QUERIES'], @args);
  $self->queries($queries);
 return $self;
}

sub run {
  my ($self) = @_;

  # get_input files - CM, results file
  my $cm_path = path( $self->analysis->db_file );
  my $seeds = path($self->analysis->gff_feature);
  my $output = path($self->datadir);
  my $program = $self->program;

  my $slice = $self->queries; #$slice_adaptor->fetch_by_seq_region_id($slice_id);
  my $slice_id = $slice->name;
  my $results_path = path($output . "/" . $slice_id . ".results");
  my $tblout_path = path($output . "/" . $slice_id . ".out");
  
  # run cmsearch
  my $filename = $output . "/" . $slice_id . ".fa";
  my $outseqfile = path($filename);
  $outseqfile->spew(">".$slice_id."\n".$slice->seq());

  my $options = " --rfam --cpu 4 --nohmmonly --cut_ga --tblout $tblout_path ";

  $options =~ s/\-\-toponly//;

  my $command = "$program $options $cm_path $filename > $results_path ";

  system($command);

  print STDOUT "(I) Finished running infernal " . $command."\n" if $verbose;

  # extract and filter rfam model mappings
  my $cm = $cm_path->slurp;
  my $results = $tblout_path->slurp; 
  my %descriptions = $self->get_descriptions($seeds);

  my %cm_models = $self->extract_metrics($cm); print STDOUT "(I) Finished extracting CM models and metrics \n" if $verbose;
  my @temp_map = $self->extract_rfam_mappings_tblout($results); print STDOUT "(I) Finished extracting Rfam results \n" if $verbose;
  
  my @res_map = $self->remove_overlapped_structures(\@temp_map, \@temp_map); print STDOUT "(I) Finished removing overlapped Rfams  \n" if $verbose;
  
  my @res_map_filtered = $self->filter_results(\%cm_models, \@res_map); print STDOUT "(I) Finished filtering by score threshold \n" if $verbose;

  # build and save sncRNAs to DB
  foreach my $result (@res_map_filtered){
    my $rfam = $result->{'query'};
    my $accession = $result->{'accession'};

    if (exists($cm_models{$rfam})){
      my $cm_model = $cm_models{$rfam};
      my $description = exists($descriptions{$accession}) ? $descriptions{$accession} : undef;
      my $gene = $self->make_gene($cm_model, $result, $description);
      $self->output($gene) if $gene;
    }

  }

  # clean up
  $command = "rm $tblout_path $results_path $filename ";
  system($command);

}


sub filter_results{
	my ($self, $cm_models, $results) = @_;

	my @filtered;
	foreach my $result (@$results){
		my $rfam = $result->{'query'};
		if (exists($cm_models->{$rfam})){
			my $threshold = $cm_models->{$rfam}->{-threshold};
      my $min_length = $cm_models->{$rfam}->{-length} - 5;
      my $max_length = $cm_models->{$rfam}->{-maxlength};
      my $mapping_length = abs($result->{'end'} - $result->{'start'});

      # although not included in RefSeq filters, additional filters that consider sizes and score_to_size ratios can be applied
      # in future work to further exclude FPs
      #
      # my $is_valid_size = $mapping_length > $min_length && $mapping_length < $max_length ? 1 : 0;
      # my $score_size_ratio = $result->{'score'} / $mapping_length;

      if ($rfam eq 'LSU_rRNA_eukarya'){
        $threshold = 1700;
      }
      if ($rfam eq 'LSU_rRNA_archaea'){ # remove from CM ?
        $threshold = 2000;
      }
      if ($rfam eq 'LSU_rRNA_bacteria'){ # remove from CM ?
        $threshold = 2000;
      }
      if ($rfam eq 'SSU_rRNA_eukarya'){
        $threshold = 1600;
      }
      if ($rfam eq '5_8S_rRNA'){
        $threshold = 85;
      }
      if ($rfam eq '5S_rRNA'){
        $threshold = 75;
      }
      
			if ($result->{'score'} > $threshold){
				push @filtered, $result;
			}
		}
	}
	return @filtered;
}

# parses input CM file to extract required metrics for each rfam model
sub extract_metrics{
	my ($self, $cm) = @_;
	my @cm_models = split(/\/\/\n/, $cm);

	my %cm_models = {};

	foreach my $cm_model (@cm_models){
		my @temp = split(/\n/, $cm_model);
		if( $cm_model =~ m/NAME\s+(\S+)/ ) {
			my $rfam = $1;
			foreach my $kv (@temp){
				if( $kv =~ m/^NAME\s+(\S+)/ ) {
				  $cm_models{$rfam}->{-name} = $1;
				}
				if( $kv =~ m/^DESC\s+(\S+)/ ) {
				   $cm_models{$rfam}->{-description} = $1;
				}
				if( $kv =~ m/^CLEN\s+(\d+)/ ) {
	     		  $cm_models{$rfam}->{-length} = $1;
	    	}
        if($kv =~ m/^W\s+(\d+)/ ){
          $cm_models{$rfam}->{-maxlength} = $1;
        }
			  if( $kv =~ m/^ACC\s+(\S+)/ ) {
			      $cm_models{$rfam}->{-accession} = $1; #RF number
			  }
			  if( $kv =~ m/^GA\s+(\d+)/ ) {
			      $cm_models{$rfam}->{-threshold} = $1;
			  }
			}
		}
	}
	return %cm_models;
}

# extracts results from the output tables post-CMSearch run
sub extract_rfam_mappings_tblout{
  my ($self, $results) = @_;
  my @res = split(/\n/, $results);
  my %res_map = {};

  my @processed_results;

  foreach my $result (@res){
    if ($result =~ m/^#.+/){
      next;
    }
    my @hit = split(/\s+/, $result);
    my $accession = $hit[3];
    my $target_name = $hit[0];
    my $query_name = $hit[2];
    
    my ($start,$end,$hstart,$hend,$score,$evalue,$strand) =0;
    $hstart = $hit[5]; 
    $hend = $hit[6];

    $start = $hit[7];
    $end = $hit[8];
    $strand = $hit[9] eq "+" ? 1 : -1;
    $evalue = $hit[15];
    $score = $hit[14];
    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
	    -slice          => $self->queries,
	    -start          => $strand == 1 ? $start : $end,
	    -end            => $strand == 1 ? $end : $start,
	    -strand         => $strand,
	    -hstart         => $hstart,
	    -hend           => $hend,
	    -hstrand        => $strand,
	    -score          => $score,
	    -hseqname       => length($target_name) > 39 ? substr($target_name, 0, 39) : $target_name,,
	    -p_value	=> $evalue,
      -align_type => 'ensembl',
      -cigar_string  => abs($hend - $hstart) . "M",
	 );
	
   push @processed_results, {
      score => $score,
      start => $start,
      end   => $end,
      strand=> $strand,
    	query => $query_name,
    	accession => $accession,
    	daf => $daf };

  }

  return @processed_results;
}

# remove low scoring duplicate structures; some snoRNAs are labelled differently but completely overlap after mapping
sub remove_overlapped_structures{
  my ($self, $rfams_x, $rfams_y) = @_;

  my @filtered;
  my $counter = 0;
  my %seen;
  foreach my $rfam_x (@$rfams_x){
    my $structure = $rfam_x;
    foreach my $rfam_y (@$rfams_y){
      $counter++;
      my $is_overlapped = ($rfam_x->{'start'} <= $rfam_y->{'end'} ) && ($rfam_x->{'end'} >= $rfam_y->{'start'}) ? 1 : 0;
      if ($is_overlapped){
        $structure = $structure->{'daf'}->{'score'} >= $rfam_y->{'daf'}->{'score'} ? $structure : $rfam_y;
      }
    }
    splice(@$rfams_y, $counter, 1);
    $counter = 0;
    $seen{$structure->{'accession'} ."_".$structure->{'start'}} = $structure;
  }
  @filtered = values %seen;
  return @filtered;
}


# method not used anymore; parsing results from output tables instead
sub extract_rfam_mappings{ 
	my ($self, $results) = @_;
	my @res = split(/\/\/\n/, $results);

	my %res_map = {};

	my @processed_results;

	foreach my $result (@res){
		if ($result =~ m/^Query:\s+(\S+)\s+/ ){
			my $rfam = $1;

			$result =~ m/(RF\d+)/ ;
			my $accession = $1;
			my @hits = split(/rank/, $result);
      my $counter = 0;
			foreach my $hit (@hits){
        $counter++;
        next if $counter lt 3;
				my $align = {};
				my $line = -1;
				my @lines = split(/\n/, $hit);
        my $daf;
				my ($start,$end,$hstart,$hend,$score,$evalue,$strand,$str) =0;
        my $index = 0;
				foreach my $kv (@lines){
          $index++;
					if ($kv =~ /NC$/){
						$line = -1;
					}
					$line++;
          my $hit_name;
          if ($kv =~ /\s+\S+\s+\d+\s+[AGCUagcu.]+\s+\d+\s+$/){
            $kv =~ m/\s+(\S+)\s+\d+\s+([AGCUagcu.]+)\s+\d+\s*$/;
            $align->{'name'} = $1;
            push @{$align->{'target'}},  split//, $2 ;
            $line = 10;
          }
          if ($kv =~ /\d+$/){
            $kv =~ m/([AGCUagcu-]+)\s+\d+$/;
            push @{$align->{'query'}} ,  split//, $1 ;
          }
          if ($kv =~ /CS$/){
            $kv =~ m/([:()~{}<_>,-.]+)\s+CS$/ ;
            push @{$align->{'str'}},  split//, $1 ;
          }
          if ($line ==  11){
            push @{$align->{'match'}} , split//,substr($kv,45,65) ;
          }

						my @temp = split(/ +/, $lines[2]);
						$evalue = $temp[3];
						$score = $temp[4];
						$hstart = $temp[7];
						$hend = $temp[8];
						$start = $temp[10];
						$end = $temp[11];
						$strand = $end < $start ? -1 : 1;
					$daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
				        -slice          => $self->queries,
				        -start          => $strand == 1 ? $start : $end,
				        -end            => $strand == 1 ? $end : $start,
				        -strand         => $strand,
				        -hstart         => $hstart,# < $hend ? $hstart : $hend,
				        -hend           => $hend,# > $hstart ? $hend : $hstart,
				        -hstrand        => $strand,
				        -score          => $score,
				        -hseqname       => $hit_name,
				        -p_value			=> $evalue,
                -align_type => 'ensembl',
                -cigar_string  => abs($end - $start) . "M",
				      );
				}
				$str = $self->parse_structure($align);
				push @processed_results, {str => $str,
			        score => $score,
			        start => $start,
			        end   => $end,
			        strand=> $strand,
			    	query => $rfam,
			    	accession => $accession,
			    	daf => $daf };
			}

		}

	}

	return @processed_results;
}

# make_gene
sub make_gene{
  my ($self, $cm_model, $result, $description) = @_;

  my $slice = $self->queries;
  my $domain = $result->{'query'};
  my $padding =  $cm_model->{-length};
  my $type = $description->{'type'};
  my %gene_hash;
  my @attributes;
  my ($start,$end,$str,$score,$exon, $daf, $accession);

  if ($result->{'strand'} == 1){
  	$start = $result->{'start'}; # - 1;
  	$end = $result->{'end'}; # - 1;
  } else {
  	$start = $result->{'end'}; # + 1;
  	$end = $result->{'start'}; # + 1;
  }

  $str = "";#$result->{'str'};
  $score = $result->{'score'};
  $daf = $result->{'daf'};
  $accession = $result->{'accession'};

  # exons
  $exon = Bio::EnsEMBL::Exon->new
      (
       -start => $start,
       -end   => $end,
       -strand => $result->{'strand'},
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
  #last if ($exon->overlaps($daf));
  # return undef if no suitable candidates are found
  return unless ($exon->start >= 1);
  return unless ($exon->overlaps($daf));
  
  #Biotypes
  my $biotype = "misc_RNA";
  $biotype = "snRNA"  if($type =~ /^snRNA;/ );
  $biotype = "snoRNA"  if($type =~ /^snRNA; snoRNA;/);
  $biotype = "scaRNA" if($type =~ /^snRNA; snoRNA; scaRNA;/);
  $biotype = "rRNA"   if($type =~ /rRNA;/);
  $biotype = $domain   if($domain =~ /RNaseP/);
  $biotype = "Vault_RNA"   if($domain =~ /Vault/);
  $biotype = "Y_RNA"   if($domain =~ /Y_RNA/);
  $biotype = "antisense"   if($type =~ /antisense;/);
  $biotype = "antitoxin"   if($type =~ /antitoxin;/);
  $biotype = "ribozyme"    if($type =~ /ribozyme;/);
  
  my @temp_id = split(":", $slice->name);
  my $rel_start = $temp_id[3] + $start;
  my $rel_end = $rel_start + abs($end - $start);
  my $chrom = "chr$temp_id[2]"; 
  my $seq = Bio::PrimarySeq->new(
    #-display_id => $biotype . "_" . $temp_id[2] ."_". $rel_start ."_" . $rel_end . "_" . $start,
    -display_id => $domain . "_" . $rel_end,
    -seq => $daf->seq,
   );

  my $RNAfold = Bio::EnsEMBL::Analysis::Runnable::RNAFold->new
  (
    -analysis  => $self->analysis,
    -sequence  => $seq,
  );
  $RNAfold->run;
  return unless $RNAfold->structure;

  my $new_daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
    -slice          => $self->queries,
    -start          => $start,
    -end            => $end,
    -strand         => $daf->strand,
    -hstart         => $daf->hstart,# < $hend ? $hstart : $hend,
    -hend           => $daf->hend,# > $hstart ? $hend : $hstart,
    -hstrand        => $daf->strand,
    -score          => $daf->score,
    -hseqname       => length($daf->{'hit_name'}) > 39 ? substr($daf->{'hit_name'}, 0, 39) : $daf->{'hit_name'},
    -p_value                 => $daf->p_value,
    -align_type => 'ensembl',
    -cigar_string  => abs($end - $start) . "M",
    -hcoverage    => $RNAfold->score,
    );
  
  my $hit_name_id = $accession;

  $daf = $new_daf;
  $daf->analysis($self->analysis);
  $daf->hseqname($hit_name_id);
  $exon->add_supporting_features($daf);
  
  # transcripts
  my $transcript = Bio::EnsEMBL::Transcript->new;
  $transcript->add_Exon($exon);
  $transcript->start_Exon($exon);
  $transcript->end_Exon($exon);
  $transcript->source("ensembl");
  my $gene = Bio::EnsEMBL::Gene->new;
  $gene->biotype($biotype);

  $gene->description($description->{'description'} ." [Source: RFAM;Acc:$accession]");
  print STDERR "Rfam_id $accession ".$description."\n"if $verbose;;
  $gene->analysis($self->analysis);
  $gene->add_Transcript($transcript);
  $transcript->biotype($gene->biotype);
  $gene->source("ensembl");
  
  # XREFS
  my $xref = Bio::EnsEMBL::DBEntry->new
    (
     -primary_id => $accession,
     -display_id => $domain,
     -dbname => 'RFAM',
     -version => 1,
     -description => $description." [Source: RFAM;Acc:$accession]",
    );
  
  my @final_str = @{$RNAfold->encoded_str};
  foreach my $str (@final_str){
    # add the transcript attribute to the gene hash
    my $attribute = Bio::EnsEMBL::Attribute->new
      (-CODE => 'ncRNA',
       -NAME => 'Structure',
       -DESCRIPTION => 'RNA secondary structure line',
       -VALUE => $RNAfold->structure, #$str
      );
    push @attributes,$attribute;
  }

  my $daf_path = $self->datadir . "/daf_metrics.tsv";
  open(FH, '>>', $daf_path) or die "Could not write to $daf_path";
  my $strand = $gene->strand > 0 ? "+" : "-";

  print FH  $chrom . "\t" . 
            $rel_start . "\t" . 
            $rel_end . "\t" . 
            $accession . "\t" . 
            abs($start - $end) . "\t" .
            $strand . "\t" . 
            $daf->{'score'} . "\t" . 
            $daf->{'p_value'} . "\t" . 
            $RNAfold->score . "\t" . 
            $daf->{'score'} / abs($start - $end) . "\t" .
            $gene->biotype . "\n";
  
  close(FH);


  # add the final structure to the gene as a transcript attribute
  $gene_hash{'attrib'} = \@attributes;
  $gene_hash{'gene'} = $gene;
  $gene_hash{'xref'} = $xref;
  print "Chosen hit and structure constraint : $start $end " . $gene->strand ." $description $domain\n$str\n";
  return \%gene_hash;
}

# methods not used in new approach
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

sub get_descriptions{
  my($self, $seed_path)= @_;
  my %descriptions;
  my $domain;
  my $name;

  # read descriptions file
  return undef unless -e $seed_path;
  open( T,$seed_path) or
    $self->throw("can't file the ".$seed_path . " file");
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
  return %descriptions if scalar(keys %descriptions) > 0;
  $self->throw("Unable to find descriptions");
  return undef;
}

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

