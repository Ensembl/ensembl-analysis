#!/usr/bin/env perl

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
# Usage examples:

# This script finds Ensembl readthrough and misaligned transcripts as follows:

# if an Ensembl transcript has exons overlapping more than
# READTHROUGH_MAX_GENE_OVERLAP transcripts belonging to different genes
# (including the gene the Ensembl transcript belongs to) on the same strand,
# the Ensembl transcript is reported as "readthrough"

# if an Ensembl transcript has exons overlapping (on any strand) more than
# MISALIGNED_MAX_GENE_OVERLAP transcripts belonging to different genes
# (including the gene the Ensembl transcript belongs to)
# and the percentage of exonic overlap of the other transcript is greater or equal than
# MISALIGNMENT_EXONIC_PERCENTAGE for at least one of the other transcripts, the Ensembl transcript is reported as "misaligned"

# if an Ensembl transcript overlaps any transcript on the opposite strand OR
# it overlaps a transcript on the same strand but there is not any exon overlap,
# and there is not another havana or ensembl_havana transcript within
# the gene the Ensembl transcript belongs to that has a similar* start and end,
# the Ensembl transcript is reported as "other overlap" if there are more than OTHER_MAX_GENE_OVERLAP overlap genes involved.
# *Error defined by the following fraction:
# the start of the transcript lies between start-length_of_overlap_gene/MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION and start+length_of_overlap_gene/MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION
# and the end of the transcript lies between end-length_of_overlap_gene/MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION and end+length_of_overlap_gene/MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION 

# if an Ensembl transcript has introns overlapping more than
# INTRON_MAX_GENE_OVERLAP transcripts
#   whose exons overlap the introns (belonging to different genes (including the gene the Ensembl transcript belongs to)) on the same strand,
#   AND this Havana transcript does not overlap any other Havana transcript that also overlaps both this overlap gene transcript (exonic overlap) and the Ensembl transcript (with MISALIGNMENT_EXONIC_PERCENTAGE),
#   (# for these exonic overlap cases, the Ensembl transcript is reported as "ambiguous" separatedly)
#   AND this Havana transcript does not overlap any other Havana transcript that also overlaps both this overlap gene transcript (coding exonic overlap) and the Ensembl transcript (with MISALIGNMENT_EXONIC_PERCENTAGE),
# the Ensembl transcript is reported as "intron"

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(exon_overlap coding_exon_overlap features_overlap overlap_length);
use List::Util qw(max min);

$| = 1;

my $READTHROUGH_MAX_GENE_OVERLAP = 1;

my $MISALIGNMENT_EXONIC_PERCENTAGE = 80;
my $MISALIGNED_MAX_GENE_OVERLAP = 2;

my $MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION = 10;
my $OTHER_MAX_GENE_OVERLAP = 2;

my $INTRON_MAX_GENE_OVERLAP = 1; # 1 gene whose exons overlap the transcript introns is allowed to cover the same gene the transcript belongs to

my $dbname;
my $dbhost;
my $dbport = '3306';
my $dbuser = 'ensro';
my $dbpass = undef;
my $chromosome = undef;
my $coordsystem;

GetOptions('dbname:s'       => \$dbname,
           'dbhost:s'       => \$dbhost,
           'dbport:s'       => \$dbport,
           'dbuser:s'       => \$dbuser,
           'dbpass:s'       => \$dbpass,
           'chr:s'          => \$chromosome,
           'coord_system:s' => \$coordsystem);

if (!$dbuser || !$dbhost || !$dbname) {
  throw("DB connection parameters missing. Can't connect.");
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbname => $dbname,
                                            -host => $dbhost,
                                            -user => $dbuser,
                                            -pass => $dbpass,
                                            -port => $dbport);

my $sa = $db->get_SliceAdaptor();

my @slices = ();
my $all_slices = $sa->fetch_all('toplevel',$coordsystem,1,undef);
if ($chromosome) { # if chr was defined as parameter, choose the corresponding chromosome slices
  foreach my $sl (@{$all_slices}) {
    if ($sl->seq_region_name() eq $chromosome) {
      push(@slices,$sl);
    }
  }
} else { # if chr was not defined, choose all toplevel slices
  @slices = @{$all_slices};
}

print("Working on ".scalar(@slices)." toplevel slices.\n");

foreach my $slice (@slices) {
  foreach my $gene (@{$slice->get_all_Genes()}) {
  	my $g_source = $gene->source();
  	
  	# find ensembl t within ensembl and ensembl_havana g
  	if (($g_source eq 'ensembl') or ($g_source eq 'ensembl_havana')) {
  	  TRANSCRIPT: foreach my $transcript (@{$gene->get_all_Transcripts()}) {
  	  	if ($transcript->source() eq 'ensembl') {
  	  	  my $transcript_strand = $transcript->seq_region_strand();
  	  	  my $transcript_slice = $transcript->feature_Slice();

          my @intron_exon_overlap_genes = ();
          my $intron_exon_overlap = 0;

          my @same_strand_exon_overlap_genes = ();
          my $same_strand_exon_overlap = 0;
          
          my @misalignment_exon_overlap_genes = ();
          my $misalignment_exon_overlap = 0;
          
          my @other_overlap_genes = ();
          my $other_overlap = 0;
          my $other_overlap_gene_ignore = 0;
          
          # find all g that overlap the current ensembl t
          SLICE_GENE: foreach my $overlap_gene (@{$transcript_slice->get_all_Genes()}) {

          	my $same_strand_exon_overlap_gene_pushed = 0;
          	my $misalignment_exon_overlap_gene_pushed = 0;
          	my $other_overlap_gene_pushed = 0;
          	my $intron_exon_overlap_gene_pushed = 0;
          	
          	my $overlap_gene_source = $overlap_gene->source();
          	
          	# look at genes on both strands
          	# (including the gene which the transcript belongs to
          	# because there might be other transcripts we also want to look at)
          	# and the target gene sources are 'havana' and 'ensembl_havana'
          	# because they contain 'havana' and 'ensembl_havana' transcripts
            if (($overlap_gene_source eq 'havana') or ($overlap_gene_source eq 'ensembl_havana')) {
          	    	
              foreach my $slice_transcript (@{$overlap_gene->get_all_Transcripts()}) {
              	my $slice_transcript_source = $slice_transcript->source();

                # only look at havana and ensembl_havana transcripts
              	if (($slice_transcript->dbID() != $transcript->dbID()) and
              	    (($slice_transcript_source eq 'havana') or ($slice_transcript_source eq 'ensembl_havana'))
              	   ) {

                  # look for any overlap between the introns and the slice transcript exons on the same strand
                  if ( (!($intron_exon_overlap_gene_pushed)) 
                   and ($overlap_gene->seq_region_strand() == $transcript_strand)
                   and ($overlap_gene->seq_region_name() eq $transcript->seq_region_name())
                   ) {
                   
                    INTRON: foreach my $intron (@{$transcript->get_all_Introns()}) {

                      foreach my $slice_exon (@{$slice_transcript->get_all_Exons()}) {

                        if (features_overlap($intron,$slice_exon)) {

                          # ignore the overlap if there is another Havana transcript that also overlaps both this overlap_gene transcript (coding exonic overlap) and the Ensembl transcript (with MISALIGNMENT_EXONIC_PERCENTAGE %)
                          my $overlap_slice = $overlap_gene->feature_Slice();
                          my $overlap_slice_fixed;                       
                          
                          # make sure that the slice start and end lie within the assembly exception start and end.
                          # Otherwise, the API returns truncated transcripts in cases where a transcript crosses an assembly exception boundary.
                          if (!($overlap_slice->is_reference())) {
                            my @overlap_slice_unique_list = @{$sa->fetch_by_region_unique($overlap_slice->coord_system_name(),$overlap_slice->seq_region_name())};
                            my $overlap_slice_unique = shift(@overlap_slice_unique_list);

                            my @overlap_slice_fixed_list = @{$sa->fetch_by_region_unique($overlap_slice_unique->coord_system_name(),
                                                                                         $overlap_slice_unique->seq_region_name(),
                                                                                         max($overlap_slice->seq_region_start(),$overlap_slice_unique->seq_region_start()),
                                                                                         min($overlap_slice->seq_region_end(),$overlap_slice_unique->seq_region_end()))};
                            $overlap_slice_fixed = shift(@overlap_slice_fixed_list);
                          } else {
                            $overlap_slice_fixed = $overlap_slice;
                          }
                          
                          my @transcripts_overlapping_the_overlap_gene = @{$overlap_slice_fixed->get_all_Transcripts_by_source('havana',0)};
                          push(@transcripts_overlapping_the_overlap_gene,@{$overlap_slice_fixed->get_all_Transcripts_by_source('ensembl_havana',0)});
                          
                          foreach my $transcript_overlapping_the_overlap_gene (@transcripts_overlapping_the_overlap_gene) {

                            if ( ($transcript_overlapping_the_overlap_gene->dbID() != $slice_transcript->dbID())
                             and ($transcript_overlapping_the_overlap_gene->seq_region_strand() == $overlap_gene->seq_region_strand())
                             and ($transcript_overlapping_the_overlap_gene->seq_region_name() eq $overlap_gene->seq_region_name())
                               ) {
                              my $overlap_length_totog = coding_exon_overlap($transcript_overlapping_the_overlap_gene,$transcript);

                              # report the cases where there is enough overlap between an Ens and an Hav exonic structure as ambiguous for manual review
                              my $length = $transcript_overlapping_the_overlap_gene->length();
                              my $overlap_length = exon_overlap($transcript_overlapping_the_overlap_gene,$transcript);
                              if (($overlap_length*100/$length >= $MISALIGNMENT_EXONIC_PERCENTAGE) and
                                    ($transcript_overlapping_the_overlap_gene->seq_region_strand() == $transcript_strand) and
                                    ($transcript_overlapping_the_overlap_gene->seq_region_name() eq $transcript->seq_region_name())) {
                                  print_overlap_genes("ambiguous",0,$transcript);
                                  last INTRON;
                              }

                              my $translateable_length = length($transcript_overlapping_the_overlap_gene->translateable_seq());
                              if ($translateable_length > 0) {

                                if (($overlap_length_totog*100/$translateable_length >= $MISALIGNMENT_EXONIC_PERCENTAGE) and
                                    ($transcript_overlapping_the_overlap_gene->seq_region_strand() == $transcript_strand) and
                                    ($transcript_overlapping_the_overlap_gene->seq_region_name() eq $transcript->seq_region_name())) {
                                  last INTRON;
                                }
                              }
                            }
                          }
 
                          # if the Ensembl transcript does not have enough overlap with any other Havana transcript that overlaps the overlap_gene, then flag this overlap_gene
                          push(@intron_exon_overlap_genes,$overlap_gene);
                          $intron_exon_overlap++;
                          $intron_exon_overlap_gene_pushed = 1;
                          last INTRON;
                        }
                      }
                    }
                  }

              	  my $exonic_overlap_length = exon_overlap($slice_transcript,$transcript);

              	  if (($exonic_overlap_length*100/$slice_transcript->length() >= $MISALIGNMENT_EXONIC_PERCENTAGE) and
              	      (!$misalignment_exon_overlap_gene_pushed)) {
              	    # take overlap gene into account for potential labelling of misalignment
              	    push(@misalignment_exon_overlap_genes,$overlap_gene);
                    $misalignment_exon_overlap++;
                    $misalignment_exon_overlap_gene_pushed = 1;
              	  }

                  if (($overlap_gene->seq_region_strand() == $transcript_strand) and ($overlap_gene->seq_region_name() eq $transcript->seq_region_name()) and 
              	  #if (($overlap_gene->strand() == $transcript_strand) and 
              	      $exonic_overlap_length) {
                    # take overlap gene into account for potential labelling of readthrough
              	      if (!$same_strand_exon_overlap_gene_pushed) {
              	        push(@same_strand_exon_overlap_genes,$overlap_gene);
              	        $same_strand_exon_overlap++;
              	        $same_strand_exon_overlap_gene_pushed = 1;
              	      } else {
              	      	;
              	      }
  
                      # if a transcript within the same gene having similar start and end is found
                      # do not consider the current transcript for the other overlap category
              	      if (!$other_overlap_gene_ignore and ($overlap_gene->dbID() == $gene->dbID())) {
              	        if (same_start_end($transcript,
                                           $overlap_gene,
                                           (abs($overlap_gene->seq_region_end()-$overlap_gene->seq_region_start)/$MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION))) {
                          @other_overlap_genes = ();
                          $other_overlap_gene_pushed = 0;
                          $other_overlap_gene_ignore = 1;
                        } else {
                          ;	
                        }
                      } else {
                      	;
                      }
              	      
              	  } elsif (!$other_overlap_gene_pushed and !$other_overlap_gene_ignore) {
              	  	
              	  	# do not report the transcript if another havana or ensembl_havana transcript within
              	  	# the same gene has a similar start and end
              	  	if (!same_start_end($transcript,
              	  	                    $overlap_gene,
              	  	                    (abs($overlap_gene->seq_region_end()-$overlap_gene->seq_region_start)/$MAX_SAME_START_END_ERROR_OTHER_OVERLAP_FRACTION))) {
              	  	# take overlap gene into account for potential labelling of other overlap
              	  	  push(@other_overlap_genes,$overlap_gene);
              	  	  $other_overlap++;
              	  	  $other_overlap_gene_pushed = 1;
              	  	} else {
              	  	  # if a transcript within the same gene having similar start and end is found
              	  	  # do not consider the current transcript again for the other overlap category
              	  	  @other_overlap_genes = ();
              	  	  $other_overlap_gene_pushed = 0;
              	  	  $other_overlap_gene_ignore = 1;
              	  	}
                  } # } elsif (!$other_overlap_gene_pushed
              	  
              	} # if ($slice_transcript->dbID()

              } # foreach my $slice_transcript
            } # if (($overlap_gene_source eq 'havana')
          } # SLICE_GENE:

          # skip if the transcript matches any CCDS
          if (scalar(@{$transcript->get_all_Attributes('ccds_transcript')})) {
           
            print($transcript->stable_id()." matches a CCDS transcript (attribute), skipped.");
          
          } else {

            # still try to get if there is any CCDS tsf as the ccds_transcript attribute was only introduced in release e87
            my @tsfs = @{$transcript->get_all_supporting_features()};
            foreach my $tsf (@tsfs) {
              if (substr($tsf->hseqname(),0,4) eq 'CCDS') {
                print($transcript->stable_id()." matches a CCDS transcript (tsf), skipped.");
                next TRANSCRIPT;
              }
            }
            print_overlap_genes("readthrough",$READTHROUGH_MAX_GENE_OVERLAP,$transcript,@same_strand_exon_overlap_genes);
            print_overlap_genes("misaligned",$MISALIGNED_MAX_GENE_OVERLAP,$transcript,@misalignment_exon_overlap_genes);
            print_overlap_genes("other",$OTHER_MAX_GENE_OVERLAP,$transcript,@other_overlap_genes);
            print_overlap_genes("intron",$INTRON_MAX_GENE_OVERLAP,$transcript,@intron_exon_overlap_genes);
          }
          
  	  	} # if $transcript->source()
  	  } # foreach my $transcript
  	} # if $g_source
  } # foreach my $gene
} # foreach my $slice

sub print_overlap_genes {
  my $overlap_type = shift;
  my $max_genes = shift;
  my $transcript = shift;
  my @overlap_genes = @_;
  
  my $transcript_slice = $transcript->feature_Slice();
  
  if (scalar(@overlap_genes) > $max_genes) {
    print("Ensembl ".$overlap_type." ".$transcript->biotype()." ".$transcript->stable_id()." ".$transcript_slice->name()." overlaps the following ".scalar(@overlap_genes)." genes containing havana or ensembl_havana transcripts:");

    foreach my $overlap_gene (@overlap_genes) {
      print(" ".$overlap_gene->biotype()." ".$overlap_gene->stable_id()." , ");
    }
    print("\n");
  }
}

sub same_start_end {
  my ($transcript,$gene,$max_error) = @_;
  
  foreach my $t (@{$gene->get_all_Transcripts()}) {
  	
  	if ((abs($transcript->seq_region_start()-$t->seq_region_start()) <= $max_error) and
  	    (abs($transcript->seq_region_end()-$t->seq_region_end()) <= $max_error) and
  	    (($t->source() eq 'havana') or ($t->source() eq 'ensembl_havana'))) {
      return 1;
  	}
  }
  return 0;

}

1;
