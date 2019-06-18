# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

use warnings ;
use strict;

use Getopt::Long;
use Bio::SeqIO;


my (@good_tax_strings,
    @bad_tax_strings,
    @good_keywords,
    @bad_keywords,
    $kill_list, %kill_list,
    $gene_name_in_id,
    $translate,
    $keep_pseudos,
    $keep_partials,
    $keep_stops,
    $min_length);

&GetOptions(
            'taxstring=s@'       => \@good_tax_strings,
            'rejecttaxstring=s@' => \@bad_tax_strings, 
            'keyword=s@'         => \@good_keywords,
            'rejectkeyword=s@'   => \@bad_keywords,
            'genenameinid'       => \$gene_name_in_id,
            'translate'          => \$translate,
            'pseudokeep'         => \$keep_pseudos,
            'partialkeep'        => \$keep_partials,
            'stopkeep'           => \$keep_stops,
            'minlength=s'        => \$min_length,
            'kill_list=s'        => \$kill_list,
            );

# the following tag set should cover mammalian structures; 
# the situation is a little more complex in fish; will need
# another specialised script for that

my @tags = qw(
              CH1 CH2 CH3 CH4 CH-S M M1 M2
              H H1 H2 H3 H4
              CL
              EX1 EX2 EX3 EX4 EX2A EX2B EX2C EX2R EX2T
              );

my %tags = map { $_ => 1 } @tags;


my (%good_keywords, $bad_keywords);

@good_tax_strings = map { s/\_/ /g; lc($_) } @good_tax_strings;
@bad_tax_strings  = map { s/\_/ /g; lc($_) } @bad_tax_strings;
@good_keywords    = map { s/\_/ /g; lc($_) } @good_keywords;
@bad_keywords     = map { s/\_/ /g; lc($_) } @bad_keywords;

if (defined $kill_list) {
  open KILL, $kill_list or die "Could not open kill list file '$kill_list'";
  while(<KILL>) {
    /^(\S+)/ and $kill_list{$1} = 1;
  }
  close(KILL);
}

my $outio = Bio::SeqIO->new(-format => 'fasta',
                            -fh     => \*STDOUT);

my $seqio = Bio::SeqIO->new(-format => 'IMGT_embl',
                            -fh     => \*ARGV);
  
SEQ: while(my $seq = $seqio->next_seq) {
  next if exists $kill_list{$seq->id};
  next if $seq->data_class =~ /keyword/;
  next if $seq->molecule !~ /DNA/i;

  if (@good_tax_strings) {
    next SEQ if not defined $seq->species;
    my @taxonomy = map { lc($_) } reverse $seq->species->classification;
    
    pop @taxonomy; push @taxonomy, lc($seq->species->species);
    my $taxstring = join(" ", @taxonomy);
    next SEQ if not grep { index($taxstring, $_) >= 0 } @good_tax_strings;
  }
  
  my @keywords = split(/\s*;\s*/, $seq->description);
  push @keywords, $seq->keywords;
  @keywords = map { $_ =~ s/\.\s*$//; lc($_) } @keywords;

  my $entry_is_pseudo = 0;
  #if (grep { $_ eq 'pseudogene' or $_ eq 'functionality pseudogene' } @keywords) {
  #  $entry_is_pseudo = 1;
  #} 
  
  if (@good_keywords) {
    my $all_found = 1;
    foreach my $ukw (@good_keywords) {
      if (not grep { $ukw eq $_ } @keywords) {
        $all_found = 0; last;
      }
    }
    next SEQ if not $all_found;
  }
  push @bad_keywords, 'functionality unproductive';
  my $one_found = 0;
  foreach my $ukw (@bad_keywords) {
    if (grep { $ukw eq $_ } @keywords) {
      $one_found = 1; last;
    }
  }
  next SEQ if $one_found;
  
  my (@segments, %seg_groups, @gene_regions, @stops, @segment_lists);
  
  foreach my $val ($seq->get_SeqFeatures) {
    my $location = $val->location;
    
    if ($val->primary_tag eq 'C-GENE') {
      my $gene_seg = {
        location => $location,
      };
      if ($val->annotation->get_Annotations("gene")) {
        my ($tag_val) = $val->annotation->get_Annotations("gene");
        $gene_seg->{gene} = $tag_val->value if $tag_val->value !~ /\s+/;
      }
      if ($val->annotation->get_Annotations("pseudo")) {
        $gene_seg->{pseudo} = 1;
      }
      if ($val->annotation->get_Annotations("IMGT_note")) {
        foreach my $tag_val ($val->annotation->get_Annotations("IMGT_note")) {
          if ($tag_val->value =~ /^Pseudo/i) {
            $gene_seg->{pseudo} = 1;
          }
        }
      }
      if ($val->annotation->get_Annotations("note")) {
        foreach my $tag_val ($val->annotation->get_Annotations("note")) {
          if ($tag_val->value =~ /^Pseudo/i) {
            $gene_seg->{pseudo} = 1;
          }
        }
      }

      push @gene_regions, $gene_seg;
      next;
    }

    if ($val->primary_tag eq 'STOP-CODON') {
      push @stops, $location;
      next;
    }
    
    if (exists $tags{$val->primary_tag}) {                
      next if $location->isa("Bio::Location::Fuzzy");
      
      my $seg =  {
        code         => $val->primary_tag,
        location     => $location,
        pseudo       => $entry_is_pseudo,
        partial      => 0,
        codon_start  => 0,
      };

      if ($val->annotation->get_Annotations("gene")) {
        my ($tag_val) = $val->annotation->get_Annotations("gene");
        $seg->{gene} = $tag_val->value if $tag_val->value !~ /\s+/;
      }
      if ($val->annotation->get_Annotations("partial")) {
        $seg->{partial} = 1;
      }
      if ($val->annotation->get_Annotations("codon_start")) {
        my ($tag_val) = $val->annotation->get_Annotations("codon_start");
        $seg->{codon_start} = $tag_val->value - 1;            
      }
      if ($val->annotation->get_Annotations("translation")) {
        my ($tag_val) = $val->annotation->get_Annotations("translation");
        $seg->{translation} = $tag_val->value;
        $seg->{translation} =~ s/\s+//g;
        $seg->{translation} =~ s/[^\w]/X/g;
      }
      if ($val->annotation->get_Annotations("pseudo")) {
        $seg->{pseudo} = 1;
      }
      if ($val->annotation->get_Annotations("IMGT_note")) {
        foreach my $tag_val ($val->annotation->get_Annotations("IMGT_note")) {
          if ($tag_val->value =~ /^Pseudo/i) {
            $seg->{pseudo} = 1;
          }
        }
      }
      if ($val->annotation->get_Annotations("note")) {
        foreach my $tag_val ($val->annotation->get_Annotations("note")) {
          if ($tag_val->value =~ /^Pseudo/i) {
            $seg->{pseudo} = 1;
          }
          if ($tag_val->value =~ /untranslated/i) {
            $seg->{untranslated} = 1;
          }
        }
      }
      
      push @segments, $seg;      
    }
  }
    
  @gene_regions = sort { 
    $a->{location}->start <=> $b->{location}->start; 
  } @gene_regions;
          
  my @geneless_segs;
  
  foreach my $seg (@segments) {
    # identify gene region that it belongs to
    my @matched_gregs;
    foreach my $greg (@gene_regions) {
      if ($seg->{location}->start <= $greg->{location}->end and
          $seg->{location}->end   >= $greg->{location}->start) {
        push @matched_gregs, $greg;
      } elsif ($greg->{location}->start > $seg->{location}->start) {
        last;
      }
    }

    map { $seg->{parent_genes}->{$_} = $_ } @matched_gregs;

    if (scalar(@matched_gregs) == 1) {
      push @{$seg_groups{$matched_gregs[0]}}, $seg;
    }
    #else {
      # this segment cannot be matched to exactly one gene
      # region, so it is impossible to group segments into
      # genes at all for this sequence; skip it
      #push @geneless_segs, $seg;
    #}
  }

  if (@geneless_segs and not keys %seg_groups) {
    # if all segs have a different code, we will assume that this
    # is a single-gene entry
    my %counts;
    map { $counts{$_->{code}}++ } @segments;
    if (not grep { $counts{$_} > 1 } keys %counts) {
      map { push @{$seg_groups{"single_gene"}}, $_ } @segments;
    }
  }

  GENE: foreach my $gid (keys %seg_groups) {
    # check completeness
    my @segs = sort { 
      $a->{location}->start <=> $b->{location}->start; } 
    @{$seg_groups{$gid}};
              
    push @segment_lists, \@segs;
  }
  
  # reject segment lists that have any 
  # - partial components
  # - pseudo components
  my @complete_seg_lists;
  SEGLIST: foreach my $seglist (@segment_lists) {
    my $is_pseudo = 0;
    my $is_partial = 0;

    foreach my $seg (@$seglist) {
      if ($seg->{pseudo}) {
        $is_pseudo = 1;
      } elsif (exists $seg->{parent_genes}) {
        foreach my $greg (values %{$seg->{parent_genes}}) {
          if ($greg->{pseudo}) {
            $is_pseudo = 1;
          }
        }
      }
      if ($seg->{partial}) {
        $is_partial = 1;
      } else {
        foreach my $loc ($seg->{location}->each_Location) {
          if ($loc->start_pos_type ne 'EXACT' or
              $loc->end_pos_type ne 'EXACT') {
            $is_partial = 1;
          }
        }
      }
    }

    next if $is_partial and not $keep_partials;
    next if $is_pseudo != $keep_pseudos;


    # definitions of completness:
    # EX1..EX3 
    # EX1..EX4 
    # CH1..CH-S 
    # CH1..M2
    # CH1..M
    # CL

    my $complete = 0;
    foreach my $end_pair ([$seglist->[0], $seglist->[-1]],
                          [$seglist->[-1], $seglist->[0]]) {
      my ($f, $l) = @$end_pair;
      if ($f->{code} eq 'CL' or
          ($f->{code} eq 'EX1' and $l->{code} eq 'EX3') or
          ($f->{code} eq 'EX1' and $l->{code} eq 'EX4') or
          ($f->{code} eq 'CH1' and $l->{code} eq 'CH-S') or
          ($f->{code} eq 'CH1' and $l->{code} eq 'M') or
          ($f->{code} eq 'CH1' and $l->{code} eq 'M2')) {
        $complete = 1;
        last;
      }
    }

    if ($complete or $keep_partials) {
      push @complete_seg_lists, $seglist;
    }
  }

  my (%used_names);

  
  SEGLIST: foreach my $orig_seglist (@complete_seg_lists) {

    my $strand = &determine_strand($seq, $orig_seglist);

    my @nov_seglists = &deal_with_overlapping_segments($orig_seglist);

    #if (not @nov_seglists) {
    #  printf(STDERR "Complex overlap for %s (%d) : %s\n", 
    #         $seq->id, 
    #         $orig_seglist->[0]->{location}->start,
    #         join(" ", map { $_->{code} } @$orig_seglist));
    #}
    
    foreach my $seglist (@nov_seglists) {
      &adjust_for_bad_codon_boundaries($strand, $seglist);
      
      my @locs = sort {
        $a->start <=> $b->start;
      } map { $_->{location}->each_Location } @$seglist;
      
      my $contains_utr = 0;
      if ($translate) {
        if (grep  { exists($_->{untranslated}) } @$seglist) {
          $contains_utr = 1;
        } else {
          foreach my $loc (@locs) {
            foreach my $stop (@stops) {
              if ($stop->start >= $loc->start and 
                  $stop->end   <= $loc->end) {
                $contains_utr = 1;
                last;
              }
            }
          }
        }
      }
    

      my $codestring = join(":", map { $_->{code} } @$seglist); 
      my $seqstring = join("", map { $seq->subseq($_->start, $_->end) } @locs);
      
      my (%gnames);
      foreach my $seg (@$seglist) {
        if (exists $seg->{gene}) {
          $gnames{$seg->{gene}} = $seg->{location}->start;
        }
        if (exists $seg->{parent_genes}) {
          foreach my $pg (values %{$seg->{parent_genes}}) {
            if (exists $pg->{gene}) {
              $gnames{$pg->{gene}} = $seg->{location}->start;
            }
          }
        }
      }
      
      my $id = $seq->primary_id 
          . "." . $locs[0]->start . "-" . $locs[-1]->end;
      if (exists $used_names{$id}) {
        $used_names{$id}++;
        $id .= "-" . $used_names{$id};
      } else {
        $used_names{$id} = 1;
      }
      
      my @gnames = sort {$gnames{$a}<=>$gnames{$b}} keys %gnames;
      if (@gnames and $gene_name_in_id) {
        $id .= "." . join(".", @gnames);
      }
      
      my @this_desc = ("tag=$codestring", map {"gene=$_"} @gnames);
      if (grep { $_->{pseudo} } @$seglist) {
        push @this_desc, "Pseudogene"; 
      }
      
      my $newseq = Bio::PrimarySeq->new(-seq => $seqstring,
                                        -id  => $id,
                                        -description => join("; ", @this_desc));
      
      my $frame = 0;
      if ($strand < 0) {
        $frame = $seglist->[-1]->{codon_start};
      } elsif ($strand > 0)  {
        $frame = $seglist->[0]->{codon_start};
      }

      if ($translate) {        
        my $pepseq = &translate($newseq, $frame, $contains_utr);
        if ($pepseq->seq !~ /\*/ or $keep_stops) {
          $outio->write_seq($pepseq) 
              if not defined $min_length or $min_length <= $pepseq->length;
        }
      } else {
        $outio->write_seq($newseq)
            if not defined $min_length or $min_length <= $newseq->length;
      }
    }
  }
}


sub determine_strand {
  my ($seq, $segs) = @_;

  my @trans_segs = grep { exists($_->{translation}) } @$segs;

  # translation, we can try to guess what the strand should
  # be based on that translation
  
  my $forward_count = 0;
  my $reverse_count = 0;
  my $forward_stop_count = 0;
  my $reverse_stop_count = 0;

  foreach my $seg (@trans_segs) {
    my $tr = $seg->{translation};
    
    my $forward = $seq->subseq($seg->{location}->start, $seg->{location}->end);
    my $reverse = Bio::PrimarySeq->new(-seq => $forward )->revcom->seq;

    if ($seg->{codon_start}) {
      $forward = substr($forward, $seg->{codon_start});
      $reverse = substr($reverse, $seg->{codon_start});
    }
    if (length($forward) % 3) {
      $forward = substr($forward, 0, length($forward) - (length($forward) % 3));
      $reverse = substr($reverse, 0, length($reverse) - (length($reverse) % 3));
    }

    my $forpep = Bio::PrimarySeq->new(-seq => $forward)->translate->seq;
    my $revpep = Bio::PrimarySeq->new(-seq => $reverse)->translate->seq;

    if (index($tr, $forpep) >= 0) {
      $forward_count++;
    } elsif (index($tr, $revpep) >= 0) {
      $reverse_count++;
    }

    $forpep =~ s/\*$//; $revpep =~ s/\*$//;
    $forward_stop_count += $forpep =~ tr/\*/\*/;
    $reverse_stop_count += $revpep =~ tr/\*/\*/;
  }

  if ($forward_count and not $reverse_count) {
    map { $_->{location}->strand(1) } @$segs;
    return 1;
  } elsif ($reverse_count and not $forward_count) {
    map { $_->{location}->strand(-1) } @$segs;
    return  -1;
  } elsif ($forward_stop_count and not $reverse_stop_count) {
    map { $_->{location}->strand(-1) } @$segs;
    return  -1;
  } elsif ($reverse_stop_count and not $forward_stop_count) {
    map { $_->{location}->strand(1) } @$segs;
    return 1;
  } else {
    map { $_->{location}->strand(0) } @$segs;
    return 0;
  }
}


sub adjust_for_bad_codon_boundaries {
  my ($strand, $segs) = @_;
  
  my @segs;
  if ($strand < 0) {
    @segs = sort { $b->{location}->start <=> $a->{location}->start } @$segs;
  } else {
    @segs = sort { $a->{location}->start <=> $b->{location}->start } @$segs;
  }

  my $prev_over;
  for(my $i=0; $i < @segs; $i++) {
    my $seg = $segs[$i];

    if (defined $prev_over) {
      my $over5 = $seg->{codon_start};
      
      if (($prev_over + $over5) % 3) { 
        # bad; just trim the excess off each segment; 
        if ($strand < 0) {
          $seg->{location}->end($seg->{location}->end - $over5);
          $segs[$i-1]->{location}->start($segs[$i-1]->{location}->start + $prev_over);
        } else {
          $seg->{location}->start($seg->{location}->start + $over5);
          $segs[$i-1]->{location}->end($segs[$i-1]->{location}->end - $prev_over);

          #$segs[$i-1]->{location}->end($segs[$i-1]->{location}->end + 1);
        }
        $seg->{codon_start} = 0;
      }
    }
    my $len = $seg->{location}->end - $seg->{location}->start + 1;
    $len -= $seg->{codon_start};
    $prev_over = $len % 3;
  }
  
  return 1;
}


sub translate {
  my ($seq, $frame, $contains_utr) = @_;
  
  my $trans_seq = $seq->translate(-frame => $frame);
  
  if ($trans_seq->seq =~ /\*$/) {
    my $pep = $trans_seq->seq;
    $pep =~ s/\*$//;
    
    $trans_seq->seq($pep);
  }
  if ($contains_utr and $trans_seq->seq =~ /\*/) {
    my ($pep) = ($trans_seq->seq =~ /^([^\*]+)\*/);
    $trans_seq->seq($pep);
  }
  
  return $trans_seq;
}


sub deal_with_overlapping_segments {
  my ($segs) = @_;

  my @segs = sort { 
    $a->{location}->start <=> $b->{location}->start 
        or $a->{location}->end <=> $b->{location}->end;
  } @$segs;

  my @seg_groups;

  my $strand = $segs[0]->{location}->strand;

  SEG: foreach my $seg (@segs) {
    foreach my $grp (@seg_groups) {
      foreach my $os (@$grp) {
        if ($os->{location}->start <= $seg->{location}->end and
            $os->{location}->end >= $seg->{location}->start) {
          push @$grp, $seg;
          next SEG;
        }
      }
    }
    push @seg_groups, [$seg];
  }

  # need to resolve overlaps
  foreach my $grp (@seg_groups) {
    if (scalar(@$grp) > 2) {
      # can't deal with this situation
      return ();
    } elsif (scalar(@$grp) == 2) {
      # check that the smaller segment is flush to the
      # boundary of the second segment

      if ($grp->[0]->{location}->start == $grp->[1]->{location}->start) {
        $grp->[1]->{location}->start($grp->[0]->{location}->end + 1);

        if ($strand > 0) {
          # only need to adjust codon_start if on forward strand
          my $bps_lost = $grp->[0]->{location}->end - $grp->[0]->{location}->start + 1;
          $bps_lost -= $grp->[1]->{codon_start};
          $grp->[1]->{codon_start} = (3 - ($bps_lost % 3)) % 3;
        }

      } elsif ($grp->[0]->{location}->end == $grp->[1]->{location}->end) {
        $grp->[0]->{location}->end($grp->[1]->{location}->start - 1);

        if ($strand < 0) {
          my $bps_lost = $grp->[1]->{location}->end - $grp->[1]->{location}->start + 1;
          $bps_lost -= $grp->[0]->{codon_start};
          $grp->[0]->{codon_start} = (3 - ($bps_lost % 3)) % 3;
        }
      } else {
        return ();
      }           
    }
  }

  @segs =  map { @$_ } @seg_groups;
  @segs = sort { $a->{location}->start <=> $b->{location}->start } @segs; 

  # finally, deal with the special case of the secreted and
  # membrane-bound alternative splice forms
  if (grep { $_->{code} eq 'CH-S' } @segs and
      grep { $_->{code} =~ /^M(1|2)?$/ } @segs) {
    my @core = grep { $_->{code} ne 'CH-S' and $_->{code} !~ /^M(1|2)?$/ } @segs;
    my @secreted = grep { $_->{code} eq 'CH-S' } @segs;
    my @membrane = grep { $_->{code} =~ /^M(1|2)?$/ } @segs;

    return ([@core, @secreted], [@core, @membrane]);

  } else {
    return (\@segs);
  }
 
}


=pod

=head1 NAME

constant_regions_from_imgt.pl

=head1 SYNOPSIS

This script extracts constant-region gene segments from genomic DNA entries
in the IMGT database flat file. The results are written in FASTA format to stdout

=head1 OPTIONS

    (Note: when giving an option value that contains maore thanone word, 
     e.g. homo sapiens, spaces should be replaced with "_", i.e. homo_sapiens)

    -taxstring       : Restrict to entries with taxonomy containing this pattern

    -rejecttaxstring : Ignore entries with taxonomy containing this pattern

    -keyword         : Restrict to entries containing these keywords
    
    -rejectkeyword   : Ignore entries including these keywords
                       
    -translate       : Translate the obtained segments to peptide

    -pseudokeep      : Only output segments annotated as pseudogene (default:
                       do not output segments annotated as pseudogene)

    -stopkeep        : When used with -translate, retaine segments that contain
                       in-frame stops (default: reject them)

    -minlength       : Reject segments shorter than this length


    -kill_list       : Path to a file containing names of IMGT entries that should
                       be ignored


=head1 EXAMPLES

The following example extracts human C segments, ignoring entries with the "Orphon" 
keyword (i.e. notpart of the core Ig/Tcr sets), and the "Ig-surrogate" keyword 
(i.e. Ig/TcR related gene, but not a gene segment), translating the entries.

constant_regions_from_imgt.pl -tax Homo_sapiens -rejectkeyword orphon -rejectkeyword Ig-surrogate -rejectkeyword Ig-surrogate-vpreb -translate -minlength 150 -kill /path/to/kill_list.txt /path/to/imgt.dat

=cut
