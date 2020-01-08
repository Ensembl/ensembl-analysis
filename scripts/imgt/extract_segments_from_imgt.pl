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
use warnings ;
use strict;

use Getopt::Long;
use Bio::SeqIO;

my (@good_tax_strings,
    @bad_tax_strings,
    @moltype,
    @dataclass,
    @good_keywords,
    @bad_keywords,
    @tagstrings,
    $kill_list, %kill_list,
    $gene_name_in_id,
    $translate,
    $keep_pseudos,
    $keep_partials,
    $keep_stops,
    $guess_frame,
    $n_stop_allow,
    $c_stop_allow,
    $min_length,
    $flanking_seq);


$min_length = 0;
$n_stop_allow = 0;
$c_stop_allow = 0;

&GetOptions(
            'taxstring=s@'       => \@good_tax_strings,
            'rejecttaxstring=s@' => \@bad_tax_strings, 
            'moltype=s@'         => \@moltype,
            'class=s@'           => \@dataclass,
            'tagstring=s@'       => \@tagstrings,
            'keyword=s@'         => \@good_keywords,
            'rejectkeyword=s@'   => \@bad_keywords,
            'genenameinid'       => \$gene_name_in_id,
            'translate'          => \$translate,
            'partialkeep'        => \$keep_partials,
            'pseudokeep'         => \$keep_pseudos,
            'stopkeep'           => \$keep_stops,
            'guessframe'         => \$guess_frame,
            'Nstopallow=s'       => \$n_stop_allow,
            'Cstopallow=s'       => \$c_stop_allow,
            'minlength=s'        => \$min_length,
            'flank=s'            => \$flanking_seq,
            'kill_list=s'        => \$kill_list,
            );

die "You must spply a tag string\n" if not @tagstrings;
die "You cannot use -flank with -translate\n" if $flanking_seq and $translate;

my (%moltype, %dataclass, %good_keywords, $bad_keywords);

@good_tax_strings = map { s/\_/ /g; lc($_) } @good_tax_strings;
@bad_tax_strings  = map { s/\_/ /g; lc($_) } @bad_tax_strings;
@good_keywords    = map { s/\_/ /g; lc($_) } @good_keywords;
@bad_keywords     = map { s/\_/ /g; lc($_) } @bad_keywords;

map { $moltype{$_}   = 1 } map { s/_/ /g; lc($_) } @moltype;
map { $dataclass{$_} = 1 } map { s/_/ /g; lc($_) } @dataclass;

my (%all_tags);
foreach my $tag (@tagstrings) {
  if ($tag =~ /:/) {
    map { $all_tags{$_} = 1 } split(/:/, $tag);
  } else {
    $all_tags{$tag} = 1;
  }
}

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

  if (@good_tax_strings) {
    next SEQ if not defined $seq->species;
    my @taxonomy = map { lc($_) } reverse $seq->species->classification;
    
    pop @taxonomy; push @taxonomy, lc($seq->species->species);
    my $taxstring = join(" ", @taxonomy);
    next SEQ if not grep { index($taxstring, $_) >= 0 } @good_tax_strings;
  }
  
  if (@moltype) {
    next SEQ if not exists $moltype{lc($seq->molecule)};
  }
  if (@dataclass) {
    next SEQ if not exists $dataclass{lc($seq->data_class)};
  }
  
  my @keywords = split(/\s*;\s*/, $seq->description);
  #push @keywords, $seq->get_keywords;
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
  
  my (%segments, @gene_regions, @segment_lists);
  
  foreach my $val ($seq->get_SeqFeatures) {
    my $location = $val->location;
    
    if ($val->primary_tag =~ /^(gene)|(\S+\-GENE)$/) {
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
    
    if (exists $all_tags{$val->primary_tag}) {                
      my $seg =  {
        code         => $val->primary_tag,
        location     => $location,
        codon_start  => 1,
        pseudo       => $entry_is_pseudo,
        partial      => 0,
      };

      if ($val->annotation->get_Annotations("gene")) {
        my ($tag_val) = $val->annotation->get_Annotations("gene");
        $seg->{gene} = $tag_val->value if $tag_val->value !~ /\s+/;
      }
      if ($val->annotation->get_Annotations("partial") or
          $location->isa("Bio::Location::Fuzzy")) {
        $seg->{partial} = 1;
      }
      if ($val->annotation->get_Annotations("codon_start")) {
        my ($tag_val) = $val->annotation->get_Annotations("codon_start");
        $seg->{codon_start} = $tag_val->value;            
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
        }
      }      

      push @{$segments{$val->primary_tag}}, $seg;      
    }
  }
  
  @gene_regions = sort { 
    $a->{location}->start <=> $b->{location}->start; 
  } @gene_regions;
    
  TAG: foreach my $tagstring (@tagstrings) {    
    if (exists $segments{$tagstring}) {
      foreach my $seg (@{$segments{$tagstring}}) {
        my @matching_gregs;
        foreach my $greg (@gene_regions) {
          if ($seg->{location}->start >= $greg->{location}->start and
              $seg->{location}->end <= $greg->{location}->end) {
            push @matching_gregs, $greg;
          } 
        }
        map { $seg->{parent_genes}->{$_} = $_ } @matching_gregs; 
        
        push @segment_lists, [$seg];
      }        
    } elsif ($tagstring =~ /:/) {
      my @bits = split(/:/, $tagstring);
      my %bits = map { $_, 1 } @bits;
      
      if (not (grep { not exists $segments{$_} } @bits)) {
        # all components of the composite are represented...
        my (%seg_groups);
        
        if (not @gene_regions) {
          # assume that this is a "single entry" sequence
          # check that we can identify exactly one segment 
          # for each element of our composite
          my $relevant = 0;
          map { $relevant += scalar(@{$segments{$_}}) } keys %bits;
          if ($relevant == scalar(@bits)) {    
            foreach my $bit (keys %bits) {
              push @{$seg_groups{'single_gene'}}, @{$segments{$bit}};
            }
          }
          
        } else {
          # need to segregate segments according to gene locus coords
          
          foreach my $bit (keys %bits) {
            if (exists $segments{$bit}) {
              
              foreach my $seg (@{$segments{$bit}}) {
                # identify gene region that it belongs to

                my @matched_gregs;
                foreach my $greg (@gene_regions) {
                  if ($seg->{location}->start >= $greg->{location}->start and
                      $seg->{location}->end <= $greg->{location}->end) {
                    push @matched_gregs, $greg;
                  } elsif ($greg->{location}->start > $seg->{location}->end) {
                    last;
                  }
                }

                map { $seg->{parent_genes}->{$_} = $_ } @matched_gregs; 

                if (scalar(@matched_gregs) == 1) {
                  push @{$seg_groups{$matched_gregs[0]}}, $seg;
                } else {
                  # this segment cannot be matched to exactly one gene
                  # region, so skip it
                  next;
                }
              }
            }
          }
        }
        
        
        GENE: foreach my $gid (keys %seg_groups) {
          # check completeness
          my @segs = sort { 
            $a->{location}->start <=> $b->{location}->start; } 
          @{$seg_groups{$gid}};
          
          # we should have exactly one segment for each component of the tagstring
          if (scalar(@segs) != scalar(@bits)) {
            next GENE;
          }

          push @segment_lists, \@segs;
        }
      }
    }
  }
  
  # reject segment lists that have any partial components
  my @complete_segment_lists;
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
    
    push @complete_segment_lists, $seglist;
  }
  
  # Make segment list non-redundant, preferring lists with higher
  # total length
  
  @segment_lists = sort {
    my $alen = 0;
    my $blen = 0;
    map { $alen += $_->{location}->end - $_->{location}->start + 1 } @$a;
    map { $blen += $_->{location}->end - $_->{location}->start + 1 } @$b;

    $blen <=> $alen;
  } @complete_segment_lists;
  
  my @nr_seg_lists;
  foreach my $list (@segment_lists) {
    my @locs;
    foreach my $seg (@$list) {
      foreach my $loc ($seg->{location}->each_Location) {
        push @locs, $loc;
      }
    }
    @locs = sort { $a->start <=> $b->start } @locs;

    my $overlaps = 0;
    CHECK: foreach my $olist (@nr_seg_lists) {
      my @olocs;
      foreach my $seg (@$olist) {
        foreach my $loc ($seg->{location}->each_Location) {
          push @olocs, $loc;
        }
      }
      @olocs = sort { $a->start <=> $b->start } @olocs;

      foreach my $loc (@locs) {
        foreach my $oloc (@olocs) {
          if ($loc->start <= $oloc->end and
              $loc->end   >= $oloc->start) {
            $overlaps = 1;
            last CHECK;
          } elsif ($oloc->start > $loc->end) {
            last;
          }
        }
      }
    }
    if (not $overlaps) {
      push @nr_seg_lists, $list;
    }
  }
  
  @nr_seg_lists = sort { 
    $a->[0]->{location}->start <=> $b->[0]->{location}->start ;
  } @nr_seg_lists;
  
  # Finally, output the gene segments

  my (%used_names);
  
  SEGLIST: foreach my $seglist (@nr_seg_lists) {
    # many reverse-strand entries in the IMGT flat file are
    # not designated as such. If a translation has been 
    # supplied in one or more of the segments, we can try to 
    # guess the strand.     

    my @segs = sort { $a->{location}->start <=> $b->{location}->start } @$seglist;
    my $strand = $segs[0]->{location}->strand;

    if ($strand >= 0) {
      $strand = &determine_strand($seq, $seglist);
    }
    
    if ($translate) {
      &adjust_for_bad_codon_boundaries($strand, $seglist);
    }

    my $codestring = join(":", map { $_->{code} } @segs); 

    my $seqstring = ""; 
    my (%gnames, $tr_start, $tr_end);

    foreach my $seg (@segs) {      
      my $position = length($seqstring) + 1;
      foreach my $loc (sort { $a->start <=> $b->start } $seg->{location}->each_Location) {
        $seqstring .= $seq->subseq($loc->start, $loc->end);
      }

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
    }

    my $frame = 0;
    if ($strand > 0 and exists $seglist->[0]->{codon_start}) {
      $frame = $seglist->[0]->{codon_start} - 1;
    } elsif ($strand < 0 and exists $seglist->[-1]->{codon_start}) {
      $frame = $seglist->[-1]->{codon_start} - 1;
    }

    if ($flanking_seq) {
      my $left_en = $segs[0]->{location}->start - 1;
      my $left_st = $left_en - $flanking_seq + 1;
      next SEGLIST if $left_st < 1;
      my $left_seq = $seq->subseq($left_st, $left_en);

      my $right_st = $segs[-1]->{location}->end + 1;
      my $right_en = $right_st + $flanking_seq - 1;
      next SEGLIST if $right_en > $seq->length;
      my $right_seq = $seq->subseq($right_st, $right_en);

      $tr_start = $flanking_seq + 1;
      $tr_end   = $tr_start + length($seqstring) - 1;
      $seqstring = ($left_seq . $seqstring . $right_seq);

      # trim to flush codon boundaries
      $tr_start += $frame;
      $tr_end -= (($tr_end - $tr_start + 1) % 3);
    }
    
    my $newseq = Bio::PrimarySeq->new(-seq => $seqstring);

    if ($strand < 0) {
      $newseq = $newseq->revcom;
      if (defined $tr_start and defined $tr_end) {
        ($tr_start, $tr_end) = ($newseq->length - $tr_end + 1, $newseq->length - $tr_start + 1);
      }
    }

    ####
    # set unique name

    my $id = $seq->primary_id 
        . "." . $segs[0]->{location}->start . "-" . $segs[-1]->{location}->end;
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
    $newseq->id($id);
      
    ###
    # set meaningful description

    my @this_desc = ("tag=$codestring", map {"gene=$_"} @gnames);
    if (grep { $_->{pseudo} } @$seglist) {
      push @this_desc, "Pseudogene"; 
    }
    if (defined $tr_start and defined $tr_end and not $translate) {
      push @this_desc, "cds=$tr_start-$tr_end";
    }
    $newseq->description(join("; ", @this_desc));
    if ($translate) {
		  unless ($newseq->seq =~ /^[AGCTagct]+$/){
		  next SEQ;
			}
			#print "SEQUENCE: ",$newseq->seq,"\n";
      my $pepseq = &translate($newseq, $frame);
      # we optionally allow stop-containing translations if the
      # stops are within the given distance of the terminus; these
      # stops can be "spliced out" during the rearrangement reaction      
      if ($pepseq->seq !~ /\*/ or 
          $keep_stops or 
          $pepseq->seq !~ /^.{$n_stop_allow}.*\*.*.{$c_stop_allow}$/) {
        $outio->write_seq($pepseq) 
            if $min_length <= $pepseq->length;
      } elsif ($guess_frame) {
        my %frame; 
        map { $frame{$_} = 1 } (0,1,2); 
        delete $frame{$frame};
        my @good_seqs;

        foreach my $fr (keys %frame) {
          my $pep = &translate($newseq, $fr);

          if ($pep->seq !~ /\*/ or 
              $keep_stops or 
              $pep->seq !~ /^.{$n_stop_allow}.*\*.*.{$c_stop_allow}$/) {
            push @good_seqs, $pep;
          }
        }
        if (scalar(@good_seqs) == 1) {
          $outio->write_seq($good_seqs[0]) 
            if $min_length <= $good_seqs[0]->length;
        }
      }
    } else {
      $outio->write_seq($newseq)
          if $min_length <= $newseq->length;
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

    if (exists $seg->{codon_start}) {
      $forward = substr($forward, $seg->{codon_start} - 1);
      $reverse = substr($reverse, $seg->{codon_start} - 1);
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




sub translate {
  my ($seq, $frame) = @_;

  my $trans_seq = $seq->translate();

  if ($trans_seq->seq =~ /\*$/) {
    my $pep = $trans_seq->seq;
    $pep =~ s/\*$//;

    $trans_seq->seq($pep);
  }

  return $trans_seq;
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
      my $over5 = (exists($seg->{codon_start})) ? $seg->{codon_start} - 1 : 0;
      
      if (($prev_over + $over5) % 3) { 
        # bad; just trim the excess off each segment; 
        if ($strand < 0) {
          $seg->{location}->end($seg->{location}->end - $over5);
          $segs[$i-1]->{location}->start($segs[$i-1]->{location}->start + $prev_over);
        } else {
          $seg->{location}->start($seg->{location}->start + $over5);
          $segs[$i-1]->{location}->end($segs[$i-1]->{location}->end - $prev_over);
        }
        delete $seg->{codon_start};
      }
    }
    my $len = $seg->{location}->end - $seg->{location}->start + 1;
    $len -= ($seg->{codon_start} - 1) if exists($seg->{codon_start});
    $prev_over = $len % 3;
  }
  
  return 1;
}

=pod

=head1 NAME

extract_segments_from_imgt.pl

=head1 SYNOPSIS

This script extracts sequence segments from the given IMGT database 
flat-file that satisfy the specified conditions. The results are written
in FASTA format to stdout

=head1 OPTIONS

    (Note: when giving an option value that contains maore thanone word, 
     e.g. homo sapiens, spaces should be replaced with "_", i.e. homo_sapiens)

    -tagstring       : Extract segments with these feature table (FT) keys. Single
                       segments formed by concatenation of different features can
                       be achieved by combinging their keys with a ":"

    -taxstring       : Restrict to entries with taxonomy containing this pattern

    -rejecttaxstring : Ignore entries with taxonomy containing this pattern

    -moltype         : Restrict to entries with this moledule type

    -class           : Restrict to entries with this entry class (where IMGT classes, 
                       from least-annotated to most-aanotated, are 'keyword_level', 
                       'automatic' and 'by_annotators')

    -keyword         : Restrict to entries containing these keywords
    
    -rejectkeyword   : Ignore entries including these keywords
                       
    -translate       : Translate the obtained segments to peptide

    -partialkeep     : Retain partial segments (default: reject them)
    
    -pseudokeep      : Only output segments annotated as pseudogene (default:
                       do not output segments annotated as pseudogene)

    -stopkeep        : When used with -translate, retaine segments that contain
                       in-frame stops (default: reject them)

    -Nstopallow      : When used with -translate, is a stop occurs within this 
                       many residues of the N-terminus, retain the segment

    -Cstopallow      : When used with -translate, is a stop occurs within this 
                       many residues of the C-terminus, retain the segment

    -minlength       : Reject segments shorter than this length

    -flank           : Include this much flanking sequence with the segment; the 
                       location of the segments in the result sequences will included
                       in their header descriptions

    -kill_list       : Path to a file containing names of IMGT entries that should
                       be ignored

    -genenameinid    : Include in the header description a gene name (or names), 
                       if any can be inferred from the entry

=head1 EXAMPLES

The following example extracts human L-V segments, scanning only the automatic
and by_annotators entries, ignoring entries with the "Orphon" keyword (i.e. not
part of the core Ig/Tcr sets), and the "Ig-surrogate" keyword (i.e. Ig/TcR related
gene, but not a gene segment), translating the entries, allowing in-frame stops
within 2 residues of the C-terminus of the translation:

extract_segments_from_imgt.pl -tag L-PART1:L-PART2:V-REGION -tag L-REGION:V-REGION -tag L-V-REGION -tax Homo_sapiens -rejectkeyword orphon -rejectkeyword Ig-surrogate -rejectkeyword Ig-surrogate-vpreb -class automatic -class by_annotators -translate -minlength 30 -Cstopallow 2 -kill /path/to/kill_list.txt /path/to/imgt.dat

=cut
