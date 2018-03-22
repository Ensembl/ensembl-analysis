=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::HaplotypeMapper - 

=head1 SYNOPSIS

# This is the main analysis database

    my $genebuilder = new Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection
      (       
       '-hap_slice' => $self->query,
       '-slice'   => $self->target,
       '-input_id' => $self->input_id,
      );



=head1 DESCRIPTION

This module aligned the genomic sequence of a haplotype region and the corresponding
reference chromosome region and projects the gene annotations

=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::HaplotypeMapper;

use warnings ;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Runnable::Blastz;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster;

use Bio::EnsEMBL::Analysis::Config::HaplotypeProjection qw (
                                                 
                                                    );

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($hap_slice,$slice,$input_id) = $self->_rearrange([qw(HAP_SLICE SLICE INPUT_ID)],
					      @args);

    $self->throw("Must input a reference slice and a haplotype slice to HaplotypeProjection") unless (defined($slice) && defined($hap_slice));
    $self->{_final_genes} = [];
    $self->{_gene_types}  = [];

    $self->query($hap_slice);
    $self->target($slice);

    #$self->gene_types($GB_ENSEMBL_INPUT_GENETYPE);
  
    $self->input_id($input_id);

    return $self;
}

############################################################

=head2 input_id

 Function: get/set for input id
 Returns : string
 Args    : string (it expects a string of the format chr_name.start_coord-end_coord

=cut
  
sub input_id {
  my ($self,$id) = @_;
  
  if (defined($id)) {
    $self->{_input_id} = $id;
  }
  return $self->{_input_id};
}
############################################################


=head2 create_alignment

 Description: Generates Blastz alignment between the ref and the hap sequences 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub create_alignment {
  my ($self) = @_;

  my $target_slice = $self->target;
  my $query_slice =  $self->query;

  my $target_top_slice = $self->ensembl_db->get_SliceAdaptor->
      fetch_by_region($target_slice->coord_system->name,
                      $target_slice->seq_region_name,
                      undef,
                      undef,
                      1,
                      $target_slice->coord_system->version);
  
  my $query_top_slice = $self->ensembl_db->get_SliceAdaptor->
      fetch_by_region($query_slice->coord_system->name,
                      $query_slice->seq_region_name,
                      undef,
                      undef,
                      1,
                      $query_slice->coord_system->version);
  
#my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic);
#if (not defined $ana) {
  my $ana = Bio::EnsEMBL::Analysis->new(-logic_name => 'Hap_proj_blastz');
  $self->ensembl_db->get_AnalysisAdaptor->store($ana);
#}
  print "Starting the Blastz alignment. This may take a while so sit and relax.\n";

  my $run = Bio::EnsEMBL::Analysis::Runnable::Blastz->new(
                                                          -analysis => $ana,
                                                          -query => $query_slice, 
                                                          -options => "C=2 B=0",
                                                          -database => [Bio::PrimarySeq->
                                                                        new(
                                                                            -id  => $target_slice->seq_region_name, 
                                                                            -seq => $target_slice->seq)]);
  
  $run->run;
  
  my @hits;
  foreach my $f (@{$run->output}) {
    
    printf("PRELIM: %d %d %d %d %d %d %d\n", 
           $f->start,
           $f->end,
           $f->strand,
           $f->hstart,
           $f->hend,
           $f->hstrand,
           $f->score);
    
    $f->analysis($ana);
    $f->seqname($query_slice->seq_region_name);
    $f->hseqname($target_slice->seq_region_name);
    
    $f->invert($target_slice);
    $f = $f->transfer($target_top_slice);
    $f->invert($query_slice);
    $f = $f->transfer($query_top_slice);

    printf("HIT: %s %d %d %d %s %d %d %d (%d, %.2f)\n", 
           $f->seqname,
           $f->start,
           $f->end,
           $f->strand,
           $f->hseqname,
           $f->hstart,
           $f->hend,
           $f->hstrand, 
           $f->score,
           $f->percent_id);

    # BEWARE this is a partial solution that may not really fix the problem
    # May be better just to skip this cases and not add them to the @hits
    
      push @hits, $f;
      
  }
  
  @hits = sort { $a->start <=> $b->start } @hits;
  
  $self->raw_hits(@hits);

}

############################################################


=head2 filter_alignment

 Description: Filters the alignments between . 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub filter_alignment {
  my ($self) = @_;

  my @hits = $self->raw_hits;

  my @keep = sort { $a->start <=> $b->start } @hits;
    
  @keep = sort { $b->score <=> $a->score } @hits;

  my @retained;
  
  foreach my $kept(@keep){ 
    push @retained, $self->trim_hit($kept, @retained);
  }

  @retained = sort { $a->start <=> $b->start } @retained;

  #Check integrity of retained hit. Check that end is bigger than start
  my @final_hits;

  foreach my $retained (@retained){

    unless($retained->hend < $retained->hstart && $retained->end < $retained->start){
      push (@final_hits,$retained);
    }
    $self->hits(@final_hits);
  }
}

sub trim_hit {
  my ($self, $h, @others) = @_;

  @others = sort { $a->start <=> $b->start } @others;

  # get list of "unique" parts of $h
  my (@ov_regs, @cut_feats);

  foreach my $oh (@others) {
    if (($oh->end >= $h->start and $oh->start <= $h->end) or
        ($oh->hend >= $h->hstart and $oh->hstart <= $h->hend))  { 
      
      my $ov_start = $h->start;
      $ov_start = $oh->start if $oh->start > $ov_start;
      
      my $ov_end = $h->end;
      $ov_end = $oh->end if $oh->end < $ov_end;
      
      my $ov_hstart = $h->hstart;
      $ov_hstart = $oh->hstart if $oh->hstart > $ov_hstart;
      
      my $ov_hend = $h->hend;
      $ov_hend = $oh->hend if $oh->hend < $ov_hend;

      my $f1 = $h->restrict_between_positions($ov_start,
                                              $ov_end,
                                              'SEQ');
      my $f2 = $h->restrict_between_positions($ov_hstart,
                                              $ov_hend,
                                              'HSEQ');
      if (defined $f2) {
        $ov_start = $f2->start if $f2->start < $ov_start;
        $ov_end   = $f2->end   if $f2->end   > $ov_end;
      }
      if (defined $f1) {
        $ov_hstart = $f1->hstart if $f1->hstart < $ov_hstart;
        $ov_hend   = $f1->hend if $f1->hend > $ov_hend;
      }
      
      push @ov_regs, {
        start => $ov_start,
        end   => $ov_end,
        hstart => $ov_hstart,
        hend   => $ov_hend,
      };
    }
  }
  
  foreach my $reg (@ov_regs) {

  }


  if (@ov_regs) {
    if ($ov_regs[0]->{start} > $h->{start}) {
      my $new_f = $h->restrict_between_positions($h->start,
                                                 $ov_regs[0]->{start} - 1,
                                                 'SEQ');
      if (defined $new_f) {
        push @cut_feats, $new_f;
      }
    }
    for(my $i=1; $i < @ov_regs; $i++) {      
      my $new_f =  $h->restrict_between_positions($ov_regs[$i-1]->{end} + 1,
                                                  $ov_regs[$i]->{start} - 1,
                                                  'SEQ');
      if (defined $new_f) { 
        push @cut_feats, $new_f;
      }
    }
    if ($ov_regs[-1]->{end} < $h->{end}) {
      my $new_f =  $h->restrict_between_positions($ov_regs[-1]->{end} + 1,
                                                  $h->end,
                                                  'SEQ');
      if (defined $new_f) {
        push @cut_feats, $new_f;
      }
    }
  } else {
    @cut_feats = ($h);
  }

  my @all = sort { $a->start <=> $b->start } (@others, @cut_feats);
  my $consistent = 1;
  for(my $i=1; $i < @all; $i++) {
    if ($all[$i]->hstart <= $all[$i-1]->hend) {
      $consistent = 0;
      last;
    }
  }



  if ($consistent) {
    return @cut_feats;
  } else {
    return ();
  }
}

 

############################################################


=head2 make_map_regions

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub make_map_regions {
  my ($self) = @_;
  
  my @mappings;
  my @assembly_writes;
  my $cutoff = 50.0;
  
  my $target_slice =  $self->target;
  
  my $tl_target_slice = $self->ensembl_db->get_SliceAdaptor->fetch_by_region($target_slice->coord_system->name,
                                                          $target_slice->seq_region_name);
  
  my $target_dbid = $self->ensembl_db->get_SliceAdaptor->get_seq_region_id($tl_target_slice);
  
  my @genes = @{$target_slice->get_all_Genes};
  @genes = map { $_->transfer($tl_target_slice) } @genes;
  
  my @ex_regs;
  foreach my $g (@genes) {
    foreach my $t (@{$g->get_all_Transcripts}) {
      foreach my $e (@{$t->get_all_Exons}) {
        push @ex_regs, {
          start => $e->start,
          end   => $e->end,
        }
      }
    }
  }
  @ex_regs = sort { $a->{start} <=> $b->{start} } @ex_regs;
  my @nr_ex_regs;
  foreach my $e (@ex_regs) {
    if (not @nr_ex_regs or
        $nr_ex_regs[-1]->{end} < $e->{start} - 1) {
      push @nr_ex_regs, $e;
    } else {
      if ($e->{end} > $nr_ex_regs[-1]->{end}) {
        $nr_ex_regs[-1]->{end} = $e->{end};
      }
    }
  }
  
  my $query_slice = $self->query;
  
  my $query_dbid = $self->ensembl_db->get_SliceAdaptor->get_seq_region_id($query_slice);
  
  my @f = $self->hits;
  
  print STDERR "Fetched ", scalar(@f), " features\n";
  
  @f = grep { $_->percent_id >= $cutoff } @f;
  @f = sort { $a->start <=> $b->start } @f;
  
  my @seq_bits;
  
  foreach my $bit (@{$query_slice->project('seqlevel')}) {
    if (not @seq_bits or $seq_bits[-1]->{end} < $bit->from_start - 1) {
      push @seq_bits, {
        start => $bit->from_start,
        end   => $bit->from_end,
      };
    } else {
      $seq_bits[-1]->{end} = $bit->from_end;
    }
  }
  
  my @ungapped;
  
  foreach my $gf (@f) { 
    #print "YOUR GF is: ",$gf,"\n";
    print "start: ",$gf->start,"\n";
    print "end: ",$gf->end,"\n";
    print "Strand: ",$gf->strand,"\n";
    print "hit_start: ",$gf->hstart,"\n";
    print "hit_end: ",$gf->hend,"\n";
    print "hit_Strand: ",$gf->hstrand,"\n";

    foreach my $f ($gf->ungapped_features) {
      push @ungapped, $f;
    }
  }
  
  # merge consistent ungapped;
  my @merged_ungapped;
  foreach my $f (@ungapped) {
    if (not @merged_ungapped or
        $f->start - $merged_ungapped[-1]->end !=
        $f->hstart - $merged_ungapped[-1]->hend) {
      push @merged_ungapped, $f;
    } else {
      # need to make sure that be merging, we dont extend
      # across a gap
      my $bridge_start = $merged_ungapped[-1]->end + 1;
      my $bridge_end   = $f->start - 1;
      
      my $inside_seq = 0;
      foreach my $b (@seq_bits) {
        if ($bridge_start >= $b->{start} and
            $bridge_end   <= $b->{end}) {
          $inside_seq = 0;
          last;
        }
      }
      
      if ($inside_seq) {
        $merged_ungapped[-1]->end($f->end);
        $merged_ungapped[-1]->hend($f->hend);
      } else {
        push @merged_ungapped, $f;
      }
    }
  }
  
  foreach my $f (@merged_ungapped) {
    push @assembly_writes, [$query_dbid,
                            $target_dbid, 
                            $f->start, 
                            $f->end,
                            $f->hstart, 
                            $f->hend, 
                            1];
  }
  
  my $mapping_path = $query_slice->coord_system->name 
      . ":" 
      . $query_slice->seq_region_name 
      . "#" 
      . $target_slice->coord_system->name 
      . ":" 
      . $target_slice->coord_system->version;
  
  # Update the meta table with the new mapping between the reference chromosome and the haplotype
  my $st = "INSERT into meta(meta_key, meta_value) VALUES(\"assembly.mapping\", \"".$mapping_path."\")"; 
  
  push (@mappings, $st);
  
  my $contig_mapping_path = $query_slice->coord_system->name 
      . ":" 
      . $query_slice->seq_region_name 
      . "#contig"; 

  my $cst = "INSERT into meta(meta_key, meta_value) VALUES(\"assembly.mapping\", \"".$contig_mapping_path."\")"; 
  
  push (@mappings, $cst);


  # calculate coverage of query sequence
  
  # temporary add an entry in the coord_system table;
      
  my $rq = $self->ensembl_db->dbc->prepare("select MAX(rank) from coord_system");
  #print "select MAX(rank) from coord_system\n";
  $rq->execute();
  my @ranks  = $rq->fetchrow_array;
  #print "RANKS: ",@ranks,"\n";
  
  my $rank_val = $ranks[0]+1;
  
  # Create a new coordinate system where the Haplotype name looks like a version of the reference chromosome
  my $cs = "INSERT into coord_system(name,version,rank,attrib) VALUES(\"chromosome\", \"".$query_slice->seq_region_name."\", ".$rank_val.",\" \")";
  push (@mappings, $cs);

  my $total_len = 0;
  my $align_len = 0;
  
  map { $total_len += ($_->{end} - $_->{start} + 1) } @seq_bits;
  map { $align_len += ($_->end - $_->start + 1) } @f;
  
  printf("# Total seq length = %d, total align length = %d, = %.2f\%\n", 
         $total_len,
         $align_len,
         100 * ($align_len / $total_len));
  
  
  # find bits of @seq_bits that are not covered by alignment
  foreach my $f (@f) {
    @seq_bits = $self->remove_from_list($f, \@seq_bits);
  }
  
  foreach my $bit (@seq_bits) {
    printf("# Missing: %d %d (%d)\n", 
           $bit->{start}, 
           $bit->{end}, $bit->{end} - $bit->{start} + 1);
  }
  
  my $total_ex_len = 0;
  my $total_ex_non_aln_len = 0;
  
  map { $total_ex_len += ($_->{end} - $_->{start} + 1) } @nr_ex_regs;
  foreach my $f (@f) {
    @nr_ex_regs = $self->remove_from_list($f, \@nr_ex_regs, 1);
  }
  map { $total_ex_non_aln_len += ($_->{end} - $_->{start} + 1) } @nr_ex_regs;
  printf ("total length = %d \n",$total_ex_len);
  printf ("covered length = %d \n",$total_ex_len - $total_ex_non_aln_len);
  printf ("percent length = %.2f \n",100 * (($total_ex_len - $total_ex_non_aln_len) / $total_ex_len));

  printf("# Total exon length = %d, covered by alignment = %d, = %.2f\n", 
         $total_ex_len, 
         $total_ex_len - $total_ex_non_aln_len,
         100 * (($total_ex_len - $total_ex_non_aln_len) / $total_ex_len));

  #Dirty way to remove duplicates
  my %unique_aw;
  foreach my $all_vals(@assembly_writes){
    #print "All vals: ",@$all_vals,"\n";
    $unique_aw{join(":",@$all_vals)}=1;
  }    
  
  my @unique_assembly_writes;
  foreach my $uniq(keys %unique_aw){
    push (@unique_assembly_writes, $uniq); 
    
  }

  # Insert the new assembly mapping between the reference chromosome and the Haplotype 
  foreach my $el (@unique_assembly_writes) {
    my @real_vals = split(/:/,$el);
    my $stl = "INSERT into assembly values(".$real_vals[0].", ".$real_vals[1].", ".$real_vals[2].", ".$real_vals[3].", ".$real_vals[4].", ".$real_vals[5].", ".$real_vals[6].")";
  push (@mappings, $stl);

  } 

  #Get the new coord region for the Haplotype sequence
  #my $cr = $self->ensembl_db->dbc->prepare("SELECT coord_system_id from coord_system where name ='chromosome' and version =?");
  #$cr->execute($self->query->seq_region_name);
  #my @coord_ids = $cr->fetchrow_array;
  #print "This is what it returns: ",@coord_ids,"\n";
  #$cr->finish;

  # Update the Haplotaype sequence coord_system_id in the seq_region table to be able to perform the projection
  my $sr = "UPDATE seq_region sr, coord_system cs set sr.coord_system_id =cs.coord_system_id where sr.name = cs.version and sr.name =\"".$self->query->seq_region_name."\"";
  push (@mappings, $sr);

  print "TEST REACH POINT 1\n";

  return @mappings;

}

sub remove_from_list {
  my ($self, $f, $bits, $hit_coords) = @_;

  # $f will overlap at most one of @b

  my $start = $f->start;
  my $end   = $f->end;
  if (defined $hit_coords and $hit_coords) {
    $start = $f->hstart;
    $end   = $f->hend;
  }

  my $overlap; 
  my @keep;

  foreach my $bit (@$bits) {
    if ($start <= $bit->{end} and
        $end   >= $bit->{start}) {
      $overlap = $bit;
    } else {
      push @keep, $bit;
    }
  }

  if (defined $overlap) {
    if ($start > $overlap->{start} and 
        $end < $overlap->{end}) {
      my $new_b = {
        start => $end + 1,
        end   => $overlap->{end},
      };      
      push @keep, $new_b;
      $overlap->{end} = $start - 1;
      push @keep, $overlap;
    } elsif ($start > $overlap->{start}) {
      $overlap->{end} = $start - 1;
      push @keep, $overlap;
    } elsif ($end < $overlap->{end}) {
      $overlap->{start} = $end + 1;
      push @keep, $overlap;
    }
  }

  @keep = sort { $a->{start} <=> $b->{start} } @keep;

  return @keep;
}

############################################################


sub check_transcript_in_discarded_db{
  my ($self, $tran) = @_;
 
  my @exons = @{$tran->get_all_Exons};

  my $discardedslice = $self->discarded_db->get_SliceAdaptor->fetch_by_region('toplevel',$tran->slice->seq_region_name,$tran->seq_region_start,$tran->seq_region_end);
  #print STDERR "Fetching discarded genes\n"; 
  #print "NUMBER OF DISCARDED GENES: ",scalar(@{$discardedslice->get_all_Genes}),"\n"; 
  DGENE: 
  foreach my $dgene (@{$discardedslice->get_all_Genes}){
    DTRANS:foreach my $dtran (@{$dgene->get_all_Transcripts}){
      my @dexons = @{$dtran->get_all_Exons};
      if(scalar(@exons) == scalar(@dexons)){
        #print "Number of exons: ",scalar(@exons),"\n";
        for (my $i=0; $i < scalar(@exons); $i++){

          if ($exons[$i]->seq_region_start   != $dexons[$i]->seq_region_start ||
              $exons[$i]->strand  != $dexons[$i]->strand ||
              $exons[$i]->seq_region_end     != $dexons[$i]->seq_region_end){
            # if you enter here means that these two transcripts are not the same
            #print "transcript exon coordinates are different\n";
            next DTRANS;
          }
        }
        # If you are here means that both transcripts are the same and $trans must be discarded
        print "transcript found in discarded db\n";
        return 0;
      }else{
      # if you enter here means that these two transcripts are not the same
        #print "transcript number of exons is different\n";
        next DGENE;
      }
    }
  }
  #If we reach here means that no transcript in the discarded db is the same as our transcript so we keep it
  return 1;
}

sub load_tables {
  my ($self, @mappings) = @_;

  print "REACH TEST 3\n";  
  print "MAppings: ",@mappings,"\n";
  foreach my $map(@mappings){
    print "YOUR SQL: ",$map,"\n";
    my $loading = $self->ensembl_db->dbc->prepare($map);
    $loading->execute();
    $loading->finish;
  }
  
  return 1;
}

############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes,
# havana genes and discarded genes
# this reference is set in Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeProjection

sub ensembl_db{
 my ($self,$ensembl_db) = @_;
 if ( $ensembl_db ){
   $self->{_ensembl_db} = $ensembl_db;
 }
 
 return $self->{_ensembl_db};
}

sub output_db{
 my ($self,$output_db) = @_;

 if ( $output_db ){
   $self->{_output_db} = $output_db;
 }
 
 return $self->{_output_db};
}


sub discarded_db{
  my ($self, $discarded_db) = @_;

  if ( $discarded_db ){
    $self->{_discarded_db} = $discarded_db;;
  }

  return $self->{_discarded_db};
}


############################################################

=head2 gene_types

 Description: get/set for the type(s) of genes (usually TGE_gw, similarity_genewise and combined_e2g genes) 
              to be used in the genebuilder they get set in new()
              Does not include the ab inition predictions
=cut

sub gene_types {
  my ($self,$type) = @_;

  if (defined($type)) {
     push(@{$self->{_gene_types}},$type);
  }

  return @{$self->{_gene_types}};
}

sub features {
  my ($self,@features) = @_;
  
  if (!defined($self->{_feature})) {
    $self->{_feature} = [];
  }
  if ( scalar @features ) {
    push(@{$self->{_feature}},@features);
  }
  return @{$self->{_feature}};
}

############################################################

sub query {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_query} = $slice;
  }
  return $self->{_query};
}

sub target {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_target} = $slice;
  }
  return $self->{_target};
}

sub raw_hits {
  my ($self,@raw_hits) = @_;
  
  if (@raw_hits) {
    push (@{$self->{_raw_hits}}, @raw_hits);
  }

  return @{$self->{_raw_hits}};
}

sub hits {
  my ($self,@hits) = @_;

  if (@hits) {
    push (@{$self->{_hits}}, @hits);
  }

  return @{$self->{_hits}};
}


sub clean_queries {
  my ($self,$sql_q) = @_;
  
  if ($sql_q) {
    push (@{$self->{_clean_queries}}, $sql_q);
  }

  return $self->{_clean_queries};
}

#fetches sequence from appropriate database

sub fetch_sequence{
  my ($self, $name, $db) = @_;

  my $sa = $db->get_SliceAdaptor; 

  my $slice = $sa->fetch_by_name($name);

  return $slice;
}

1;
