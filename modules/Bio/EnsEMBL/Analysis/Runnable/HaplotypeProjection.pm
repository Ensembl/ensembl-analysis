#
#
# BioPerl module for GeneBuilder
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection

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

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection;

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


=head2 project_genes

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub project_genes {
  my($self) = @_;
# Call get genes on target slice to get the genes to project.

  my @projected_genes;
  my $target_slice =  $self->target;
  my $query_slice =  $self->query;
  
  my $full_query_slice = $self->ensembl_db->get_SliceAdaptor->fetch_by_region($query_slice->coord_system->name,
                                                                              $query_slice->seq_region_name,
                                                                              undef,
                                                                              undef,
                                                                              undef,
                                                                              $query_slice->seq_region_name);
  
  
  my $tl_slice = $self->ensembl_db->get_SliceAdaptor->fetch_by_region($target_slice->coord_system->name,
                                                                      $target_slice->seq_region_name);

  print "About to start fetching genes in reference sequence\n";
  my @genes = @{$target_slice->get_all_Genes(undef,undef,1)};
  @genes = map { $_->transfer($tl_slice) } @genes;
  @genes = sort { $a->start <=> $b->start } @genes;

  
# fully load genes
  foreach my $g (@genes) {
    #print "Your genes has: ",scalar(@{$g->get_all_Transcripts})," transcripts\n";
    foreach my $t (@{$g->get_all_Transcripts}) {
      $t->translation;
      #$t->get_all_supporting_features;
      foreach my $e (@{$t->get_all_Exons}) {
        ; #$e->get_all_supporting_features;
      }
    }
  }

  my $projected = 0;
  my $not_projected = 0;

  my $query_version = $full_query_slice->coord_system->version;
  my $query_name = $full_query_slice->seq_region_name;

  print "Processing $query_name...\n";
  
  foreach my $g (@genes) {
    my $tg = Bio::EnsEMBL::Gene->new;
    $tg->analysis($g->analysis);
    $tg->biotype($g->biotype);
    $tg->status($g->status);
    #$tg->stable_id($g->stable_id);

    printf("Transferring gene %s %s %d %d\n", 
           $g->stable_id,
           $g->biotype,
           $g->start,
           $g->end);

    my @proj_trans;

    #print "Number of transcripts before projection: ",scalar(@{$g->get_all_Transcripts}),"\n";
    my $before = scalar(@{$g->get_all_Transcripts});
    foreach my $t (@{$g->get_all_Transcripts}) {
      my $complete_t = $self->get_complete_transcript($t);
      my $tt = $complete_t->transform("chromosome", $query_name);      

      if (defined $tt) {
 
        push @proj_trans, $tt;

        printf(" %s %s %d %d transferred successfully:\n", 
               $t->stable_id,
               $t->biotype,
               $t->start,
               $t->end);
        #$tt->status("CLEAN_TRANSFER");
        #$self->log_compare_transcripts($t, $tt);    
      } else {
        if ($self->transcript_is_missed($t, $query_name)) {
         } else {
          if ($complete_t->translation) {
            my $cds_t = $self->get_coding_part_of_transcript($complete_t);
            
            my $tt = $cds_t->transform("chromosome", $query_name);
            
            if (defined $tt) {
              
            #  $tt->status("CLEAN_CDS_ONLY_TRANSFER");
              push @proj_trans, $tt;
             
            #  $self->log_compare_transcripts($cds_t, $tt);
            } else {

              my $new_t = $self->project_transcript_the_hard_way($complete_t, $query_name, $full_query_slice);
              if ($new_t){
                push @proj_trans, $new_t;
              }
              #  $new_t->status("COMPLEX_CODING_WITHOUT_CDS");
               # $self->log_summarise_projection_by_exon($cds_t, $query_name);
            }
          } else {
            # can project the transcript ignoring CDS
            my $new_t = $self->project_transcript_the_hard_way($complete_t, $query_name, $full_query_slice);    
            if ($new_t){
              push @proj_trans, $new_t;
              
            }
            
            #$new_t->status("COMPLEX_NON_CODING");
            # $self->log_summarise_projection_by_exon($t, $query_name);
          }
        }
      }
    }

    if (@proj_trans) {
      foreach my $proj_t (@proj_trans){
        if ($proj_t->translation){
          $self->check_translation($proj_t);
        }
      }

      map { $tg->add_Transcript($_) } @proj_trans;
      #print "Number of transcripts after projection: ",scalar(@{$tg->get_all_Transcripts}),"\n";
      if($before > scalar(@{$tg->get_all_Transcripts})){
        print "MISSING TRANSCRIPTS\n";
      }
      push(@projected_genes, $tg);
    }

  #  printf("Summary for %s : %d out of %d transcripts transferred\n", 
  #         $g->stable_id, 
  #         scalar(@proj_trans),
  #         scalar(@{$g->get_all_Transcripts}));
  }

  print "Gene projection Finished\n";
  return @projected_genes;

}

sub clean_tables {

  my ($self) = @_;

  my $target_slice =  $self->target;
  my $query_slice =  $self->query;

  print "Cleaning meta table\n";
  my $meta_q = "DELETE from meta where meta_key=\"assembly.mapping\" and meta_value like \"%".$query_slice->seq_region_name."\%\"";

  print "CHECK $meta_q\n";
  my $clean_meta = $self->ensembl_db->dbc->prepare($meta_q);
  $clean_meta->execute();
  $clean_meta->finish;
  
  print "Cleaning meta_coord table\n";
  my $meta_coord = $self->ensembl_db->dbc->prepare("DELETE mt.* from meta_coord mt, coord_system cs where mt.coord_system_id = cs.coord_system_id and cs.version = ?"); 

  $meta_coord->execute($query_slice->seq_region_name);
  $meta_coord->finish;

  print "Cleaning coord_system table\n";
  my $clean_coord_system = $self->ensembl_db->dbc->prepare("DELETE from coord_system where version = ?");
  $clean_coord_system->execute($query_slice->seq_region_name);
  $clean_coord_system->finish;

  my $sr_id = $self->ensembl_db->dbc->prepare("SELECT seq_region_id from seq_region where name =?");
  $sr_id->execute($query_slice->seq_region_name);
  my @hap_sr = $sr_id->fetchrow_array;

  $sr_id->execute($target_slice->seq_region_name);
  my @chr_sr = $sr_id->fetchrow_array;

  $sr_id->finish;

  print "Cleaning assembly table\n";
  my $clean_assembly = $self->ensembl_db->dbc->prepare("DELETE from assembly where asm_seq_region_id = ? and cmp_seq_region_id = ?");
  $clean_assembly->execute($hap_sr[0],$chr_sr[0]);
  $clean_assembly->finish;

  my $old_sr = $self->ensembl_db->dbc->prepare("SELECT coord_system_id from seq_region where name =?");
  $old_sr->execute($target_slice->seq_region_name);
  my @old_coord = $old_sr->fetchrow_array;

  $old_sr->finish;

  print "Cleaning seq_region table\n";
  my $clean_seq_req = $self->ensembl_db->dbc->prepare("UPDATE seq_region set coord_system_id = ? where name = ?");

  $clean_seq_req->execute($old_coord[0],$query_slice->seq_region_name);
  $clean_seq_req->finish;

  print "Table clean up done\n";
}

sub transcript_is_missed {
  my ($self, $t, $hapname) = @_;

  my $found_part = 0;
  foreach my $e (@{$t->get_all_Exons}) {
    my @bits = @{$e->project("chromosome", $hapname)};

    if (@bits) {
      $found_part = 1;
      last;
    }
  }
  
  return not $found_part;
}

sub get_complete_transcript {
  my ($self, $t) = @_;

  my @exons = @{$t->get_all_Exons};
  
  my $full_t = Bio::EnsEMBL::Transcript->new();
  $full_t->add_supporting_features(@{$t->get_all_supporting_features});

  for(my $i=0; $i < @exons; $i++) {
    $exons[$i]->get_all_supporting_features();
    $full_t->add_Exon($exons[$i]);
  }

  $full_t->analysis($t->analysis);
  $full_t->biotype($t->biotype);
  $full_t->slice($t->slice);
  $full_t->status($t->status);

  if($t->translation){
    my @tr_exons = @{$t->get_all_translateable_Exons};
            
    my $tr = Bio::EnsEMBL::Translation->new;
    $tr->start($t->translation->start);
    $tr->end($t->translation->end);
    
    for(my $i=0; $i < @exons; $i++) {
      
      if ($tr_exons[0]->dbID eq $exons[$i]->dbID){
        $tr->start_Exon($exons[$i]);
      }
      if ($tr_exons[-1]->dbID eq $exons[$i]->dbID){
        $tr->end_Exon($exons[$i]);
      }
    }
       
    $full_t->translation($tr);
  }

  for(my $i=0; $i < @exons; $i++) {
    $exons[$i]->{'adaptor'} = '';
    $exons[$i]->{'dbID'} = '';
  }
  
  my $attribute = Bio::EnsEMBL::Attribute->new
      (-CODE => 'SourceTran',
       -NAME => 'source transcript',
       -DESCRIPTION => 'source transcript',
       -VALUE => $t->stable_id);
  
  $full_t->add_Attributes($attribute);
  
  return $full_t;
}


sub get_coding_part_of_transcript {
  my ($self, $t) = @_;

  my @cds = @{$t->get_all_translateable_Exons};
  
  my $cds_t = Bio::EnsEMBL::Transcript->new();
  for(my $i=0; $i < @cds; $i++) {
    $cds[$i]->get_all_supporting_features();
    $cds_t->add_Exon($cds[$i]);
  }
  
  $cds_t->analysis($t->analysis);
  $cds_t->biotype($t->biotype);
  $cds_t->slice($t->slice);
  $cds_t->status($t->status);
  my $tr = Bio::EnsEMBL::Translation->new;
  $tr->start_Exon($cds[0]);
  $tr->start(1);
  $tr->end_Exon($cds[-1]);
  $tr->end($tr->end_Exon->length);
  $cds_t->translation($tr);
  
  return $cds_t;
}


sub project_transcript_the_hard_way {
  my ($self, $t, $hapname, $hap_slice) = @_;

  my $new_t = Bio::EnsEMBL::Transcript->new();
  $new_t->analysis($t->analysis);
  $new_t->biotype($t->biotype);
  $new_t->status($t->status);
            
  my @new_e;

  my @exons;
  if($t->translation){
    @exons = @{$t->get_all_translateable_Exons};
  }else{
    @exons = @{$t->get_all_Exons};
  }

  foreach my $e (@exons) {
    $e->get_all_supporting_features();
    my $te = $e->transform("chromosome", $hapname);
    
    if (defined $te) {
      push @new_e, $te;
    } else {
      my @bits = @{$e->project("chromosome",
                               $hapname)};
      if (@bits) {
        # need to make a new exon from the bits
        my ($ex_st, $ex_en, $strand);
        foreach my $bit (@bits) {
          if (not defined $ex_st or $bit->to_Slice->start < $ex_st) {
            $ex_st = $bit->to_Slice->start;
          }
          if (not defined $ex_en or $bit->to_Slice->end > $ex_en) {
            $ex_en = $bit->to_Slice->end;
          }
          if (not defined $strand) {
            $strand = $bit->to_Slice->strand;
          }
        }
        if (defined $ex_st and defined $ex_en) {
          my $new_e = Bio::EnsEMBL::Exon->new;
          $new_e->start($ex_st);
          $new_e->end($ex_en);
          $new_e->strand($strand);
          $new_e->phase(-1);
          $new_e->end_phase(-1);
          $new_e->slice($hap_slice);
          push @new_e, $new_e;
        }
      }
    }
  }          
  if($t->translation && scalar(@new_e) > 0){
    
    my $tr = Bio::EnsEMBL::Translation->new;
    $tr->start(1);
    $tr->end($new_e[-1]->seq->length);
    $tr->start_Exon($new_e[0]);
    $tr->end_Exon($new_e[-1]);
    
    #print "START: ",$tr->start ," END: ",$tr->end , " EXON1: ",$tr->start_Exon," EXON2: ",$tr->end_Exon,"\n";
    
    $new_t->translation($tr);
    
  }
    
  if (scalar(@new_e) > 0){
    
    map { $new_t->add_Exon($_) } @new_e;
    
    if($new_t->translation){
      $self->fix_coding_phases($new_t);
    }
    return $new_t;
    
  }else{
    
    return 0;
  }
}

sub check_translation {
  my ($self, $trans) = @_;

  my $tseq = $trans->translate();
  if ( $tseq->seq =~ /\*/ ) {
    $trans->{'translation'} = undef;
    $trans->biotype("processed_transcript");
  }
}

sub fix_coding_phases {
  my ($self, $trans) = @_;

  #print "TRANSCRIPT: ", %{$trans},"\n";
  
  if (defined($trans->translation)) {
    my @exons = @{$trans->get_all_translateable_Exons};
    
    my $prev_exon = $exons[0];
    my $start_phase = $exons[0]->phase;
    
    if ($start_phase == -1) {
      $start_phase = 0;
    }

    my @calc_phases;
    my $phase = $start_phase;
    for (my $i=0; $i<scalar(@exons);$i++) {
      push @calc_phases, $phase;
      $phase = (($exons[$i]->length + $phase) % 3);
    }

    for (my $i=0; $i<scalar(@exons);$i++) {
      my $exon=$exons[$i];

      if ($exon->phase != $calc_phases[$i] && $i != 0) {
        print "Had to fix coding exon phase in exon \n";# . $exon->dbID . " in trans " . $trans->dbID . "\n";

        $exon->phase($calc_phases[$i]);

        # Note this condition is here because get_all_translateable_Exons may have made a temporary start or end exon 
        if ($i == $#exons) {
          $trans->translation->end_Exon->phase($calc_phases[$i]);
        }
      }

      my $calc_endphase = (($exon->length + $calc_phases[$i]) % 3);
      if ($exon->end_phase != $calc_endphase && $i != $#exons) {
        print "Had to fix coding exon end_phase in exon\n ";# . $exon->dbID . " in trans " . $trans->dbID . "\n";

        $exon->end_phase($calc_endphase);

        # Note this condition is here because get_all_translateable_Exons may have made a temporary start or end exon 
        if ($i == 0) {
          $trans->translation->start_Exon->end_phase($calc_endphase);
        }
      }

      # Now should be consistent with one another - check
      if ($exon->phase != $prev_exon->end_phase && $i != 0) {
        die "EndPhase/Phase mismatch \n";#(" . $prev_exon->end_phase . " and " . $exon->phase .
                         #") between exons " .  $prev_exon->dbID . " and " . $exon->dbID . "\n";
      }

      $prev_exon = $exon;
    }
  }
}




sub log_summarise_projection_by_exon {
  my ($self, $t, $hapname) = @_;

  my $exons_total = 0;
  my $exons_transferred = 0;
  my $exons_missed = 0;
  my $exons_part_missed = 0;
  my $exons_over_indel = 0;
  
  foreach my $e (@{$t->get_all_Exons}) {
    $exons_total++;

    my $te = $e->transform("chromosome",
                           $hapname);
    
    if (defined $te) {
      $exons_transferred++;
    } else {
      my @bits = @{$e->project("chromosome",
                               $hapname)};

      if (not @bits) {
        $exons_missed++;
      } elsif ($bits[0]->from_start != 1 or
               $bits[-1]->from_end != $e->length) {
        $exons_part_missed++;
      } else {
        $exons_over_indel++;
      }
    }
  }
  
  printf(" %d exons; %d transferred, %d missed, %d part missed, %d over indel\n",
         $exons_total, 
         $exons_transferred, 
         $exons_missed,
         $exons_part_missed,
         $exons_over_indel);
}
  
sub log_compare_transcripts {
  my ($self, $old, $new) = @_;

  my $cdnaseq_old = $old->spliced_seq;
  my $cdnaseq_new = $new->spliced_seq;

  my @exons_old = @{$old->get_all_Exons};
  my @exons_new = @{$new->get_all_Exons};

  my ($pepseq_old, $pepseq_new);
  if ($old->translation) {
    $pepseq_old = $old->translate;
  }
  if ($new->translation) {
    $pepseq_new = $new->translate;
  }

  if ($cdnaseq_old eq $cdnaseq_new) {
    print "  transfer to identical DNA\n";
  } else {
    my $cdna_diffs = $self->compare_seqs($cdnaseq_old, $cdnaseq_new);
    print "  transfer to $cdna_diffs cDNA diffs";
    
    #print $logfh "\nOLD:", $cdnaseq_old, "\n";
    #print $logfh "NEW:", $cdnaseq_new, "\n";
    if (defined $pepseq_old and defined $pepseq_new) {
      #print $logfh "OLD:", $pepseq_old->seq, "\n";
      #print $logfh "NEW:", $pepseq_new->seq, "\n";
      my $prot_diffs = $self->compare_seqs($pepseq_old->seq, $pepseq_new->seq);
      print " and $prot_diffs protein diffs";
      my $stop_count = $pepseq_new->seq =~ tr/\*/\*/;
      if ($stop_count) {
        print " ($stop_count STOPS)";
      }
    }
    print "\n";
  }

}

sub compare_seqs {
  my ($self, $seq1, $seq2) = @_;

  my @seq1 = split(//, $seq1);
  my @seq2 = split(//, $seq2);

  my $total = 0;
  my $diffs = 0;

  for(my $i=0; $i < @seq1; $i++) {
    if ($seq1[$i] ne $seq2[$i]) {
      $diffs++;
    }
    $total++;
  }

  return $diffs;
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

sub hits {
  my ($self,@hits) = @_;
  
  #print "My hits reach this: ",@hits,"\n";

  #if (!defined($self->{_hits})){
  #  $self->{_hits} = [];
  #}
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
