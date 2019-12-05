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

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ProjectedTranscriptEvidence

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::RunnableDB::ProjectedTranscriptEvidence->new(
    );

 $runnable->run;
 my @results = $runnable->output;
 
=head1 DESCRIPTION
Designed to be used after the genes have been projected to the assembly patches and once the final stable IDs are set. 

This aligns the two transcripts (original and projected), creating a 
dna_align_feature, which is then moved from transcript to genomic coords and added
as an extra piece of supporting evidence for the projected transcript.

The original transcripts are fetched via the parent_exon_key transcript attribute added at the time of projection.

=head1 CONTACT

http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProjectedTranscriptEvidence;

use vars qw(@ISA);
use strict;
use warnings;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ProjectedTranscriptEvidence;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(parse_config create_file_name write_seqfile);
use Bio::EnsEMBL::Analysis::Config::General qw (ANALYSIS_WORK_DIR);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);
$| = 1;

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->db->dbc->disconnect_when_inactive(1);

  my( $program, $options, $genedb, $outgenedb) =
    rearrange(['PROGRAM','OPTIONS','GENEDB','OUTGENEDB'], @args);

  # config
  # config precedence is default, param hash, constructor and finally logic_name config
  # read default config entries and do checks
  parse_config($self, $PROJECTED_TRANSCRIPT_EVIDENCE_CONFIG_BY_LOGIC, 'DEFAULT');

  # Defaults are over-ridden by parameters given in analysis table...
  my $ph = $self->parameters_hash;
  $self->PROGRAM($ph->{-program})     if $ph->{-program};
  $self->OPTIONS($ph->{-options})     if $ph->{-options};
  $self->OUTGENEDB($ph->{-outgenedb}) if $ph->{-outgenedb};
  $self->GENEDB($ph->{-genedb})       if $ph->{-genedb};

  # ...which are over-ridden by constructor arguments.
  $self->PROGRAM($program);
  $self->OPTIONS($options);
  $self->OUTGENEDB($outgenedb);
  $self->GENEDB($genedb);

  #Finally, analysis specific config
  #use uc as parse_config call above switches logic name to upper case
  $self->PROGRAM(${$PROJECTED_TRANSCRIPT_EVIDENCE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{PROGRAM});
  $self->OPTIONS(${$PROJECTED_TRANSCRIPT_EVIDENCE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{OPTIONS});
  $self->OUTGENEDB(${$PROJECTED_TRANSCRIPT_EVIDENCE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{OUTGENEDB});
  $self->GENEDB(${$PROJECTED_TRANSCRIPT_EVIDENCE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{GENEDB});

  return $self;
}

# fetch input
sub fetch_input {
  my ($self) = @_;

  my $gene_db = $self->get_dbadaptor($self->GENEDB);
  my $outgene_db = $self->get_dbadaptor($self->OUTGENEDB);

  my $out_sa = $outgene_db->get_SliceAdaptor;
  my $out_ga = $outgene_db->get_GeneAdaptor;
  my $out_ta = $outgene_db->get_TranscriptAdaptor;
  my $ga = $gene_db->get_GeneAdaptor;
  my $ta = $gene_db->get_TranscriptAdaptor;
  $self->outta($out_ta);

  my $patch_slice = $out_sa->fetch_by_name($self->input_id);
  my @out_genes = @{$out_ga->fetch_all_by_Slice($patch_slice)};
  my @out_transcripts;

  #make sure get transcripts outside the patch
  foreach my $out_gene (@out_genes) {
    my @out_trans = @{$out_gene->get_all_Transcripts};
    foreach my $out_t (@out_trans) {
      if ($out_t->analysis->logic_name() =~ m/^proj/) {
        push @out_transcripts, $out_t;
      }
    }
  }

  #get reference slice
  my $ref_slice = get_ref_slice($patch_slice);
  if (!$ref_slice) {
    throw("Cannot find reference slice for patch slice ".$patch_slice->name());
  }

  #create reference transcripts hash
  my @ref_genes = @{$ga->fetch_all_by_Slice($ref_slice)};
  my @ref_transcripts;
  my %ref_trans_hash = ();
  foreach my $ref_gene (@ref_genes){ #make sure get transcripts outside the patch
    my @ref_trans = @{$ref_gene->get_all_Transcripts};
    foreach my $t (@ref_trans) {
      $ref_trans_hash{get_transcript_exon_key($t)} = $t->stable_id;
      push @ref_transcripts, $t;
    }
  }

  $self->out_transcripts(\@out_transcripts);
  print "Fetched ".scalar(@out_transcripts)." projected transcripts\n";
  print "Fetched ".scalar(@ref_transcripts)." reference transcripts\n";

  foreach my $out_transcript (@out_transcripts){
    #get the transcripts (reference and output from projection)

    my @parent_exon_key_attribs = @{$out_transcript->get_all_Attributes('parent_exon_key')};
    my $parent_exon_key = "";

    if (scalar(@parent_exon_key_attribs) > 1) {
      warning("Projected transcript ".$out_transcript->stable_id." has more than 1 parent_exon_key attribute.");
    } elsif (scalar(@parent_exon_key_attribs) > 0) {
      $parent_exon_key = $parent_exon_key_attribs[0]->value();
    }

    my @parent_sid_attribs = @{$out_transcript->get_all_Attributes('parent_sid')};
    my $parent_sid = "";
    if (scalar(@parent_sid_attribs) > 1) {
      warning("Projected transcript ".$out_transcript->stable_id." has more than 1 parent_sid attribute.");
    } elsif (scalar(@parent_sid_attribs) > 0) {
      $parent_sid = $parent_sid_attribs[0]->value();
    }

    if ($parent_exon_key) {

      if (exists($ref_trans_hash{$parent_exon_key})) {

        my $ref_stable_id = $ref_trans_hash{$parent_exon_key};

        # print comparison between previous way of getting the stable id and the current one
        if (!($ref_stable_id eq $parent_sid)) {
          print("Parent stable ID $ref_stable_id for projected transcript ".$out_transcript->stable_id." is different from the parent stable ID set in the old way $parent_sid. This is only for statistical purposes.\n");
        }

        my $ref_transcript = $ta->fetch_by_stable_id($ref_stable_id);

        if ($ref_transcript) {
          my $ref_seq = $ref_transcript->seq;
          my $out_seq = $out_transcript->seq;

          my $ref_seq_file = create_file_name( "refseq_", "fa", $self->ANALYSIS_WORK_DIR );
          my $out_seq_file = create_file_name( "outseq_", "fa", $self->ANALYSIS_WORK_DIR );

          write_seqfile($ref_seq, $ref_seq_file, 'fasta');
          write_seqfile($out_seq, $out_seq_file, 'fasta');

          my %parameters;

          if (not exists( $parameters{-options} )
              and defined $self->OPTIONS) {
            $parameters{-options} = $self->OPTIONS;
          }
          #set up the runnable and add the files to delete
          my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature->new(
            -analysis           => $self->analysis,
            -program            => $self->PROGRAM,
            -query_type         => 'dna',
            -query_file         => $ref_seq_file,
            -query_chunk_number => undef,
            -query_chunk_total  => undef,
            -target_file        => $out_seq_file,
            %parameters,
          );

          $runnable->files_to_delete($ref_seq_file);
          $runnable->files_to_delete($out_seq_file);
          $self->runnable($runnable);
        } else {
          warning("Parent transcript $ref_stable_id could not be fetched for projected transcript ".$out_transcript->stable_id);
        }
      } else {
        warning("Projected transcript ".$out_transcript->stable_id." does not have a parent.");
      }
    } else {
      throw("Projected transcript ".$out_transcript->stable_id." does not have a parent_exon_key attribute. It should have been added during the projection stage!");
    }
  }
}

sub run {
  my ($self) = @_;
  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my @p_transcripts = @{$self->out_transcripts};
  my $pt_count = 0;

  foreach my $runnable (@{$self->runnable}){
    $runnable->run;

    my $features = $runnable->output;
    my @feats = @{$features};
    my $p_transcript = $p_transcripts[$pt_count];
    my $genomic_features =  $self->process_features($features, $p_transcript);

    $self->output($genomic_features);
    $pt_count++;
  }
}


=head2 write_output

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Overrides the write_output method from the superclass.
              Adds the align feature as a TranscriptSupportingFeature
              of the projected transcript.              
  Returntype: 1
  Exceptions: Throws if there isn't one AlignFeature (there should be 
              one and only one, representing the alignment of the two
              Transcripts).

=cut

sub write_output {
  my ( $self ) = @_;
  my @output = @{$self->output};
  my @p_transcripts = @{$self->out_transcripts};
  print "Got " .  scalar(@output) ." genomic features \n";
  warning("Should only be one dna align feature per projected transcript.\n") unless  scalar(@output) == scalar(@p_transcripts);

  my $out_count = 0;

  foreach my $transcript (@p_transcripts){
    my $tsfa = $self->get_dbadaptor($self->OUTGENEDB)->get_TranscriptSupportingFeatureAdaptor;
    my $t_id = $transcript->dbID;
    my $gf = $output[$out_count];
    $tsfa->store($t_id, [$gf]) if ($gf);
    $out_count++;
  }
  print "out_count is: ".$out_count."\n";
}




=head2 process_features

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Uses cdna2genomic to convert the ungapped align 
              features into genomic coords
              The strand shold be 1, if it is not, it does not
              create a Bio::EnsEMBL::DnaDnaAlignFeature object
  Returntype: Arrayref of Bio::EnsEMBL::DnaDnaAlignFeature

=cut

sub process_features {
  my ( $self, $flist, $trans ) = @_;

  # first do all the standard processing, adding a slice and analysis etc.

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my @dafs;
  my $count = 0;
  FEATURE: foreach my $f (@$flist) {
    $count++;
    #my $trans = $self->out_transcript;
    my @mapper_objs;
    my @features;
    my $start = 1;
    my $end = $f->length;
    my $out_slice = $slice_adaptor->fetch_by_name($trans->slice->name);
    # get as ungapped features
    foreach my $ugf ( $f->ungapped_features ) {
      # Project onto the genome
	    foreach my $obj ($trans->cdna2genomic($ugf->start, $ugf->end)){
	      if( $obj->isa("Bio::EnsEMBL::Mapper::Coordinate")){
	        # make into feature pairs
          # NOTE: this is an unusual case. The original and projected transcripts should always align to each other 
          # with both strands = 1
          # This means the hit is always forward and the target gets updated to the strand of the projected transcript.
          if(!($f->hstrand == 1) or !($f->strand == 1)){
            warning("Feature strands are not as expected\n");
            next FEATURE;
          }
	        my $strand = $trans->strand;
	        my $hstrand = $f->hstrand;
	        my $fp;
	        $fp = Bio::EnsEMBL::FeaturePair->new
	          (-start     => $obj->start,
	          -end        => $obj->end,
	          -strand     => $strand,
	          -slice      => $trans->slice,
	          -hstart     => 1,
	          -hend       => $obj->length,
	          -hstrand    => $hstrand,
	          -percent_id => $f->percent_id,
	          -score      => $f->score,
	          -hseqname   => $f->hseqname,
	          -hcoverage  => $f->hcoverage,
	          -p_value    => $f->p_value,
	          );
	        push @features, $fp->transfer($out_slice);
	      }
	    }
    }
    @features = sort { $a->start <=> $b->start } @features;
    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);
    # corect for hstart end bug
    $feat->hstart($f->hstart);
    $feat->hend($f->hend);
    $feat->analysis($self->analysis);
    # transfer the original sequence of the read
    $feat->{"_feature_seq"} = $f->{"_feature_seq"};
    push @dafs,$feat;
  } 
  return \@dafs;
}



###########################################################
# containers


sub out_transcript {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_out_transcript'} = $value;
  }

  if (exists($self->{'_out_transcript'})) {
    return $self->{'_out_transcript'};
  } else {
    return undef;
  }
}

sub out_transcripts {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_out_transcripts'} = $value;
  }

  if (exists($self->{'_out_transcripts'})) {
    return $self->{'_out_transcripts'};
  } else {
    return undef;
  }
}

sub runnables {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_runnabless'} = $value;
  }

  if (exists($self->{'_runnables'})) {
    return $self->{'_runnables'};
  } else {
    return undef;
  }
}


#transcript adaptor to output gene db
sub outta {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_outta'} = $value;
  }

  if (exists($self->{'_outta'})) {
    return $self->{'_outta'};
  } else {
    return undef;
  }
}

sub PROGRAM {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_PROGRAM'} = $value;
  }

  if ( exists( $self->{'_CONFIG_PROGRAM'} ) ) {
    return $self->{'_CONFIG_PROGRAM'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OPTIONS'} ) ) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}

sub ANALYSIS_WORK_DIR {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_ana_dir} = $val;
  }
  return $self->{_ana_dir};
}

sub GENEDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_GENEDB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_GENEDB'})) {
    return $self->{'_CONFIG_GENEDB'};
  } else {
    return undef;
  }
}

sub OUTGENEDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUTGENEDB'} = $value;
  }

  if (exists($self->{'_CONFIG_OUTGENEDB'})) {
    return $self->{'_CONFIG_OUTGENEDB'};
  } else {
    return undef;
  }
}

sub get_transcript_exon_key {
  my $transcript = shift;
  my $string = $transcript->slice->seq_region_name.":".$transcript->biotype.":".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->seq_region_strand.":";

  my $exons = sort_by_start_end_pos($transcript->get_all_Exons);
  foreach my $exon (@{$exons}) {
    $string .= ":".$exon->seq_region_start.":".$exon->seq_region_end;
  } 

  return $string;
}

sub sort_by_start_end_pos {
  my ($unsorted) = @_;

  my @sorted = sort { if ($a->seq_region_start < $b->seq_region_start) {
        return -1;
    } elsif ($a->seq_region_start == $b->seq_region_start) {
      if ($a->seq_region_end < $b->seq_region_end) {
        return-1;
      } elsif ($a->seq_region_end == $b->seq_region_end) {
        return 0;
      } elsif ($a->seq_region_end > $b->seq_region_end) {
        return 1;
      }
        return 0;
    } elsif ($a->seq_region_start > $b->seq_region_start) {
        return 1;
    }
  } @$unsorted;

  return \@sorted;
}

sub get_ref_slice {
  my $patch_slice = shift;
  my $ref_slice;

  print "patch slice is: ".$patch_slice->name."\n";

  my @excs = $patch_slice->get_all_AssemblyExceptionFeatures();

  if (@excs) {
    foreach my $exc (@excs) {
      if (@$exc[0]) {
        if (@$exc[0]->type() !~ m/REF/) {
            $ref_slice = @$exc[0]->alternate_slice();
        }
      }
    }
  }
  print "reference slice is: ".$ref_slice->name."\n";
  return $ref_slice;
}

1;
