
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
Designed to be used after genes have been projected to assembly patches. 

Takes a transcript stable_ID as an input_ID (which immediately after projection 
should be the same in the original reference gene and the copy on the patch).

This aligns the two copies of the transcript (original and projected), creating a 
dna_align_feature, which is then moved from transcript to genomic coords and added
as an extra piece of supporting evidence for the projected transcript.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

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
  $self->db->disconnect_when_inactive(1);

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

  my $ta = $gene_db->get_TranscriptAdaptor;
  my $outta = $outgene_db->get_TranscriptAdaptor;
  $self->outta($outta);


  #get the transcripts (reference and output from projection)
  my $ref_transcript = $ta->fetch_by_stable_id($self->input_id);
  my $out_transcript = $outta->fetch_by_stable_id($self->input_id);
  $self->out_transcript($out_transcript);

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

}

sub run {
  my ($self) = @_;
  throw("Can't run - no runnable objects") unless ( $self->runnable );
  my ($runnable) = @{$self->runnable};
  $runnable->run;

  my $features = $runnable->output;
  my @feats = @{$features};
  
  my $genomic_features = $self->process_features($features);

  $self->output($genomic_features);
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
  print "Got " .  scalar(@output) ." genomic features \n";
  throw("Should only be one dna align feature now.\n") unless  scalar(@output) == 1;

  my $transcript = $self->out_transcript;
  my $tsfa = $self->get_dbadaptor($self->OUTGENEDB)->get_TranscriptSupportingFeatureAdaptor;

  my $t_id = $transcript->dbID;
  $tsfa->store($t_id, \@output);
}




=head2 process_features

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Uses cdna2genomic to convert the ungapped align 
              features into genomic coords
  Returntype: 1

=cut

sub process_features {
  my ( $self, $flist ) = @_;

  # first do all the standard processing, adding a slice and analysis etc.

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my @dafs;
  my $count = 0;
  FEATURE: foreach my $f (@$flist) {
    $count++;
    my $trans = $self->out_transcript;
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
            throw("Feature strands are not as expected\n");
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

1;
