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

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaLocalAlignment - 

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaLocalAlignment->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->fetch_input();
$runnableDB->run();
$runnableDB->write_output(); #writes to DB

=head1 DESCRIPTION

Extends Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa to allow
reads to be realigned against a small piece of genomic with high sensitivity

=head1 METHODS


=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaLocalAlignment;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateSolexaLocalAlignment;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexaLocalAlignment;
use Bio::EnsEMBL::FeaturePair;
use vars qw(@ISA);

@ISA =  ("Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa");


sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($EXONERATE_SOLEXA_LOCAL_ALIGNMENT_CONFIG_BY_LOGIC);
  return $self;
}

=head2 fetch_input

  Title   :   fetch_input
  Usage   :   $self->fetch_input
  Function:   Fetches Trancripts and exons for the rough model stable_id specified in the input_id
              Checks to see if it leads to split refined models
	      Fetches reads from database to realign
  Returns :   none
  Args    :   gene stable id
=cut

sub fetch_input {
  my ($self) = @_;
  # do all the normal stuff
  # then get all the transcripts and exons
  my $gene_db = $self->get_dbadaptor($self->TRANSDB);
  my $gene_adaptor = $gene_db->get_GeneAdaptor;
  my %gene_by_id;
  my $biotype = $self->BIOTYPE;
  my @gene;
  my %keep_analyses;
  my $score = $self->SCORE;
  # fetch genes, transcripts
  if ( $biotype ){
    push @gene , @{$gene_adaptor->generic_fetch("biotype = \"$biotype\"")};
  } else {
    push @gene , @{$gene_adaptor->fetch_all(undef,undef,undef)};
  }
  
  foreach my $gene ( @gene ) {
    $gene_by_id{$gene->stable_id} = $gene;
  }
 
  # fetch the rough model from the inout id
  my $model = $gene_by_id{$self->input_id};
  print "Found model " . $model->stable_id . " ". 
    $model->seq_region_name ." " . 
       $model->start ." " . 
	 $model->end ." " .
	   $model->strand . "\n";
  # fetch the reads overlapping this model that are only partly aligned
  my $genomic_db = $self->get_dbadaptor($self->GENOMICDB);
  my $daf_adaptor = $genomic_db->get_DnaAlignFeatureAdaptor;
  my $slice_adaptor = $genomic_db->get_SliceAdaptor;
  $self->slice_adaptor($slice_adaptor);
  my $slice = $slice_adaptor->fetch_by_region('toplevel',$model->seq_region_name,$model->start,$model->end);
  $self->query($slice);

  # figure out where the break in the models are that we want to fix
  # get all the models with the appropriate logic name and look for possible breaks
  # these are the slices from which to look for partial alignments
  my $refined_db = $self->get_dbadaptor($self->REFINED_DB);
  my $refined_ga = $refined_db->get_GeneAdaptor;
  my @refined_plus;
  my @refined_minus;
  my @tmp;
   if ( $self->REFINED_LN ){
    push @tmp , @{$refined_ga->fetch_all_by_Slice($slice,$self->REFINED_LN)};
  } else {
    push @tmp , @{$refined_ga->fetch_all_by_Slice($slice,undef,1)};
  }
  # filter   
  my $iid = $self->input_id ;
  foreach my $g ( @tmp ) {
    next unless $g->stable_id =~ /^$iid.+/;
    if ( $g->strand == 1 ) {
      push @refined_plus,$g;
    } else {
      push @refined_minus,$g;
    }
  }
  print "Got " . scalar ( @refined_plus ) . " forward and " .  scalar ( @refined_minus ) . " reverse refined genes \n";
  my @daf_slices;
  # now we need to order them
  @refined_plus = sort { $a->start <=> $b->start }   @refined_plus ;
  @refined_minus = sort { $a->start <=> $b->start }   @refined_minus ;
  # look for places where we have a gap between the genes on the same strand that could be a split gene
  # look in that gap for partially aligned reads you could realign in a more sentitive way to find missing
  # introns do we want a limit to the size of the gap? probably
  if ( scalar  @refined_plus > 1 ) {
    for ( my $i = 1 ; $i < scalar(@refined_plus) ; $i++ ) {
      my $left_gene  = $refined_plus[$i-1];
      my $right_gene = $refined_plus[$i];
      my $gap = $right_gene->start - $left_gene->end +1;
      if ( $gap > 0  && $gap < $self->MAX_GAP) {
	my @right_exons = sort  { $a->start <=> $b->start } @{$right_gene->get_all_Exons};
	my @left_exons = sort  { $a->start <=> $b->start } @{$left_gene->get_all_Exons};
	my $new_slice = $slice->sub_Slice($left_exons[-1]->start,$right_exons[0]->end,1);
	print $new_slice->name ."\n";
	push @daf_slices,$new_slice;
      }
    }
  }

  my %logic_names;
  my @dafs;
  my @all_dafs;
  foreach my $daf_slice ( @daf_slices ) {
    # single reads
    @all_dafs = @{$daf_adaptor->fetch_all_by_Slice_constraint($daf_slice,"score < $score")};
    push  @all_dafs , @{$daf_adaptor->fetch_all_by_Slice_constraint($daf_slice,"score < ( $score *2  ) AND right(hit_name,1) not in ('a','b')")};
  }
  if ( $self->LOGIC_NAMES ) {
    foreach my $ln ( @{$self->LOGIC_NAMES} ) {
      $logic_names{$ln} = 1;
    }
    foreach my $daf ( @all_dafs ) {
      push @dafs, $daf if $logic_names{$daf->analysis->logic_name};
    }
  } else {
    @dafs = @all_dafs;
  }   

  $self->reads(\@dafs);
  # store all the analyses so we can pair the reads up with the correct
  # alignments post exonerate
  foreach my $daf ( @all_dafs ) {
     my $name = $daf->hseqname;
     if ( $name =~ /(\S+):(a$|b$|a3p$|b3p|3p)/ ) {
	$name = $1;
      } 
     #   print "NAME $name\n";
      $keep_analyses{$name} = $daf->analysis;
  }
  $self->analysis_store(\%keep_analyses);
  # paired reads
  print "Got " . scalar(@dafs) . " dafs returned  in total \n"; 
  if ( scalar(@dafs) == 0 ) {
    print "No reads found - bailing out\n";
    return;
  }
  # set up the runnable 
  my %parameters = %{ $self->parameters_hash };
  
  if (not exists( $parameters{-options} )
      and defined $self->OPTIONS) {
    $parameters{-options} = $self->OPTIONS;
  }
  
  #print STDERR "PROGRAM FILE: ".$self->analysis->program_file."\n";
  my $program = $self->PROGRAM;
  $program = $self->analysis->program_file if not defined $program;
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateSolexaLocalAlignment->new
    (
     -analysis   => $self->analysis,
     -program    => $program,
     -query_type => $self->QUERYTYPE,
     -reads      => \@dafs,
     -genomic    => $slice,
     %parameters,
    );
  
  $self->runnable($runnable);
  
  return ;
}

=head2 run

  Args       : none
  description: Runs ExonerateSolexa with the reads found from the indexes
  Returntype : none

=cut 

sub run {
  my ($self) = @_;
  return unless  scalar(@{$self->reads} > 0 ) ;
  my ($runnable) = @{$self->runnable};
  # do all the normal stuff
#  $self->SUPER::run();
  $runnable->run();
  my $features = $runnable->output;
  # clear the output arrray
  $self->{'output'} = [];
  my $intron_features = $self->process_features($features);
  $self->output($intron_features);
}

=head2 process_features

  Args       : List of dna_align_features
  description: Checks splice sites
  Returntype : none

=cut 


sub process_features {
  my ( $self, $flist ) = @_;
  my @dafs;
  my $analyses = $self->analysis_store;
  my $chr = $self->slice_adaptor->fetch_by_region('toplevel',$self->query->seq_region_name);
  print $chr->display_id ."\n";
 FEATURE:  foreach my $f (@$flist) {
     # move it from the slice to the chromosome
    $f->slice($self->query);
    $f = $f->transfer($chr);

    my $canonical = 1;
    my $accept = 0;
    $accept = 1 unless   $self->INTRON_MODELS ;
    # restrict just to splice models where read is within an exon
    if ( $self->INTRON_MODELS ) {
      $accept = 1 if $f->{"_intron"};
    }
    next FEATURE unless $accept;
    my @features = $f->ungapped_features;
    @features = sort { $a->start <=> $b->start } @features;
    # if we have a spliced alignment check to see if it's non-canonical
    # if so tag it so we can tell later on
    if ( $f->{'_intron'} ) {
      if ( scalar(@features) == 2 ) {
	my $left_splice = $self->slice_adaptor->fetch_by_region('toplevel',
								$self->query->seq_region_name,
								$features[0]->end+1,
								$features[0]->end+2,
								$features[0]->strand
							       );
	my $right_splice = $self->slice_adaptor->fetch_by_region('toplevel',
								 $self->query->seq_region_name,
								 $features[1]->start-2,
								 $features[1]->start-1,
								 $features[0]->strand
								);	
	if ( $left_splice->seq eq 'NN' && $right_splice->seq eq 'NN' ) {
	  warn("Cannot find dna sequence for " . $f->display_id .
	       " this is used in detetcting non canonical splices\n");
	} else {
	  
	  # is it canonical
	  if ( $features[0]->strand  == 1 ) {
	    print "Splice type " . $left_splice->seq ."- ".  $right_splice->seq ."\n";
	    # is it GTAG?
	    unless ( $left_splice->seq eq 'GT' && $right_splice->seq eq 'AG' ) {
	      $canonical = 0;
	    }
	  } else {
	    print "Splice type " . $right_splice->seq ."- ".  $left_splice->seq ."\n";
	    # is it GTAG?
	    unless ( $right_splice->seq eq 'GT' && $left_splice->seq eq 'AG' ) {
	      $canonical = 0;
	    }
	  }
	}
      }
    }
    print "Making feat " . $f->hseqname . "\n";
    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);
    # put the analysis back
     my $name = $feat->hseqname;
     if ( $name =~ /(\S+):(a:aa$|b:bb$|a$|b$)/ ) {
	$name = $1;
      }
   if ( $analyses->{$name} ) {
      $feat->analysis($analyses->{$name});
      } else {
      	$self->throw("Read not found " . $feat->hseqname ." using key $name \n");
      }
    unless ( $canonical ) {
      print "Non canonical \n";
      # mark it as non canonical splice
      $feat->hseqname($feat->hseqname.":NC");
    } else {
      print "Canonical \n";
    }
    push @dafs,$feat;
  }
  return \@dafs;
}


##########################################################################

sub reads {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_reads'} = $value;
  }

  if (exists($self->{'_reads'})) {
    return $self->{'_reads'};
  } else {
    return undef;
  }
}

sub slice_adaptor {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_slice_adaptor'} = $value;
  }

  if (exists($self->{'_slice_adaptor'})) {
    return $self->{'_slice_adaptor'};
  } else {
    return undef;
  }
}


sub analysis_store {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_analysis_store'} = $value;
  }

  if (exists($self->{'_analysis_store'})) {
    return $self->{'_analysis_store'};
  } else {
    return undef;
  }
}




sub GENOMICDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_GENOMICDB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_GENOMICDB'})) {
    return $self->{'_CONFIG_GENOMICDB'};
  } else {
    return undef;
  }
}


sub LOGIC_NAMES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_LOGIC_NAMES'} = $value;
  }
  
  if (exists($self->{'_CONFIG_LOGIC_NAMES'})) {
    return $self->{'_CONFIG_LOGIC_NAMES'};
  } else {
    return undef;
  }
}

sub REFINED_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_REFINED_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_REFINED_DB'})) {
    return $self->{'_CONFIG_REFINED_DB'};
  } else {
    return undef;
  }
}

sub REFINED_LN {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_REFINED_LN'} = $value;
  }
  
  if (exists($self->{'_CONFIG_REFINED_LN'})) {
    return $self->{'_CONFIG_REFINED_LN'};
  } else {
    return undef;
  }
}

sub MAX_GAP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_GAP'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_GAP'})) {
    return $self->{'_CONFIG_MAX_GAP'};
  } else {
    return undef;
  }
}

sub SCORE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SCORE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_SCORE'})) {
    return $self->{'_CONFIG_SCORE'};
  } else {
    return undef;
  }
}
1;
