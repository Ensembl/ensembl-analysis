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
=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BlastRNASeqPep

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBAdaptor->new($locator);
my $btpep     = Bio::EnsEMBL::Analysis::RunnableDB::BlastRNASeqPep->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$btpep->fetch_input();
$btpep->run();
$btpep->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Analysis::Runnable::Blast on peptides
obtained by translating a representative transcript from each gene
in the region. The resulting blast hits are written back as
DnaPepAlignFeature's.

The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters.

=head1 CONTACT


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastRNASeqPep;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptPep;
use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;
use Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastRNASeqPep;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptPep Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($BLASTRNASEQPEP_CONFIG_BY_LOGIC);
  return $self;
}

=head2 fetch_input

  Args       : none
  Example    : $runnable->fetch_input
  Description: Fetches input data for BlastRNASeqPep and makes runnable
  Returntype : none
  Exceptions : $self->input_id is not defined
  Caller     : run_RunnableDB, Bio::EnsEMBL::Pipeline::Job

=cut

sub fetch_input {
  my($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);
  
  my $slice = $self->fetch_sequence;
  $self->query($slice);

  my %blast = %{$self->BLAST_PARAMS};
  my $parser = $self->make_parser;
  my $filter;
  my %store_genes;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  my $ga = $self->get_dbadaptor($self->MODEL_DB)->get_GeneAdaptor;
  my $genes;
  if ( $self->LOGICNAME ) {
    $genes = $ga->fetch_all_by_Slice($self->query,$self->LOGICNAME,1);
  } else {
    $genes =  $ga->fetch_all_by_Slice($self->query,undef,1);
  }
  print "Found " .  scalar(@$genes) . " genes ";
  if (  $self->LOGICNAME ) {
    print " of type " . $self->LOGICNAME . "\n";
  } else {
    print "\n";
  }
  foreach my $gene (@$genes) {
    fully_load_Gene($gene);
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      $store_genes{$tran->dbID} = $gene;
      if ($tran->translation) {
	foreach my $db (split ',', ($self->analysis->db_file)) {
	  $self->runnable(Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
			  ->new(
				-transcript     => $tran,
				-query          => $self->query,
				-program        => $self->analysis->program_file,
				-parser         => $parser,
				-filter         => $filter,
				-database       => $db,
				-analysis       => $self->analysis,
				%blast,
			       ));
	}
      }
    }
  }
  $self->genes_by_tran_id(\%store_genes);
  return 1;
}

sub run {
  my ($self) = @_;
  my @runnables = @{$self->runnable};
  foreach my $runnable(@runnables){
    eval{
      $runnable->run;
    };
    if(my $err = $@){
      chomp $err;

      # only match '"ABC_DEFGH"' and not all possible throws
      if ($err =~ /^\"([A-Z_]{1,40})\"$/i) {
        my $code = $1;
        # treat VOID errors in a special way; they are really just
        # BLASTs way of saying "won't bother searching because
        # won't find anything"

        if ($code ne 'VOID') {
          $self->failing_job_status($1);          
          throw("Blast::run failed $@");
        }
      } elsif ($err) {
        throw("Blast Runnable returned unrecognised error string: $err");
      }
    }
    $self->add_supporting_features($runnable->transcript,$runnable->output);
  }
  1;
}

sub add_supporting_features {
  my ($self,$transcript,$features) = @_;
  my %gene_store = %{$self->genes_by_tran_id};
  my @output;
  my %feature_hash;
  my %name_hash;
  my @best;
  my $gene = $gene_store{$transcript->dbID};
   $self->throw("Gene not found using id " .$transcript->dbID ."\n")
    unless $gene;
  print "Got gene " . $gene->stable_id ."\n";
  foreach my $tran ( @{$gene->get_all_Transcripts} ) {
    my $sa = $self->get_dbadaptor($self->OUTPUT_DB)->get_SliceAdaptor;
    my $chr_slice = $sa->fetch_by_region('toplevel',$gene->seq_region_name);
    $tran = $tran->transfer($chr_slice);
    # order them by name
    foreach my $f ( @$features ) {
      $f = $f->transfer($chr_slice);
      push @{$name_hash{$f->hseqname}},$f;
    }
    # get the lowest eval and highest score for each protein
    foreach my $name (  keys %name_hash ) {
      my $protein = $name_hash{$name};
      my @p_val = sort { $a->p_value <=> $b->p_value } @{$protein};
      my @pid   = sort { $b->percent_id <=> $a->percent_id } @{$protein};
      my @score = sort { $b->score   <=> $a->score   } @{$protein};
      print STDERR "$name ";
      print STDERR " highest scoring HSP " . $p_val[0]->p_value . " " . $score[0]->score ."\n";
      # store all this data protein by highest score - lowest eval
      push @{$feature_hash{$p_val[0]->p_value}{$score[0]->score}{$name}},@{$protein};
    }
    
    # sort the alignments first by p_val then by score
    my @sorted_alignments;
    print "Sorted hits\n";
    foreach my $key ( sort { $a <=> $b }  keys %feature_hash ) {
      foreach my $score ( sort { $b <=> $a } keys %{$feature_hash{$key}} ) {
	foreach my $name (keys %{$feature_hash{$key}{$score}} ) {
	  print "$name $score $key\n"; 
	  push @sorted_alignments, $feature_hash{$key}{$score}{$name};
	}
      }
    }
    
  ALIGN:  while ( scalar(@sorted_alignments) > 0 ) {
      my @tmp;
      my @hsps = @{shift(@sorted_alignments)};
      # make the HSPS non overlapping unless they are full length
      my @contigs;
    HSP:   foreach my $hsp ( sort  { $b->score   <=> $a->score   } @hsps ) {
        print "HSP " . $hsp->hseqname ." "  . $hsp->hstart . " " . $hsp->hend ."\n";
	# does it overlap with another HSP
	foreach my $contig ( @contigs ) {
	  if ( $hsp->start <= $contig->end && 
	       $hsp->end >= $contig->start ) {
	    # reject it
	    print "REJECT\n";
	    next HSP;
	  }
	  if ( $hsp->hstart <= $contig->hend && 
	       $hsp->hend >= $contig->hstart ) {
	    # reject it
	    print "REJECT\n";
	    next HSP;
	  }
	}
	# not overlapping keep it
	push @contigs, $hsp;
      }

      @best =   sort { $a->hstart <=> $b->hstart } @contigs ;
      # we need to know they are all in order
      if ( $self->check_order(\@best) ) {
	# order got messed up, try the next best hit
	# empty the array
	@best = [];
	pop(@best);
      } else {
	last;
      }
    }

    if ( scalar(@best > 0 )  ) {
      # we want to make a single transcript supporting feature
      # using the highest scoring alignment and computing the peptide coverage
      print STDERR "Best hit\n";
      my $hlen= 0;
      my @filtered_best;      
      # set the p_val score and %ID to be the same for all hsps or
      # it will throw when you try to make a protien_align_feature
      my $pval = 10000000 ;
      my $pid = 0;
      my $score = 0;
      foreach my $f (@best ) {     
        $score = $f->score if $f->score > $score;
        $pid   = $f->percent_id if $f->percent_id > $pid;	
        $pval   = $f->p_value if $f->p_value < $pval;	
      } 
      foreach my $gf (@best ) {
	#need to break them back into ungapped features to get them to work consistantly
	foreach my $f ($gf->ungapped_features) {
	  print  $f->seq_region_name ." " .
	    $f->start ." " .
	      $f->end ." " .
		$f->strand ." " .
		  $f->hstart ." " .
		    $f->hend ." " .
		      $f->score ." " .
			$f->p_value ." " .
			  $f->percent_id . " " .
			      $f->hseqname ."\n";	
	  my $flen = $f->hend - $f->hstart + 1 ;
	  $hlen +=  $flen ;
	  $f->score($score);
	  $f->percent_id($pid);
	  $f->p_value($pval);
      	  push @filtered_best,$f;
	}
      }
      # make a transcript supporting feature
      my $coverage = sprintf("%.3f",( $hlen / length($tran->translation->seq) ) * 100);
      
      push my @tsf , Bio::EnsEMBL::DnaPepAlignFeature->new
	(-features => \@filtered_best);
      # hijack percent_id - not going to use it 
      $tsf[0]->percent_id(  $coverage ) ;
      $tsf[0]->analysis($tran->analysis);
      print STDERR "coverage $coverage\n";
      # calculate the hcoverage
      my $seqfetcher =Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher->new
	( -db => [$self->INDEX],
	  -format => 'fasta'
	);
      my $seq = $seqfetcher->get_Seq_by_acc($tsf[0]->hseqname);
      my $hcoverage =  sprintf("%.3f",( $hlen / length($seq->seq) ) * 100);
      $tsf[0]->hcoverage(  $hcoverage ) ;
      print STDERR "hcoverage $hcoverage\n";
      $tran->add_supporting_features(@tsf);
    }
  }

  push @output,$gene;
  $self->output(\@output);
}

sub write_output{
  my ($self) = @_;
  
  my $outdb = $self->get_dbadaptor($self->OUTPUT_DB);
  my $gene_adaptor = $outdb->get_GeneAdaptor; 
  
  my @output = @{$self->output};
  
  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output){
    #$gene->analysis($self->analysis);
    empty_Gene($gene);
    eval {
      $gene_adaptor->store($gene);
    };    
    if ($@){
      $self->warning("Unable to store gene!!\n");
      print STDERR "$@\n";
      $fails++;
    }
    $total++;
  }
  if ($fails > 0) {
    throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }
}

sub check_order {
  my ( $self, $hsps ) = @_;
  # should be the same order if sorted by hstart or start dependant on strand
  my $hstring;
  my $tstring;
  my @hsps =   sort { $a->hstart <=> $b->hstart } @$hsps ;
  foreach my $h ( @hsps ) {
    $hstring .= $h->start.":";
  }
  if ( $hsps[0]->strand == 1 ) {
    @hsps =   sort { $a->start <=> $b->start } @$hsps ;
  } else {
    @hsps =   sort { $b->start <=> $a->start } @$hsps ;
  }
  foreach my $h ( @hsps ) {
    $tstring .= $h->start.":";
  }
  print STDERR "\nCheck all hsps are contiguous\n$hstring\n$tstring\n\n";
  if ( $hstring eq $tstring ) {
    return;
  } else {
    print STDERR "HSP ORDER DOESNT MATCH EXON ORDER GETTING NEXT BEST HIT $hstring\n$tstring\n";
    return 1;
  }
}


sub genes_by_tran_id {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_gbtid'} = $value;
  }
  
  if (exists($self->{'_gbtid'})) {
    return $self->{'_gbtid'};
  } else {
    return undef;
  }
}



sub OUTPUT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUTPUT_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUTPUT_DB'})) {
    return $self->{'_CONFIG_OUTPUT_DB'};
  } else {
    return undef;
  }
}


sub MODEL_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MODEL_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MODEL_DB'})) {
    return $self->{'_CONFIG_MODEL_DB'};
  } else {
    return undef;
  }
}


sub LOGICNAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_LOGICNAME'} = $value;
  }
  
  if (exists($self->{'_CONFIG_LOGICNAME'})) {
    return $self->{'_CONFIG_LOGICNAME'};
  } else {
    return undef;
  }
}


sub INDEX {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INDEX'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INDEX'})) {
    return $self->{'_CONFIG_INDEX'};
  } else {
    return undef;
  }
}

1;
