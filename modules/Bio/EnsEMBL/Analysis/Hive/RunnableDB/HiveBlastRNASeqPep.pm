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
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BlastRNASeqPep

=head1 SYNOPSIS

my $db    = Bio::EnsEMBL::DBAdaptor->new($locator);
my $btpep = Bio::EnsEMBL::Analysis::RunnableDB::BlastRNASeqPep->new (
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;
use Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Analysis_to_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(align_proteins);
use Bio::EnsEMBL::IO::Parser::Fasta;
use Data::Dumper;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast');


=head2 fetch_input

 Arg [1]    : None
 Description: Fetches genes from the input database, either all of them or the one specified by a logic_name.
 Returntype : None
 Exceptions : Throw if input_id is not defined

=cut

sub fetch_input {
  my($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);

  $self->create_analysis(1);
  $self->analysis->parameters($self->param('commandline_params')) if ($self->param_is_defined('commandline_params'));


  my $input_dba = $self->hrdb_get_dba($self->param('input_db'));
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));

  $self->hrdb_set_con($input_dba,'input_db');
  $self->hrdb_set_con($output_dba,'output_db');

  my $input_sa = $input_dba->get_SliceAdaptor();
  my $input_slice = $input_sa->fetch_by_name($self->input_id);
  $self->query($input_slice);

  my $output_sa = $output_dba->get_SliceAdaptor;
  my $output_slice = $output_sa->fetch_by_region('toplevel',$input_slice->seq_region_name);
  unless($output_slice) {
    $self->throw("Could not fetch an output slice from the output db. Slice name: ".$input_slice->seq_region_name);
  } else {
    say "Fetched output slice: ".$output_slice->name();
  }
  $self->param('toplevel_slice',$output_slice);

  my $parser = $self->make_parser;
  my $filter;
  my %store_genes;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  my $logicname = $self->param('logic_name');

  my $genes = $input_slice->get_all_Genes($logicname);
  unless(scalar(@$genes)) {
    $self->warning("Found no input genes for ".$self->input_id." of type ".$self->LOGICNAME.
                   "\nAnalysis dbID: ".$self->analysis->dbID.
                   "\nSlice db name: ".$input_slice->adaptor->dbc->dbname);
    $self->complete_early('No genes to process');
  }

  print "Found ".scalar(@$genes)." genes";
  if ($logicname) {
    print " of type ".$logicname."\n";
  } else {
    print "\n";
  }

  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  $input_dba->dnadb($dna_dba);
  $output_dba->dnadb($dna_dba);
  $self->hrdb_set_con($dna_dba,'dna_db');

  if ($self->param_is_defined('create_translation') and $self->param('create_translation')) {
    foreach my $gene (@$genes) {
      foreach my $tran (@{$gene->get_all_Transcripts}) {
        compute_translation($tran);
      }
    }
  }

  foreach my $gene (@$genes) {
      foreach my $tran (@{$gene->get_all_Transcripts}) {
          $store_genes{$tran->dbID} = $gene;
          if ($tran->translation) {
              $tran->translation->seq;
#              foreach my $db (split ',', ($self->analysis->db_file)) {
              foreach my $db (@{$self->param('uniprot_index')}) {
                  $self->runnable(Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep->new(
                              -transcript     => $tran,
                              -query          => $self->query,
                              -program        => $self->param('blast_program'),
                              -parser         => $parser,
                              -filter         => $filter,
                              -database       => $db,
                              -analysis       => $self->analysis,
                              %{$self->BLAST_PARAMS},
                              ));
              }
          } else {
            $self->warning("No translation found for refine gene with dbID: ".$gene->dbID);
	  }
      }
  }
  $self->genes_by_tran_id(\%store_genes);
  return 1;
}


=head2 run

 Arg [1]    : None
 Description: Run blast against UniProt to assess the coding potential of the transcripts
 Returntype : None
 Exceptions : None

=cut

sub run {
    my ($self) = @_;
    if ($self->param('disconnect_jobs')) {
      $self->dbc->disconnect_if_idle;
      $self->hrdb_get_con('dna_db')->dbc->disconnect_if_idle;
    }
    foreach my $runnable (@{$self->runnable}){
        eval{
            $runnable->run;
        };
        if(my $err = $@){
            chomp $err;

            $self->warning($err);
            # only match '"ABC_DEFGH"' and not all possible throws
            if ($err =~ /^\"([A-Z_]{1,40})\"$/i) {
              my $code = $1;
              # treat VOID errors in a special way; they are really just
              # BLASTs way of saying "won't bother searching because
              # won't find anything"

              if ($code ne 'VOID') {
                $self->throw("Blast::run failed $@");
              }
              my $failed_transcript_id = $runnable->transcript()->dbID;
              my $failed_gene = $self->genes_by_tran_id->{$failed_transcript_id};
	      unless($failed_gene) {
                $self->throw("Could not fetch failed gene from the initial gene hash. Something has gone wrong. Hash transcript dbID checked: ".$failed_transcript_id);
	      }
              $self->failed_genes($failed_gene);
            }
            elsif ($err) {
                $self->throw("Blast Runnable returned unrecognised error string: $err");
            }
        }
        $self->add_supporting_features($runnable->transcript,$runnable->output);
#        if($self->param('calculate_coverage_and_pid')) {
          $self->realign_results($self->output);
#        }
    }
    return 1;
}


=head2 add_supporting_features

 Arg [1]    : Bio::EnsEMBL::Transcript
 Arg [2]    : Hashref
 Description: Add supporting evidence to a transcript
 Returntype : None
 Exceptions : Throws if it cannot find the associated gene

=cut

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
  my $chr_slice = $self->param('toplevel_slice');
  my $transcripts = $gene->get_all_Transcripts;
  $gene->flush_Transcripts;
  foreach my $tran ( @$transcripts ) {
    $tran = $tran->transfer($chr_slice);
    # order them by name
    foreach my $f ( @$features ) {
      push(@{$name_hash{$f->hseqname}}, $f->transfer($chr_slice));
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
	  print  join(' ', $f->seq_region_name, $f->start, $f->end, $f->strand, $f->hstart, $f->hend, $f->score, $f->p_value, $f->percent_id, $f->hseqname)."\n";
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

      my $tsf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@filtered_best, -align_type => 'ensembl');
      # hijack percent_id - not going to use it
      $tsf->percent_id(  $coverage ) ;
      $tsf->analysis($tran->analysis);
      print STDERR "coverage $coverage\n";
      # calculate the hcoverage
      my $seqfetcher = Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher->new(
        -db => [$self->INDEX],
        -format => 'fasta'
      );
      my $seq = $seqfetcher->get_Seq_by_acc($tsf->hseqname);
      my $hcoverage =  sprintf("%.3f",( $hlen / length($seq->seq) ) * 100);
      $tsf->hcoverage(  $hcoverage ) ;
      print STDERR "hcoverage $hcoverage\n";
      $tran->add_supporting_features($tsf);
    }
    $gene->add_Transcript($tran);
  }

  push @output,$gene;
  $self->output(\@output);
}


=head2 write_output

 Arg [1]    : None
 Description: Write the genes in the output database
 Returntype : None
 Exceptions : Throws if it fails to write at least one gene

=cut

sub write_output {
  my ($self) = @_;

  my $outdb = $self->hrdb_get_con('output_db');
  my $gene_adaptor = $outdb->get_GeneAdaptor;

  my $fails = 0;
  my $total = 0;
  my $analysis = $self->analysis;
  foreach my $gene (@{$self->output}){
    attach_Slice_to_Gene($gene,$self->param('toplevel_slice'));
    attach_Analysis_to_Gene($gene,$analysis);
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

  my $runnable_failed_genes = $self->failed_genes;
  foreach my $gene (@{$runnable_failed_genes}) {
    my $failed_biotype = $gene->biotype()."_failed";
    my $failed_source = $gene->source()."_failed";
    $gene->biotype($failed_biotype);
    $gene->source($failed_source);

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
    $self->throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }

  my $total_input_genes = scalar(keys(%{$self->genes_by_tran_id()}));
  my $success_percent = $total / $total_input_genes * 100;
  if($total_input_genes >= 10 && $success_percent <= 90) {
    $self->throw("The number of output genes differs too much from the number of input genes\n".
                 "Total input: ".$total_input_genes."\nTotal output: ".$total);
  }
}


=head2 check_order

 Arg [1]    : Arrayref of Bio::EnsEMBL::FeaturePair
 Description: Sort the hits by start and check that they don't overlap
 Returntype : Boolean, 1 if order is wrong
 Exceptions : None

=cut

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


=head2 genes_by_tran_id

 Arg [1]    : (optional) Hashref of Bio::EnsEMBL::Gene, the key is the transcript dbID
 Description: Getter/setter for a gene based on a transcript id
 Returntype : Hashref of Bio::EnsEMBL::Gene
 Exceptions : None

=cut

sub genes_by_tran_id {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->param('_gbtid', $value);
    }

    return $self->param('_gbtid');
}


=head2 OUTPUT_DB

 Arg [1]    : None
 Description: Getter for the name of the output database, default is 'output_db'
 Returntype : String 'output_db'
 Exceptions : None

=cut

sub OUTPUT_DB {
    my ($self) = @_;

    return 'output_db';
}


=head2 MODEL_DB

 Arg [1]    : None
 Description: Getter for the name of the output database, default is 'input_db'
 Returntype : String 'input_db'
 Exceptions : None

=cut

sub MODEL_DB {
    my ($self) = @_;

    return 'input_db';
}


=head2 LOGICNAME

 Arg [1]    : None
 Description: Getter for the logic name of the genes to fetch using param 'logic_name'
 Returntype : String
 Exceptions : None

=cut

sub LOGICNAME {
    my ($self) = @_;

    if ($self->param_is_defined('logic_name')) {
      return $self->param('logic_name');
    }
    else {
      return;
    }
}


=head2 INDEX

 Arg [1]    : None
 Description: Getter for the directory to the Indicate index used to fetch protein sequences
              using param 'indicate_index'
 Returntype : String
 Exceptions : None

=cut

sub INDEX {
    my ($self,$value) = @_;

    return $self->param('indicate_index');
}


sub failed_genes {
  my ($self,$value) = @_;

  unless($self->param('_failed_genes')) {
    $self->param('_failed_genes',[])
  }

  if($value) {
    push(@{$self->param('_failed_genes')},$value)
  }

  return $self->param('_failed_genes');
}

1;



sub realign_results {
  my ($self,$genes) = @_;

  foreach my $gene (@{$genes}) {
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      my $sfs = $transcript->get_all_supporting_features;
      unless(scalar(@$sfs)) {
        next;
      }

      my $sf = ${$sfs}[0];
      my $hit_name = $sf->hseqname;
      my $original_seq;

      my $uniprot_db = ${$self->param('uniprot_index')}[0];
      my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($uniprot_db);

      while($parser->next()) {
        my ($accession) = $parser->getHeader =~ /^(\S+)/;
        if($accession eq $hit_name) {
          $original_seq = $parser->getSequence;
          last;
        }
      }

      unless($original_seq) {
        $self->throw("Was not able to find the original sequence in the blast db for realignment.".
                     "\nOriginal seq: ".$hit_name."\nBlast db location: ".$uniprot_db);
      }

      my ($coverage,$percent_id) = align_proteins($original_seq,$transcript->translation->seq);
      foreach my $sf (@{$sfs}) {
        $sf->hcoverage($coverage);
        $sf->percent_id($percent_id);
      }
    }
  }
}
