# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::Runnable::RepeatMasker->
  new(
      -query => 'slice',
      -program => 'repeatmasker',
      -options => '-low'
      -analysis => $analysis,
     );
  $repeat_masker->run;
  my @repeats = @{$repeat_masker->output};

=head1 DESCRIPTION

This module is base class for our Runnables. Runnables are there to
provide modules which can run different analyses and then parse the
results into core api objects

This module provides some base functionatily

The constructor can take 9 different arguments. The analysis object is
obligatory and must be passed in. The next 3 arguments, query, program and
options are the most important as it is with these the Runnable knows what
to run and on what sequences with what command line options. The next 4
are directory paths which can be determined from the config file
Bio::EnsEMBL::Analysis::Config::General but arguments are placed here so
they can be overidden if desired

The other base functionality includes some container methods
an output method aswell as methods for finding files and executables
and writing sequence to fasta files

All Runnables are expected to have 2 methods, run and output
run is the method which should run the analysis and output is where
the results should be stored

Generic versions of these methods are provided here. The run
method expects the runnables program to fit the commandline model used
by run_analysis, program options queryfile > resultsfile or to implement
its own run_analysis. If the run method is used the child Runnable
must implement a parse_results method as each analysis general has its
own output format and as such it cant be genericized

The output method provided simple holds an array of results and can
be given an arrayref to push onto that array

For more details about the specification look at Runnable.spec
in the ensembl-doc cvs module

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::InternalStopFix;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(contains_internal_stops);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : Bio::EnsEMBL::Slice
  Arg [3]   : string, name/path of program
  Arg [4]   : string commandline options for the program
  Arg [5]   : string path to working dir
  Arg [6]   : string, path to bin dir
  Arg [7]   : string, path to libary dir
  Arg [8]   : string, path to data dir
  Arg [9]   : Bio::EnsEMBL::Analysis;
  Function  : create a new Bio::EnsEMBL::Analysis::Runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable
  Exceptions: throws if not passed an analysis object
  Example   : $runnable = Bio::EnsEMBL::Analysis::Runnable::RepeatMasker
  ->new
  (
   -query => $self->query,
   -program => $self->analysis->program_file,
   $self->parameters_hash,
  );

=cut


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($genes, $edited_biotype, $stop_codon_biotype) = rearrange([qw( GENES EDITED_BIOTYPE STOP_CODON_BIOTYPE ) ], @args);

    throw("Should pass an array of genesinstead of ".ref($genes)."\n") if (ref($genes) ne 'ARRAY');
    $self->genes($genes);
    $self->edited_biotype($edited_biotype) if ($edited_biotype);
    $self->stop_codon_biotype($stop_codon_biotype) if ($stop_codon_biotype);

    return $self;
}

sub run {
    my $self = shift;

    my @new_genes;
    my $genes = $self->genes;
    foreach my $gene (@$genes) {
        my $contains_stop_translations = 0 ;
        my @transcripts;
        foreach my $tr ( @{ $gene->get_all_Transcripts } ) {
            if (!$tr->translation or $tr->translation->length < 1) {
                warning('Skipping transcript '.$tr->dbID."because it has no translation.\n");
            }
            else {
                my $nb_stops = contains_internal_stops($tr);
                if ($nb_stops == 1) {
                    my $new_tr = replace_stops_with_introns($tr);
                    if (!$new_tr) {
                        warning('Replace stops returned 0, meaning int stop was beside an intron - no model written '.$tr->display_id."\n");
                    } else {
                        push(@transcripts, $new_tr);
                        $contains_stop_translations++;
                    }
                }
                elsif ($nb_stops > 1) {
                    warning('Transcript '.$tr->display_id.' from gene '.$gene->display_id."has $nb_stops stop codons. New biotype: ".$self->stop_codon_biotype."\n");
                    $contains_stop_translations++;
                }
            }
        }
        if ($contains_stop_translations){
            if (@transcripts) {
                my $new_gene = Bio::EnsEMBL::Gene->new();
                foreach my $transcript (@transcripts) {
                    $transcript->biotype($self->edited_biotype);
                    $new_gene->add_Transcript($transcript);
                    $new_gene->analysis($self->analysis);
                    $new_gene->biotype($self->edited_biotype);
                }
                push(@new_genes, $new_gene);
            }
            else {
              $gene->biotype($self->stop_codon_biotype);
              foreach my $t (@{$gene->get_all_Transcripts}) {
                $t->biotype($self->stop_codon_biotype);
              }
              push(@new_genes, $gene);
            }
        }
    }
    $self->output(\@new_genes);
}

sub edited_biotype {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_edited_biotype} = $val;
  }
  return $self->{_edited_biotype};
}

sub stop_codon_biotype {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_stop_codon_biotype} = $val;
  }
  return $self->{_stop_codon_biotype};
}

sub genes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_genes} = $val;
  }
  return $self->{_genes};
}

1;
