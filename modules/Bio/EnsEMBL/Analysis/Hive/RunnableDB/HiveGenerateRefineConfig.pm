# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Bam2Genes

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::Bam2Genes->new (
          -db      => $db,
          -input_id   => $input_id
          -analysis   => $analysis );
  $refine_genes->fetch_input();
  $refine_genes->run();
  $refine_genes->write_output(); #writes to DB

=head1 DESCRIPTION

The module creates "proto-transcripts" based on the alignments of short reads.
It will first create blocks from overlapping reads which represent possible exons and
it will link these blocks by using pairing information if it is available or by
using a predefined "proto-transcript" length when it uses single reads.
The "proto-transcripts" will be stored in an Ensembl database.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig;

use strict;
use warnings;
use File::Spec;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use constant DEFAULT_SEPARATOR => ';';
use constant DEFAULT_TAB => '  ';

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _intron_overlap_threshold => 4, # This is a experimental value but it worked well
  }
}


sub fetch_input {
  my ($self) = @_;

    my @output_ids;
    my $other_isoforms = '';
    my $bad_models = '';
    my $default_iot = $self->param('_intron_overlap_threshold');
    my $merged_iot = $default_iot;
    if ($self->param('single_tissue')) {
        my $tissue_count = 0;
        my $table_adaptor = $self->db->get_NakedTableAdaptor;
        $table_adaptor->table_name($self->param('csvfile_table'));
        my %tissue_hash;
        my $results = $table_adaptor->fetch_all();
        foreach my $result (@$results) {
            $tissue_hash{lc($result->{$self->param('sample_column')})}->{$result->{$self->param('sample_id_column')}} = 1;
        }
        foreach my $key (keys %tissue_hash) {
          my %analysis_hash = (BEST_SCORE => "best_$key", SINGLE_EXON_MODEL => "single_$key", INTRON_OVERLAP_THRESHOLD => $default_iot);
          $analysis_hash{OTHER_ISOFORMS} = $self->param('other_isoforms').'_'.$key if ($self->param_is_defined('other_isoforms'));
          $analysis_hash{BAD_MODELS} = $self->param('bad_models').'_'.$key if ($self->param_is_defined('bad_models'));
          push(@output_ids, [File::Spec->catfile($self->param('wide_output_dir'), $self->param('wide_species').'_'.$key.'.conf'), [{FILE => $self->param('wide_intron_bam_file').'.bam', GROUPNAME => [keys %{$tissue_hash{$key}}], DEPTH => 0, MIXED_BAM => 0}], $self->param('wide_species').'_'.$key.'_rnaseq', \%analysis_hash]);
          $tissue_count++;
        }
# This feature is still experimental, maybe it can be lowered.
# Anyway people can specify the INTRON_OVERLAP_THRESHOLD for any samples
# in the config file
         $merged_iot = $tissue_count*$default_iot;
    }
    my %analysis_hash = (BEST_SCORE => 'best', SINGLE_EXON_MODEL => 'single', INTRON_OVERLAP_THRESHOLD => $merged_iot);
    $analysis_hash{OTHER_ISOFORMS} = $self->param('other_isoforms').'_merged' if ($self->param_is_defined('other_isoforms'));
    $analysis_hash{BAD_MODELS} = $self->param('bad_models').'_merged' if ($self->param_is_defined('bad_models'));
    push(@output_ids, [File::Spec->catfile($self->param('wide_output_dir'), $self->param('wide_species').'_merged.conf'), [{FILE => $self->param('wide_intron_bam_file').'.bam', GROUPNAME => [], DEPTH => 0, MIXED_BAM => 0}], $self->param('wide_species').'_merged_rnaseq', \%analysis_hash]);
  $self->param('analyses', \@output_ids);
  $self->param('database_file', File::Spec->catfile($self->param('wide_output_dir'), $self->param('wide_species').'_database.conf'));
}

sub run {
  my ($self) = @_;

  $self->output([$self->generate_databases_file, $self->generate_config_file]);
}

sub write_output {
  my ($self) = @_;

  my $database_file = $self->param('database_file');
  my $output = $self->output;
  open(DBH, ">$database_file") || $self->throw("Could not open $database_file for writing");
  print DBH @{$output->[0]};
  close(DBH) || $self->throw("Could not close $database_file");
  foreach my $config (@{$output->[1]}) {
    my $config_file = $config->[0];
    open(CH, ">$config_file") || $self->throw("Could not open $config_file for writing");
    print CH @{$config->[1]};
    close(CH) || $self->throw("Could not close $config_file");
  }
  my @output_ids;
  foreach my $analysis (@{$self->param('analyses')}) {
    push(@output_ids, {config_file => $analysis->[0], logic_name => $analysis->[2]});
  }
  $self->dataflow_output_id(\@output_ids, 2);
}

sub generate_databases_file {
  my ($self) = @_;

  my %config_file = (
    DATABASES => {
      REFERENCE_DB => $self->param('dna_db'),
      ROUGH_DB => $self->param('input_db'),
      REFINED_DB => $self->param('output_db'),
    },
    DNA_DBNAME => "REFERENCE_DB",
  );
  foreach my $db ('REFERENCE_DB', 'ROUGH_DB', 'REFINED_DB') {
    $config_file{DATABASES}->{$db}->{-pass} = "" unless (exists $config_file{DATABASES}->{$db}->{-pass});
  }

  return print_hash('Config', \%config_file, undef, '');
}

sub generate_config_file {
  my ($self) = @_;

  my @results;
  foreach my $analysis (@{$self->param('analyses')}) {
    my %config_file = (
        DATABASES_FILE => $self->param('database_file'),
        REFINESOLEXAGENES_CONFIG_BY_LOGIC => {
          DEFAULT => {
            OUTPUT_DB => "REFINED_DB",
            INTRON_DB => "",
            MODEL_DB  => "ROUGH_DB",
            INTRON_BAM_FILES => $analysis->[1],
            WRITE_INTRONS => 1,
            MAX_RECURSIONS => 100000,
            LOGICNAME => [],
            MODEL_LN  => "",
            RETAINED_INTRON_PENALTY => "2.0",
            MIN_INTRON_SIZE  => 30,
            MAX_INTRON_SIZE  => 200000,
            SINGLE_EXON_MODEL => "single_exon",
            MIN_SINGLE_EXON => 1000,
            SINGLE_EXON_CDS => "66.0",
            STRICT_INTERNAL_SPLICE_SITES => 1,
            STRICT_INTERNAL_END_EXON_SPLICE_SITES => 1,
            INTRON_OVERLAP_THRESHOLD => 4,
            BEST_SCORE => "best",
            OTHER_ISOFORMS => "",
            OTHER_NUM      => 10,
            MAX_NUM      => 1000,
            BAD_MODELS     => "",
            TRIM_UTR => 1,
            MAX_3PRIME_EXONS => 2,
            MAX_3PRIME_LENGTH => 5000,
            MAX_5PRIME_EXONS => 3,
            MAX_5PRIME_LENGTH => 1000,
            REJECT_INTRON_CUTOFF => "5.0",
            TYPE_PREFIX    => "ccode",
            CONSLIMS    => [ "1.0" ],
            NONCONSLIMS => [ "1.0" ],
            RESTART_NONCONSLIM => "-1.0",
          },
          $analysis->[2] => $analysis->[3],
        }
    );
    foreach my $key (keys %{$config_file{REFINESOLEXAGENES_CONFIG_BY_LOGIC}->{DEFAULT}}) {
      $config_file{REFINESOLEXAGENES_CONFIG_BY_LOGIC}->{$analysis->[2]}->{$key} = $self->param($key)
        if ($self->param_is_defined($key));
    }
    push(@results, [$analysis->[0], print_hash('Config', \%config_file, undef, '')]);
  }

  return \@results;
}

sub print_hash {
  my ($name, $hashref, $separator, $tab) = @_;

  $separator = DEFAULT_SEPARATOR unless ($separator);
  my @lines;
  if ($name) {
    push(@lines, $tab, $name, ' : {', "\n");
  }
  else {
    push(@lines, $tab, '{', "\n");
  }

  while( my ($key, $value) = each %$hashref) {
    $key =~ s/^-//;
    if (ref($value) eq 'ARRAY') {
      push(@lines, @{print_array($key, $value, $separator, $tab)});
    }
    elsif (ref($value) eq 'HASH') {
      push(@lines, @{print_hash($key, $value, ',', DEFAULT_TAB.$tab)});
    }
    else {
      push(@lines, DEFAULT_TAB, $tab, $key , ' = ', $value =~ /^-?[0-9.]+$/ ? $value : "\"$value\"", DEFAULT_SEPARATOR, "\n");
    }
  }
  push(@lines, $tab, '}', "$separator\n");

  return \@lines;
}

sub print_array {
  my ($name, $arrayref, $separator, $tab) = @_;

  $separator = DEFAULT_SEPARATOR unless ($separator);
  my @lines = (DEFAULT_TAB, $tab, $name, ' = ', '[ ');
  foreach my $value (@$arrayref) {
    if (ref($value) eq 'ARRAY') {
      push(@lines, @{print_array($name, $value, $separator, DEFAULT_TAB.$tab)});
    }
    elsif (ref($value) eq 'HASH') {
      $lines[-1] = "(\n";
      push(@lines, @{print_hash(undef, $value, $separator, DEFAULT_TAB.$tab)});
      $lines[-1] = "\n";
      push(@lines, $tab, ')'.DEFAULT_SEPARATOR);
    }
    else {
      push(@lines, $value =~ /^-?[0-9.]+$/ ? $value : "\"$value\"",  ',');
    }
  }
  if ($lines[-1] eq ',') {
    $lines[-1] = ' ]';
  }
  elsif ($lines[-1] eq '[ ') {
    $lines[-1] = '[]'.DEFAULT_SEPARATOR;
  }
  elsif ($lines[-1] ne ');') {
    push(@lines, ' ]');
  }
  push(@lines, "\n");
  return \@lines;
}

1;
