=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig

=head1 SYNOPSIS

{
  -logic_name => 'create_ccode_config',
  -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig',
  -parameters => {
    single_tissue => $self->o('single_tissue'),
    sample_column => $self->o('read_group_tag'),
    sample_id_column => $self->o('read_id_tag'),
    csvfile_table => $self->o('summary_csv_table'),
    input_db => $self->o('rough_db'),
    dna_db => $self->o('dna_db'),
    output_db => $self->o('refine_db'),
    # write the intron features into the OUTPUT_DB along with the models
    write_introns => 1,
    # maximum number of times to loop when building all possible paths through the transcript
    max_recursions => 10000000000000,
    # analysis logic_name for the dna_align_features to fetch from the INTRON_DB
    # If left blank all features will be fetched
    logicname => [],
    # logic name of the gene models to fetch
    model_ln  => '',
    # penalty for removing a retined intron
    retained_intron_penalty => 2,
    #Remove introns that overlap X introns
    filter_on_overlap => 0,
    # minimum size for an intron
    min_intron_size  => 30,
    max_intron_size  => $self->o('maxintron'),
    # biotype to give to single exon models if left blank single exons are ignored
    # minimum single exon size (bp)
    min_single_exon => 1000,
    # minimum percentage of single exon length that is coding
    single_exon_cds => 66,
    # Intron with most support determines the splice sites for an internal exon
    # lower scoring introns with different splice sites are rejected
    strict_internal_splice_sites => 1,
    # In some species alternate splice sites for end exons seem to be common
    strict_internal_end_exon_splice_sites => 1,
    # biotypes to give gene models if left blank these models will not get written to the output database
    # best score - model with most supporting intron features
    # all other possible models
    # max number of other models to make - blank = all
    other_num      => '10',
    # max number of other models to process - blank = all
    max_num      => '1000',
    other_isoforms => $self->o('other_isoforms'),
    # biotype to label bad models ( otherwise they are not written )
    # do you want to trim UTR
    trim_utr => 1,
    # config for trimming UTR
    max_3prime_exons => 2,
    max_3prime_length => 5000,
    max_5prime_exons => 3,
    max_5prime_length => 1000,
    # % of average intron score that a UTR intron must have
    reject_intron_cutoff => 5,
  },
  -rc_name          => '1GB',
  -flow_into => {
    2 => ['create_ccode_input_ids'],
  },
},

=head1 DESCRIPTION

This module generate the configuration file needed by the RefineSolexaGenes
software which is available at https://github.com/Ensembl/ensc-core

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig;

use strict;
use warnings;

use File::Spec;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use constant DEFAULT_SEPARATOR => ';';
use constant DEFAULT_TAB => '  ';


=head2 param_defaults

 Arg [1]    : None
 Description: Returns the deafult parameters
               _intron_overlap_threshold => 4, An intron overlapping other introns
                must have more than _intron_overlap_threshold to be kept
               _merged_name => 'merged', Name for the merged set can be used to
                when you have only one tissue sample
               _gene_suffix => 'gene', Suffix to create the logic_name for the genes
               _intron_suffix => 'daf', Suffix to create the logic_name for the introns
               _ise_suffix => 'ise', Suffix to create the logic_name for the introns
                supporting evidences
 Returntype : Hash ref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _intron_overlap_threshold => 4, # This is a experimental value but it worked well
    _merged_name => 'merged',
    _gene_suffix => 'gene',
    _intron_suffix => 'daf',
    _ise_suffix => 'ise',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Prepare an hashref based on 'single_tissue' and other parameters like
              'sample_id_column', 'sample_id_column', 'other_isoforms', 'bad_models',
              'introns_ln_suffix', 'ise_logic_name' to generate the analysis specific
              parameters
              if 'target_db' is defined, check that the analysis exists in order to
              generate the config. This does NOT apply for the merge analysis
              if there is only one tissue sample, we do not create the merge analysis
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  my @output_ids;
  my $default_iot = $self->param('_intron_overlap_threshold');
  my $merged_iot = $default_iot;
  my $merged_name = $self->param('_merged_name');
  my $gene_suffix = $self->param('_gene_suffix');
  my $intron_suffix = $self->param('_intron_suffix');
  my $ise_suffix = $self->param('_ise_suffix');
  my $analysis_adaptor = $self->get_database_by_name('target_db')->get_AnalysisAdaptor;
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
        my $base_logic_name = $self->param('wide_species').'_'.$key.'_rnaseq';
        next unless ($analysis_adaptor->fetch_by_logic_name($base_logic_name."_$gene_suffix"));
        my %analysis_hash = (BEST_SCORE => "best_$key", SINGLE_EXON_MODEL => "single_$key", INTRON_OVERLAP_THRESHOLD => $default_iot);
        $analysis_hash{OTHER_ISOFORMS} = $self->param('other_isoforms').'_'.$key if ($self->param_is_defined('other_isoforms'));
        $analysis_hash{BAD_MODELS} = $self->param('bad_models').'_'.$key if ($self->param_is_defined('bad_models'));
        $analysis_hash{INTRONS_LOGIC_NAME} = $base_logic_name."_$intron_suffix" if ($intron_suffix);
        $analysis_hash{ISE_LOGIC_NAME} = $base_logic_name."_$ise_suffix" if ($ise_suffix);
        push(@output_ids, [
          File::Spec->catfile($self->param('wide_output_dir'), $self->param('wide_species').'_'.$key.'.conf'),
          [
            {FILE => $self->param('wide_intron_bam_file').'.bam',
            GROUPNAME => [keys %{$tissue_hash{$key}}],
            DEPTH => 0,
            MIXED_BAM => 0}
          ],
          $base_logic_name."_$gene_suffix",
          \%analysis_hash
        ]);
        $tissue_count++;
      }
# This feature is still experimental, maybe it can be lowered.
# Anyway people can specify the INTRON_OVERLAP_THRESHOLD for any samples
# in the config file
       $merged_iot = $tissue_count*$default_iot;
  }
  if (@output_ids > 1) {
    my $merged_logic_name = $self->param('wide_species').'_'.$merged_name.'_rnaseq';
    my %analysis_hash = (BEST_SCORE => 'best', SINGLE_EXON_MODEL => 'single', INTRON_OVERLAP_THRESHOLD => $merged_iot);
    $analysis_hash{OTHER_ISOFORMS} = $self->param('other_isoforms').'_'.$merged_name if ($self->param_is_defined('other_isoforms'));
    $analysis_hash{BAD_MODELS} = $self->param('bad_models').'_merged' if ($self->param_is_defined('bad_models'));
    $analysis_hash{INTRONS_LOGIC_NAME} = $merged_logic_name."_$intron_suffix" if ($intron_suffix);
    $analysis_hash{ISE_LOGIC_NAME} = $merged_logic_name."_$ise_suffix" if ($ise_suffix);
    push(@output_ids, [
      File::Spec->catfile($self->param('wide_output_dir'), $self->param('wide_species')."_$merged_name.conf"),
      [
        {FILE => $self->param('wide_intron_bam_file').'.bam',
        GROUPNAME => [],
        DEPTH => 0,
        MIXED_BAM => 0}
      ],
      $merged_logic_name."_$gene_suffix",
      \%analysis_hash
    ]);
  }
  elsif (@output_ids == 0) {
    $self->complete_early('Not enough tissue samples to work on');
    $self->input_job->autoflow(0);
  }
  $self->param('analyses', \@output_ids);
  $self->param('database_file', File::Spec->catfile($self->param('wide_output_dir'), $self->param('wide_species').'_database.conf'));
}


=head2 run

 Arg [1]    : None
 Description: Generate the configuration and store it in output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  $self->output([$self->_generate_databases_file, $self->_generate_config_file]);
}


=head2 write_output

 Arg [1]    : None
 Description: Write the analysis configuration file and the database configuration file
 Returntype : None
 Exceptions : Throws if it can not open or close a file to write

=cut

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


=head2 _generate_databases_file

 Arg [1]    : None
 Description: Generate the database configuration file. It will not set the driver
 Returntype : Arrayref
 Exceptions : None

=cut

sub _generate_databases_file {
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
    delete $config_file{DATABASES}->{$db}->{-driver} if (exists $config_file{DATABASES}->{$db}->{-driver});
  }

  return __print_hash('Config', \%config_file, undef, '');
}

=head2 _generate_config_file

 Arg [1]    : None
 Description: Generate the analysis configuration file with a default
              section and an analysis section
 Returntype : Arrayref
 Exceptions : None

=cut

sub _generate_config_file {
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
            INTRONS_LOGIC_NAME => $analysis->[2],
            ISE_LOGIC_NAME => $analysis->[2],
          },
          $analysis->[2] => $analysis->[3],
        }
    );
    foreach my $key (keys %{$config_file{REFINESOLEXAGENES_CONFIG_BY_LOGIC}->{DEFAULT}}) {
      $config_file{REFINESOLEXAGENES_CONFIG_BY_LOGIC}->{$analysis->[2]}->{$key} = $self->param($key)
        if ($self->param_is_defined($key));
    }
    push(@results, [$analysis->[0], __print_hash('Config', \%config_file, undef, '')]);
  }

  return \@results;
}


=head2 __print_hash

 Arg [1]    : Hashref
 Arg [2]    : String $separator
 Arg [3]    : String $tab
 Description: Generate the list of string that will represent the hash Arg[1] for the config
              and applying the separator Arg[2] and the tabulation Arg[3]
 Returntype : Arrayref
 Exceptions : None

=cut

sub __print_hash {
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
      push(@lines, @{__print_array($key, $value, $separator, $tab)});
    }
    elsif (ref($value) eq 'HASH') {
      push(@lines, @{__print_hash($key, $value, ',', DEFAULT_TAB.$tab)});
    }
    else {
      push(@lines, DEFAULT_TAB, $tab, $key , ' = ', $value =~ /^-?[0-9.]+$/ ? $value : "\"$value\"", DEFAULT_SEPARATOR, "\n");
    }
  }
  push(@lines, $tab, '}', "$separator\n");

  return \@lines;
}


=head2 __print_array

 Arg [1]    : Arrayref
 Arg [2]    : String $separator
 Arg [3]    : String $tab
 Description: Generate the list of string that will represent the array Arg[1] for the config
              and applying the separator Arg[2] and the tabulation Arg[3]
 Returntype : Arrayref
 Exceptions : None

=cut

sub __print_array {
  my ($name, $arrayref, $separator, $tab) = @_;

  $separator = DEFAULT_SEPARATOR unless ($separator);
  my @lines = (DEFAULT_TAB, $tab, $name, ' = ', '[ ');
  foreach my $value (@$arrayref) {
    if (ref($value) eq 'ARRAY') {
      push(@lines, @{__print_array($name, $value, $separator, DEFAULT_TAB.$tab)});
    }
    elsif (ref($value) eq 'HASH') {
      $lines[-1] = "(\n";
      push(@lines, @{__print_hash(undef, $value, $separator, DEFAULT_TAB.$tab)});
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
