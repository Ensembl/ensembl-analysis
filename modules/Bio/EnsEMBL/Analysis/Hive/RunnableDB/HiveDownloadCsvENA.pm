=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA;

use strict;
use warnings;
use feature 'say';

use JSON::PP;
use LWP::UserAgent;
use File::Spec::Functions qw(splitpath);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters for the analysis
               ena_base_url => 'http://www.ebi.ac.uk/ena/data/warehouse/search?display=report',
               files_domain => 'domain=read&result=read_run',
               files_fields => 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title',
               sample_domain => 'domain=sample&result=sample',
               sample_fields => 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,submitted_sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
               download_method => 'ftp',
               separator => '\t',
               _read_length => 1, # This is a default that should not exist. Some data do not have the read count and base count
               _centre_name => 'ENA',
               print_all_info => 0, # It will print all the sample information in the description instead of a selection, good to use when checking the CSV file
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ena_base_url => 'http://www.ebi.ac.uk/ena/data/warehouse/search?display=report',
    files_domain => 'domain=read&result=read_run',
    files_fields => 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title',
    sample_domain => 'domain=sample&result=sample',
    sample_fields => 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,submitted_sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
    download_method => 'ftp',
    separator => '\t',
    _read_length => 1,
    _centre_name => 'ENA',
    print_all_info => 0,
    paired_end_only => 1, #by default, module will only add paired-end data to the csv, add "paired_end_only => 0" to pipeline config to include single end data
    read_type => 'short_read',
    instrument_platform => 'ILLUMINA',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Create the query to be based on 'study_accession' or 'taxon_id'.
              An array of values can be given.
              If the 'inputfile' already exists, it will complete early
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  if($self->param_required('read_type') eq 'isoseq') {
    $self->param('paired_end_only',0);
    $self->param('instrument_platform','PACBIO_SMRT'),
  }
  $self->param_required('inputfile');
  if (-e $self->param('inputfile')) {
    $self->complete_early("'inputfile' exists so I will use that");
  } elsif ($self->param_is_defined('study_accession') and $self->param('study_accession')) {
    $self->_populate_query($self->param('study_accession'), 'study_accession=%s');
  } elsif ($self->param_is_defined('taxon_id') and $self->param('taxon_id')) {
    $self->_populate_query($self->param('taxon_id'), 'tax_tree(%s) AND instrument_platform='.$self->param('instrument_platform').' AND library_source=TRANSCRIPTOMIC');
  } else {
    $self->throw('"inputfile" does not exist and neither "study_accession" nor "taxon_id" were defined');
  }
}


=head2 _populate_query

 Arg [1]    : String or Arrayref of Strings, it should be a study accession or a taxon id or an array of them
 Arg [2]    : String $format, the format string for the param
 Description: It populates the 'query' parameters with an arrayref of string base on the format Arg[2] and the
              parameter(s) in Arg[1]. You cannot mix taxon_ids and study_accessions.
 Returntype : None
 Exceptions : None

=cut

sub _populate_query {
  my ($self, $params, $format) = @_;

  if (ref($params) eq 'ARRAY') {
    my @queries;
    foreach my $param (@$params) {
      push(@queries, sprintf($format, $param));
    }
    $self->param('query', \@queries);
  }
  else {
    $self->param('query', [sprintf($format, $params)]);
  }
}


=head2 run

 Arg [1]    : None
 Description: Query ENA to find all possibles project that could be used for the RNASeq pipeline
              It will avoid samples which are for non coding RNA analyses or if the individual was
              infected, immunised,...
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my %csv_data;
  my %samples;
  foreach my $query (@{$self->param('query')}) {
    my $ua = LWP::UserAgent->new;
    $ua->env_proxy;
    my $url = join('&', $self->param('ena_base_url'), 'query="'.$query.'"', $self->param('files_domain'), 'fields='.$self->param('files_fields'));
    $self->say_with_header($url);
    my $response = $ua->get($url);
    my $fastq_file = 'fastq_'.$self->param('download_method');
    if ($response->is_success) {
      my $content = $response->decoded_content(ref => 1);
      if ($content) {
        my %fields_index;
        while ($$content =~ /^(\w+.*)$/mgc) {
          my $line = $1;
          if ($line =~ /^[a-z]/) {
            my $index = 0;
            %fields_index = map { $_ => $index++} split('\t', $line);
          }
          else {
            # if these two checks below are removed, more time might be needed to prepare the CSV file
            next if ($line =~ / infected | [iIu]mmune| challenge |tomi[zs]ed/); # I do not want to do that but I don't think we have a choice
            next if ($line =~ /[Mm]i\w{0,3}RNA|lncRNA|circRNA|small RNA/); # I do not want to do that but I don't think we have a choice

            my @row = split("\t", $line);
            my $read_length = $self->param('_read_length');
            my $nominal_length = 0;
            my $calculated_length = 0;

            if ($self->param('paired_end_only')){
	      my $third_file = "ftp[^_]*\.fastq\.gz\;ftp";
	      $row[$fields_index{$fastq_file}] =~ s/$third_file/ftp/; # sometimes a third combined fastq file exists (it has single end naming format) - discard it

              next if ($row[$fields_index{library_layout}] eq 'SINGLE'); # don't include single end reads
	      next if ($row[$fields_index{$fastq_file}] !~ m/.*_1\.fastq\.gz.*_2\.fastq\.gz/);# this will throw out the paired end data that is stored in a single file, i.e. it looks like single end data to a regex
            }

            if ($row[$fields_index{nominal_length}]) {
              $nominal_length = $row[$fields_index{nominal_length}];
              $read_length = $nominal_length;
            }
            if ($row[$fields_index{base_count}] and $row[$fields_index{read_count}]) {
              $calculated_length = $row[$fields_index{base_count}]/$row[$fields_index{read_count}];
              $read_length = $calculated_length;
            }
            if ($nominal_length != $calculated_length and $nominal_length != 0 and $calculated_length != 0) {
              $self->warning("${row[$fields_index{run_accession}]} NOMINAL $nominal_length CALC $calculated_length");
            }
            if ($row[$fields_index{library_layout}] eq 'PAIRED') {
              $read_length /= 2;
            }
            next if ($read_length < 75);
            my %line = (
              run_accession => $row[$fields_index{run_accession}],
              instrument_model => $row[$fields_index{instrument_model}],
              instrument_platform => $row[$fields_index{instrument_platform}],
              library_layout => $row[$fields_index{library_layout}],
              fastq_file => $row[$fields_index{$fastq_file}],
              fastq_md5 => $row[$fields_index{fastq_md5}],
              study_title => $row[$fields_index{study_title}],
              experiment_title => $row[$fields_index{experiment_title}],
              instrument_model => $row[$fields_index{instrument_model}],
              read_length => $read_length,
              center_name => $row[$fields_index{center_name}],
            );
            $samples{$row[$fields_index{sample_accession}]} = $row[$fields_index{secondary_sample_accession}];
            push(@{$csv_data{$row[$fields_index{study_accession}]}->{$row[$fields_index{sample_accession}]}}, \%line);
          }
        }
        my @sample_names = keys %samples;
        my $header;
        SAMPLE: foreach my $sample (@sample_names) {
          $url = join('&', $self->param('ena_base_url'), 'query="accession='.$sample.'"', $self->param('sample_domain'), 'fields='.$self->param('sample_fields'));
          $response = $ua->get($url);
          if ($response->is_success) {
            $content = $response->decoded_content();
            if ($content) {
              while ($content =~ /^(\w+.*)$/mgc) {
                my $line = $1;
                if ($line =~ /^[a-z]/) {
                  my $index = 0;
                  %fields_index = map { $_ => $index++} split('\t', $line);
                  $header = $line;
                }
                else {
                  $self->say_with_header($line);
                  my @row = split("\t", $line);
                  my %line = (
                    center_name => $row[$fields_index{center_name}],
                    cell_line => $row[$fields_index{cell_line}],
                    cell_type => $row[$fields_index{cell_type}],
                    dev_stage => $row[$fields_index{dev_stage}],
                    sex => $row[$fields_index{sex}],
                    strain => $row[$fields_index{strain}],
                    sub_species => $row[$fields_index{sub_species}],
                    sub_strain => $row[$fields_index{sub_strain}],
                    tissue_lib => $row[$fields_index{tissue_lib}],
                    tissue_type => $row[$fields_index{tissue_type}],
                    variety => $row[$fields_index{variety}],
                    tax_id => $row[$fields_index{tax_id}],
                    description => $row[$fields_index{description}],
                    sample_collection => $row[$fields_index{sample_collection}],
                    sequencing_method => $row[$fields_index{sequencing_method}],
                    sample_alias => $row[$fields_index{sample_alias}],
                  );
                  my $dh = $ua->default_headers;
                  $ua->default_header('Content-Type' => 'application/json');
                  my $biosd = $ua->get('http://www.ebi.ac.uk/biosamples/samples/'.$sample);
                  if ($biosd->is_success) {
                    $content = $biosd->decoded_content();
                    my $json = JSON::PP->new();
                    my $data = $json->decode($content);
                    if (exists $data->{characteristics}->{immunization}) {
                      delete $samples{$sample};
                      $self->warning("Removed $sample from the set as it has immunization value: ".$data->{characteristics}->{immunization}->[0]->{text});
                      next SAMPLE;
                    }
                    $line{dev_stage} = $data->{characteristics}->{'development stage'}->[0]->{text}
                      if (exists $data->{characteristics}->{'development stage'});
                    $line{status} = $data->{characteristics}->{healthStatusAtCollection}->[0]->{text}
                      if (exists $data->{characteristics}->{healthStatusAtCollection});
                    if (exists $data->{characteristics}->{age}) {
                      $line{age} = $data->{characteristics}->{age}->[0]->{text};
                      if (exists $data->{characteristics}->{age}->[0]->{unit}) {
                        $line{age} .= ' '.$data->{characteristics}->{age}->[0]->{unit};
                      }
                    }
                    if (exists $data->{characteristics}->{tissue}) {
                      $line{organismPart} = $data->{characteristics}->{tissue}->[0]->{text};
                      if (exists $data->{characteristics}->{tissue}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{tissue}->[0]->{ontologyTerms}->[-1];
                      }
                    }
                    elsif (exists $data->{characteristics}->{cellType}) {
                      $line{cellType} = $data->{characteristics}->{cellType}->[0]->{text};
                      if (exists $data->{characteristics}->{cellType}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{cellType}->[0]->{ontologyTerms}->[-1];
                      }
                    }
                  }
                  else {
                    $self->warning("Could not connect to BioSample with $sample");
                  }

                  $ua->default_headers($dh);
                  $samples{$sample} = \%line;
                }
              }
            }
          }
        }
      }
      else {
        $self->throw("There was a problem with '$url'");
      }
    }
    else {
      $self->warning('No results with "'.$url);
    }
  }
  if (keys %csv_data) {
    foreach my $project (keys %csv_data) {
      my %dev_stages;
      my %celltypes;
      foreach my $sample (keys %{$csv_data{$project}}) {
        unless($sample =~ /^SAM/) {
          next;
	}
        next unless (exists $samples{$sample});
        if (exists $samples{$sample}->{dev_stage} and $samples{$sample}->{dev_stage}) {
#        if (exists $samples{$sample}->{dev_stage}) {
          next if ($samples{$sample}->{dev_stage} eq 'sexually immature stage');
          $dev_stages{$samples{$sample}->{dev_stage}} = 1;
        }
      }
      if (scalar(keys(%dev_stages)) > 1) {
        foreach my $sample (keys %{$csv_data{$project}}) {
          unless($sample =~ /^SAM/) {
            next;
  	  }
          next unless (exists $samples{$sample});
          if (exists $samples{$sample}->{dev_stage} and $samples{$sample}->{dev_stage}) {
            $samples{$sample}->{sample_name} = $samples{$sample}->{dev_stage};
          }
          else {
            #instead of breaking pipeline, show warning where sample has no development stage
            #use unknown for such cases
            $self->warning('No dev stages for '.$sample.' "'.join('", "', keys %dev_stages).'"');
            $samples{$sample}->{sample_name} = 'unknown';
          }

          if ($samples{$sample}->{sex}) {
            $samples{$sample}->{sample_name} = $samples{$sample}->{sex}.'_'.$samples{$sample}->{sample_name};
          }
        }
      }
      else {
        foreach my $sample (keys %{$csv_data{$project}}) {
          unless($sample =~ /^SAM/) {
            next;
  	  }
          next unless (exists $samples{$sample});
          $samples{$sample}->{sample_name} = $samples{$sample}->{cellType} || $samples{$sample}->{organismPart} || $samples{$sample}->{sample_alias} || $samples{$sample}->{description};
          if (!$samples{$sample}->{sample_name}) {
            $self->throw('No sample name for '.$sample);
          }
          if ($samples{$sample}->{sex}) {
            $samples{$sample}->{sample_name} = $samples{$sample}->{sex}.'_'.$samples{$sample}->{sample_name};
          }
        }
      }
    }
    $self->output([\%csv_data, \%samples]);
  }
  else {
    $self->complete_early('Could not find any data for this job');
  }
}


=head2 write_output

 Arg [1]    : None
 Description: It writes a csv file named 'inputfile' which can be process by the RNA-seq
              pipeline. The header should be
              SM\tID\tis_paired\tfilename\tis_mate_1\tread_length\tis_stranded\tCN\tPL\tDS
              is_mate_1 will always be -1 as we get the data from ENA so the filename will
                be informative of the mates
              is_stranded is 0 as we don't have a correct way of getting this information yet
              DS is made of "study_accession, sample_accession, study_title, experiment_title
                and cell_type if the sample is a cell_type
              It will also return the list of file to download on channel '_branch_to_flow_to',
              usually #2
 Returntype : None
 Exceptions : Throws if it cannot open or close the file 'inputfile'

=cut

sub write_output {
  my ($self) = @_;

  open(FH, '>'.$self->param('inputfile')) || $self->throw('Could not open '.$self->param('inputfile'));
  my $data = $self->output;
  my $samples = $data->[1];
  my $download_method = $self->param('download_method');
  my @output_ids;
  foreach my $study_accession (keys %{$data->[0]}) {
    my $study = $data->[0]->{$study_accession};
    foreach my $sample (keys %{$study}) {
      unless($sample =~ /^SAM/) {
        next;
      }
      next unless (exists $samples->{$sample});
      foreach my $experiment (@{$study->{$sample}}) {
        my @files = split(';', $experiment->{fastq_file});
        my @checksums = split(';', $experiment->{fastq_md5});
        my $index = 0;
        foreach my $file (@files) {
          my (undef, undef, $filename) = splitpath($file);
          my $sample_name = $samples->{$sample}->{sample_name};
          $sample_name =~ s/\s+-\s+\w+:\w+$//;
          $sample_name =~ s/[[:space:][:punct:]]+/_/g;
          $sample_name =~ s/_{2,}/_/g;
          $sample_name =~ s/^_|_$//g;
          my $description = sprintf("%s, %s%s",
            $experiment->{study_title},
            $experiment->{experiment_title},
            $samples->{$sample}->{cell_type} ? ', '.$samples->{$sample}->{cell_type} : '', );
          if ($self->param('print_all_info')) {
            foreach my $field (values %{$samples->{$sample}}) {
              $description .= ';'.$field if ($field);
            }
          }
          $description =~ tr/:\t/ /;
          if($self->param('read_type') eq 'short_read') {
            print FH sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s, %s, %s\n",
              lc($sample_name),
              $experiment->{run_accession},
              $experiment->{library_layout} eq 'PAIRED' ? 1 : 0,
              $filename,
              -1,
              $experiment->{read_length},
              0,
              $experiment->{center_name} || $self->param('_centre_name'),
              $experiment->{instrument_platform},
              $study_accession,
              $sample,
              $description,
            );
	  } elsif($self->param('read_type') eq 'isoseq') {
            print FH sprintf("%s\t%s\t%s\n",
              lc($sample_name),
              $filename,
              $description,
            );
	  } else {
            $self->throw('Read type unknown: '.$self->param('read_type'));
	  }
          push(@output_ids, {url => $file, download_method => $download_method, checksum => $checksums[$index++]});
        }
      }
    }
  }
  close(FH) || $self->throw('Could not close '.$self->param('inputfile'));
  $self->dataflow_output_id(\@output_ids, $self->param('_branch_to_flow_to'));
}

1;
