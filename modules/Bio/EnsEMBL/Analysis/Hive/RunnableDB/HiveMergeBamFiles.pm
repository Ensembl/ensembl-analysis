=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::MergeBamFiles - Merge BAM files

=head1 SYNOPSIS

{
  -logic_name => 'merged_tissue_file',
  -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
  -parameters => {
    java       => 'java',
    java_options  => '-Xmx2g',
    # If 0, do not use multithreading, faster but can use more memory.
    # If > 0, tells how many cpu to use for samtools or just to use multiple cpus for picard
    use_threading => $self->o('use_threads'),
    # Path to MergeSamFiles.jar
    picard_lib    => $self->o('picard_lib_jar'),
    # Use this default options for Picard: 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'
    # You will need to change the options if you want to use samtools for merging
    options       => 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT',
    # target_db is the database where we will write the files in the data_file table
    # You can use store_datafile => 0, if you don't want to store the output file
    target_db => $self->o('blast_db'),
    disconnect_jobs => 1,
  },
  -rc_name    => '3GB_multithread',
  -flow_into => {
    1 => [ ':////accu?filename=[]' ],
  },
},

=head1 DESCRIPTION

Merge BAM files using either picard or samtools. Picard is the prefered choice.
It uses samtools to check if it was successfull. All the input files are retrieved
from the accu table of the Hive database and it stores the output file name in
the accu table on branch 1.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles;

use warnings ;
use strict;

use File::Copy;
use File::Spec;
use File::Basename qw(basename);

use Bio::EnsEMBL::DataFile;
use Bio::EnsEMBL::Analysis::Runnable::Samtools;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default values for this module are:
               bam_version => 1,
               rename_file => 0,
               sample_name => 'merged',
               store_datafile => 1,
               _index_ext => 'bai',
               _file_ext => 'bam',
               _logic_name_ext => 'bam',
               samtools => 'samtools',
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    bam_version => 1,
    rename_file => 0,
    sample_name => 'merged',
    store_datafile => 1,
    _index_ext => 'bai',
    _file_ext => 'bam',
    _file_type => 'BAMCOV',
    _logic_name_ext => 'bam',
    samtools => 'samtools',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: It will fetch all the BAM files. It will create a logic_name based on 'species',
              'sample_name' and '_logic_name_ext'. It will create a filename based on 'assembly_name',
              'rnaseq_data_provider', 'sample_name', 'bam_version' and '_file_ext'.
              If there is only one file, the file can either be renamed by settting 'rename_file' to 1
              or left as it is. The index will be generated if not present and the corect data will be
              inserted in the data_file table. The filename would be send to the '_branch_to_flow_to'
              as 'alignment_bam_file'. 'bam_file' will be set for the accu table to work properly.
              If 'picard_lib' is defined it will use Picard to merge the files otherwise it will
              use samtools
 Returntype : None
 Exceptions : Throws if 'filename' is empty
              Throws if the filename is not an absolute path
              Throws if it cannot rename files when there is only one file and rename_file is set

=cut

sub fetch_input {
    my ($self) = @_;

    my $input_files = $self->param('filename');
    if (@$input_files) {
      $self->create_analysis;
      if (!$self->param_is_defined('logic_name')) {
        $self->analysis->logic_name(join('_', $self->param('species'), $self->param('sample_name'), 'rnaseq', $self->param('_logic_name_ext')));
      }
      my $alignment_bam_file = File::Spec->catfile($self->param_required('output_dir'),
        join('.', $self->param_required('assembly_name'), $self->param_required('rnaseq_data_provider'), $self->param('sample_name'), $self->param('bam_version'), $self->param('_file_ext')));
      # We need 'bam_file' to be set for the accu table to work
      $self->param('bam_file', $alignment_bam_file);

      # We are creating the connection to make sure that the details are correct
      if ($self->param('store_datafile')) {
        my $db = $self->get_database_by_name('target_db');
        $db->dbc->disconnect_when_inactive(1) if ($self->param('disconnect_jobs'));
        $self->hrdb_set_con($db, 'target_db');
      }
      if (@$input_files == 1) {
        # In other cases I just want to push the filename but I don't need to run the BAM merge
        # First pushing the filename
        my $abs_filename = $input_files->[0];
        if (-e File::Spec->catfile($self->param('input_dir'), $abs_filename)) {
          $abs_filename = File::Spec->catfile($self->param('input_dir'), $abs_filename);
        }
        $self->throw($abs_filename.' is not an absolute path!') unless (File::Spec->file_name_is_absolute($abs_filename));
        if ($self->param('rename_file')) {
          my $index_ext = '.'.$self->param('_index_ext');
          if (!move($abs_filename, $alignment_bam_file)) {
            $self->throw("Could not rename file $abs_filename to $alignment_bam_file  Error is $!");
          }
          if (-e $abs_filename.$index_ext) {
            if (!move($abs_filename.$index_ext, $alignment_bam_file.$index_ext)) {
              $self->throw("Could not rename file $abs_filename$index_ext to $alignment_bam_file$index_ext  Error is $!");
            }
          }
        }
        else {
          # We need 'bam_file' to be set for the accu table to work
          $self->param('bam_file', $abs_filename);
          $alignment_bam_file = $abs_filename;
        }
        if (!-e $alignment_bam_file.'.'.$self->param('_index_ext')) {
          my $samtools = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(
              -program => $self->param('samtools'),
              -use_threading => $self->param('use_threading')
              );

          $samtools->index($alignment_bam_file);
        }
        # save file to datafile table. 
        if ($self->param('store_datafile')) {
          $self->store_filename_into_datafile;
        }        
        # Finally tell Hive that we've finished processing
        $self->dataflow_output_id({filename => $alignment_bam_file}, $self->param('_branch_to_flow_to'));
        $self->complete_early('There is only one file to process');
      }
      if ($self->param_is_defined('picard_lib')) {
          $self->require_module('Bio::EnsEMBL::Analysis::Runnable::PicardMerge');
          $self->runnable(Bio::EnsEMBL::Analysis::Runnable::PicardMerge->new(
              -program => $self->param('java') || 'java',
              -java_options => $self->param('java_options'),
              -lib => $self->param('picard_lib'),
              -options => $self->param('options'),
              -analysis => $self->analysis,
              -output_file => $alignment_bam_file,
              -input_files => $input_files,
              -use_threading => $self->param('use_threading'),
              -samtools => $self->param('samtools'),
              ));
      }
      else {
          $self->require_module('Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge');
          $self->runnable(Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge->new(
              -program => $self->param('samtools'),
              -options => $self->param('options'),
              -analysis => $self->analysis,
              -output_file => $alignment_bam_file,
              -input_files => $input_files,
              -use_threading => $self->param('use_threading'),
              ));
      }
    }
    else {
      $self->throw('You did not have input files.');
    }

}


=head2 run

 Arg [1]    : None
 Description: Merge the files and check that the merge was successful
 Returntype : None
 Exceptions : None

=cut

sub run {
    my ($self) = @_;

    $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
    foreach my $runnable (@{$self->runnable}) {
        $runnable->run;
        $runnable->check_output_file;
        $self->output([$runnable->output_file]);
    }
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the merged BAM file on branch '_branch_to_flow_to' as 'filename'
              In the pipeline it expects 'bam_file' to be set so the accu table can be used
 Returntype : None
 Exceptions : None

=cut

sub write_output {
    my ($self) = @_;

    if ($self->param('store_datafile')) {
      $self->store_filename_into_datafile;
    }
    $self->dataflow_output_id({filename => $self->output->[0]}, $self->param('_branch_to_flow_to')); 
}


=head2 store_filename_into_datafile

 Arg [1]    : None
 Description: It checks if the filename is already present. When absent, it uses
              the filename from 'bam_file' to generate the correct name and it
              adds the logic_name created in fetch_input if the logic_name is not
              present
 Returntype : None
 Exceptions : None

=cut

sub store_filename_into_datafile {
  my ($self) = @_;

  my $db = $self->hrdb_get_con('target_db');
  $db->dbc->disconnect_when_inactive(0);
  my $datafile_adaptor = $db->get_DataFileAdaptor;
  my $datafile_name = basename($self->param('bam_file'), '.'.$self->param('_file_ext'));
  if (!$datafile_adaptor->fetch_by_name_and_type($datafile_name, $self->param('_file_type'))) {
    my $analysis_adaptor = $db->get_AnalysisAdaptor;
    my $analysis = $analysis_adaptor->fetch_by_logic_name($self->analysis->logic_name);
    if (!$analysis) {
      $analysis_adaptor->store($self->analysis);
      $analysis = $self->analysis;
    }
    my @coord_systems = sort {$a->rank <=> $b->rank} @{$db->get_CoordSystemAdaptor->fetch_all_by_version($self->param('assembly_name'))};
    my $datafile = Bio::EnsEMBL::DataFile->new();
    $datafile->analysis($analysis);
    $datafile->name($datafile_name);
    $datafile->file_type($self->param('_file_type'));
    $datafile->version_lock(0);
    $datafile->absolute(0);
    $datafile->coord_system($coord_systems[0]);
    $db->get_DataFileAdaptor->store($datafile);
  }
}

1;
