=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::MergBamFiles - 

=head1 SYNOPSIS

  [mergebamfiles]
  module=MergeBamFiles
  db_version=10 #You will have 10 bam files to merge, its for the checking
  db_file=/path/to/my/bam/files #Can be comma seperated
  parameters=use_threads=1, samtools=/path/to/samtools, outfile=/path/to/file/merge.bam, min_mapped=50, min_paired=50
  program=java
  program_file=/path/to/MergeSamFiles.jar

=head1 DESCRIPTION

  Merge the bam files to create the "super" bam file needed by
  Bam2Genes, the rough model step.
  It uses samtools to check but it uses picard for merging,
  no config file as everything is fetched from the analysis table
  The module is looking for files named: *_sorted.bam

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::MergeBamFiles;

use warnings ;
use strict;

use File::Copy;
use File::Find;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config;
  return $self;
}

=head2 read_and_check_config

 Example    : $self->read_and_check_config;
 Description: Read the parameters from the analysis table
 Returntype : Void
 Exceptions : Void

=cut

sub read_and_check_config {
    my ($self) = @_;

    logging_info( "Checking directories...\n");
    my @params = split('\s*,\s*', $self->analysis->params);
    $self->_use_threads(0);
    foreach my $param (@params) {
        my ($key, $value) = $param =~ /(\w+)\s*=\s*(\S+)/;
        if ($key eq 'use_threads') {
            throw("Need an integer for $key") unless ($value =~ /\d/);
            $self->_use_threads($value);
        }
        elsif ($key eq 'samtools') {
            throw("$value is not executable") unless (-x $value);
            $self->_samtools($value);
        }
        elsif($key eq 'outfile') {
            my ($dir) = $value =~ /(.*)\/[^\/]+$/;
            throw("$dir does exists") unless (-d $dir);
            $self->outfile($value);
        }
    }
    foreach my $dir (split('\s*,\s*', $self->analysis->db_file)) {
        throw("$dir does not exists") unless (-d $dir);
        $self->_bam_dir($dir);
    }
}

sub fetch_input {
    my ($self) = @_;

    foreach my $dir (@{$self->_bam_dir}) {
        opendir(my $dh, $dir) || throw("Could not open $dir");
        while (readdir $dh) {
            $self->_bam_files("$dir/$_") if ($_ =~ /_sorted.bam$/ and -f "$dir/$_");
        }
        closedir($dh) || throw("Could not close directory $dir");
    }
    my $num_bam_files = scalar(@{$self->_bam_files});
    throw("You're missing some bam files, you have $num_bam_files bam files instead of ".$self->analysis->db_version) unless ($self->analysis->db_version == $num_bam_files);
    if ($num_bam_files == 1) {
        warning("There is only one file, no need to merge...");
        $self->input_is_void(1);
    }
}

sub run {
    my ($self) = @_;

    my $picard_cmd_files = '';
    my @stat_checks;
    foreach my $bam_file (@{$self->_bam_files}) {
        my $cmd = $self->_samtools.' flagstat '.$bam_file;
        my $index = 0;
        open(PR, "$cmd > 2&1 | ") || throw("Could not open: $cmd");
        while(<PR>) {
            throw("File $bam_file is truncated") if (/truncated/);
            if (/mapped\s+\(([0-9.]+)/) {
                if ($1 < $self->_min_mapped) {
                    warning("Bam file $bam_file has less than ".$self->_min_mapped."% reads aligned");
                }
                $index = 0;

            }
            elsif (/properly paired \(([0-9.]+)/) {
                if ($1 < $self->_min_paired) {
                    warning("Bam file $bam_file has less than ".$self->_min_paired."% reads properly paired");
                }
            }
            if (/^\s*(\d+)/) {
                $stat_checks[$index++] += $1;
            }
        }
        close(PR) || throw("Could not close: $cmd");
        $picard_cmd_files .= ' INPUT='.$bam_file;
    }
    $self->_stat_check(\@stat_checks);
    my $options = ' MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT';
    $options .= ' USE_THREADING=true' if ($self->_use_threads);
    my $picard_cmd = $self->analysis->program.' -Xmx2g -jar '.$self->analysis->program_file.$options.' OUTPUT='.$self->outfile.$picard_cmd_files;
    throw("Picard cmd: $picard_cmd FAILED") unless (system($picard_cmd));
}

sub write_output {
    my ($self) = @_;

    my $outfile = $self->outfile;
    my $cmd = $self->_samtools.' flagstat '.$outfile;
    my $index = 0;
    my @checks;
    open(PR, "$cmd > 2&1 | ") || throw("Could not open: $cmd");
    while(<PR>) {
        throw("File ".$self->outfile." is truncated") if (/truncated/);
        if (/^\s*(\d+)/) {
            $checks[$index++] += $1;
        }
    }
    close(PR) || throw("Could not close: $cmd");
    my $stat_check = $self->_stat_check;
    for (my $i = 0; $i < $#checks; $i++) {
        throw("Merge failed for $outfile: ".$checks[$i].' <=> '.$stat_check->[$i]) unless ($checks[$i] == $stat_check->[$i]);
    }
    my $index_file = $outfile;
    $index_file =~ s/bam$/bai/;
    # Return 1 on success
    throw("Failed to move index file $index_file to $outfile.bai") unless (move($index_file, $outfile.'.bai'));
    # Return 0 on success
    throw("Failed to update the timestamp for $outfile.bai") if (system("touch $outfile.bai"));
}

###
# Containers
###

=head2 _stat_check

 Arg [1]    : Arrayref, (optional)
 Example    : $self->_stat_check(\@checks);
 Description: Getter/Setter for the stat from samtools flagstat
 Returntype : Arrayref
 Exceptions : None

=cut

sub _stat_check {
    my ($self, $value) = shift;

    if ($value) {
        $self->{'__stat_check'} = $value;
    }
    return $self->{'__stat_check'};
}

=head2 _use_threads

 Arg [1]    : Boolean, 0 or 1
 Example    : $sef->_use_threads(1);
 Description: Getter/Setter, set this to 1 if you want picard to use threads, is ~20% faster but use ~20% more memory
 Returntype : Boolean
 Exceptions : None

=cut

sub _use_threads {
    my ($self, $value) = shift;

    if ($value) {
        $self->{'__use_threads'} = $value;
    }
    return $self->{'__use_threads'};
}

=head2 _samtools

 Arg [1]    : String, '/usr/bin/samtools'
 Example    : $self->_samtools('/usr/bin/samtools');
 Description: Getter/Setter for Samtools path
 Returntype : String, samtools path
 Exceptions : None

=cut

sub _samtools {
    my ($self, $value) = shift;

    if ($value) {
        $self->{'__samtools'} = $value;
    }
    return $self->{'__samtools'};
}

=head2 outfile

 Arg [1]    : String, '/my/data/dir/merged.bam'
 Example    : $self->outfile('/my/data/dir/merged.bam');
 Description: Getter/Setter for the merged bam file
 Returntype : String, merged bam file path
 Exceptions : None

=cut

sub outfile {
    my ($self, $value) = shift;

    if ($value) {
        $self->{'_outfile'} = $value;
    }
    return $self->{'_outfile'};
}

=head2 _min_mapped

 Arg [1]    : Int, any value between 0 and 100
 Example    : $self->_min_mapped(50);
 Description: Getter/Setter for the minimum number of mapped read before a warning
              If you want to know which data align less of XX% of your reads
 Returntype : Int
 Exceptions : None
 
 
=cut

sub _min_mapped { 
    my ($self, $value) = shift;

    if ($value) {
        $self->{'__min_mapped'} = $value;
    }
    return $self->{'__min_mapped'};
}

=head2 _min_paired

 Arg [1]    : Int, between 0 and 100
 Example    : $self->_min_paired;
 Description: Getter/Setter for the minimum nnumber of properly paired reads before a warning
              If you want to know which data align less of XX% of your reads
 Returntype : Int
 Exceptions : None

=cut

sub _min_paired {
    my ($self, $value) = shift;

    if ($value) {
        $self->{'__min_paired'} = $value;
    }
    return $self->{'__min_paired'};
}

=head2 _bam_dir

 Arg [1]    : String, a directory containing bam files
 Example    : $self->_bam_dir('/my/data/dir/tobam');
 Description: Getter/Setter to add one or multiple directories containing bam files to be merged
 Returntype : Arrayref, list of string
 Exceptions : None

=cut

sub _bam_dir {
    my ($self, $value) = shift;

    if ($value) {
        push(@{$self->{'__bam_dir'}}, $value);
    }
    return $self->{'__bam_dir'};
}

=head2 _bam_files

 Arg [1]    : String, path to a bam files
 Example    : $self->_bam_files('/my/data/files/tobam/myfile_sorted.bam');
 Description: Getter/Setter to add one or multiple bam files to be merged
 Returntype : Arrayref, list of string
 Exceptions : None

=cut

sub _bam_files {
    my ($self, $value) = shift;

    if ($value) {
        push(@{$self->{'__bam_files'}}, $value);
    }
    return $self->{'__bam_files'};
}

1;
