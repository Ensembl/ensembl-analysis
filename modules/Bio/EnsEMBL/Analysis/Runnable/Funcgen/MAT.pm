# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Funcgen::MAT

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::MAT->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @features = @{$runnable->output};

=head1 DESCRIPTION

MAT expects to run the program MAT (Model-based Analysis of Tiling-array, Johnson 
et al., PMID: 16895995) and predicts features which can be stored in the 
annotated_feature table in the eFG database


=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::MAT;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);

=head2 query

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : integer
  Function  : container for the chip number
  Returntype: interger
  Exceptions: throws if passed an object which isnt a slice
  Example   : 

=cut

sub query {
    my ($self, $chip) = @_;
    if ($chip){
        throw("Must pass RunnableDB::Funcgen::MAT::query an integer ".
              "specifying the chip to br processed not ".$chip) 
            unless($chip =~ m/^\d+$/);
        $self->{'chip'} = $chip;
    }
    return $self->{'chip'};
}

=head2 write_infile

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::MAT
    Arg [2]     : filename
    Description : 
    Returntype  : 
    Exceptions  : 
    Example     : 

=cut

sub write_infile {
	
	my ($self, $filename) = @_;

    my $chipno = sprintf("%02d", $self->query);

    if (exists $ENV{INFILE_EXISTS} && $ENV{INFILE_EXISTS}) {

        opendir(DIR, $self->workdir())
            or throw("Can't open dir ".$self->workdir());
        my @files = grep { /.$chipno.tag$/ } readdir DIR;

        throw("More than one tag file fir chip $chipno in ".$self->workdir)
            if (! @files);
        $filename = $self->workdir.'/'.$files[0];

        closedir DIR;

    }

	if (! $filename) {
		($filename = $self->infile()) =~ s/\.dat$/.$chipno.tag/;
	}

    (my $logfile = $filename) =~ s/\.tag$/.log/;
    
    $ENV{BPMAPFILE} =~ s/<CHIPNO>/$chipno/;
    $ENV{CELFILES} =~ s/<CHIPNO>/$chipno/g;

    open(TAG, ">$filename")
        or throw("Can't open .tag-file $filename.");
    print TAG <<EOT;

[data]
BpmapFolder = $ENV{BPMAPFOLDER}
CelFolder = $ENV{CELFOLDER}
GenomeGrp = $ENV{GENOMEGRP}
RepLib = $ENV{REPLIB}
Group = $ENV{GROUP}
Pair = $ENV{PAIR}

[bpmap]
1 = $ENV{BPMAPFILE}

[cel]
1 = $ENV{CELFILES}

[intensity analysis]
BandWidth = $ENV{BANDWIDTH}
MaxGap = $ENV{MAXGAP}
MinProbe  = $ENV{MINPROBE}

[interval analysis]
Matscore =
Pvalue = 1e-5
FDR =
Extend =

[output]
Log = $logfile

EOT
    
    close TAG;

    $self->infile($filename);

    (my $bedfile = $filename) =~ s/\.tag$/.bed/;
    $self->resultsfile($bedfile);

    # set columns (fields) for output
    my @fields = (1..2,4); # bed

    $self->output_fields(\@fields);

    return $filename;

}

=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::MAT
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

sub run_analysis {

    my ($self, $program) = @_;
    
    if(!$program){
        $program = $self->program;
    }
    throw($program." is not executable") 
        unless($program && -x $program);

    if (! $ENV{INFILE_EXISTS}) {
        
        my $command = $self->program . ' ' . $self->infile;
        
        warn("Running analysis " . $command . "\n");
        
        #eval { system($command) };
        throw("FAILED to run $command: ", $@) if ($@);
    }
}

=head2 infile

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::MAT
  Arg [2]     : filename (string)
  Description : will hold a given filename or if one is requested but none
  defined it will use the create_filename method to create a filename
  Returntype  : string, filename
  Exceptions  : none
  Example     : 

=cut

sub infile{

    my ($self, $filename) = @_;
    
    if($filename){
        $self->{'infile'} = $filename;
    }
    if(!$self->{'infile'}){
        $self->{'infile'} = $self->create_filename($self->analysis->logic_name, 'dat');
    }
    
    return $self->{'infile'};
    
}

sub config_file {
    my $self = shift;
    $self->{'config_file'} = shift if(@_);
    return $self->{'config_file'};
}

1;
