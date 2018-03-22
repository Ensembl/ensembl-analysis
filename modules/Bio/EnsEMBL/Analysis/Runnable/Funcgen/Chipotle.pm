# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @annotated_features = @{$runnable->output};

=head1 DESCRIPTION

Chipotle expects to run the program ChIPOTle (PMID: 16277752) and 
predicts features which can be stored in the predicted_feature table 
in the eFG database

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);

=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

sub run_analysis {
        
    print "Analysis::Runnable::Funcgen::Chipotle::run_analysis\n";
    my ($self, $program) = @_;
    
    if(!$program){
        $program = $self->program;
    }
    throw($program." is not executable Chipotle::run_analysis ") 
        unless($program && -x $program);


	# parse parameters
	my $param = $self->get_parameters();

    (my $resultsfile = $self->infile) =~ s/\.dat$/.peaks/;
    $self->resultsfile($resultsfile);
    $self->files_to_delete($resultsfile);

    my @fields = (1..3);
    $self->output_fields(\@fields);

    my $command = join(' ', $self->program, 'get_peaks', '0', '0', 
					   $self->infile(), $self->resultsfile(), 
					   $param->{alpha}, $param->{stepSize},
					   $param->{windowSize}, $param->{adjustPvalue}, 'NORM');
    
    warn("Running analysis " . $command . "\n");

    eval { system($command) };
    throw("FAILED to run $command: ", $@) if ($@);

    # Check if adjustments via --transform need to be applied
	### Version 1 logfile processing disabled ###
    
    #(my $logfile = $self->resultsfile()) =~ s/_peaks.tsv/_STDOUT.txt/;
    #open(LOG, "$logfile")
    #    or throw("Unable to open logfile $logfile.");
    #my $transform;
    #while (<LOG>) {
    #    print;
    #    next until (/^Re-run chipotle using \'(--transform -?\d+\.\d+)\'/);
    #    $transform = $1;
    #    #print "transform: ", $transform, "\n";
    #}
    #close LOG;

    #if (defined $transform) {
    #    $command .= " $transform";
    #    warn("Re-running analysis " . $command . "\n");
    #    
    #    eval { system($command) };
    #    if ($@) {
    #        throw("FAILED to run ".$command);
    #    }
    #}
}

=head2 write_infile

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
    Arg [2]     : filename
    Description : 
    Returntype  : 
    Exceptions  : 
    Example     : 

=cut

sub write_infile {

    my ($self, $filename) = @_;
    print "Analysis::Runnable::Funcgen::Chipotle::write_infile\n";

    if (! $filename) {
        $filename = $self->infile();
    }
    
    # merge infiles if more than one replicate have been selected
    my $replicates = scalar(keys %{$self->result_features});
    #print Dumper $replicates;

    if ($replicates > 1) {

        my @tmpfiles = ();
        foreach my $rset (sort keys %{$self->result_features}) {
            
            my $i = 1;
            
            (my $tmpfile = $self->infile) =~ s/\.dat$/.$rset.dat/;
            print Dumper $tmpfile;
            open(DAT, ">$tmpfile") 
                or throw ("Can't open file $tmpfile");

            map {
                printf DAT "%d\t%s\t%d\t%d\t%f\n", 
                $i++, (@$_)[0..3];
            } @{$self->result_features->{$rset}};

            close DAT 
                or throw ("Can't close file $tmpfile");
            
            push(@tmpfiles, $tmpfile);
            

        }
        
        my $command = join(' ', $self->program, 'merge_files', '0',
                           $replicates, @tmpfiles, $filename);

        warn("Merging infiles " . $command . "\n");
        
        eval { system($command) };
        throw("FAILED to run $command: ", $@) if ($@);

        map {$self->files_to_delete($_)} @tmpfiles;


    } else {

        open(F, ">".$filename)
            or throw("Can't open file $filename.");
        
        #print Dumper $self->result_features;
        foreach (values %{$self->result_features}) {
            foreach my $ft (@$_) {
                #print join("\t", (@$ft)[4,0..3]), "\n";
                print F join("\t", (@$ft)[4,0..3]), "\n";
            }
        }
        close F;
    }
        
    return $filename;

}

1;
