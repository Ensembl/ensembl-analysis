=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::FirstEF - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::FirstEF->new
  (
   -query => $slice,
   -program => 'cpg',
  );
  $runnable->run;
  my @simple_features = @{$runnable->output};

=head1 DESCRIPTION

FirstEF expects to run the program FirstEF and produces SimpleFeature which
can be stored in the simple_feature table in the core database


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::FirstEF;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, path to parameters directory
  Arg [3]   : string, path to parse script
  Function  : create a Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Returntype: Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Exceptions: 
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    

  my($param_dir, $parse_script) = rearrange(["PARAM_DIR", "PARSE_SCRIPT"], @args);

  ##################
  #SETTING DEFAULTS#
  ##################
  $self->program('firstef') if(!$self->program);
 ##################

  $self->param_dir($param_dir) if($param_dir);
 
  $self->parse_script($parse_script) if($parse_script);

  return $self;
}


##container methods




=head2 parsed_output

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, file path
  Function  : container for filename, if none is passed in one is generated
  the first time it is requested
  Returntype: string
  Exceptions: 
  Example   : 

=cut


sub parsed_output{
  my ($self, $file) = @_;
  if($file){
    $self->files_to_delete($file);
    $self->{'parsed_output'} = $file;
  }
  if(!$self->{'parsed_output'}){
    my $file = $self->create_filename('first_parse');
    $self->files_to_delete($file);
    $self->{'parsed_output'} = $file;
  }
  return $self->{'parsed_output'};
}


=head2 listfile

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, filename
  Function  : this is a container for the listfile required by firstef
  If no filename is passed in or the file doesnt already exist the
  file is prepared
  Returntype: string, filename
  Exceptions: 
  Example   : 

=cut



sub listfile{
  my ($self, $listfile) = @_;
  if($listfile){
    $self->{'listfile'} = $listfile;
    $self->files_to_delete($listfile);
  }
  if(!$self->{'listfile'} || ! -e $self->{'listfile'}){
    my $file = $self->prepare_listfile($self->{'listfile'});
    $self->{'listfile'} = $file;
  }
  return $self->{'listfile'};
}


=head2 parse_script

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, path to script
  Function  : container for path to script
  Returntype: string
  Exceptions: throws if script doesnt exist
  Example   : 

=cut


sub parse_script{
  my ($self, $script) = @_;
  if($script){
    throw($script." does not exist can't use FirstEF:parse_script")
      unless(-e $script);
    $self->{'parse_script'} = $script;
  }
  return $self->{'parse_script'};
}


=head2 param_dir

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, directory path
  Function  : stores the parameters directory and checks all the files
  which should exist in it do
  Returntype: string
  Exceptions: throws if any of the specified files dont exist
  Example   : 

=cut


sub param_dir {
  my $self = shift;

  if (@_) {
    $self->{'param_dir'} = shift;

    $self->{'param_dir'} .= '/' unless $self->{'param_dir'} =~ /\/$/;

    if($self->{'param_dir'}){
      my @known_param_files = ('donor.3mer_wts_GChighDown',
                               'donor.3mer_wts_GClowDown',
                               'donor.6mer_wts_GChighDown',
                               'donor.6mer_wts_GChighUp',
                               'donor.6mer_wts_GClowDown',
                               'donor.6mer_wts_GClowUp',
                               'donor.decisiontree',
                               'donor.decisiontree.orig',
                               'donor.qdamodel.GChigh',
                               'donor.qdamodel.GClow',
                               'exon.qdamodel.CpGpoor_GChigh',
                               'exon.qdamodel.CpGpoor_GClow',
                               'exon.qdamodel.CpGrich_GChigh',
                               'exon.qdamodel.CpGrich_GClow',
                               'promoter.5mer_wts_CpGpoor_430.510',
                               'promoter.5mer_wts_CpGpoor_490.570',
                               'promoter.5mer_wts_CpGrich_430.510',
                               'promoter.5mer_wts_CpGrich_490.570',
                               'promoter.6mer_wts_CpGpoor_1.250',
                               'promoter.6mer_wts_CpGpoor_1.450',
                               'promoter.6mer_wts_CpGpoor_200.450',
                               'promoter.6mer_wts_CpGrich_1.250',
                               'promoter.6mer_wts_CpGrich_1.450',
                               'promoter.6mer_wts_CpGrich_200.450',
                               'promoter.qdamodel.CpGpoor',
                               'promoter.qdamodel.CpGrich');
      my @missing_files;
      foreach my $param_file (@known_param_files) {
        unless (-e $self->{'param_dir'} . "/$param_file"){
          push (@missing_files, $self->{'param_dir'}."/$param_file");
        }
      }
      if(@missing_files > 0){
        print STDERR join("\n", @missing_files);
        throw("The above parameter files are missing.");
      }
    }
  }

  return $self->{'param_dir'}
}


=head2 prepare_listfile

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, filename
  Function  : createa a filename or uses the one passed in
  opens this file and write a string, queryfilename -1500
  and then closes the file and adds it plus a couple of other names
  to the files to delete list
  Returntype: string
  Exceptions: throws if fails to open or close file
  Example   : 

=cut


sub prepare_listfile{
  my ($self, $listfile) = @_;
  if(!$listfile){
    $listfile = $self->create_filename('firstef_listfile_', '', 
                                       $self->workdir);
    open(LISTFILE, ">$listfile") or throw("FAILED to open $listfile ".
                                          "FirstEF:prepare_listfile");

    print LISTFILE $self->queryfile."  -1500\n";

    close(LISTFILE) or throw("FAILED to close $listfile ".
                             "FirstEF:prepare_listfile");
    $self->files_to_delete($listfile);
    $self->files_to_delete($listfile."_domain");
    $self->files_to_delete($listfile."_domain_comp");
    $self->resultsfile($listfile."_out");
    $self->files_to_delete($self->resultsfile);
  }
  return $listfile;
}


=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, program name
  Function  : constructs commandline and runs commandline
  Returntype: 
  Exceptions: throws if system doesnt return a 0
  Example   : 

=cut


sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  my $command = $program." 1500 ".$self->listfile." ".$self->param_dir.
    " 0 0.4 0.4 0.5";
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, results filename
  Function  : parses results, initially using an external script which
  also does some statistics and then using a simple regex to produce
  simple features
  Returntype: none
  Exceptions: throws if fails to open or close first parses output file
  Example   : 

=cut


sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $ff = $self->feature_factory;
  my @output;
  my $output = $self->first_parse($results);
  throw("FAILED to run ".$self->parse_script." ".$output." doesn't exist")
    unless(-e $output);
  open(FH, $output) or throw("FAILED to open ".$output);
  my $strand;
 LINE:while(<FH>){
    chomp;
    $strand = 1 if (/direct strand/);
    $strand = -1 if (/complementary strand/);
    if (/\d+\s+\S+\s+\S+([^\.]+)\.\.(\d+)\s+(\S+)\s+\S+\s+\S+\s+(\d+)/) {
      my ($start, $end) = sort {$a <=> $b} ($1 * 1 , $2 * 1);
      my $sf = $ff->create_simple_feature($start, $end, $strand, $3,
                                          "rank = $4", '', $self->query);
      push(@output, $sf);
    }
  }
  $self->output(\@output);
}


=head2 first_parse

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::FirstEF
  Arg [2]   : string, a results filename
  Function  : runs the parsing script across the initial firstef output
  Returntype: string, filename
  Exceptions: throws if script doesnt return 0
  Example   : 

=cut



sub first_parse{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $command = $self->parse_script." ".$results." ".$self->parsed_output;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
  return $self->parsed_output;
}


1;
