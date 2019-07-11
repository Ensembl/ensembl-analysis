=head1 LICENSE

# Copyright [2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE->new
  (
   -query => $slice,
   -program => 'trnascan-SE',
  );
  $runnable->run;
  my @simple_features = @{$runnable->output};

=head1 DESCRIPTION

tRNAscan_SE expects to run the program tRNAscan-SE and produces 
SimpleFeature which can be stored in the simple_feature table in the 
core database

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE
  Function  : produce a tRNAscan_SE runnable and set the options
  and program if undefined
  Returntype: Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  ######################
  #SETTING THE DEFAULTS#
  ######################
  if(!$self->options){
    $self->options('-q');
  }
  if(!$self->program){
    $self->program('tRNAscan-SE');
  }

  ######################
  return $self;
}


=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string, directory
  Function  : a generic run method. This checks the directory specifed
  to run it, write the query sequence to file, marks the query sequence
  file and results file for deletion, runs the analysis parses the
  results and deletes any files
  Returntype: 1
  Exceptions: throws if no query sequence is specified
  Example   :

=cut

sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);

  unless($self->query) {
    throw("Can't run ".$self." without a query sequence")
  }

  $self->checkdir();

#  my $program = $self->program;
  my $program = "/hps/nobackup2/production/ensembl/fergal/coding/trnascan/install_dir/bin/tRNAscan-SE";
  # options = -H -q --detail
  my $query_filename = $self->write_seq_file();
  my $initial_results_filename = $self->create_filename.".out";
  my $secondary_structure_filename = $initial_results_filename.".ss";
  my $filtered_results_prefix =  $initial_results_filename.".filt";
  $self->files_to_delete($initial_results_filename);
  $self->files_to_delete($secondary_structure_filename);
  $self->files_to_delete($filtered_results_prefix.".out");
  $self->files_to_delete($filtered_results_prefix.".ss");

  say "FERGAL QUERY: ".$query_filename;
  say "FERGAL INITIAL: ".$initial_results_filename;
  say "FERGAL SS: ".$secondary_structure_filename;

  my $command = $program." ";
  $command .= $query_filename;
  $command .= " -o ".$initial_results_filename." -f ".$secondary_structure_filename;
  $command .= " -H -q --detail";#.$self->options." " if($self->options);

  system($command) == 0 or throw("FAILED to run ".$command);

  # This is slightly annoying because of how the filtering module requires an output dir to be specified
  # and the file names already have tmp on them. So in this case we need to remove tmp from the name and
  # specify it as the output dir
  my $filtered_prefix = $initial_results_filename;
  $filtered_prefix =~ s/\/tmp\///;
  my $high_confidence_command = "/hps/nobackup2/production/ensembl/fergal/coding/trnascan/install_dir/bin/EukHighConfidenceFilter".
                                " --result ".$initial_results_filename.
                                " --ss ".$secondary_structure_filename.
                                " --output /tmp/ ".
                                " --prefix ".$filtered_prefix.".filt";
  system($high_confidence_command) == 0 or throw("FAILED to run ".$high_confidence_command);
  $self->parse_results($filtered_results_prefix);
  $self->delete_files;
  return 1;
}


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE
  Arg [2]   : string, filename
  Function  : to open and parse the results file
  Returntype: none
  Exceptions: throws on failure to open or close the results file
  or if the results file doesnt exist
  Example   :

=cut

sub parse_results{
  my ($self, $filtered_results_prefix) = @_;

  my $filtered_file = $filtered_results_prefix.".out";
  my $structure_file = $filtered_results_prefix.".ss";
  unless(-e $filtered_file) {
    throw("Could not find the filtered output file. Filename used: ".$filtered_file);
  }

  unless(-e $structure_file) {
    throw("Could not find the structure output file. Filename used: ".$structure_file);
  }

  my $output = [];
  open(IN,$filtered_file) or throw("FAILED to open ".$filtered_file);
  while(<IN>){
    my $line = $_;
    if($line =~ /high confidence set/) {
      my $biotype = 'tRNA';
      my $gene = $self->create_gene($line,$biotype,$structure_file);
      push(@$output,$gene);
    } elsif($line =~ /\tpseudo\,/) {
      my $simple_feature = $self->create_simple_feature($line);
      push(@$output,$simple_feature);
    } else {
      next;
    }
  }
  close IN;
  $self->output($output);
}


sub create_gene {
  my ($self,$line,$biotype,$structure_file) = @_;

  my @element = split (/\s+/, $line);
  my ($name, $start, $end, $trna_type, $score) = @element[0, 2, 3, 4, 8];
  my $strand = 1;
  if($start > $end){
    $strand = -1;
    my $temp_end = $start;
    $start = $end;
    $end = $temp_end;
  }

  my $slice = $self->slice();

  my $gene = Bio::EnsEMBL::Gene->new();
  $gene->analysis($self->analysis);
  $gene->biotype($biotype);
  $gene->slice($self->slice);
  $gene->start($start);
  $gene->end($end);
  $gene->strand($strand);
  $gene->description($trna_type);

  my $transcript = Bio::EnsEMBL::Transcript->new();
  $transcript->analysis($self->analysis);
  $transcript->biotype($biotype);
  $transcript->slice($self->slice);
  $transcript->start($start);
  $transcript->end($end);
  $transcript->strand($strand);
  $transcript->description($trna_type);

  my $exon = Bio::EnsEMBL::Exon->new();
  $exon->slice($slice);
  $exon->start($start);
  $exon->end($end);
  $exon->strand($strand);
  $exon->phase(-1);
  $exon->end_phase(-1);

  $transcript->add_Exon($exon);
  $gene->add_Transcript($transcript);

  return($gene);
}

sub create_simple_feature {
  my ($self,$line) = @_;

  my @element = split (/\s+/, $line);
  my ($name, $start, $end, $display_label, $score) = @element[0, 2, 3, 4, 8];
  my $strand = 1;
  if($start > $end){
    $strand = -1;
    my $temp_end = $start;
    $start = $end;
    $end = $temp_end;
  }

  unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
    warning("The start and end were not integers, something went wrong. Skipping");
    next;
  }

  my $ff = $self->feature_factory;
  my $simple_feature = $ff->create_simple_feature($start,$end,$strand,$score,$display_label,$name,$self->query);
  return($simple_feature);
}

1;
