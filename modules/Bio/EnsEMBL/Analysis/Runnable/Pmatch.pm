=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Pmatch - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Runnable::Pmatch;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);
use Bio::EnsEMBL::Analysis::Tools::Pmatch::First_PMF;
use Bio::SeqIO;

use parent qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($protein_file, $protein_length, $max_intron_length, $min_coverage) = 
    rearrange(['PROTEIN_FILE', 'PROTEIN_LENGTHS', 'MAX_INTRON_LENGTH', 
               'MIN_COVERAGE'], @args);

  ###SETTING DEFAULTS###
  $self->program('pmatch') if(!$self->program);
  $self->max_intron_length(50000);
  $self->min_coverage(25);
  #######################

  $self->protein_file($protein_file);
  throw("Need to pass in a defined readable protein file not ".$protein_file) 
    unless($self->protein_file && (-s $self->protein_file));
  $self->protein_lengths($protein_length);
  $self->max_intron_length($max_intron_length);
  $self->min_coverage($min_coverage);
  return $self;
}


sub protein_file{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'protein_file'} = $arg;
  }
  return $self->{'protein_file'};
}

sub protein_lengths{
  my ($self, $arg) = @_;
  if($arg){
    throw("Runnable::Pmatch::protein_lengths must be passed a hashref not a ".$arg)
      unless(ref($arg) eq 'HASH');
    $self->{'protein_lengths'} = $arg;
  }
  if(!$self->{'protein_lengths'}){
    my $hash = $self->make_protein_length_hash;
    $self->{'protein_lengths'} = $hash;
  }
  return $self->{'protein_lengths'};
}

sub max_intron_length{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'max_intron_length'} = $arg;
  }
  return $self->{'max_intron_length'};
}

sub min_coverage{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'min_coverage'} = $arg;
  }
  return $self->{'min_coverage'};
}


sub make_protein_length_hash{
  my ($self) = @_;
  my %hash;
  my $seqio = Bio::SeqIO->new(-format => 'Fasta',
                              -file => $self->protein_file);
  while(my $seq = $seqio->next_seq){
    $hash{$seq->id} = $seq->length;
  }
  return \%hash;
}

sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Pmatch::run_analysis ") 
    unless($program && -x $program); 
  my $command = $self->program." -D ".$self->options." ".$self->protein_file." ".
     $self->queryfile." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  logger_info("Running analysis ".$command);
  execute_with_wait($command);
}

sub parse_results {
  my ($self) = @_;

  my (%prot_ids, @all_merged_hits);

  open RES, $self->resultsfile or $self->throw("Could not open " . $self->resultsfile . " results file\n");
  while(<RES>) {
    my @cols = split;
    $prot_ids{$cols[5]}++;
  }
  close(RES);

  my @idlist = sort keys %prot_ids;
  my @lists;
  # chop list into group ensuring that the total number of features
  # in each group does not exceed 1M
  foreach my $id (sort keys %prot_ids) {
    my $size = $prot_ids{$id};
    if (not @lists or $lists[-1]->{size} + $size > 1000000) {
      push @lists, {
        size => 0,
        ids  => [],
      };
    }
    push @{$lists[-1]->{ids}}, $id;
    $lists[-1]->{size} += $size;
  }

  foreach my $list (@lists) {
    my %these_ids = map { $_ => 1 } @{$list->{ids}};
    my (%hits);

    open RES, $self->resultsfile or $self->throw("Could not re-open " . $self->resultsfile . " results file\n");
    while(<RES>) {
      my @cols = split;
      if (exists $these_ids{$cols[5]}) {
        push @{$hits{$cols[5]}}, [$cols[2], $_];
      }
    }
    close(RES);

    foreach my $id (keys %hits) {
      my $pmf = new Bio::EnsEMBL::Analysis::Tools::Pmatch::First_PMF(
                                                                     -plengths => $self->protein_lengths,
                                                                     -maxintronlen => $self->max_intron_length,
                                                                     -min_coverage => $self->min_coverage,
                                                                     );
      foreach my $el (sort { $a->[0] <=> $b->[0] } @{$hits{$id}}) {
        $_ = $el->[1];
        $pmf->make_coord_pair($_);
      }
      my @merged_hits = @{$pmf->merge_hits};

      push @all_merged_hits, @merged_hits;
    }
  }

  my $ff = $self->feature_factory;

  foreach my $hit (@all_merged_hits) {
    my $start = $hit->qstart;
    my $end = $hit->qend;
    if($start > $end){
      ($start, $end) = ($end, $start);
    }
    my $hstart = 1;
    my $hlength = ($end - $start +1)/3;
    my $hend = sprintf("%.0f", $hlength);
    my $fp = $ff->create_feature_pair($start, 
                                      $end, 
                                      $hit->strand,
                                      $hit->coverage, 
                                      $hstart, 
                                      $hend,
                                      1, 
                                      $hit->target, 
                                      100, 
                                      0, 
                                      $hit->query, 
                                      $self->query, 
                                      $self->analysis);
    my $paf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => [$fp], -align_type => 'ensembl');
    $self->output([$paf]);
  }
}

1;
