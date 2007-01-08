package Bio::EnsEMBL::Analysis::Runnable::Pmatch;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::Pmatch::First_PMF;
use Bio::SeqIO;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


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

sub pm_filters{
  my ($self, $arg) = @_;
  
  if(!$self->{'_pm_filters'}){
    $self->{'_pm_filters'} = [];
  }
  
  if($arg){
    push(@{$self->{'_pm_filters'}}, $arg);
  }
  return $self->{'_pm_filters'};
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
  system($command) == 0 or throw("FAILED to run ".$command);
}


sub parse_results{
  my ($self) = @_;

  my $command = "sort -k6,6 -k3,3n " . $self->resultsfile;
  open(PM, "$command |") or throw("Runnable::Pmatch::parse_results error sorting ".
                                  "the results file" . $self->resultsfile." $!");
  my ($current_pmf, $prot_id);
  while(<PM>){
    #print;
    my @cols = split;
    if(!$prot_id || $cols[5] ne $prot_id){
      $self->pm_filters($current_pmf) if($current_pmf);
   
      $current_pmf = Bio::EnsEMBL::Analysis::Tools::Pmatch::First_PMF
        ->new(
              -plengths => $self->protein_lengths,
              -maxintronlen => $self->max_intron_length,
              -min_coverage => $self->min_coverage,
             );
      $prot_id = $cols[5];
      $current_pmf->make_coord_pair($_);
    }else{
      $current_pmf->make_coord_pair($_);
    }
  }
  $self->pm_filters($current_pmf) if($current_pmf);
  close(PM) or throw("Runnable::Pmatch::parse_results error closing".
                     "the results file" . $self->resultsfile." $! ");
  my @filters = @{$self->pm_filters};
  my $ff = $self->feature_factory;
  my %unique;
  foreach my $f(@filters){
    my @hits = @{$f->merge_hits};
    #print "Have ".@hits." hits\n";
    foreach my $hit(@hits){
      my $start = $hit->qstart;
      my $end = $hit->qend;
      if($start > $end){
        my $temp = $start;
        $start = $end;
        $end = $temp;
      }
      my $hstart = 1;
      my $hlength = ($end - $start +1)/3;
      my $hend = sprintf("%.0f", $hlength);
      my $fp = $ff->create_feature_pair($start, $end, $hit->strand,
                                        $hit->coverage, $hstart, $hend,
                                        1, $hit->target, 100, 0, $hit->query, 
                                        $self->query, $self->analysis);
      my $paf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => [$fp]);
      $self->output([$paf]);
    } 
  }
}
