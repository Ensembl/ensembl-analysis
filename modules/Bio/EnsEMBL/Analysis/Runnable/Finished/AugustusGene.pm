
package Bio::EnsEMBL::Analysis::Runnable::Finished::AugustusGene;
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

Bio::EnsEMBL::Analysis::Runnable::Finished::AugustusGene

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus->new
  (
   -query => $slice,
   -program => 'augustus',
   -species => 'human',
   -analysis => $analysis,
  );
$runnable->run;
my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the augustus gene predictor and then parse the results
into prediction transcripts


=head1 CONTACT



=cut

use strict;
use warnings;

use File::Spec::Functions qw(splitpath catfile);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_timer);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;

use parent('Bio::EnsEMBL::Analysis::Runnable');

=head2 new

  Returntype : Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus
=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($species) = rearrange(['SPECIES'], @args);
  $self->species($species);
  $self->program('/homes/thibaut/src/Augustus/bin/augustus') if(!$self->program);
  return $self;
}
 


sub species {
  my ($self,$species) = @_;
  if($species){
    $self->{'species'} = $species;
  }
  return $self->{'species'};
}


sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  throw("Can't run ".$self." without a query sequence") 
    unless($self->query);
  $self->checkdir();
  my $filename = $self->write_seq_file();
  $self->files_to_delete($filename);
  $self->files_to_delete($self->resultsfile);
  $self->run_analysis();
  $self->parse_results;
  $self->delete_files;
  return 1;
}

sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Augustus::run_analysis ")
    unless($program && -x $program);

  my $command = $program." --species=".$self->species." ";
  $command .= $self->options." " if($self->options);
  $command .= $self->queryfile." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);

  $self->resultsfile()

}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus
  Arg [2]   : string, filename
  Function  : parse the output from Augustus into prediction transcripts
  Returntype: none
  Exceptions: throws if cannot open or close results file
  Example   :

=cut

sub parse_results{
  my ($self, $results) = @_;

 if(!$results){
    $results = $self->resultsfile;
  }

  open (OUT, "<".$results) or throw("FAILED to open ".$results."\nAugustus:parse_results");
  my (%transcripts, @transcripts);

  while(<OUT>) {
  chomp;
  next if /^#/;
  my $verbose = 0;

    if(/CDS/){
      my @elements = split"\t",$_;

      if(scalar @elements != 9) {
        throw("Can't parse ".$_." splits into wrong number of elements "."Augustus:parse_results");
      }

      my ($chromosome, $type, $start, $end, $score, $strand, $other) =  @elements[0, 2, 3, 4, 5, 6];
#      my @temp = split' ',$chromosome;
#      my $new_chrom = shift(@temp);
#      my $slice = $slice_adaptor->fetch_by_name($new_chrom);

      my ($transcript_id,$gene_id) = $elements[8] =~ /transcript_id "(.*)"; gene_id "(.*)";/;
      my $group = $transcript_id."_".$gene_id;

      my $exon = Bio::EnsEMBL::Exon->new(-start => $start,
                                         -end   => $end,
                                         -strand => $strand eq '-' ? -1 : 1,
                                         -analysis => $self->analysis,
                                         -slice => $self->query);

      push @{$transcripts{$group}->{exons}}, $exon;
    }
  }

  close(OUT) or throw("FAILED to close ".$results."Augustus:parse_results");

  foreach my $tid (keys %transcripts) {

    my @exons = sort { $a->start <=> $b->start } @{$transcripts{$tid}->{exons}};


    my $tran = Bio::EnsEMBL::Transcript->new(-exons => \@exons,
                                             -analysis => $self->analysis,
                                             -stable_id => $transcripts{$tid}->{hitname},
                                             -slice => $transcripts{$tid}->{slice});

#    $tran->{'accession'} = $transcripts{$tid}->{hitname};
#    $tran->{'pid'} = $transcripts{$tid}->{pid};
#    $tran->{'cov'} = $transcripts{$tid}->{cov};
#    $tran->{'rank'} = $transcripts{$tid}->{rank};
    $tran->{'augustus_id'} = $tid;
    $tran->{'slice_name'} = $transcripts{$tid}->{slice_name};

    # Reverse the exons for negative strand to calc the translations
    my $strand = $exons[0]->strand;
    if($strand == -1) {
      @exons = sort { $b->start <=> $a->start } @{$transcripts{$tid}->{exons}};
    }

    my $start_exon = $exons[0];
    my $end_exon = $exons[scalar(@exons)-1];

    my $translation = Bio::EnsEMBL::Translation->new();

    $translation->start_Exon($start_exon);
    $translation->start(1);
    $translation->end_Exon($end_exon);
    $translation->end($end_exon->length());
    $tran->translation($translation);

    # Set the phases
    calculate_exon_phases($tran, 0);

    # Set the exon and transcript supporting features
#    $self->set_supporting_features($tran);

    push @transcripts, $tran;

#    my $pep = $tran->translate->seq;
#    if ($pep =~ /\*/) {
#      printf STDERR "Bad translation for $tid : $pep\n";
#    }
  }

  $self->clean_output;
  $self->output(\@transcripts);
}

sub database_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_database_adaptor} = $val;
  }

  return $self->{_database_adaptor};
}


=head2 query

    Title   :   query
    Usage   :   $self->query($seq)
    Function:   Get/set method for query.  If set with a Bio::Seq object it
                will get written to the local tmp directory
    Returns :   filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub query {
  my ($self, $val) = @_;

  if (defined $val) {
    if (not ref($val)) {
      throw("[$val] : file does not exist\n") unless -e $val;
    } elsif (not $val->isa("Bio::PrimarySeqI")) {
      throw("[$val] is neither a Bio::Seq not a file");
    }
    $self->{_query} = $val;
  }

  return $self->{_query}
}

sub database {
  my ($self, $val) = @_;

  if (defined $val) {
    if (not ref($val)) {
      $self->warning("[$val] : file does not exist\n") unless -e $val;
    } else {
      if (ref($val) eq 'ARRAY') {
        foreach my $el (@$val) {
          throw("All elements of given database array should be Bio::PrimarySeqs")
        }
      } elsif (not $val->isa("Bio::PrimarySeq")) {
        throw("[$val] is neither a file nor array of Bio::Seq");
      } else {
        $val = [$val];
      }
    }
    $self->{_database} = $val;
  }

  return $self->{_database};
}

1;
