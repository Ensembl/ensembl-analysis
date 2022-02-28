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

Bio::EnsEMBL::Analysis::Runnable::PSILC - 

=head1 SYNOPSIS

  my $PSILC = Bio::EnsEMBL::Analysis::Runnable::PSILC->new 
    (
     '-trans'     => $transcript,
     '-homologs'  => $homologs,
     '-analysis'  => $analysis,
     '-domain'    => $domains,
     '-input_id'  => $self->input_id,
    );

$runnabledb->fetch_input();
$runnabledb->run();
my @array = @{$runnabledb->output};
$runnabledb->write_output();

=head1 DESCRIPTION

Runnable for PSILC

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::PSILC;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::AlignIO;
use Bio::SimpleAlign;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

  Arg [1]    : none
  Description: Creates runnable
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_trans'} = {};	# transcript to test;
  $self->{'_homologs'} = ();	# array of holmologs to cluster with transcript;
  $self->{'_domains'} = ();	# array ref of protein feature ids;

  my($trans,$domains,$homologs,$input_id, $psilc_work_dir) = $self->_rearrange([qw(
							TRANS
							DOMAIN
							HOMOLOGS
							INPUT_ID
                            PSILC_WORK_DIR
						       )], @args);
  if ($trans) {
    $self->trans($trans);
  }
  if ($homologs) {
    $self->homologs($homologs);
  }
  if ($domains) {
    $self->domains($domains);
  }
  $self->output_dir($input_id);
  $self->program('/nfs/acari/sw4/Pseudogenes/PSILC/psilc1.21/psilc.jar');
  return $self;
}

=head2 run

  Arg [1]    : none
  Description: Makes alignment out of transcripts and runs PSILC
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run{
  my ($self) = @_;
  $self->checkdir();
  $self->write_pfam_ids;
  my $filename = $self->create_filename('PSILC','fasta');
  $self->make_codon_based_alignment;
  $self->run_analysis();
  $self->files_to_be_deleted($filename);
  $self->parse_results;
#  $self->delete_files;
}

=head2 write_pfam_ids

  Arg [1]    : none
  Description: Writes identifiers for PFAM domains that PSILC will search for
in the protein alignment
  Returntype : none
  Exceptions : throws if it cannot open the file to write to
  Caller     : general

=cut

sub write_pfam_ids{
  my ($self)=@_;
  my $dir =  $self->workdir;
  my $domains = $self->domains;
  print STDERR "writing pfam ids to $dir/pfamA\n";
  open (IDS,">$dir/pfamA") or $self->throw("Cannot open file for writing pfam ids $dir/pfamA\n");
  print IDS "$domains";
  close IDS;
  return 1;
}

=head2 run_analysis

  Arg [1]    : scalar - program name 
  Description: runs the PSILC executable
  Returntype : none
  Exceptions : throws if program is not executable
  Caller     : general

=cut

sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Runnable::run_analysis ") 
    unless($program && -x $program);
  my $file = $self->queryfile;
  my $dir =  $self->workdir;
  $file =~ s/^$dir\///;
  my $command = "java -Xmx400m -jar $program";
  $command .= " --align $file.aln";
  $command .= " --seq ".$self->id;
  $command .= " --restrict 1";
 # $command .= " --alignMethod 1";
  $command .= " --max_nodes 10:6";
  $command .= " --repository /ecs2/work2/sw4/PFAM";
  $command .= " > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or  $self->throw("FAILED to run ".$command);
}

=head2 parse_results

  Arg [1]    : scalar - filename 
  Description: parses PSILC results
  Returntype : none
  Exceptions : warns if the results are unreadable
  Caller     : general

=cut

sub parse_results{
  my ($self, $filename) =@_;
  my ($nuc_dom,$prot_dom,$postPMax,$postNmax) = 0;
  my $dir =  $self->workdir;
  my $trans = $self->trans;
  my $id    = $self->id;
  open (PSILC,"$dir/PSILC_WAG_HKY/summary") or $self->warn("Cannot open results file $dir/PSILC_WAG_HKY/summary\n$@\n");
  while(<PSILC>){
    chomp;
    next unless $_ =~ /^$id\/\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $nuc_dom =  $1;
    $prot_dom = $2;
    $postPMax = $3;
    $postNmax = $4;
    last;
  }
  close PSILC;
  unless  ($nuc_dom){
    $self->warn("ERROR unable to parse results file $dir/PSILC_WAG_HKY/summary\n");
    $self->warn("Failed for $id"); 
    return 0;
  }
  my %results = (
		 'id'       => $id,
		 'prot_dom' => $prot_dom,
		 'nuc_dom'  => $nuc_dom,
		 'postPMax' => $postPMax,
		 'postNMax' => $postNmax,
		);
  $self->output(\%results);
  return 1;
}

=head2 files_to_be_deleted

  Arg [1]    : scalar - filename 
  Description: marks the PSILC output files for deletion
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub files_to_be_deleted{
  my ($self,$filename) = @_;
  my $dir =  $self->workdir;
  my $domains = $self->domains;
  $domains =~ s/PF.+\t//;
  chomp $domains;
  # Break filename into its component parts
  my @fnc = split /\./,$filename;
#  $self->files_to_delete($filename.".aln");
  $self->files_to_delete($filename.".out");
  $self->files_to_delete($fnc[0].".nhx");
  $self->files_to_delete("$dir/pfamA");
  $self->files_to_delete("$dir/PSILC_WAG_HKY/domdom");
  $self->files_to_delete("$dir/PSILC_WAG_HKY/posteriorN");
  $self->files_to_delete("$dir/PSILC_WAG_HKY/posteriorP");
  $self->files_to_delete("$dir/PSILC_WAG_HKY/psilcN");
  $self->files_to_delete("$dir/PSILC_WAG_HKY/psilcP");
#  $self->files_to_delete("$dir/PSILC_WAG_HKY/summary");  
  $self->files_to_delete("$dir/PSILC_WAG_HKY/$domains");
  return 1;
}

=head2 output_dir

  Arg [1]    : scalar - input id of the anlysis
  Description: Creates the output directory for the PSILC files
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub output_dir{
  my ($self,$input_id) = @_;
  my $output_dir;
  $output_dir = $self->id;
  if ($self->psilc_work_dir){
    system("mkdir $self->psilc_work_dir/$$input_id/$output_dir");
    $self->workdir($self->psilc_work_dir."/$$input_id/$output_dir");
  }
  else{
    $self->throw("Cannot make output directory\n");
  }
  return 1;
}

sub id {
  my ($self) = @_;
  if ($self->trans->stable_id){
    return $self->trans->stable_id;
  }
  return $self->trans->dbID;
}

=head2 make_codon_based_alignment

  Arg [1]    : none
  Description: Makes a dna alignment based on a protein alignment using 
Bio::Align::Utilities DNA alignment is for the translateable part of the 
transcript only
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub make_codon_based_alignment{
  my ($self) = @_;
  my @aaseqs;
  my $aaseq;
  my %dnaseqs;
  # make alignment file for output
  my $filename = $self->queryfile.".aln";
  my $alignIO = Bio::AlignIO->new
    (
     -file => ">$filename",
     -format => "clustalw",
    );
  my $pep_alignIO = Bio::AlignIO->new
    (
     -file => ">$filename.pep",
     -format => "clustalw",
    );
  # run clustal
  my $clustalw = Bio::Tools::Run::Alignment::Clustalw->new(); 
  # prepare sequences, aaseqs has the proteins used for the alignment
  # dnaseqs is a hash of the tranlateable nucleotide sequence
  # Both sets need to use transcript ids which as also used as the hash keys
  $aaseq = $self->trans->translate;
  # change id of translation to match transcript
  $aaseq->display_id($self->id);
  push @aaseqs,$aaseq;
  $dnaseqs{$self->id} = Bio::Seq->new
    (
     -display_id => $self->id,
     -seq        => $self->trans->translateable_seq,
    );
  foreach my $homolog(@{$self->homologs}){
    next if ($homolog->stable_id eq $self->id);
    $aaseq = $homolog->translate;  
    # change id of translation to match transcript
    $aaseq->display_id($homolog->stable_id);
    push @aaseqs,$aaseq;
    $dnaseqs{$homolog->stable_id} = Bio::Seq->new
      (
       -display_id => $homolog->stable_id,
       -seq        => $homolog->translateable_seq,
      );
  }
  foreach my $seq (@aaseqs){
    print $seq->display_id."\n";
  }

  my $alignment = $clustalw->align(\@aaseqs);
  my $dna_aln = aa_to_dna_aln($alignment,\%dnaseqs);
  foreach my $seq ($dna_aln->each_seq){
    print $seq->display_id."\n";
  }  
  $pep_alignIO->write_aln($alignment);
  $alignIO->write_aln($dna_aln);
  return 1;
}

#######################################
# Containers

=head2 trans

Arg [1]    : array ref
  Description: get/set trans set to run over
  Returntype : array ref to Bio::EnsEMBL::Transcript objects
  Exceptions : throws if not a Bio::EnsEMBL::Transcript
  Caller     : general

=cut

sub trans {
  my ($self, $trans) = @_;
  if ($trans) {
    unless  ($trans->isa("Bio::EnsEMBL::Transcript")){
      $self->throw("Input isn't a Bio::EnsEMBL::Transcript, it is a $trans\n$@");
    }
    $self->{'_trans'} = $trans;
  }
  return $self->{'_trans'};
}

=head2 homologs

Arg [1]    : array ref
  Description: get/set gene set to run over
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub homologs {
  my ($self, $homologs) = @_;
  if ($homologs) {
  foreach my $hom (@{$homologs}){
    unless  ($hom->isa("Bio::EnsEMBL::Transcript")){
      $self->throw("Input isn't a Bio::EnsEMBL::Transcript, it is a $hom\n$@");
      }
    }
    $self->{'_homologs'} = $homologs;
  }
  return $self->{'_homologs'};
}

=head2 domains

Arg [1]    : array ref
  Description: get/set pfam domain identifiers
  Returntype : string of pfam identifiers
  Exceptions : none
  Caller     : general

=cut

sub domains {
  my ($self, $domains) = @_;
  if ($domains) {
    $self->{'_domains'} = $$domains;
  }
  return $self->{'_domains'};
}

=head2 output

  Arg [1]    : hash_ref
  Description: overrides output array
  Returntype : hash_ref
  Exceptions : none
  Caller     : general

=cut

sub output {
  my ($self, $hash_ref) = @_;
  if ($hash_ref) {
    $self->{'_output'} = $hash_ref;
  }
  return $self->{'_output'};
}

=head2 psilc_work_dir

  Arg [1]    : String
  Description: path to work dir
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub psilc_work_dir {
  my ($self, $path) = @_;
  if ($path) {
    $self->{'_psilc_work_dir'} = $path;
  }
  return $self->{'_psilc_work_dir'};
}

1;

