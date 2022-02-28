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
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::IPRScan;

use warnings ;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->options('-cli -appl blastprodom -appl fprintscan -appl hmmpfam -appl superfamily -appl hmmpir -appl scanregexp -format raw') if(!$self->options);
  #$self->options('-cli -appl blastprodom -appl fprintscan -appl hmmpfam -appl hmmpir -appl hmmpanther -appl hmmtigr -appl hmmsmart -appl superfamily  -format raw') if(!$self->options);
  #self->options('-cli -appl blastprodom -appl fprintscan -appl hmmpfam -appl hmmpir -appl hmmpanther -appl hmmtigr -appl hmmsmart -appl superfamily -appl gene3d -format raw')
  ######################

  return $self;
}


#testing the different programs 
#  -appl <name>      Application(s) to run (optional), default is all.
#                    Possible values (dependent on set-up):
# *                            blastprodom
# *                            fprintscan
# 8                            hmmpfam
# 8                            hmmpir
# *                            hmmpanther
# 8                            hmmtigr
# 7                            hmmsmart
# 8                            superfamily
# *                            gene3d
# 8                            scanregexp
# 8                            profilescan
# 8                            seg
# 8                            coils
# 8                            [tmhmm]
# 8                            [signalp]


#################


sub run_analysis{
  my ($self) = @_;
  #example command
  #/nfs/acari/lec/code/interpro/interpro_ftp/iprscan/bin/iprscan -cli -format raw 
  #-i /nfs/acari/lec/code/interpro/interpro_ftp/iprscan/test.seq
  print "RUNNING ANALYSIS\n";
  my $command = $self->program." ".$self->options." -i ".$self->queryfile.
    " >& ".$self->resultsfile;
  print "Running ".$command."\n";
  throw ("Error running ".$command)unless ((system ($command)) == 0);
}




sub parse_results{
  my ($self) =  @_;
  print "PARSING RESULTS\n";
  open(OUT, "<".$self->resultsfile) or throw("FAILED to open ".
                                             $self->resultsfile.
                                             " IPRSCAN:parse_results");
  my @out;
  while(<OUT>){
    print;
    chomp;
    if(/SUBMITTED\s+(\S+)/){
      my $dir_name = $1;
      #print "GETTING filename HAVE ".$dir_name."\n";
      #iprscan-20071115-10500980
      my ($base_dir) = $dir_name =~ /iprscan\-(\d+)\-\d+/;
      my $full_path = $self->workdir."/".$base_dir."/".$dir_name;
      $self->files_to_delete($full_path);
      #print $full_path."\n";
      next;
    }
    next if((!$_) || (!$_ =~ /\d+/));
    #print $_."\n";
    #19897	4AB50DD3FDE96DB6	468	ProfileScan	PS50157	ZINC_FINGER_C2H2_2	243	270	16.685	T	14-Nov-2007
    my ($translation_id, $crd, $aa_length, $method, $hit_name, $desc, $start, $end, $evalue, $status, $date) = split /\t/, $_;
    if(!$start || !$end){
      next;
    }
    if($evalue && !($evalue =~ /\d+/)){
      $evalue = undef;
    }
    if($start > $end){
      my $temp = $start;
      $start = $end;
      $end = $temp;
    }
    my $pf = $self->create_protein_feature($start, $end, 0, $translation_id,
                                           0, 0, $hit_name, $self->analysis,
                                           $evalue, 0);
    throw("Protein feature has no translation id") if(!$pf->seqname);
    push(@out, $pf);
  }
  $self->output(\@out);
  close(OUT);
}

#NF00181542	is the id of the input sequence.
#27A9BBAC0587AB84	is the crc64 (checksum) of the protein sequence (supposed to be unique).
#272	is the length of the sequence (in AA).
#HMMPIR	is the anaysis method launched.
#PIRSF001424	is the database members entry for this match.
#Prephenate dehydratase	is the database member description for the entry.
#1	is the start of the domain match.
#270	is the end of the domain match.
#6.5e-141	is the evalue of the match (reported by member database anayling method).
#T	is the status of the match (T: true, M: marginal).
#06-Aug-2005	is the date of the run.
