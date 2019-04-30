=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::BaseExonerate - 

=head1 SYNOPSIS

  Do NOT instantiate this class directly: must be instantiated
  from a subclass (see ExonerateTranscript, for instance).

=head1 DESCRIPTION

This is an abstract superclass to handle the common functionality for 
Exonerate runnables: namely it provides
- a consistent external interface to drive exonerate regardless of
  what features youre finally producing, and
- a common process to stop people duplicating function (eg how to
  arrange command-line arguments, what ryo-string to use etc).

It does NOT provide the parser to convert the exonerate output
into Transcripts or AffyFeatures etc. That is the job of the
subclasses, which MUST implement the parse_results method.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile parse_timer execute_with_timer);

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [QUERY_TYPE]: String
 Arg [QUERY_SEQS]: Arrayref of Bio::Seq
 Arg [QUERY_FILE]: String
 Arg [QUERY_CHUNK_NUMBER]: Integer
 Arg [QUERY_CHUNK_TOTAL]: Integer
 Arg [TARGET_SEQS]: Arrayref of Bio::Seq
 Arg [TARGET_FILE]: String
 Arg [TARGET_CHUNK_NUMBER]: Integer
 Arg [TARGET_CHUNK_TOTAL]: Integer
 Arg [ANNOTATION_FEATURES]: Hashref
 Arg [ANNOTATION_FILE]: String
 Arg [VERBOSE]: Integer
 Arg [BASIC_OPTIONS]: String
 Description: Creates a new Bio::EnsEMBL::Analysis::Runnable::BaseExonerate object
 Returntype : Bio::EnsEMBL::Analysis::Runnable::BaseExonerate
 Exceptions : Throws if QUERY_SEQS or TARGET_SEQS is defined and not an arrayref
              Throws if ANNOTATION_FILE, QUERY_FILE or TARGET_FILE is defined and does not exists
              Throws if ANNOTATION_FEATURES is not a hashref

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my 
    (
     $query_type, $query_seqs, $query_file, $q_chunk_num, $q_chunk_total,
     $target_seqs, $target_file, $t_chunk_num, $t_chunk_total, 
     $annotation_features, $annotation_file, $verbose, $basic_options
    ) =
    rearrange(
      [
        qw(
          QUERY_TYPE
          QUERY_SEQS
          QUERY_FILE
          QUERY_CHUNK_NUMBER
          QUERY_CHUNK_TOTAL          
          TARGET_SEQS
          TARGET_FILE
          TARGET_CHUNK_NUMBER
          TARGET_CHUNK_TOTAL
          ANNOTATION_FEATURES
          ANNOTATION_FILE
          VERBOSE
		  BASIC_OPTIONS
        )
      ], 
      @args
    );

  $self->_verbose($verbose) if $verbose;

  if (defined($query_seqs)) {
    if(ref($query_seqs) ne "ARRAY"){
      throw("You must supply an array reference with -query_seqs");
    }
    $self->query_seqs($query_seqs);
  } elsif (defined $query_file) {
    throw("The given query file ".$query_file." does not exist") 
      if ! -e "$query_file";
    $self->query_file($query_file);
    
  }

  if ($query_type){
    $self->query_type($query_type);
  } else{
    # default to DNA for backwards compatibilty
    $self->query_type('dna');
  }
  
  if (defined $target_seqs) {
    if (ref($target_seqs) ne "ARRAY") {
      throw("You must supply an array reference with -target_seqs");
    }
    $self->target_seqs($target_seqs);
  } elsif (defined $target_file) {
    if ($target_file !~ /^[\w\-\.]+:\d+$/) {
      throw("The given database ($target_file) does not exist") if ! -e "$target_file";
    }
    $self->target_file($target_file);
  }

  if (defined $annotation_features) {
    if (ref($annotation_features) ne "HASH") {
      throw("You must supply a hash reference with -annotation_features");
    }
    $self->annotation_features($annotation_features);
  } elsif (defined $annotation_file) {
    throw("The given annotation file does not exist") if ! -e "$annotation_file";
    $self->annotation_file($annotation_file);
  }

  if (not $self->program) {
    $self->program('exonerate-0.9.0');
  }

 
  #
  # These are what drives how we gather up the output
  $basic_options ||= "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" ";
  

  if (defined $q_chunk_num and defined $q_chunk_total) {
    $basic_options .= "--querychunkid $q_chunk_num --querychunktotal $q_chunk_total ";
  }

  if (defined $t_chunk_num and defined $t_chunk_total) {
    $basic_options .= "--targetchunkid $t_chunk_num --targetchunktotal $t_chunk_total ";
  }

  if ($self->options){
    $basic_options .= $self->options;
  }
  
  $self->options($basic_options);

  return $self;
}



############################################################
#
# Analysis methods
#
############################################################

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs exonerate script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output

=cut

sub run {
  my ($self) = @_;

  my $output_file = $self->create_filename();
  $self->files_to_delete($output_file);

  my $write_results = $self->write_to_file();

  if ($self->annotation_features) {
    my $annot_file = $self->create_filename("exonerate_a");
    $self->annotation_file($annot_file);
    open F, ">$annot_file" or 
        throw "Could not open temp $annot_file for writing";
    foreach my $id (keys %{$self->annotation_features}) {
      my $f = $self->annotation_features->{$id};
      printf(F "%s %s %d %d\n", 
             $f->seqname, 
             $f->strand < 0 ? "-" : "+",
             $f->start,
             $f->length);
    }
    close(F);
    $self->files_to_delete($annot_file);
    close($annot_file) || throw('Could not close annotation file for writing');
  } elsif ($self->annotation_file) {
    my %feats;
    open F, $self->annotation_file or 
        throw("Could not open supplied annotation file for reading");
    while(<F>) {
      /^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/ and do {
        $feats{$1} = Bio::EnsEMBL::Feature->new(-seqname => $1, 
                                                -strand  => $2 eq "-" ? -1 : 1,
                                                -start   => $3,
                                                -end     => $3 + $4 - 1); 
      };
    }
    close(F) || throw('Could not close supplied annotation file for reading');
    $self->annotation_features(\%feats);
  }


  if ($self->query_seqs) {
    # Write query sequences to file if necessary
    my $query_file = write_seqfile($self->query_seqs);
    $self->query_file($query_file);
  }

  if ($self->target_seqs) {
    # Write query sequences to file if necessary
    my $target_file = write_seqfile($self->target_seqs);
    $self->target_file($target_file);
  }

  # Build exonerate command
  my $command =
    $self->program . " " .$self->options .
    " --querytype "  . $self->query_type .
    " --targettype " . $self->target_type .
    ' --query "'  . $self->query_file .
    '" --target "' . $self->target_file.'"';
  $command .= " --annotation " . $self->annotation_file if $self->annotation_features;
  $command .= " > ".$output_file if $write_results;

  # Execute command and parse results
  print STDERR "Exonerate command : $command\n";

  my $timer = $self->timer();
  unless($timer) {
    warning("No timer set, defaulting to 5h");
    $timer = "5h";
  }

  if($write_results) {
    execute_with_timer($command, $self->timer);
    $self->output($self->parse_results($output_file,$write_results));
  } else {
    my $exo_fh;
    open( $exo_fh, "$command |" ) or throw("Error opening exonerate command: $? : $!");
    $self->output($self->parse_results( $exo_fh ));

    if (!close( $exo_fh )) {
      sleep 30;
      throw ("Error closing exonerate command: $? : $!");
    }
  } # end else
  return 1;
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseExonerate
  Arg [2]   : pointer to file-handle
  Function  : This method MUST be coded by the base class to return an array-ref of features.
              arguments - is passed a pointer to a filehandle which is the output
              of exonerate.
              Exonerate's basic options - it's always run with - are:
              --showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\"
              so this tells you what the output file will look like: you have
              to code the parser accordingly.
  Returntype: Listref of <things>
  Example   : 
    my ( $self, $fh ) = @_;
    while (<$fh>){
      next unless /^RESULT:/;
      chomp;
      my (
        $tag, $q_id, $q_start, $q_end, $q_strand, 
        $t_id, $t_start, $t_end, $t_strand, $score, 
        $perc_id, $q_length, $t_length, $gene_orientation,
        @vulgar_blocks
      ) = split;
      ...now do something with the match information and / or vulgar blocks
    }
=cut
sub parse_results {
  throw ("This method must be provided by a subclass and not invoked directly! \n"); 
}

############################################################
#
# get/set methods
#
############################################################


=head2 annotation_features

 Arg [1]    : (optional) Hashref of Bio::EnsEMBL::Feature
 Description: Getter/setter for annotation features
 Returntype : Hashref of Bio::EnsEMBL::Feature
 Exceptions : Throws if any of the feature is not a Bio::EnsEMBL::Feature

=cut

sub annotation_features {
  my ($self, $feats) = @_;
  
  if ($feats){
    foreach my $k (keys %$feats) {
      my $f = $feats->{$k};
      unless ($f->isa("Bio::EnsEMBL::Feature")) {
        throw("annotation features must be Bio::EnsEMBL::Features");
      }
    }
    $self->{_annot_feats} = $feats;
  }
  return $self->{_annot_feats};
}
  

############################################################

=head2 annotation_file

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub annotation_file {
  my ($self, $file) = @_;

  if (defined $file) {
    $self->{_annot_file} = $file;
  }
  return $self->{_annot_file};
}

############################################################

=head2 query_type

 Arg [1]    : (optional) String 'dna' or 'protein'
 Description: Getter/setter
 Returntype : String
 Exceptions : Throws if Arg[1] is not 'dna' or 'protein'

=cut

sub query_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'protein' ){
      throw("not the right query type: $type");
    }
    $self->{_query_type} = $type;
  }
  return $self->{_query_type};
}

############################################################

=head2 query_seqs

 Arg [1]    : (optional) Bio::seq
 Description: getter/setter
 Returntype : Bio::Seq
 Exceptions : Throws if Arg[1] is not a Bio::PrimarySeqI or a Bio::SeqI

=cut

sub query_seqs {
  my ($self, $seqs) = @_;
  if ($seqs){
    unless ($seqs->[0]->isa("Bio::PrimarySeqI") || $seqs->[0]->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_query_seqs} = $seqs;
  }
  return $self->{_query_seqs};
}


############################################################

=head2 queryfile

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub query_file {
  my ($self, $file) = @_;
  
  if ($file) {
    $self->{_query_file} = $file;
  }
  return $self->{_query_file};
}

############################################################

=head2 query_type

 Arg [1]    : None
 Description: Getter
 Returntype : String, 'dna'
 Exceptions : None

=cut

sub target_type {
  my ($self) = @_;

  # the target type has to be DNA, because we are making transcripts

  return 'dna';
}

############################################################

=head2 target_seqs

 Arg [1]    : (optional) Bio::seq
 Description: getter/setter
 Returntype : Bio::Seq
 Exceptions : Throws if Arg[1] is not a Bio::PrimarySeqI or a Bio::SeqI

=cut

sub target_seqs {
  my ($self, $seqs) = @_;
  if ($seqs){
    unless ($seqs->[0]->isa("Bio::PrimarySeqI") || $seqs->[0]->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_target_seqs} = $seqs;
  }
  return $self->{_target_seqs};
}


############################################################

=head2 target_file

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub target_file {
  my ($self, $file) = @_;
  
  if ($file) {
    $self->{_target_file} = $file;
  }
  return $self->{_target_file};
}

############################################################

=head2 _verbose

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub _verbose {
  my ($self, $val) = @_;
  
  if ($val){
    $self->{_verbose} = $val;
  }
  
  return $self->{_verbose};
}


1;

