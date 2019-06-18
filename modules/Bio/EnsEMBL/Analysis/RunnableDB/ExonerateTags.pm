=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateTags - 

=head1 SYNOPSIS

my $analysis_obj = $db->get_AnalysisAdaptor->fetch_by_logic_name("Exonerate_Tags");
my $runnabledb =
    Bio::EnsEMBL::Analysis::RunnableDB::ExonerateTags->new( -db          => $db,
                                                            -analysis    => $analysis_obj};
$runnabledb->fetch_input();
$runnabledb->run();
$runnabledb->write_output();

=head1 DESCRIPTION

RunnnableDB module for the Ditag analysis. Fetches Ditag sequences from database,
tries to align them to the genome with the Runnable and stores DitagFeatures
after filtering the hits.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateTags;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTags;
use Bio::SeqIO;
use Bio::Seq;
use Bio::EnsEMBL::Analysis::Config::ExonerateTags;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);


=head2 new

  Arg       : various
  Function  : create a runnableDB for the Ditag analysis
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::ExonerateTags
  Caller    : general

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($DITAG_CONFIG);

  return $self;
}


=head2 fetch_input

  Args      : none
  Function  : fetch data from the database and create runnable
  Returntype: none
  Caller    : general

=cut

sub fetch_input{
  my ($self) = @_;

  #define target
  my $target = $self->GENOMICSEQS;
  if ( -e $target ){
    if(!(-d $target)and !(-s $target)) {
      throw("'$target' isn't a file or a directory?");
    }
  } else {
    throw("'$target' could not be found");
  }

  #define query: batch of ditags, split in two
  my @fasta_entries = ();
  my $tmp_dir       = $self->TMPDIR;
  #expected format of input ids: ditag.<TYPE>.<NUMBER>
  $self->input_id   =~ /ditag\.(\S+)\.([0-9]+)/;
  my $tagtype       = $1;
  my $input_index   = $2;
  throw("Couldn t get index number from input-id ".$self->input_id)
      unless($input_index);
  throw("Couldn t get tag-type from input-id ".$self->input_id)
      unless($tagtype);

  #fetch all desired ditags
  my $batch_size = $self->BATCHSIZE();
  #my $tagtype    = $self->TAGTYPE();
  my $ditags     = $self->db->get_ditagAdaptor->fetch_all_with_limit(
                   $tagtype, $batch_size, ($batch_size * ($input_index - 1)) );
  my $chop_first = $self->CHOPFIRST();

  my $chop_last  = $self->CHOPLAST();

  foreach my $ditag (@$ditags) {
    #avoid polyA sequences
    $ditag = $self->_check_seq($ditag);
    if($ditag) {

      my $ditag_fasta_entry;
      if($self->SPLITSEQS()) {
        #split the seq and store as fasta line, remove specific bases if necessary
        $ditag_fasta_entry = $self->_split_seq($ditag, $chop_first, $chop_last);
      } else {
        my $tagseq = $ditag->sequence;
        #remove leading bases? (e.g. G from GIS tags)
        if($chop_first){
          $tagseq = $self->_chop_first($tagseq, $chop_first);
        }
        #remove trailing bases? (e.g. AA)
        if($chop_last){
          $tagseq = $self->_chop_last($tagseq, $chop_last);
        }
        $ditag_fasta_entry = ">".$ditag->dbID."_F"."\n".$tagseq."\n";
      }
      push(@fasta_entries, $ditag_fasta_entry);
    }
  }

  #write fasta lines to temp-file for exonerate
  my $tmp_file = $tmp_dir."/".$input_index.".tmp";
  open(OUTFILE, ">$tmp_file") or throw "couldnt create fasta chunk file $tmp_file";
  foreach my $ditag_fasta_entry (@fasta_entries) {
    print OUTFILE $ditag_fasta_entry;
  }
  close(OUTFILE) or throw "couldnt close fasta chunk file $tmp_file";

  #define other parameters
  my %parameters = %{ $self->parameters_hash };

  if (
    not exists( $parameters{-options} )
    and defined $self->OPTIONS
  ){
    $parameters{-options} = $self->OPTIONS;
  }

  if($self->analysis->db_file){
    $parameters{'-STS_FILE'} = $self->analysis->db_file 
      unless($parameters{'-STS_FILE'});
  }

  my $program = ($self->PROGRAM()) ? $self->PROGRAM() : $self->analysis->program_file;
  if(!$program){ throw "Cant decide what program to use for analysis!\n" }

  #create runnable
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTags->new
    (
     -program         => $program,
     -analysis        => $self->analysis,
     -target_file     => $target,
     -query_file      => $tmp_file,
     -deletequeryfile => 0,
     -maxmismatch     => $self->MAXMISMATCH(),
     -specoptions     => $self->SPECOPTIONS(),
     -splitseqs       => $self->SPLITSEQS(),
     -maxdistance     => $self->MAXDISTANCE(),
     -keeporder       => $self->KEEPORDER(),
     %parameters
    );
  $self->runnable($runnable);

}


=head2 _chop_first

  Arg [1]   : sequence string of tag
  Arg [2]   : sequence string of bases to remove
  Function  : Remove undesired sequences at the start of a tag,
              eg. G from GIS tags if this specified for the tag lib.
  Returntype: sequence string
  Caller    : private

=cut

sub _chop_first {
  my ($self, $ditagseq, $chop_first) = @_;

  if($ditagseq =~ /^$chop_first.+$/){
    $ditagseq =  substr($ditagseq, length($chop_first), length($ditagseq)-length($chop_first));
  }

  return $ditagseq;
}


=head2 _chop_last

  Arg [1]   : sequence string of tag
  Arg [2]   : sequence string of bases to remove
  Function  : Remove undesired sequences from the end of a tag,
              eg. AA from GIS tags.
  Returntype: sequence string
  Caller    : private

=cut

sub _chop_last {
  my ($self, $ditagseq, $chop_last) = @_;

  if($ditagseq =~ /^$chop_last.+$/){
    $ditagseq =  substr($ditagseq, 0, length($ditagseq)-length($chop_last));
  }

  return $ditagseq;
}


=head2 _check_seq

  Arg [1]   : Bio::EnsEMBL::Map::Ditag
  Function  : Remove undesired sequences, eg. repetative seqs.
  Returntype: Bio::EnsEMBL::Map::Ditag or undef
  Caller    : private

=cut

sub _check_seq {
  my ($self, $ditag) = @_;

  my $repeat_number = $self->REPEATNUMBER();
  #remove repetetive sequences
  if( ($ditag->sequence =~ /A{$repeat_number}/) or ($ditag->sequence =~ /T{$repeat_number}/) or
      ($ditag->sequence =~ /C{$repeat_number}/) or ($ditag->sequence =~ /G{$repeat_number}/) ){
#  if( $ditag->sequence =~ /A{$repeat_number} | T{$repeat_number} | C{$repeat_number} | G{$repeat_number}/ ){
    print "removing repetitive ditag ".($ditag->dbID)." ".($ditag->tag_count)."\n";
    $ditag = undef;
  }
  #remove too short or too long sequences
  elsif(defined($self->MINSEQLENGTH) and (length($ditag->sequence) < $self->MINSEQLENGTH)){
    print "removing short ditag ".($ditag->dbID)." ".($ditag->tag_count)."\n";
    $ditag = undef;
  }
  elsif(defined($self->MAXSEQLENGTH) and (length($ditag->sequence) > $self->MAXSEQLENGTH)){
    print "removing long ditag ".($ditag->dbID)." ".($ditag->tag_count)."\n";
    $ditag = undef;
  }
  #remove sequences with Ns
  elsif($ditag->sequence =~ /N/){
    print "removing ditag with 'N' ".($ditag->dbID)." ".($ditag->tag_count)."\n";
    $ditag = undef;
  }

  return $ditag;
}


=head2 _split_seq

  Arg [1]   : Bio::EnsEMBL::Map::Ditag
  Function  : split a ditag into two parts to save to file
              FOR NOW THESE ARE JUST 2 HALVES
  Returntype: string with sequence header (_L & _R) and sequence 
              for the parts (FASTA format)
  Caller    : private

=cut

sub _split_seq {
  my ($self, $ditag, $chop_first, $chop_last) = @_;
  my $addone = 0;

  my $motherseq = $ditag->sequence;
  my $id = $ditag->dbID;

  #remove leading bases? (e.g. G from GIS tags)
  if($chop_first){
    $motherseq = $self->_chop_first($motherseq, $chop_first);
  }
  #remove trailing bases? (e.g. AA)
  if($chop_last){
    $motherseq = $self->_chop_last($motherseq, $chop_last);
  }

  #split seq in two parts, TODO: adjust for other types
  #define border
  my $cutoffpoint = (length($motherseq))/2;
  if($cutoffpoint =~ /\./){
    $addone = 1;
  }
  $cutoffpoint = int($cutoffpoint);

  #first part of seq
  my $start_A    = 0;
  my $end_A      = $cutoffpoint;
  if($addone){ $end_A++; }
  my $seqpart_A  = ">".$id."_L_".$start_A."\n";
  $seqpart_A    .= (substr($motherseq, 0, $end_A))."\n";

  #second part of seq
  my $start_B    = $cutoffpoint;
  my $seqpart_B  = ">".$id."_R_".($start_B)."\n";
  $seqpart_B    .= (substr($motherseq, $start_B))."\n";

  return($seqpart_A.$seqpart_B);
}


=head2 run

  Args      : none
  Function  : run the Ditag runnable
  Returntype: none
  Caller    : general

=cut

sub run {
  my ($self) = @_;
  my @out_features;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my $runnable = @{ $self->runnable }[0];

  $runnable->run;

  @out_features = @{$runnable->output};
  $self->output(\@out_features);
}


=head2 write_output

  Arg [1]   : array ref of Bio::EnsEMBL::Map::DitagFeatures
  Function  : get ditag feature adaptor with new features,
              clean and store them.
  Returntype: none
  Caller    : general

=cut

sub write_output {
  my ( $self, $ditagfeatures ) = @_;

  my $ditagfeature_adaptor;
  if($self->OUTDB) {
    my $outdb             = $self->create_output_db;
    $ditagfeature_adaptor = $outdb->get_DitagFeatureAdaptor;
  } else {
    $ditagfeature_adaptor = $self->db->get_DitagFeatureAdaptor;
  }

  if(!$ditagfeatures or !(scalar @$ditagfeatures)){
    $ditagfeatures = $self->output();
  }
  print "\nhave ".(scalar @$ditagfeatures)." features to store.\n";
  if(!$ditagfeatures or !(scalar @$ditagfeatures)){
    warn "no features to store.";
  }
  else {
    $ditagfeatures = $self->clean_features($ditagfeatures);
    #eval { $ditagfeature_adaptor->store($ditagfeatures) };
    eval { $ditagfeature_adaptor->batch_store($ditagfeatures) };
    if ($@) {
      $self->throw("Unable to store ditagFeatures!\n $@");
    }
  }
}


=head2 clean_features

  Arg [1]   : array ref of Bio::EnsEMBL::Map::DitagFeatures
  Function  : re-set analysis, attach slice
  Returntype: array ref of Bio::EnsEMBL::Map::DitagFeatures
  Caller    : general

=cut

sub clean_features {
  my ( $self, $ditagfeatures ) = @_;
  my %genome_slices;
  my @cleanfeatures;
  my $sa = $self->db->get_SliceAdaptor;

  foreach my $ditagfeature (@$ditagfeatures) {

    #(re-)set analysis
    $ditagfeature->analysis($self->analysis);

    #attach slice
    my $slice_id = $ditagfeature->slice;
    if ( not exists $genome_slices{$slice_id} ) {
      print "\nTRYING TO FETCH $slice_id.";
      $genome_slices{$slice_id} = $sa->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};
    $ditagfeature->slice($slice);
    push(@cleanfeatures, $ditagfeature);

  }

  return \@cleanfeatures;
}


=head2 create_output_db

  Args      : none
  Function  : get DBAdaptor for output db
  Returntype: Bio::EnsEMBL::DBSQL::DBAdaptor
  Caller    : general

=cut

sub create_output_db {
  my ($self) = @_;

  my $outdb;
  my $dnadb;

  print "\nCreating OUTDB.";
  if ( $self->OUTDB && $self->DNADB) {
    $dnadb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(%{ $self->OUTDB }); 

    $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
              %{ $self->OUTDB }, 
              -dnadb => $dnadb 
             );

  } elsif( $self->OUTDB) {
    $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{ $self->OUTDB }); 
  } else {
    $outdb = $self->db;
  }

  return $outdb;
}


=head2 runnable_path

  Args      : none
  Function  : GETTER functions for path to Runnable
  Returntype: string
  Caller    : general

=cut

sub runnable_path{
  my ($self);
  return "Bio::EnsEMBL::Analysis::Runnable::DitagAlign";
}


=head2 GENOMICSEQS

  Arg [1]   : (optional) new value
  Function  : GETTER/SETTER functions for Config values
  Returntype: string or undef
  Caller    : private

=cut

sub GENOMICSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_GENOMICSEQS'} = $value;
  }

  return $self->{'_GENOMICSEQS'};
}

sub QUERYSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_QUERYSEQS'} = $value;
  }

  return $self->{'_QUERYSEQS'};
}

sub TAGTYPE {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_TAGTYPE'} = $value;
  }

  return $self->{'_TAGTYPE'};
}

sub IIDREGEXP {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_IIDREGEXP'} = $value;
  }

  return $self->{'_IIDREGEXP'};
}

sub OUTDB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_OUTDB'} = $value;
  }

  return $self->{'_OUTDB'};
}

sub DNADB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_DNADB'} = $value;
  }

  return $self->{'_DNADB'};
}

sub OPTIONS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_OPTIONS'} = $value;
  }

  return $self->{'_OPTIONS'};
}

sub PROGRAM {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_PROGRAM'} = $value;
  }

  return $self->{'_PROGRAM'};
}

sub TMPDIR {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_TMPDIR'} = $value;
  }

  return $self->{'_TMPDIR'};
}

sub BATCHSIZE {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_BATCHSIZE'} = $value;
  }

  return $self->{'_BATCHSIZE'};
}

sub QUERYFILES {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_QUERYFILE'} = $value;
  }

  return $self->{'_QUERYFILE'};
}

sub SPLITSEQS {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_SPLITSEQS'} = $value;
  }

  return $self->{'_SPLITSEQS'};
}

sub CHOPFIRST {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_CHOPFIRST'} = $value;
  }

  return $self->{'_CHOPFIRST'};
}

sub CHOPLAST {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_CHOPLAST'} = $value;
  }

  return $self->{'_CHOPLAST'};
}

sub MAXDISTANCE {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_MAXDISTANCE'} = $value;
  }

  return $self->{'_MAXDISTANCE'};
}

sub MAXMISMATCH {
  my $self = shift;
  $self->{'_MAXMISMATCH'} = shift if (@_);

  return $self->{'_MAXMISMATCH'};
}

sub SPECOPTIONS {
  my $self = shift;
  $self->{'_SPECOPTIONS'} = shift if (@_);

  return $self->{'_SPECOPTIONS'};
}

sub MINSEQLENGTH {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_MINSEQLENGTH'} = $value;
  }

  return $self->{'_MINSEQLENGTH'};
}

sub MAXSEQLENGTH {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_MAXSEQLENGTH'} = $value;
  }

  return $self->{'_MAXSEQLENGTH'};
}

sub REPEATNUMBER {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_REPEATNUMBER'} = $value;
  }

  return $self->{'_REPEATNUMBER'};
}

sub KEEPORDER {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_KEEPORDER'} = $value;
  }

  return $self->{'_KEEPORDER'};
}

sub TAGCOUNT {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'_TAGCOUNT'} = $value;
  }

  return $self->{'_TAGCOUNT'};
}

1;
