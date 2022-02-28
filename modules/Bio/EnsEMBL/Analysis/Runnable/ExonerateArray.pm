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

=head1 AUTHORS

Eduardo Eyras

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateArray - 

=head1 SYNOPSIS
$database  = a full path location for the directory containing the target (genomic usually) sequence,
@sequences = a list of Bio::Seq objects,
$exonerate = a location for the binary,
$options   = a string with options ,

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateArray->new(
								       -db            =>$db,
								       -query_seqs    => \@sequences,
								       -program       => $exonerate,
								       -options       => $options,
								);

 $runnable->run; #create and fill Bio::Seq object
 my $results = $runnable->output;
 
 where $results is an arrayref of MiscFeatures.
 
=head1 DESCRIPTION

ExonerateArray takes a Bio::Seq (or Bio::PrimarySeq) object and runs Exonerate
against a set of sequences.  The resulting output file is parsed
to produce a set of features.


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateArray;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(info verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($db,$query_seqs) =
    rearrange([qw(
		  DB
		  QUERY_SEQS
		 )
	      ], @args);
  

  ###$db is needed to create $slice which is needed to create DnaDnaAlignFeatures
  $self->db($db) if $db;
  # must have a target and a query sequences
  unless( $query_seqs ){
    throw("Exonerate needs a query_seqs: $query_seqs");
  }

  our (%length);
  
  my $queryfile = $self->queryfile();
  
  foreach my $query_seq (@{$query_seqs}) {
    $length{$query_seq->display_id} = $query_seq->length;

    $self->write_seq_file ($query_seq);
  }

  my @lengths = sort {$b<=>$a} values %length;
  my $max_length = $lengths[0];
  $self->max_length($max_length);

  $self->length(\%length);

  return $self;
}

############################################################
#
# Analysis methods
#
############################################################

sub write_seq_file{
  my ($self, $seq, $filename) = @_;
  if(!$seq){
    $seq = $self->query;
  }
  if(!$filename){
    $filename = $self->queryfile;
  }
  
  my $seqout = Bio::SeqIO->new(
                               -file => ">>".$filename, ###added >>
                               -format => 'Fasta',
			      );
  eval{
    $seqout->write_seq($seq) if $seq;
  };
  
  if($@){
    throw("seq is $seq\nFAILED to write $seq to $filename Runnable:write_seq_file : $@");
  }

  return $filename;
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

  if(!$dir){
    $dir = $self->workdir;
  }
  $self->checkdir($dir);
  my $filename = $self->write_seq_file();
  $self->files_to_delete($filename);
  $self->files_to_delete($self->resultsfile);
  $self->run_analysis();
  $self->parse_results;
  $self->delete_files;
  return 1;
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::ExonerateArray
  Arg [2]   : string, filename
  Function  : open and parse the results file into misc_features
  features
  Returntype: none
  Exceptions: throws on failure to open or close output file
  Example   :

=cut


sub parse_results{
  my ($self, $results) =  @_;
  if(!$results){
    $results = $self->resultsfile;
  }

  my %length = %{$self->length()};

  open( EXO, $results ) || throw("FAILED to open ".$results." ExonerateArray::parse_results");
  
  
  ############################################################
  # store each alignment as a features
  
  my (@pro_features);
  
  ############################################################
  # parse results - avoid writing to disk the output
  
  while (<EXO>){
    
    #info ($_) ;
    
    ############################################################
    # the output is of the format:
    #
    #
    # vulgar contains 9 fields
    # ( <qy_id> <qy_start> <qy_len> <qy_strand> <tg_id> <tg_start> <tg_len> <tg_strand> <score> ), 
    # 
    # The vulgar (Verbose Useful Labelled Gapped Alignment Report) blocks are a series 
    # of <label, query_length, target_length> triplets. The label may be one of the following: 
    #
    # M    Match 
    #
    # example:
    # vulgar: probe:HG-U95Av2:1002_f_at:195:583; 0 25 + 10.1-135037215 96244887 96244912 + 125 M 25 25 100
    # match_length is 25 bs, it may not be exact match, score 125->exact match, score 116 match 24 bs
    # if it is 120 M 24 24, it means exact match for 24 bs.
    
    if (/^vulgar\:/) {
      my $h={};
      chomp;
      my ( $tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, $match, $matching_length) = split;
      
      # the VULGAR 'start' coordinates are 1 less than the actual position on the sequence
      $q_start++;
      $t_start++;
      
      my $strand;
      if ($q_strand eq $t_strand) {
	$strand = 1;
      }
      else {
	$strand = -1;
      }
      
      $h->{'q_id'} = $q_id;
      $h->{'q_start'} = $q_start;
      $h->{'q_end'} = $q_end;
      $h->{'q_strand'} = $strand;
      $h->{'t_id'} = $t_id;
      $h->{'t_start'} = $t_start;
      $h->{'t_end'} = $t_end;
      $h->{'t_strand'} = $strand;
      $h->{'score'} = $score;
      $h->{'probe_length'} = $length{$h->{'q_id'}};
      $h->{'matching_length'} = $matching_length;
        
      ###for affymetrix probe sequence, they are 25 bs long, we require at least 24 bs exact match###
      
      if ($h->{'matching_length'} == $h->{'probe_length'}-1) {
	$h->{'match_status'} = "Mismatch";
	push @pro_features, $h;
      } 
      elsif ($h->{'matching_length'} == $h->{'probe_length'}) {
	if ($h->{'score'} == 125) {
	  $h->{'match_status'} = "Fullmatch";
	}
	else {
	  $h->{'match_status'} = "Mismatch";
	}
	push @pro_features, $h;
      }
    }
  }
  
  
  close(EXO) or $self->throw("couldn't close pipe ");  
  
  $self->_make_affy_features(@pro_features);

  ############################################################
  
# remove interim files (but do not remove the database if you are using one)

}


sub _make_affy_features {
  
  
  my ($self,@h) = @_;

  my @misc_feats;

  my $feature_factory = $self->feature_factory;

  ###to make MiscFeature, we need to make slice_obj###
  ###need to watch out for target fasta sequence title for seq_region_name ###
  foreach my $h (@h) {
    my ($coord_system_name,$seq_region_name,$slice) ;
    
    if ($h->{'q_id'} =~/Zebrafish/i) {
      $coord_system_name = undef;
      if ($h->{'t_id'} =~ /^(\S+)\-1\-.*$/) {
	$seq_region_name = $1;
	$seq_region_name =~ s/\-1$//;
	#print "seq_region_name is $seq_region_name\n";
      }
      elsif ($h->{'t_id'} =~ /^(\S+)\..*$/) {
	$seq_region_name = $1;
      }
    }
    else {
      $coord_system_name = "chromosome";
      ($seq_region_name) = $h->{'t_id'} =~ /^(\S+)\..*$/;
    }
    
    $slice = $self->db->get_SliceAdaptor->fetch_by_region($coord_system_name,$seq_region_name);
    if (!$slice) {
      warning("In ExonerateArray::make_affy_features: Could not obtain slice for seq_region: $coord_system_name : $seq_region_name---try with coord_system_name = undef\n");
      $coord_system_name = undef;
      $slice = $self->db->get_SliceAdaptor->fetch_by_region($coord_system_name,$seq_region_name);
      if (!$slice) {
	warning("In ExonerateArray::make_affy_features: Could not obtain slice for seq_region: $coord_system_name : $seq_region_name");
	next;
      }
    }
    my $probe_name = $h->{'q_id'};
    $probe_name =~ s/^probe\://;
    my ($array_name,$composite_name) = split /\:/, $probe_name;
    $composite_name =~ s/\;$//;
    my $xref_name = $array_name;
    $xref_name =~ s/-/_/g;  ###this to keep name same as in code corresponds to external_db.db_name
    
    my $misc_feature = $feature_factory->create_misc_feature 
      ($h->{'t_start'},$h->{'t_end'},$h->{'t_strand'},$slice);
    
    
    $feature_factory->add_misc_feature_attribute 
      ( $misc_feature, "probeName","Probe name","the name of the probe",$probe_name);
	    
    $feature_factory->add_misc_feature_attribute
      ($misc_feature, "compositeName","Composite name", "the name of the composite", $composite_name);

    $feature_factory->add_misc_feature_attribute
      ($misc_feature,"probeLength", "Probe length", "the length of the probe", $h->{'probe_length'});

    $feature_factory->add_misc_feature_attribute
      ($misc_feature, "matchLength", "Match length", "number of bases in matched alignment", $h->{'matching_length'});

    $feature_factory->add_misc_feature_attribute
      ($misc_feature,"matchStatus", "Match status", "full_match or with_mismatch", $h->{'match_status'});
    
    #
    #  Add as many Attributes as you like
    #
    
    $feature_factory->add_misc_set
      ($misc_feature, "AFFY\_$xref_name", $array_name, "array set name", $self->max_length);

    $feature_factory->add_misc_set
      ($misc_feature, "All_Affy", "All-Array-Sets", "all array sets", $self->max_length);

    push (@misc_feats, $misc_feature);
  }
  
  $self->output(\@misc_feats);
    
}

############################################################
#
# get/set methods
#
############################################################

############################################################

sub db {
  my ($self, $db) = @_;
  if ($db) {
    $self->{'db'} = $db ;
  }
  return $self->{'db'};
}

############################################################

sub length {
  my ($self, $length) = @_;
  if ($length) {
    $self->{'length'} = $length ;
  }
  return $self->{'length'};
}

############################################################

sub max_length {
  my ($self, $max_length) = @_;
  if (defined($max_length) ){
    $self->{'max_length'} = $max_length;
  }
  return $self->{'max_length'};
}

############################################################


1;

