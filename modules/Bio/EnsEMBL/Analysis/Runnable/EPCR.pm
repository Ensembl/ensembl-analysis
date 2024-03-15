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

Bio::EnsEMBL::Analysis::Runnable::EPCR - 

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::EPCR->new(
      -query => $slice,
      -program => $self->analysis->dbfile,
      %{$self->parameters_hash};
     );
  $runnable->run;
  my @marker_features = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run EPCR and parse the results into marker features

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::EPCR;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : string, sts file
  Arg [3]   : arrayref of Bio::EnsEMBL::Map::Markers
  Arg [4]   : int, margin
  Arg [5]   : int, word_size
  Arg [6]   : int, min mismatch
  Arg [7]   : int, max mismatch
  Function  : create a Bio::EnsEMBL::Analysis::Runnable::EPCR
  Returntype: Bio::EnsEMBL::Analysis::Runnable::EPCR
  Exceptions: throws if passed both an sts file and an array of features
  Example   : 

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($sts_file, $sts_features, $margin, $word_size, $min_mismatch,
      $max_mismatch) = rearrange(['STS_FILE', 'STS_FEATURES', 'M',
                                  'W', 'NMIN', 'NMAX'], @args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->program('e-PCR') if(!$self->program);
  ######################

  if($sts_file && $sts_features){
    throw("Must pass either an STS_FILE $sts_file or an array of ".
          "STS_FEATURES $sts_features not both");
  }
  $self->sts_file($sts_file) if($sts_file);
  $self->sts_features($sts_features) if($sts_features);
  $self->margin($margin) if($margin);
  $self->word_size($word_size) if($word_size);
  $self->min_mismatch($min_mismatch) if(defined $min_mismatch);
  $self->max_mismatch($max_mismatch) if(defined $max_mismatch);

  return $self;
}


#containers


=head2 margin

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : int/string variable
  Function  : container for the specified variable. This pod
  refers the the 4 methods below, margin, word_size, min_mistmatch
  and max_mismatch
  Returntype: int/string
  Exceptions: none
  Example   : 

=cut



sub margin{
  my $self = shift;
  $self->{'margin'} = shift if(@_);
  return $self->{'margin'};
}


=head2 word_size

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : int/string variable
  Function  : container for the specified variable. This pod
  refers the the 4 methods below, margin, word_size, min_mistmatch
  and max_mismatch
  Returntype: int/string
  Exceptions: none
  Example   : 

=cut
sub word_size{
  my $self = shift;
  $self->{'word_size'} = shift if(@_);
  return $self->{'word_size'};
}

=head2 min_mismatch

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : int/string variable
  Function  : container for the specified variable. This pod
  refers the the 4 methods below, margin, word_size, min_mistmatch
  and max_mismatch
  Returntype: int/string
  Exceptions: none
  Example   : 

=cut
sub min_mismatch{
  my $self = shift;
  $self->{'min_mismatch'} = shift if(@_);
  return $self->{'min_mismatch'};
}

=head2 max_mismatch

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : int/string variable
  Function  : container for the specified variable. This pod
  refers the the 4 methods below, margin, word_size, min_mistmatch
  and max_mismatch
  Returntype: int/string
  Exceptions: none
  Example   : 

=cut
sub max_mismatch{
  my $self = shift;
  $self->{'max_mismatch'} = shift if(@_);
  return $self->{'max_mismatch'};
}


=head2 sts_file

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : string, file path
  Function  : container for sts file path, will use the find file
  method for Runnable to locate the file 
  Returntype: string
  Exceptions: none
  Example   : 

=cut



sub sts_file{
  my ($self, $file) = @_;
  if($file){
    my $found = $self->find_file($file);
    $self->{'sts_file'} = $found;
  }
  return $self->{'sts_file'};
}


=head2 sts_features

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : arrayref of Bio::EnsEMBL::Map::Markers
  Function  : container for arrayref of Bio::EnsEMBL::Map::Markers
  Returntype: arrayref
  Exceptions: throw if not passed an arrayref or if first element of
  array isnt a Bio::EnsEMBL::Map::Marker
  Example   : 

=cut


sub sts_features{
  my ($self, $features) = @_;
  if($features){
    throw("Must pass EPCR sts_features an arrayref not ".$features)
      unless(ref($features) eq 'ARRAY');
    my $test = $features->[0];
    if(!$test || !($test->isa("Bio::EnsEMBL::Map::Marker"))){
      my $err = "arrayref ".$features." must contain ".
        "Bio::EnsEMBL::Map::Marker";
      $err .= " not ".$test if($test);
      throw($err);
    }
    $self->{'sts_features'} = $features;
  }
  return $self->{'sts_features'};
}



=head2 hit_list

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : int, hit database id
  Function  : take the hit ids passed and store in a hash
  Returntype: hashref
  Exceptions: none
  Example   : 

=cut


sub hit_list{
  my ($self, $hit) = @_;
  if(!$self->{'hit_list'}){
    $self->{'hit_list'} = {};
  }
  if($hit){
    $self->{'hit_list'}->{$hit} = $hit;
  }
  return $self->{'hit_list'};
}

#utility methods

=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : string, working directory
  Function  : coordinates the running and parsing of EPCR
  EPCR is run for every mismatch value between the minimum and maximum
  Returntype: none
  Exceptions: throws if doesnt have a query sequence or if doesnt have 
  either and sts file or sts features
  Example   : 

=cut


sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  throw("Can't run ".$self." without a query sequence") 
    unless($self->query);
  $self->checkdir($dir);
  my $filename = $self->write_seq_file();
  $self->files_to_delete($filename);
  my $mismatch = $self->min_mismatch;
  my $sts_file;
  while($mismatch <= $self->max_mismatch){
    if($self->sts_features){
      $sts_file = $self->dump_sts_features($self->sts_features);
    }elsif($self->sts_file){
      $sts_file = $self->copy_sts_file($self->sts_file);
    }else{
      throw("Don't have either sts feature or a file");
    }
    my $results = $self->run_epcr($mismatch, $self->queryfile, $sts_file);
    $self->parse_results($results);
    $mismatch++;
  }
  $self->delete_files;
}


=head2 run_epcr

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : int, the mismatch 
  Arg [3]   : string, query sequence filename
  Arg [4]   : string, sts filename
  Function  : construct commandline and run epcr
  Returntype: string, filename
  Exceptions: throws if system call fails
  Example   : 

=cut


sub run_epcr{
  my ($self, $mismatch, $query_file, $sts_file) = @_;
  my $results = $self->resultsfile;
  if(-e $results){
    $results .= ".".$mismatch.".results";
  }
  $mismatch = $self->min_mismatch if(!$mismatch);
  $query_file = $self->query_file if(!$query_file);
  $sts_file = $self->sts_file if(!$sts_file);
  my $options;
  $options = " M=".$self->margin if(defined $self->margin);
  $options .= " W=".$self->word_size if(defined $self->word_size);
  $options .= " N=$mismatch " if(defined $mismatch);
  my $command = $self->program." ".$sts_file." ".$query_file." ".
    $options." > ".$results;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
  $self->files_to_delete($results);
  return $results;
}



=head2 dump_sts_features

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : arrayref of Bio::EnsEMBL::Map::Markers
  Function  : dump the markers to a file in the format expected by epcr
  markers who have already been hit or whose max_primer_dist is 0 or
  has zero left or right markers are ignored
  Returntype: 
  Exceptions: 
  Example   : 

=cut


sub dump_sts_features{
  my ($self, $sts_features, $filename) = @_;
  if(!$sts_features){
    $sts_features = $self->sts_features;
  }

  my %hit_list = %{$self->hit_list};
 
  if(!$filename){
    $filename = $self->create_filename("sts", "out");
  }
  $self->files_to_delete($filename);
  open(OUT, ">".$filename) or throw("FAILED to open $filename");
 MARKER:foreach my $m (@$sts_features){
    next MARKER if($hit_list{$m->dbID});
    next MARKER if($m->max_primer_dist == 0);
    next MARKER unless(length($m->left_primer) > 0);
    next MARKER unless(length($m->right_primer) > 0);
    my $string = $m->dbID."\t".$m->left_primer."\t".
      $m->right_primer."\t".$m->min_primer_dist."-".
        $m->max_primer_dist."\n";
    print OUT $string;
  }
  close(OUT) or throw("FAILED to close $filename");
  return $filename;
}



=head2 copy_sts_file

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : string, sts file
  Arg [3]   : string, filename for copy
  Function  : copys entries from one file into another while
  skipping the entries already on the hit list
  Returntype: filename
  Exceptions: throws if given sts file doesnt exist or if any of the
  open or close file commands fail
  Example   : 

=cut



sub copy_sts_file{
  my ($self, $sts_file, $filename) = @_;
  if(!$sts_file){
    $sts_file = $self->sts_file;
  }
  if(! -e $sts_file){
    throw("Can't copy file ".$sts_file." which doesn't exist");
  }
  if(!$filename){
    $filename = $self->create_filename("sts", "out");
  }
  $self->files_to_delete($filename);
  my %hit_list = %{$self->hit_list};
  if(keys(%hit_list) == 0){
    return $sts_file;
  }
  eval{
    open(OUT, ">".$filename) or throw("FAILED to open $filename");
    open(IN, $sts_file) or throw("FAILED to open $sts_file");
  MARKER:while(<IN>){
      my $id = (split)[0];
      next MARKER if($hit_list{$id});
      print OUT;
    }
    close(OUT) or throw("FAILED to close $filename");
    close(IN) or throw("FAILED to close $sts_file");
  };
  if($@){
    throw("FAILED to copy $sts_file $@");
  }
  return $filename;
}


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EPCR
  Arg [2]   : string, filename
  Function  : parse the file given into Bio:EnsEMBL::Map::MarkerFeatures
  Returntype: none
  Exceptions: throws if results file doesnt exist or if open and closes
  fail
  Example   : 

=cut


sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  if(!-e $results){
    throw("Can't open ".$results." as it doesn't exist");
  }
  my $ff = $self->feature_factory;
  my @output;
  open(FH, $results) or throw("FAILED to open ".$results);
  while(<FH>){
    chomp;
    #chromosome:NCBI34:1:1:920598:1 615132..615340	121028
    #chromosome:NCBI34:1:1:920598:1 622477..622687	121028
    my ($name, $start, $end, $dbid) = 
      $_ =~ m!(\S+)\s+(\d+)\.\.(\d+)\s+(\w+)!;
    my $m = $ff->create_marker($dbid);
    my $mf = $ff->create_marker_feature($start, $end, 0, $m, $name,
                                        $self->query);
    push(@output, $mf);
    $self->hit_list($dbid);
  }
  $self->output(\@output);
  close(FH) or throw("FAILED to close ".$results);
}
