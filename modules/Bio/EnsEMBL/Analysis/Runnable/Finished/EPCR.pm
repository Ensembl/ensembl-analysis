# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

EPCR.pm

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::EPCR->new(
      -query => $slice,
      -program => $self->analysis->dbfile,
      %{$self->parameters_hash};
     );
  $runnable->run;
  my @marker_features = @{$runnable->output};

=head1 DESCRIPTION

The Finished version of EPCR

=head1 CONTACT

Post general queries to B<anacode@sanger.ac.uk>

=cut

package Bio::EnsEMBL::Analysis::Runnable::Finished::EPCR;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::EPCR;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::EPCR);



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
  my $hash;
  my $name;
  
  ## Parse output file and store data in hash table
  open(FH, $results) or throw("FAILED to open ".$results);
  while(<FH>){
    chomp;
    #chromosome:NCBI34:1:1:920598:1 615132..615340	121028
    #chromosome:NCBI34:1:1:920598:1 622477..622687	121028
    my ($start, $end, $dbid);
    ($name, $start, $end, $dbid)  = $_ =~ m!(\S+)\s+(\d+)\.\.(\d+)\s+(\w+)!;
    $hash->{$dbid} = [] unless $hash->{$dbid};
    # store markers in hash table
    push @{$hash->{$dbid}}, [$start, $end];
  }
  
  ## Filter data to dismiss nested marker features
  DBID:foreach my $dbid (keys %$hash) {
  	# sort the markers by coord start
  	my @markers = sort {$a->[0] <=> $b->[0]} @{$hash->{$dbid}};
  	M1:foreach my $m1 (@markers) {
		next unless $m1;
  		my $start1 = $m1->[0];
  		my $end1 = $m1->[1];
  		#print "M1 $start1 -> $end1 dbid $dbid\n";
		M2:for(my $i = 0;$i<scalar(@markers);$i++){
			my $m2 = $markers[$i];
			next unless $m2;
			my $start2 = $m2->[0];
	  		my $end2 = $m2->[1];
	  		# ignore if marker 1 & 2 identical
	  		next M2 if($start1==$start2 && $end1==$end2);
			# stop after marker 1 end
	  		next M1 if($end1<= $start2);
	  		# delete marker 2 if nested in marker 1
	  		if($start1<=$start2 && $end1>=$end2){
  				#print "\tM2 $start2 -> $end2 DELETE\n";
	  			delete $markers[$i];
	  		}
  		}		  		
  	}
  	$hash->{$dbid} = \@markers;
  }
  
  ## Create marker feature objects
  foreach my $dbid (keys %$hash) {
  	my $m = $ff->create_marker($dbid);
  	$m->adaptor($self->analysis->adaptor);
  	foreach my $coord (@{$hash->{$dbid}}) {
  		next unless $coord;
  		my $start = $coord->[0];
  		my $end = $coord->[1];
    	my $mf = $ff->create_marker_feature($start, $end, 0, $m, $name,$self->query);
    	push(@output, $mf);
  	}
    $self->hit_list($dbid);
  	
  }
  $self->output(\@output);
  close(FH) or throw("FAILED to close ".$results);
}

1;
