# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

#!/usr/bin/env perl

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateArtificialChromosome;

use warnings;
use strict;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba);
use File::Path qw(make_path);

sub param_defaults {
  my ($self) = @_;
  return {
      %{$self->SUPER::param_defaults},
      max_chromosome_size => 50000000,
      N_padding_size => 1000,
  }
}


sub fetch_input {
  my($self) = @_;
  my $dba = hrdb_get_dba($self->param('dna_db'));
  $self->create_artificial_chromosome($dba);
}

sub create_artificial_chromosome {
  my ($self, $dba) = @_;
  my $sa = $dba->get_SliceAdaptor;

  my @slices = @{$sa->fetch_all('toplevel',undef,0,0)};
  my $chrom_seq = '';
  my $chrom_size = 0;
  my $add_N = 'N' x $self->param('N_padding_size');
  while ($chrom_size <= ($self->param('max_chromosome_size')-10000000)){
    my $random_slice = $slices[int rand@slices];
    my $seq = $random_slice->seq();
    my $num = length($seq);

    if ($num > 10000000){
      my $range = $num - 10000000;
      my $start = rand($range);
      my $piece = substr($seq, $start, ($start+10000000));
      $chrom_seq .= $add_N . $piece;
    }
    else{
        $chrom_seq .= $add_N . $seq;
    }
    $chrom_size = length($chrom_seq);
  }
  $self->param('chromosome_seqs', $chrom_seq);

}

sub write_output {
  my ($self) = @_;

  if ( !-d $self->param('tmp_dir') ) {
      make_path $self->param('tmp_dir') or die "Failed to create path: ".$self->param('tmp_dir');
    }
  my @output_ids;
  foreach my $chrom_seq ($self->param('chromosome_seqs')){
    my @chars = ("A".."Z", "a".."z", 0..9);
    my $string;
    $string .= $chars[rand @chars] for 1..8;

    my $filename = "chromosome.".$string.".fa";
    my $outfile = $self->param('tmp_dir')."/".$filename;
    my $chrom_name = (split '/', $outfile)[-1];
    $chrom_name =~ s{\.fa}{};

    open(OUT, ">", $outfile);
    say OUT ">".$chrom_name;
    say OUT $chrom_seq;
    push(@output_ids, {iid => $filename});
  }
    $self->dataflow_output_id(\@output_ids, 1);

}

1;
