# Copyright [2022] EMBL-European Bioinformatics Institute
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

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $coord_system = 'toplevel';
my $dbname = '';
my $user = '';
my $host = '';
my $port;
my $pass;
my $output_file_x = '';
my $output_file_y = '';

my $options = GetOptions ("user|dbuser|u=s" => \$user,
                          "host|dbhost|h=s" => \$host,
                          "port|dbport|P=i" => \$port,
                          "dbname|db|D=s"   => \$dbname,
                          "dbpass|pass|p=s" => \$pass,
                          "output_file_x=s" => $output_file_x,
                          "output_file_y=s" => $output_file_y);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all('toplevel');
my $gene_adaptor = $db->get_GeneAdaptor();
my $meta_adaptor = $db->get_MetaContainer();
my $production_name = $meta_adaptor->get_production_name;

my $x_slice =  $slice_adaptor->fetch_by_region('toplevel','X');
my $y_slice = $slice_adaptor->fetch_by_region('toplevel','Y');

my $x_seqs = get_slice_seqs($x_slice,2500,300,10000,10000000,'x');
my $y_seqs = get_slice_seqs($y_slice,2500,300,2500,10000000,'y');

open(OUT,">".$output_file_x);
print OUT $x_seqs;
close OUT;

open(OUT,">".$output_file_y);
print OUT $y_seqs;
close OUT;


sub get_slice_seqs {
  my ($slice,$num_regions,$region_length,$interval_length,$start_offset,$header_label) = @_;

  my $record = "";
  my $current_end = $start_offset;
  my $end_offset = 5000000;
  for(my $i=1; $i<=$num_regions && $current_end < ($slice->length()-$end_offset); $i++) {
    my $current_start = ($i * $interval_length) + $start_offset;
    my $current_end = $current_start + $region_length;
    say "Region: ".$slice->seq_region_name.":".$current_start.":".$current_end;
    my $header = ">".$header_label."_".($i-1);
    my $seq = $slice->subseq($current_start,$current_end);
    if($seq =~ /N/) {
      $num_regions++;
      next;
    }
    $record .= $header."\n".$seq."\n";
  }
  return($record);
}
