=head1 LICENSE

# Copyright [2019] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License,Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

create_apollo_adapters.pl

=head1 DESCRIPTION

  A script that builds an adapter config file for apollo based on a list of core dbs with server shorthands
  An example of the input file format would be:
  prim96_macaca_fuscata_fuscata_core_96@gb6:gb5
  prim96_pan_troglodytes_core_96@gb6:gb5
  prim96_piliocolobus_tephrosceles_core_96@gb6:gb5
  prim96_pongo_abelii_core_96@gb6:gb5
  prim96_prolemur_simus_core_96@gb6:gb5
  prim96_theropithecus_gelada_core_96@gb6:gb5

  In the above, gb6 is shorthand for the core server, while gb5 is shorthand for the output db server
  These is an array of output db names that the script will automatically search for. If these dbs are
  found on the output server then their info will be added into the adapter file. This script can be
  run standalone, but the expected usuage is to be automatically run every half hour on a shared access
  file that lists the dbs so that we can have a centralised, up to date config that everyone can use

  NOTE: For the servers, you must have env vars set for the servers and ports, i.e. $GBS1/$GBP1 for server 1 and the corresponding port

=head1 OPTIONS

  -database_file_path   The path to the file listing the core dbs and the servers
  -output_dir           The directory to write the output file to
  -output_file_name     The name of the adapter config to be create (defaults at adapterbuilder.conf, which is what Apollo looks for)

=head1 EXAMPLES

  perl create_apollo_adapters.pl -database_file_path /path/to/db_file.txt -output_dir /path/to/output/dir/

=cut


use warnings;
use strict;
use feature 'say';
use Getopt::Long;

my $database_file_path;
my $output_dir;
my $output_file_name = 'adapterbuilder.conf';

my @output_types = ('genblast',
                    'genblast_rnaseq',
                    'final',
                    'gbuild',
                    'proj_cesar',
                    'proj_coding',
                    'pseudo',
                    'layer',
                    'refine',
                    'refseq',
                    'rnaseq_blast',
                    'rnaseq_layer');

GetOptions('database_file_path:s' => \$database_file_path,
           'output_dir:s'         => \$output_dir,
           'output_file_name:s'   => \$output_file_name);

unless($database_file_path && -e $database_file_path) {
  die "The database file path was not valid:\n".$database_file_path;
}

unless($output_dir) {
  die "No output dir was specified";
}

my $hosts = {
  'gb1' => [$ENV{GBS1},$ENV{GBP1}],
  'gb2' => [$ENV{GBS2},$ENV{GBP2}],
  'gb3' => [$ENV{GBS3},$ENV{GBP3}],
  'gb4' => [$ENV{GBS4},$ENV{GBP4}],
  'gb5' => [$ENV{GBS5},$ENV{GBP5}],
  'gb6' => [$ENV{GBS6},$ENV{GBP6}],
  'gb7' => [$ENV{GBS7},$ENV{GBP7}],
};

my $default_set;
my $default_node;

my $inner_template = "";

open(IN,$database_file_path);
my @cores = <IN>;
close IN;

my $first = 1;
foreach my $core (@cores) {
  chomp $core;

  unless($core =~ /(.+)\@(.+)\:(.+)/) {
    die 'Core db didn\'t match expected format of mydatabase@host:port'."\n".$core;
  }

  my $core_name = $1;
  my $core_connection_details = $hosts->{$2};
  my $out_connection_details = $hosts->{$3};
  unless($core_connection_details && $out_connection_details) {
    die 'Could not find connection details based on the following core string: '.$core;
  }

  my $core_host = ${$core_connection_details}[0];
  $core_host =~ s/\.ebi\.ac\.uk//;

  my $core_port = ${$core_connection_details}[1];

  $core_name =~ /(.+)\_core\_/;
  my $prefix = $1;

  if($first == 1) {
    $default_set = $prefix;
    $default_node = $core_name;
    $first--;
  }

  my $out_host = ${$out_connection_details}[0];
  $out_host =~ s/\.ebi\.ac\.uk//;

  my $out_port = ${$out_connection_details}[1];


  my $seq_region_info = `mysql -uensro -h$core_host -P$core_port $core_name -NB -e 'select name,length from seq_region order by length desc limit 1'`;
  chomp($seq_region_info);

  my ($seq_region_name,$seq_region_length) = split("\t",$seq_region_info);

  unless($seq_region_name && $seq_region_length) {
    die "Failed to find seq region name/length for ".$core_name;
  }

  my $assembly_name = `mysql -uensro -h$core_host -P$core_port $core_name -NB -e 'select meta_value from meta where meta_key="assembly.default"'`;
  chomp($assembly_name);

  my $region = $assembly_name.'--primary_assembly:'.$seq_region_name.':1-'.$seq_region_length;

  my $core_template =  make_adaptor_template($core_name,'MAIN_TYPE',$core_port,'localhost',$region,'core');
  my $output_templates = "";

  foreach my $output_type (@output_types) {
    my $out_name = $core_name;
    $out_name =~ s/\_core\_/\_$output_type\_/;

    my $check_exists = `mysqlshow --user=ensro -h$out_host -P$out_port $out_name | grep -v Wildcard | grep -o $out_name`;
    if($check_exists) {
      $output_templates .= make_adaptor_template($out_name,'CHILD_TYPE',$out_port,'localhost',$region,$output_type);
    }
  }

  my $combined_template = $core_template."\n".$output_templates;
  my $set_template = '<adapter_set>'."\n".
                     '<name>'.$prefix.'</name>'."\n".
                     $combined_template.
                     '</adapter_set>';

  $inner_template .= $set_template;
}

my $outer_template = make_outer_template($default_set,$default_node);
$outer_template =~ s/SUB_IN_MAIN_TEMPLATE/$inner_template/;

unless(open(OUT, ">".$output_dir."/".$output_file_name)) {
  die "Could not open the output file for writing. Path provided:\n".$output_dir."/".$output_file_name;
}

print OUT $outer_template;
close OUT;

exit;

sub make_outer_template {
my ($default_set,$default_node) = @_;
my $outer_template = '<model>'."\n".
        '<adapter_classes>'."\n".
                '<class>apollo.dataadapter.GFFAdapter</class>'."\n".
                '<class>apollo.dataadapter.BAMAdapter</class>'."\n".
                '<class>apollo.dataadapter.ensj.EnsJAdapter</class>'."\n".
                '<class>apollo.dataadapter.das.simple.SimpleDASAdapter</class>'."\n".
        '</adapter_classes>'."\n".
        '<adapter_types>'."\n".
                '<type>MAIN_TYPE</type>'."\n".
                '<type>CHILD_TYPE</type>'."\n".
        '</adapter_types>'."\n".
        '<selected_set>'.$default_set.'</selected_set>'."\n".
        '<selected_node>'.$default_node.'</selected_node>'."\n".
'SUB_IN_MAIN_TEMPLATE'."\n".
'</model>'."\n";

return($outer_template);
}

sub make_adaptor_template {

my ($db_name,$adaptor_type,$port,$host,$region,$type) = @_;

my $adaptor_template ='<adapter>'."\n".
                      '<name>'.$db_name.'</name>'."\n".
                        '<class>apollo.dataadapter.ensj.EnsJAdapter</class>'."\n".
                        '<type>'.$adaptor_type.'</type>'."\n".
                        '<property>'."\n".
                                '<key>INCLUDE_DNA_DNA_ALIGNMENT_TYPES0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>password</key>'."\n".
                                '<value></value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.Gene</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.RepeatFeature</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.Feature</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>addSupport.Transcript</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>region</key>'."\n".
                                '<value>REGION: '.$region.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.PredictionTranscript</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.ContigFeature</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>INCLUDE_DITAG_FEATURE_TYPES0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>INCLUDE_PREDICTION_TRANSCRIPT_TYPES0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>INCLUDE_DNA_PROTEIN_ALIGNMENT_TYPES0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>resetStartAndStop.Gene</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>locationHistory1</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>port</key>'."\n".
                                '<value>'.$port.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>aggressiveNaming.Gene</key>'."\n".
                                '<value>true</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>locationHistory0</key>'."\n".
                                '<value>'.$region.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.SimplePeptideFeature</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>loggingFile</key>'."\n".
                                '<value>/Users/genebuild/Desktop/Apollo.app/Contents/MacOS/conf/logging_info_level.conf</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>database</key>'."\n".
                                '<value>'.$db_name.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>host</key>'."\n".
                                '<value>'.$host.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>addResultAsAnnotation.Gene</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.DitagFeature</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>INCLUDE_GENE_TYPES0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>typePrefix.Gene</key>'."\n".
                                '<value>'.$type.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.DnaProteinAlignment</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>stableIDHistory0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>include.DnaDnaAlignment</key>'."\n".
                                '<value>false</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>INCLUDE_FEATURE0</key>'."\n".
                                '<value>END_OF_LIST</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>name</key>'."\n".
                                '<value>'.$db_name.'</value>'."\n".
                        '</property>'."\n".
                        '<property>'."\n".
                                '<key>user</key>'."\n".
                                '<value>ensro</value>'."\n".
                        '</property>'."\n".
                '</adapter>'."\n";
  return($adaptor_template);
}
