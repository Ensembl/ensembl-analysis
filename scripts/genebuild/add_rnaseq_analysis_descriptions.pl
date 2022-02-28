#!/usr/bin/env perl
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

use strict;
use warnings;

use HTTP::Tiny;
use JSON;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use feature 'say';

my ($help, $safe_mode, $dbname, $port, $host);
my $use_datafile = 0;
my $update_analysis_description = 0;

my $dbuser = 'ensro';
my $user = $ENV{USER};
GetOptions(
	   'help|h'       => \$help,
	   'safe_mode'    => \$safe_mode,
	   'dbname=s'     => \$dbname,
	   'port=s'       => \$port,
	   'host=s'       => \$host,
     'dbuser=s'     => \$dbuser,
     'user=s'       => \$user,
     'df|Use_data_file!' => \$use_datafile,
     'update!'      => \$update_analysis_description,
);

die &helptext if ( $help );

say "Adding analysis decsriptions for logic name in database: ".$dbname;
get_content($dbname, $port, $host, $use_datafile, $update_analysis_description);


=head2 get_content

 Arg [1]    : String database name
 Arg [2]    : Int server port
 Arg [3]    : String server name
 Description: Connect to the database to retrieve all rnaseq analyses, then check
              that the logic_name does exist. If it does, no changes are needed
              otherwise, add the web data and analysis description
 Returntype : String data added to the production database
 Exceptions : Throw if it fails to load the data

=cut

sub get_content {
  my ($dbname, $port, $host, $use_datafile, $update_analysis_description) = @_;

  my $json = JSON->new;
  if ($safe_mode) {
    $json->pretty([1]);
  }
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -port    => $port,
      -user    => $dbuser,
      -host    => $host,
      -dbname  => $dbname);

  my $sth_species = $db->dbc->prepare("select meta_value from meta where meta_key='species.scientific_name';");
  $sth_species->execute();
  my $species_name = lc $sth_species->fetchrow;
  $species_name =~ s/ /_/g;

  my $sth_logic = $db->dbc->prepare("select logic_name from analysis");
  $sth_logic->execute;
  my $http_client = HTTP::Tiny->new;
  my $server = 'http://production-services.ensembl.org';
  my $ext = '/api/production_db/analysisdescription';
  my @rnaseq_logic_names = ();
  while (my $logic_name = $sth_logic->fetchrow) {
    if($logic_name =~ /\_rnaseq_gene$/ || $logic_name =~ /\_rnaseq_bam$/ || $logic_name =~ /\_rnaseq_daf$/ || $logic_name =~ /\_rnaseq_ise$/ || $logic_name =~ /_isoseq$/) {
      my $response = $http_client->get("$server$ext/$logic_name");
      if (!$update_analysis_description and $response->{success}) {
        say "$logic_name already exists" if ($safe_mode);
      }
      else {
        push(@rnaseq_logic_names,$logic_name);
      }
    }
  }# end while
  foreach my $rnaseq_logic_name (@rnaseq_logic_names){
    my $logic_type = $rnaseq_logic_name;
    if ($logic_type =~ /(rnaseq_.*)/) {
      $logic_type = $1;
    }
    elsif ($logic_type =~ /_isoseq$/) {
      $logic_type = 'isoseq';
    }

    my (undef, $sample_name) = $rnaseq_logic_name =~ /(${species_name}_)?(\w+)_(rna|iso)seq.*/;
    if ($sample_name eq $rnaseq_logic_name) {
      die("Could not retrieve the sample name from $rnaseq_logic_name");
    }
    my $provider = 'ENA';
    if ($use_datafile and $rnaseq_logic_name =~ /_rnaseq_/) {
      my $analysis = $rnaseq_logic_name;
      $analysis =~ s/_[^_]+$/_bam/;
      my $data_files = $db->get_DataFileAdaptor->fetch_all_by_logic_name($analysis);
      if (@$data_files != 1) {
        warn('You have '.scalar(@$data_files).' this should not happen, I will use the default provider');
      }
      else {
        ($provider) = $data_files->[0]->name =~ /\.([^.]+)\.$sample_name/;
      }
    }
    $sample_name =~ s/_+/ /g;

    my $values_dict = get_values($sample_name, $logic_type);

    my %json_data = (
      user => $user,
      logic_name => $rnaseq_logic_name,
      description => $values_dict->{description},
      display_label => $values_dict->{display_label},
    );
    if ($logic_type eq 'isoseq') {
      $json_data{web_data} = {
        data => {
          zmenu => $values_dict->{web_data_zmenu},
          label_key => $values_dict->{web_data_label_key},
          colour_key => $values_dict->{web_data_colour_key},
          type => 'longreads',
        },
      };
    }
    elsif ($logic_type ne "rnaseq_ise"){
      $json_data{web_data} = {
        data => {
          zmenu => $values_dict->{web_data_zmenu},
          label_key => $values_dict->{web_data_label_key},
          colour_key => $values_dict->{web_data_colour_key},
          matrix => {
            column => $values_dict->{matrix_column},
            menu => 'rnaseq',
            group => $provider,
            row => ucfirst($sample_name),
            group_order => $values_dict->{matrix_group_order},
          },
          type => 'rnaseq',
        },
      };
      if ($logic_type eq "rnaseq_daf") {
        $json_data{web_data}{data}{additional_renderers} = $values_dict->{web_data_additional_renders},
      }
    }

    if ( $safe_mode ) {
      print "RUN IN SAFE MODE\nCONTENT:\n".$sample_name." ".$logic_type."\n".$json->encode(\%json_data)."\n";
    }
    else {
      my $http_action = 'POST';
      my $error_msg = 'add';
      my $url = $server.$ext;
      if ($update_analysis_description) {
        $http_action = 'PATCH';
        $error_msg = 'update';
        $url .= "/$rnaseq_logic_name"
      }
      my $response = $http_client->request($http_action, $url, {
          headers => {
          'Content-type' => 'application/json',
          'Accept' => 'application/json'
        },
        content => $json->encode(\%json_data)
      });

      print "\nFailed to $error_msg analyses descriptions for $sample_name!\n  ".$response->{status}."\n  ".$response->{content}."\n"
        unless $response->{success};
    }

  }# end foreach logic_name


}# end get_content


sub get_values {
  my ($sample_name, $logic_type) = @_;

  my %description = (
		     'rnaseq_gene' => "Annotation generated from ".$sample_name." RNA-seq data",
		     'rnaseq_bam'  => 'Alignments of '.$sample_name.' RNA-seq data. This BAM file can be downloaded from the <a href="ftp://ftp.ensembl.org/pub/data_files/">Ensembl FTP site</a>',
		     'rnaseq_ise'  => "Spliced-read support for ".$sample_name,
		     'rnaseq_daf'  => "Spliced-read support for ".$sample_name,
         'isoseq'      => ucfirst($sample_name).' PacBio long reads from <a rel="external" href="https://www.ebi.ac.uk/ena">ENA</a> aligned to the genome using <a rel="external" href="https://doi.org/10.1093/bioinformatics/bty191">Minimap2</a>',
		    );

  my %display_label = (
    rnaseq_gene => ucfirst($sample_name)." RNA-seq gene models",
    rnaseq_bam  => ucfirst($sample_name)." RNA-seq alignments",
    rnaseq_ise  => ucfirst($sample_name)." intron-spanning reads",
    rnaseq_daf  => ucfirst($sample_name)." intron-spanning reads",
    isoseq => ucfirst($sample_name).' PacBio lond reads',
        	    );

  my %web_data_zmenu = (
		      'rnaseq_gene' => "RNASeq",
		      'rnaseq_bam'  => "RNASeq_bam",
	              'rnaseq_daf'  => "Supporting_alignment",
          isoseq => 'CLS',
		    );

  my %web_data_label_key = (
		       'rnaseq_gene' => "RNASeq [biotype]",
		       'rnaseq_bam'  => "RNASeq [biotype]",
		       'rnaseq_daf'  => "",
           isoseq => '[biotype] [display_label]',
		    );

  my %web_data_additional_renders = (
		       'rnaseq_gene' => "",
	       	       'rnaseq_bam'  => "",
       'rnaseq_daf'  => ["histogram", "Variable height"],
		    );

  my %web_data_colour_key = (
		       'rnaseq_gene' => "human_rnaseq",
		       'rnaseq_bam'  => "bam",
		       'rnaseq_daf'  => "intron_support",
           isoseq => 'long_reads_isoseq',
		    );

  my %matrix_group_order = (
      rnaseq_bam  => 1,
      rnaseq_daf  => 2,
      rnaseq_gene => 3,
      );

  my %matrix_column = (
    rnaseq_gene => 'Gene models',
    rnaseq_bam  => 'BAM files',
    rnaseq_daf  => 'Intron-spanning reads',
  );

  return {
    description => $description{$logic_type},
    display_label => $display_label{$logic_type},
    web_data_zmenu => $web_data_zmenu{$logic_type},
    web_data_label_key => $web_data_label_key{$logic_type},
    web_data_additional_renders => $web_data_additional_renders{$logic_type},
    web_data_colour_key => $web_data_colour_key{$logic_type},
    matrix_group_order => $matrix_group_order{$logic_type},
    matrix_column => $matrix_column{$logic_type},
  };

}# end get_values

sub helptext {
  my $msg = <<HELPEND;

IMPORTANT: it is strongly recommended that you run this in SAFE MODE and check the content before you run it as normal, i.e. before you POST any content

Usage: perl add_rnaseq_analysis_descriptions.pl -dbname <rnaseq_dbname> -host <host> -port <port>

Options: -safe_mode -> run the script without POSTing to the production database, i.e. print the content that would be POSTed when not run in safe mode
         -user <production services username>, default to your Unix user name

HELPEND
  return $msg;
}

