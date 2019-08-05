#!/usr/bin/env perl
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

use strict;
use warnings;

use HTTP::Tiny;
use Time::HiRes qw/sleep/;
use JSON qw/decode_json/;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use feature 'say';

my ($help, $safe_mode, $dbname, $port, $host);

GetOptions(
	   'help|h'       => \$help,
	   'safe_mode'    => \$safe_mode,
	   'dbname=s'     => \$dbname,
	   'port=s'       => \$port,
	   'host=s'       => \$host,
);

die &helptext if ( $help );

say "Adding analysis decsriptions for logic name in database: ".$dbname;
get_content($dbname, $port, $host);


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
  my ($dbname, $port, $host) = @_;
  my $content;

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -port    => $port,
      -user    => 'ensro',
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
  my $ext = '/production_db/api/analysisdescription';
  my @rnaseq_logic_names = ();
  while (my $logic_name = $sth_logic->fetchrow) {
    if($logic_name =~ /\_rnaseq_gene$/ || $logic_name =~ /\_rnaseq_bam$/ || $logic_name =~ /\_rnaseq_daf$/ || $logic_name =~ /\_rnaseq_ise$/) {
      my $response = $http_client->get("$server$ext/$logic_name");
      push(@rnaseq_logic_names,$logic_name) unless ($response->{success});
    }
  }# end while
  foreach my $rnaseq_logic_name (@rnaseq_logic_names){
    my $logic_type = $rnaseq_logic_name;
    $logic_type = $1 if $logic_type =~ /(rnaseq_.*)/;

    my $sample_name = $rnaseq_logic_name;
    $sample_name =~ s/$species_name// && $sample_name =~ s/_rnaseq_.*// && $sample_name =~ s/_// && $sample_name =~ s/_/ /g;

    my %values_dict = get_values($sample_name, $logic_type);

    if ($logic_type eq "rnaseq_ise"){
      $content = "{
                       \"logic_name\": \"".$rnaseq_logic_name."\",
                       \"description\": \"".$values_dict{'description'}."\",
                       \"display_label\": \"".$values_dict{'display_label'}."\"
                     }";
    }
    else{
      $content = "{
                      \"web_data\": {
                           \"data\": {
                               \"zmenu\": \"".$values_dict{'web_data_zmenu'}."\",
                               \"label_key\": \"".$values_dict{'web_data_label_key'}."\",
                               \"colour_key\": \"".$values_dict{'web_data_colour_key'}."\",
                               \"type\": \"rnaseq\",";
      if ($logic_type eq "rnaseq_daf") {
	$content .=           "
                               \"additional_renderers\": ".$values_dict{'web_data_additional_renders'}.""
      }
      $content .=              "
                               \"matrix\": {";
      if ($logic_type eq "rnaseq_bam") {
	$content .=                         "
                                           \"group_order\": \"".$values_dict{'matrix_group_order'}."\","
      }
      $content .=                           "
                                           \"column\": \"".$values_dict{'matrix_column'}."\",
                                           \"menu\": \"rnaseq\",
                                           \"group\": \"ENA\",
                                           \"row\": \"".$sample_name."\"
                                           }
                                     }
                                  },
                      \"logic_name\": \"".$rnaseq_logic_name."\",
                      \"description\": \"".$values_dict{'description'}."\",
                      \"display_label\": \"".$values_dict{'display_label'}."\"
                     }";
    }

    if ( $safe_mode ) {
      print "RUN IN SAFE MODE\nCONTENT:\n".$sample_name." ".$logic_type."\n".$content."\n";
    }
    else {
      my $http = HTTP::Tiny->new;
      my $response = $http->request('POST', $server.$ext, {
          headers => {
		      'Content-type' => 'application/json',
		      'Accept' => 'application/json'
		     },
          content => $content
							  });

      print "\nFailed to add analyses descriptions for ".$sample_name."!\n" unless $response->{success};
    }

  }# end foreach logic_name

  return $content;

}# end get_content


sub get_values {
  my ($sample_name, $logic_type) = @_;

  my %description = (
		     'rnaseq_gene' => "Annotation generated from ".$sample_name." RNA-seq data",
		     'rnaseq_bam'  => 'BWA alignments of '.$sample_name.' RNA-seq data. This BAM file can be downloaded from the <a href=\\"ftp://ftp.ensembl.org/pub/data_files/\\">Ensembl FTP site</a>',
		     'rnaseq_ise'  => "Spliced-read support for ".$sample_name,
		     'rnaseq_daf'  => "Spliced-read support for ".$sample_name,
		    );

  my %display_label = (
         	     'rnaseq_gene' => $sample_name." RNA-seq gene models",
		     'rnaseq_bam'  => $sample_name." RNA-seq BWA alignments",
		     'rnaseq_ise'  => $sample_name." intron-spanning reads",
	             'rnaseq_daf'  => $sample_name." intron-spanning reads",
        	    );

  my %web_data_zmenu = (
		      'rnaseq_gene' => "RNASeq",
		      'rnaseq_bam'  => "RNASeq_bam",
	              'rnaseq_daf'  => "Supporting_alignment",
		    );

  my %web_data_label_key = (
		       'rnaseq_gene' => "RNASeq [biotype]",
		       'rnaseq_bam'  => "RNASeq [biotype]",
		       'rnaseq_daf'  => "",
		    );

  my %web_data_additional_renders = (
		       'rnaseq_gene' => "",
	       	       'rnaseq_bam'  => "",
		       'rnaseq_daf'  => "\[\"histogram\", \"Variable height\"\],",
		    );

  my %web_data_colour_key = (
		       'rnaseq_gene' => "human_rnaseq",
		       'rnaseq_bam'  => "bam",
		       'rnaseq_daf'  => "intron_support",
		    );

  my %matrix_group_order = (
		        'rnaseq_gene' => "",
			'rnaseq_bam'  => "1",
			'rnaseq_daf'  => "",
		    );

  my %matrix_column = (
			'rnaseq_gene' => "Gene models",
			'rnaseq_bam'  => "BAM files",
			'rnaseq_daf'  => "Intron-spanning reads",
		   );

  my %desc_info = (
		        'description' => $description{$logic_type},
		        'display_label' => $display_label{$logic_type},
		        'web_data_zmenu' => $web_data_zmenu{$logic_type},
		        'web_data_label_key' => $web_data_label_key{$logic_type},
		        'web_data_additional_renders' => $web_data_additional_renders{$logic_type},
		        'web_data_colour_key' => $web_data_colour_key{$logic_type},
		        'matrix_group_order' => $matrix_group_order{$logic_type},
		        'matrix_column' => $matrix_column{$logic_type},
                  );

  return %desc_info;

}# end get_values

sub helptext {
  my $msg = <<HELPEND;

IMPORTANT: it is strongly recommended that you run this in SAFE MODE and check the content before you run it as normal, i.e. before you POST any content

Usage: perl add_rnaseq_analysis_descriptions.pl -dbname <rnaseq_dbname> -host <host> -port <port>

Options: -safe_mode -> run the script without POSTing to the production database, i.e. print the content that would be POSTed when not run in safe mode

HELPEND
  return $msg;
}

