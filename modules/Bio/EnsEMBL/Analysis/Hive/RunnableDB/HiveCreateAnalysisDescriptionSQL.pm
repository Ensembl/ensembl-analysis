=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateAnalysisDescriptionSQL;

use warnings;
use strict;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use HTTP::Tiny;
use Time::HiRes qw/sleep/;
use JSON qw/decode_json/;
use Data::Dumper;

=head2 fetch_input

    Description : Implements fetch_input() interface method of Bio::EnsEMBL::Hive::Process that is used to read in parameters and load data.
                  Here we have nothing to do.

=cut

sub fetch_input {
}


=head2 run

 Arg [1]    : None
 Description: Creates an sql command to update the analysis_description table with info from the ensembl_production_db
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my @output_ids;
  my $logic_name = $self->param('logic_name');
  my $input_dba = $self->hrdb_get_dba($self->param('input_db'));

  my @core_logic_names = ("ensembl","ncrna","blastmirna","rfamcmsearch","rnaseq_intron_support","cdna2genome","projected_transcript","other_protein","genscan","trnascan","cpg","eponine","trf","dust","repeatmask_repeatmodeler","repeatmask_repbase_mammals","repeatmask_repbase_aves","repeatmask_repbase_rodents","repeatmask_repbase_primates","repeatmask_repbase_teleost","repeatmask_repbase_mouse","repeatmask_repbase_human","repeatmask_repbase_human_low","repeatmask_repbase_zebrafish","repeatmask_repbase_chicken","repeatmask_repbase_sus_scrofa");

  if ( grep( /^$logic_name$/, @core_logic_names ) ) {

    my $sth = $input_dba->dbc->prepare('SELECT analysis_id from analysis where logic_name="'.$logic_name.'"');
    $sth->execute();
    my $analysis_id = $sth->fetchrow_array();

    my $http = HTTP::Tiny->new();
    my $server = 'http://production-services.ensembl.org';
    my $ext = '/api/production_db/analysisdescription/';
    my $response = $http->request('GET', $server.$ext.$logic_name, {
			 headers => {
				     'Content-type' => 'application/json',
				    },
								   });

    if ($response->{success}){
      my $hash_ref = decode_json($response->{content});
      my %hash = %$hash_ref;

      local $Data::Dumper::Terse = 1;
      local $Data::Dumper::Indent = 0;
      my $web_data = Dumper($hash{'web_data'}->{data});
      if ($web_data eq 'undef') {
	$web_data = "NULL";
      }
      else {
	$web_data = '"'.$web_data.'"';
      }
# convert the web_data back to json format
      $web_data =~ s/\=\>/:/g;
      my $desc = $hash{'description'};
      $desc =~ s/\'/\\\'/g;

      my $sql = "INSERT INTO analysis_description (analysis_id, description, display_label, displayable, web_data) VALUES ($analysis_id, '$desc', '$hash{'display_label'}', $hash{'displayable'}, $web_data);";
      my $output_hash = {};
      $output_hash->{'sql_command'} = $sql;
      $self->dataflow_output_id($output_hash,2);
    }
  }

  else{
    my $fail_output_hash = {};
    $fail_output_hash->{'logic_name'} = $logic_name;
    $self->dataflow_output_id($fail_output_hash,-2);
    $self->input_job->autoflow(0);
  }
}

sub write_output {
  my ($self) = @_;
}


1;
