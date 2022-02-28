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

use warnings;
use strict;
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long qw(:config no_ignore_case);

my $blastdbname;
my $blasthost;
my $blastport;
my $blastpass;
my $blastuser;

my $refinedbname;
my $refinehost;
my $refineport;
my $refineuser;

my $species;
my $output_path;
my $uniprot_version;

my $result = GetOptions ("blastdbuser=s"   => \$blastuser,
                         "blastdbhost=s"   => \$blasthost,
                         "blastdbport=i"   => \$blastport,
                         "blastdbname=s"   => \$blastdbname,
                         "blastdbpass=s"   => \$blastpass,
                         "refinedbuser=s"  => \$refineuser,
                         "refinedbhost=s"  => \$refinehost,
                         "refinedbport=i"  => \$refineport,
                         "refinedbname=s"  => \$refinedbname,
                         "output_path=s"   => \$output_path,
                         "species=s"       => \$species,
                         "uniprot_version" => \$uniprot_version);

unless($blastuser && $blasthost && $blastport && $blastpass && $species && $refineuser && $refinehost && $refineport) {
  die "You must pass in the following connection info: blastdbuser, blastdbhost, blastdbport,".
      " blastdbname, blastdbpass, refinedbuser, refinedbhost, refinedbport, refinedbname, species";
}

my $blast_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $blastport,
  -user    => $blastuser,
  -host    => $blasthost,
  -dbname  => $blastdbname,
  -pass    => $blastpass);

my $refine_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $refineport,
  -user    => $refineuser,
  -host    => $refinehost,
  -dbname  => $refinedbname);

#unless(-e $output_path) {
#  system("mkdir -p ".$output_path)
#}


# Truncate the dna align feature table and remove any unlinked supporting features
my $sth_clean = $blast_db->dbc->prepare('TRUNCATE dna_align_feature');
$sth_clean->execute();

$sth_clean = $blast_db->dbc->prepare('DELETE supporting_feature from supporting_feature left join dna_align_feature on feature_id=dna_align_feature_id'.
                                     ' where feature_type = "dna_align_feature" and dna_align_feature_id is NULL');
$sth_clean->execute();

$sth_clean = $blast_db->dbc->prepare('DELETE transcript_supporting_feature from transcript_supporting_feature left join dna_align_feature on feature_id=dna_align_feature_id'.
                                     ' where feature_type = "dna_align_feature" and dna_align_feature_id is NULL');
$sth_clean->execute();

# This code is very slow but will work when the tables are out of sync between dbs
# There are a bunch of ways this could be done, including parallelising it outside
# this script, but I think that the pipeline should just be fixed and then this code
# deleted and the code commented out below it should be activated
#my $refine_daf_adaptor = $refine_db->get_DnaAlignFeatureAdaptor();
#my $blast_daf_adaptor = $blast_db->get_DnaAlignFeatureAdaptor();
#my $refine_slice_adaptor = $refine_db->get_SliceAdaptor();
#my $slices = $refine_slice_adaptor->fetch_all('toplevel');
#say "Transferring DAFs from refine db to blast db...";
#foreach my $slice (@$slices) {
#  my $refine_dafs = $refine_daf_adaptor->fetch_all_by_Slice($slice);
#  foreach my $refine_daf (@$refine_dafs) {
#    $blast_daf_adaptor->store($refine_daf);
#  }
#}
#say "Completed transfer of DAFs from refine db to blast db";

# Note this code cannot be used currently as the refine and blast tables seem to be
# out of sync with one another, if the pipeline is fixed to do this then reactivate it
#my $dump_daf_table = "mysqldump -u".$refineuser.
#                              " -h".$refinehost.
#                              " -P".$refineport.
#                              " ".$refinedbname." dna_align_feature > ".$output_path."/".$species."_refine_daf.sql";

#my $cmd_result = system($dump_daf_table);
#if($cmd_result) {
#  die "Failed to dump refine db dna_align_feature table. Command line used:\n".$dump_daf_table;
#}

#my $load_daf_table = "mysql -u".$blastuser.
#                          " -h".$blasthost.
#                          " -p".$blastpass.
#                          " -P".$blastport.
#                          " -D".$blastdbname.
#                          " < ".$output_path."/".$species."_refine_daf.sql";

#$cmd_result = system($load_daf_table);
#if($cmd_result) {
#  die "Failed to laod refine db dna_align_feature table. Command line used:\n".$load_daf_table;
#}


# Fix the protein align feature table
$uniprot_version =~ s/\///g;
my $sth_paf = $blast_db->dbc->prepare('INSERT IGNORE into analysis (logic_name,db_version) values ("other_protein",?)');
$sth_paf->bind_param(1,$uniprot_version);
$sth_paf->execute();

$sth_paf = $blast_db->dbc->prepare('SELECT analysis_id from analysis where logic_name="other_protein"');
$sth_paf->execute();
my $paf_analysis_id = $sth_paf->fetchrow_array;

unless($paf_analysis_id) {
  die "Problem retrieving analysis id for other_protein logic name. This is needed to update the PAF table";
}

$sth_paf = $blast_db->dbc->prepare('UPDATE protein_align_feature set analysis_id=?');
$sth_paf->bind_param(1,$paf_analysis_id);
$sth_paf->execute();

$sth_paf = $blast_db->dbc->prepare('SELECT external_db_id from external_db where db_name="UniProtKB_all"');
$sth_paf->execute();
my $paf_external_db_id = $sth_paf->fetchrow_array;

unless($paf_external_db_id) {
  die "Problem retrieving external_db_id for UniProtKB_all. This is needed to update the PAF table";
}

$sth_paf = $blast_db->dbc->prepare('UPDATE protein_align_feature set external_db_id=?');
$sth_paf->bind_param(1,$paf_external_db_id);
$sth_paf->execute();


# Retrieve all the rnaseq logic names
my $sth_select = $blast_db->dbc->prepare('SELECT analysis_id,logic_name from analysis where logic_name like "'.$species.'_%_rnaseq"');
$sth_select->execute();

my %initial_logic_names = ();

while(my ($analysis_id,$logic_name) =  $sth_select->fetchrow_array()) {
  say "Logic name: ".$logic_name." (".$analysis_id.")";
  $initial_logic_names{$logic_name} = $analysis_id;
}

unless(scalar(keys(%initial_logic_names))) {
  die "Could not find any logic names of the format: ".$species.'_%_rnaseq';
}

# Loop through the main logic names, create and insert the new ones and assign
foreach my $logic_name (keys(%initial_logic_names)) {
  my $analysis_id =  $initial_logic_names{$logic_name};
  my $sth_insert;
  my $sth_update;

  my $gene_logic_name = $logic_name."_gene";
  my $intron_supporting_evidence_logic_name = $logic_name."_ise";
  my $data_file_logic_name = $logic_name."_bam";
  my $dna_align_feature_logic_name = $logic_name."_daf";

  my $gene_analysis_id;
  my $intron_supporting_evidence_analysis_id;
  my $data_file_analysis_id;
  my $dna_align_feature_analysis_id;

  say "For ".$logic_name." created the following to be inserted and linked:\n".$gene_logic_name."\n".
                                                                               $intron_supporting_evidence_logic_name."\n".
                                                                               $data_file_logic_name."\n".
                                                                               $dna_align_feature_logic_name;

  $sth_insert = $blast_db->dbc->prepare('INSERT into analysis (logic_name) values (?)');
  $sth_insert->bind_param(1,$gene_logic_name);
  $sth_insert->execute();
  $sth_insert->bind_param(1,$intron_supporting_evidence_logic_name);
  $sth_insert->execute();
  $sth_insert->bind_param(1,$data_file_logic_name);
  $sth_insert->execute();
  $sth_insert->bind_param(1,$dna_align_feature_logic_name);
  $sth_insert->execute();

  $sth_select = $blast_db->dbc->prepare('SELECT analysis_id from analysis where logic_name = ?');

  $sth_select->bind_param(1,$gene_logic_name);
  $sth_select->execute();
  $gene_analysis_id = $sth_select->fetchrow_array;

  $sth_select->bind_param(1,$intron_supporting_evidence_logic_name);
  $sth_select->execute();
  $intron_supporting_evidence_analysis_id = $sth_select->fetchrow_array;

  $sth_select->bind_param(1,$data_file_logic_name);
  $sth_select->execute();
  $data_file_analysis_id = $sth_select->fetchrow_array;

  $sth_select->bind_param(1,$dna_align_feature_logic_name);
  $sth_select->execute();
  $dna_align_feature_analysis_id = $sth_select->fetchrow_array;

  unless($gene_analysis_id && $intron_supporting_evidence_analysis_id && $data_file_analysis_id && $dna_align_feature_analysis_id) {
    die "Issue with retieving the analysis_ids for one or several of: ".$gene_logic_name.", ".
                                                                        $intron_supporting_evidence_logic_name.", ".
                                                                        $data_file_logic_name.", ".
                                                                        $dna_align_feature_logic_name;
  }

  # Update gene and transcript tables
  $sth_update = $blast_db->dbc->prepare('UPDATE gene join analysis using(analysis_id) set'.
                                        ' gene.analysis_id=? where gene.analysis_id=?');
  $sth_update->bind_param(1,$gene_analysis_id);
  $sth_update->bind_param(2,$analysis_id);
  $sth_update->execute();
  $sth_update = $blast_db->dbc->prepare('UPDATE transcript join gene using(gene_id) set transcript.analysis_id=gene.analysis_id');
  $sth_update->execute();

  # Update intron supporting evidence table
  $sth_update = $blast_db->dbc->prepare('UPDATE intron_supporting_evidence join analysis using(analysis_id) set'.
                                        ' intron_supporting_evidence.analysis_id=? where intron_supporting_evidence.analysis_id=?');
  $sth_update->bind_param(1,$intron_supporting_evidence_analysis_id);
  $sth_update->bind_param(2,$analysis_id);
  $sth_update->execute();

  # The the momemt the pipeline does not set the data files in a convenient manner
  # I will put the following sub in to attempt to deal with the current state. When the
  # pipeline has been updated this bit should been removed
  update_data_file($blast_db,$species,$data_file_logic_name,$data_file_analysis_id);

  # This SQL does not do anything with the current pipeline, but should work when the pipeline is updated
  # to point the approriate rnaseq analyses to the data files
#  $sth_update = $blast_db->dbc->prepare('UPDATE data_file join analysis using(analysis_id) set'.
#                                  ' data_file.analysis_id=? where data_file.analysis_id=?');
#  $sth_update->bind_param(1,$data_file_analysis_id);
#  $sth_update->bind_param(2,$analysis_id);
#  $sth_update->execute();



  # This SQL does not do anything with the current pipeline, but should work when the pipeline is updated
  # to sync the analyses tables
#  $sth_update = $blast_db->dbc->prepare('UPDATE dna_align_feature join analysis using(analysis_id) set'.
#                                  ' dna_align_feature.analysis_id=? where dna_align_feature.analysis_id=?');
#  $sth_update->bind_param(1,$dna_align_feature_analysis_id);
#  $sth_update->bind_param(2,$analysis_id);
#  $sth_update->execute();


} # end foreach my $logic_name

# DELETE CODE SECTION: Delete code below when the pipeline is updated to assign correct logic names to data file and daf!!!!!
# This will fix the data file and daf tables
my $sth_fix_merged = $blast_db->dbc->prepare('SELECT name from data_file where name != "merged"');
$sth_fix_merged->execute();
my $example_name = $sth_fix_merged->fetchrow_array;
unless($example_name =~ /^(.+)\.[^\.]+\.1/) {
  die "Temp code for parsing example name to fixed merged data file name failed. Name: ".$example_name;
}
my $fixed_merged_name = $1.".merged.1";
$sth_fix_merged = $blast_db->dbc->prepare('UPDATE data_file set name = ? where name = "merged"');
$sth_fix_merged->bind_param(1,$fixed_merged_name);
$sth_fix_merged->execute();

update_daf($blast_db,$refine_db,\%initial_logic_names);
# END DELTE CODE SECTION

# Now do a final clean up of logic_names
my $sth_clean_logic_names = $blast_db->dbc->prepare('DELETE from analysis where logic_name not like "'.$species.'_%_rnaseq_gene"'.
                                                    ' and logic_name not like "'.$species.'_%_rnaseq_ise"'.
                                                    ' and logic_name not like "'.$species.'_%_rnaseq_daf"'.
                                                    ' and logic_name not like "'.$species.'_%_rnaseq_bam"'.
                                                    ' and logic_name != "other_protein"');
$sth_clean_logic_names->execute();


# Clean the meta table
my $sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key like "%\.level"');
$sth_update_meta->execute();
$sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key like "provider\.%"');
$sth_update_meta->execute();
$sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key like "assembly.web_accession%"');
$sth_update_meta->execute();
$sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key like "removed_evidence_flag\.%"');
$sth_update_meta->execute();
$sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key like "marker\.%"');
$sth_update_meta->execute();
$sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key = "repeat.analysis"');
$sth_update_meta->execute();
$sth_update_meta = $blast_db->dbc->prepare('DELETE from meta where meta_key in ("genebuild.method","genebuild.projection_source_db","genebuild.start_date")');
$sth_update_meta->execute();

# Fix the biotypes
my $sth_update_biotypes = $blast_db->dbc->prepare('UPDATE gene set biotype="protein_coding"');
$sth_update_biotypes->execute();
$sth_update_biotypes = $blast_db->dbc->prepare('UPDATE transcript set biotype="protein_coding"');
$sth_update_biotypes->execute();

# fix stable ids
my $sth_update_ids = $blast_db->dbc->prepare('UPDATE gene set stable_id=concat("RNASEQG",LPAD(gene_id, 11, 0))');
$sth_update_ids->execute();
$sth_update_ids = $blast_db->dbc->prepare('UPDATE transcript set stable_id=concat("RNASEQT",LPAD(transcript_id, 11, 0))');
$sth_update_ids->execute();
$sth_update_ids = $blast_db->dbc->prepare('UPDATE translation set stable_id=concat("RNASEQP",LPAD(translation_id, 11, 0))');
$sth_update_ids->execute();
$sth_update_ids = $blast_db->dbc->prepare('UPDATE exon set stable_id=concat("RNASEQE",LPAD(exon_id, 11, 0))');
$sth_update_ids->execute();

# Set canonicals
my $sth_update_canonicals = $blast_db->dbc->prepare('UPDATE gene join transcript using(gene_id) set canonical_transcript_id=transcript_id');
$sth_update_canonicals->execute();
$sth_update_canonicals = $blast_db->dbc->prepare('UPDATE transcript join translation using(transcript_id) set canonical_translation_id=translation_id');
$sth_update_canonicals->execute();

exit;


sub update_data_file {
 my ($blast_db,$species,$data_file_logic_name,$data_file_analysis_id) = @_;

  # This whole sub needs to be deleted once the pipeline is updated
  my $tissue = $data_file_logic_name;
  $tissue =~ s/$species\_//;
  $tissue =~ s/\_rnaseq\_bam//;

  my $name_match;
  my $sth_data_file_update = $blast_db->dbc->prepare('UPDATE data_file join analysis using(analysis_id) set'.
                                                     ' data_file.analysis_id=? where name like ?');

  $sth_data_file_update->bind_param(1,$data_file_analysis_id);
  if($tissue eq 'merged') {
    $name_match = "merged";
  } else {
    $name_match = '%\.'.$tissue.'\.%';
  }

  $sth_data_file_update->bind_param(2,$name_match);
  $sth_data_file_update->execute();
}

sub update_daf {
  my ($blast_db,$refine_db,$initial_logic_names) = @_;

  my $shift_value = 1000000;
  my $dump_path = "/tmp/".$refine_db->dbc->dbname.".".time().".tmp";

  # Retrieve all the refine rnaseq logic names and store their analysis ids
  my $sth_select_refine = $refine_db->dbc->prepare('SELECT analysis_id,logic_name from analysis where logic_name like "'.$species.'_%_rnaseq"');
  $sth_select_refine->execute();

  my %refine_logic_names = ();
  while(my ($analysis_id,$logic_name) =  $sth_select_refine->fetchrow_array()) {
    $refine_logic_names{$logic_name} = $analysis_id;
  }

  unless(scalar(keys(%refine_logic_names))) {
    die "Could not find any logic names of the format: ".$species.'_%_rnaseq';
  }

  my %blast_logic_names = %{$initial_logic_names};

  my $dump_daf_table = "mysqldump -u".$refineuser.
                              " -h".$refinehost.
                              " -P".$refineport.
                              " ".$refinedbname." dna_align_feature > ".$dump_path;

  my $cmd_result = system($dump_daf_table);
  if($cmd_result) {
    die "Failed to dump refine db dna_align_feature table. Command line used:\n".$dump_daf_table;
  }

  my $load_daf_table = "mysql -u".$blastuser.
                            " -h".$blasthost.
                            " -p".$blastpass.
                            " -P".$blastport.
                            " -D".$blastdbname.
                            " < ".$dump_path;

  $cmd_result = system($load_daf_table);
  if($cmd_result) {
    die "Failed to laod refine db dna_align_feature table. Command line used:\n".$load_daf_table;
  }

  # Clean the dump. Heheheh.
  system("rm ".$dump_path);

  # Need to shift analysis ids in case of conflict, so will temporarily shift
  # them outside of the allowed range of analysis ids. The small int range is only 65000 or so
  my $sth_update_daf = $blast_db->dbc->prepare('ALTER TABLE dna_align_feature MODIFY analysis_id INT');
  $sth_update_daf->execute();
  $sth_update_daf = $blast_db->dbc->prepare('UPDATE dna_align_feature set analysis_id=(analysis_id + ?)');
  $sth_update_daf->bind_param(1,$shift_value);
  $sth_update_daf->execute();

  # Now all analysis ids are shifted by 1000000, can cleanly update with no conflicts
  foreach my $refine_logic_name (keys(%refine_logic_names)) {
    my $refine_analysis_id = $refine_logic_names{$refine_logic_name};
    my $blast_logic_name = $refine_logic_name."_daf";

    my $sth_select_blast_id = $blast_db->dbc->prepare('SELECT analysis_id from analysis where logic_name = ?');
    $sth_select_blast_id->bind_param(1,$blast_logic_name);
    $sth_select_blast_id->execute();
    my $blast_analysis_id = $sth_select_blast_id->fetchrow_array;
    unless($blast_analysis_id) {
      die "Could not find an analysis id for the following refine db logic_name in the blast db: ".$blast_logic_name;
    }

    $sth_update_daf = $blast_db->dbc->prepare('UPDATE dna_align_feature set analysis_id=? where analysis_id=?');
    $sth_update_daf->bind_param(1,$blast_analysis_id);
    $sth_update_daf->bind_param(2,$refine_analysis_id + $shift_value);
    $sth_update_daf->execute();
  }

  # Now alter the table back
  $sth_update_daf = $blast_db->dbc->prepare('ALTER TABLE dna_align_feature MODIFY analysis_id SMALLINT');
  $sth_update_daf->execute();

}
