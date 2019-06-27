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


use warnings;
use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

# defaults
#use Bio::EnsEMBL::ApiVersion;
my $current_release = $ENV{ENSEMBL_RELEASE};
#my $current_release = 97;
my $max_len = 100000;
my $min_len = $max_len * 0.75;

my ($help, $reg_conf, $compara_db, $core_db, $host, $port, $pass);
GetOptions(
    "help|h"         => \$help,
    "reg_conf=s"   => \$reg_conf,
    "compara_db=s" => \$compara_db,
    "max_len=i"    => \$max_len,
    "release=i"    => \$current_release,
    "core_db=s"    => \$core_db,
    "port=s"       => \$port,
    "host=s"       => \$host,
);

die &helptext if ( $help );

# if -use_compara flag is specified, the orlog info is used to find a sample gene with a good set of orthologs, including human (this can only be used to find sample locations for species handed over this release and will only work after comapara and xrefs have been completed)
if ($reg_conf) {
  print "\nYou have specified a reg_conf. I will now run get_sample_genes, using the Compara database, for all new species in release ".$current_release."\n\n";
  open(OUT, '>', "./sample_locations.sql");
  $compara_db = 'compara_curr' unless defined $compara_db;
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($reg_conf, 0, 0, 0, "throw_if_missing") if $reg_conf;
  my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->go_figure_compara_dba( $compara_db );

# first, grab all gene trees that are good candidates
  my $gene_tree_adaptor = $compara_dba->get_GeneTreeAdaptor;
  my $good_gene_trees_sql = "SELECT root_id, tree_type, member_type, clusterset_id, gene_align_id, method_link_species_set_id, species_tree_root_id, stable_id, version, ref_root_id, NULL, NULL, NULL ".
						  "FROM gene_tree_root_attr a JOIN gene_tree_root r USING(root_id) ".
						  "WHERE a.taxonomic_coverage > 0.95 AND a.aln_percent_identity > 75 AND r.member_type = 'protein' and r.clusterset_id = 'default' and r.tree_type = 'tree';";
  my $gene_tree_sth = $compara_dba->dbc->prepare($good_gene_trees_sql);
  $gene_tree_sth->execute;
  my $good_gene_trees = $gene_tree_adaptor->_objs_from_sth($gene_tree_sth);
  print "Found " . scalar @$good_gene_trees . " 'good' gene trees\n";

# next, pick out genes from any tree that includes human (i.e. has a human ortholog)
  my $genome_db_adaptor = $compara_dba->get_GenomeDBAdaptor;
  my $human_genome_db = $genome_db_adaptor->fetch_by_name_assembly('homo_sapiens');
  my $new_genome_dbs = get_new_genome_dbs($genome_db_adaptor, $current_release);
  my %candidate_genes;
  my $x = 0;
  foreach my $tree ( @$good_gene_trees ) {
    $x++;
    my $human_gene_count = $tree->Member_count_by_source_GenomeDB('ENSEMBLPEP', $human_genome_db);
    next unless $human_gene_count > 0;
    foreach my $genome_db ( @$new_genome_dbs ) {
      my $tree_genes = $tree->get_Member_by_source_GenomeDB('ENSEMBLPEP', $genome_db);
      foreach my $gene ( @$tree_genes ) {
	push @{ $candidate_genes{$genome_db->db_adaptor->dbc->dbname} }, $gene if ($gene->length > $min_len && $gene->length < $max_len);
      }
    }
  }

# now find a candidate with xrefs and good support
  foreach my $db_name ( keys %candidate_genes ) {
    print "$db_name: " . scalar @{$candidate_genes{$db_name}} . " genes\n";
    my $sample_transcript;
  TRANSCRIPT:foreach my $seq_member ( @{$candidate_genes{$db_name}} ) {
      # find a transcript that has xrefs and 'good' supporting evidence
      my $transcript = $seq_member->get_Transcript;

      my $xrefs = $transcript->get_all_xrefs();
      if ($xrefs){
	if ($transcript->external_name()){
	  my $supporting_features = $transcript->get_all_supporting_features;
	  foreach my $support (@$supporting_features){
	    if ($support->hcoverage() >= 50 && $support->percent_id() >= 75){
	      $sample_transcript=$transcript;
	      last;
	    }
	  }
	}
	else{
	  next TRANSCRIPT;
	}
      }
      else{
	next TRANSCRIPT;
      }
    }
# finally, write SQL patch for meta table
    if ($sample_transcript){
        my $sample_gene=$sample_transcript->get_Gene;
	my $sample_coord=$sample_gene->seq_region_name().':'.$sample_gene->seq_region_start().'-'.$sample_gene->seq_region_end();

	print OUT "\nUSE ".$db_name.";
UPDATE meta set meta_value='".$sample_coord."' WHERE meta_key='sample.location_param');
UPDATE meta set meta_value='".$sample_coord."' WHERE meta_key='sample.location_text');
UPDATE meta set meta_value='".$sample_gene->stable_id()."' WHERE meta_key='sample.gene_param');
UPDATE meta set meta_value='".$sample_gene->external_name()."' WHERE meta_key='sample.gene_text');
UPDATE meta set meta_value='".$sample_transcript->stable_id()."' WHERE meta_key='sample.transcript_param');
UPDATE meta set meta_value='".$sample_transcript->external_name()."' WHERE meta_key='sample.transcript_text');
UPDATE meta set meta_value='".$sample_gene->external_name()."' WHERE meta_key='sample.search_text');"
}

#locations will now exist in meta as placholder info so replaces sql statesments with above ^
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.location_param', '".$sample_coord."');
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.location_text', '".$sample_coord."');
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.gene_param', '".$sample_gene->stable_id()."');
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.gene_text', '".$sample_gene->external_name()."');
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.transcript_param', '".$sample_transcript->stable_id()."');
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.transcript_text', '".$sample_transcript->external_name()."');"

    else{
      print "No suitable transcripts found for ".$db_name;
    }
  }
}#end if reg_conf

elsif ($core_db){
  print "\nYou have specified a core db. I will now run get_sample_genes, without Compara or Xref information.\n\n";

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -port    => $port,
      -user    => 'ensro',
      -host    => $host,
      -dbname  => $core_db);

  my $sa = $db->get_SliceAdaptor;

  my $sth_longest = $db->dbc->prepare('select seq_region_id from seq_region order by length desc limit 10');
  $sth_longest->execute;

  my $sample_transcript;
  LOOP: while (my $seq_region_id = $sth_longest->fetchrow_array) {
    my $region = $sa->fetch_by_seq_region_id($seq_region_id);
    my @transcripts = @{$region->get_all_Transcripts_by_type('protein_coding')};

  TRANSCRIPT:foreach my $transcript (@transcripts){
      my $supporting_features = $transcript->get_all_supporting_features;
      foreach my $support (@$supporting_features){
        if ($support->hcoverage() >= 99 && $support->percent_id() >= 75){
          $sample_transcript=$transcript;
          last LOOP;
        }
        else{
          next TRANSCRIPT;
        }
      }
    }
  }#end while
# finally, write SQL patch for meta table
  if ($sample_transcript){
    my $db_name = $db->dbc->dbname;
    my $sample_gene=$sample_transcript->get_Gene;
    my $sample_coord=$sample_gene->seq_region_name().':'.$sample_gene->seq_region_start().'-'.$sample_gene->seq_region_end();

    open(OUT, '>', "./".$db_name."_placeholder_locations.sql");

    print OUT "USE ".$db_name.";
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.location_param', '".$sample_coord."');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.location_text', '".$sample_coord."');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.gene_param', '".$sample_gene->stable_id()."');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.gene_text', 'ensembl_gene');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.transcript_param', '".$sample_transcript->stable_id()."');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.transcript_text', 'ensembl_transcript');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.search_text', 'ensembl_gene');\n"
}
  else{
    print "No suitable transcripts found for $core_db\n";
  }
}#end if core_db

else {
  print "You have not specified a reg_conf or core db... I have no idea what to do!\n";
  die &helptext;
}

sub helptext {
	my $msg = <<HELPEND;

get_sample_genes.pl

Choose a sample gene for the "Sample Location" page on ensembl.org

Usage:

1. To run on a single core db, without Compara and Xref information, i.e. before handover, to provide a placeholder sample location:
perl get_sample_genes.pl --core_db name_of_core_db --port port --host host

2. To run with Compara and Xref information, i.e. during the release cycle, for all new core databses in this release:
perl get_sample_genes.pl --reg_conf /path/to/registry/config

Options:

-max_len N            *default=100000
      Maximum length of sample gene

-release N            *default=current Ensembl release cycle
      Ensembl release


HELPEND
	return $msg;
}

sub get_new_genome_dbs {
	my ($gdb_adap, $release_num) = @_;

	my @new_this_rel;
	my $all_this_rel = $gdb_adap->fetch_all_by_release($release_num);
	foreach my $gdb ( @$all_this_rel ) {
		push @new_this_rel, $gdb if $gdb->first_release == $release_num;
	}
	return \@new_this_rel;
}
