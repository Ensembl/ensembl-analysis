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
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;

$| = 1;

my $host         = '';
my $user         = '';
my $pass         = '';
my $port         = '';
my $dbname       = '';
my $gene_host    = '';
my $gene_user    = '';
my $gene_pass    = '';
my $gene_port    = '';
my $gene_dbname  = '';
my $out_host    = '';
my $out_user    = '';
my $out_pass    = '';
my $out_port    = '';
my $out_dbname  = '';
my $coord_sys_ver = '';
my @patch_types = ('PATCH_FIX','PATCH_NOVEL');

&GetOptions( 'mapdbhost:s'       => \$host,
             'mapdbuser:s'       => \$user,
             'mapdbname:s'       => \$dbname,
             'mapdbpass:s'       => \$pass,
             'mapdbport:n'       => \$port,
             'genedbhost:s'      => \$gene_host,
             'genedbuser:s'      => \$gene_user,
             'genedbname:s'      => \$gene_dbname,
             'genedbpass:s'      => \$gene_pass,
             'genedbport:n'      => \$gene_port,
             'outdbhost:s'       => \$out_host,
             'outdbuser:s'       => \$out_user,
             'outdbname:s'       => \$out_dbname,
             'outdbpass:s'       => \$out_pass,
             'outdbport:n'       => \$out_port,
             'coord_sys_ver:s'   => \$coord_sys_ver);

#database with mapping between coord systems
my $mapdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                             -user   => $user,
                                             -pass   => $pass,
                                             -port   => $port,
                                             -dbname => $dbname );

#gene source db
my $genedb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -dnadb => $mapdb,
                                             -host   => $gene_host,
                                             -user   => $gene_user,
                                             -pass   => $gene_pass,
                                             -port   => $gene_port,
                                             -dbname => $gene_dbname);

#output db for projected genes
my $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -dnadb => $mapdb,
                                             -host   => $out_host,
                                             -user   => $out_user,
                                             -pass   => $out_pass,
                                             -port   => $out_port,
                                             -dbname => $out_dbname);


my $ga = $genedb->get_GeneAdaptor();
my $outga =$outdb->get_GeneAdaptor();
my $map_sa = $mapdb->get_SliceAdaptor();

#get the patches
print "Getting patches...\n";
my $asm_exc_adaptor = $mapdb->get_AssemblyExceptionFeatureAdaptor();
my @exceptions = @{$asm_exc_adaptor->fetch_all()};
my @patches;
EXC: foreach my $exc (@exceptions){
  foreach my $type (@patch_types){
    if($exc->type() =~ m/$type/){
      push(@patches, $exc);
      next EXC;
    }
  }
}
#Assuming that AssemblyExceptionFeatureAdaptor's fetch_all will always return two 
#entries for each patch and that those two entries are adjacent
my $num_patches = scalar(@patches)/2;
print "Have ".$num_patches." patches.\n";
#for each patch
for (my $i = 0;$i < $num_patches;$i++){

  #get the two slices
  my $ref_slice;
  my $patch_slice;
  for(my $j = 0;$j < 2; $j++){
    my $exc = pop(@patches);
    #if this is the ref version
    if($exc->type =~ m/REF/){
      #alt is only the patch slice
      $patch_slice = $exc->alternate_slice();
    }
    else{
      #alt is replaced region of ref
      $ref_slice = $exc->alternate_slice();
    }    
  }
  if(!($patch_slice and $ref_slice)){
    throw("Something is wrong, the patch and ref slices were not set correctly.\n");
  }


  #get the genes
  my @ref_genes = @{$ga->fetch_all_by_Slice($ref_slice, undef, 1)};
  print $patch_slice->name." slice ".$ref_slice->name." has ".scalar(@ref_genes)." genes\n";

  #project the genes
  foreach my $ref_gene (@ref_genes){
    $ref_gene->load;
    my $coord_system_name = $patch_slice->coord_system->name;
    my $alt_gene = $ref_gene->transform($coord_system_name, $coord_sys_ver);
    if($alt_gene){
      print "Transformed ".$ref_gene->stable_id."\n";

      # add parent exon key attribute and parent stable id attribute to the transcripts
      foreach my $transcript (@{ $alt_gene->get_all_Transcripts }) {
        print "Looking at transcript ".$transcript->stable_id." dbID: ".$transcript->dbID." region: ".$transcript->seq_region_name."\n";
        #my @trans_attributes = @{$transcript->get_all_Attributes()};
        #$transcript->{'attributes'} = []; # flush attributes
        warn ref($transcript->{attributes});
        print $transcript->{attributes};

        my $parent_stable_id = $transcript->stable_id; # at this stage just after the projection, the stable IDs are duplicated so the transcript stable ID is the same as the parent transcript stable ID
        my $parent_exon_key = "";

        # looking for $transcript's parent transcript
        foreach my $ref_transcript (@{ $ref_gene->get_all_Transcripts }) {
          if ($ref_transcript->stable_id eq $parent_stable_id) {
            $parent_exon_key = get_transcript_exon_key($ref_transcript);
            last;
          }
        }

        my $parent_sid_att = Bio::EnsEMBL::Attribute->new
          (-CODE => "parent_sid",
           -NAME => "parent_sid",
           -DESCRIPTION => "The parent stable ID to identify a projected transcript's parent transcript. For internal statistics use only since this method does not work in all cases.",
           -VALUE => "$parent_stable_id");

        my $exon_key_att = Bio::EnsEMBL::Attribute->new
          (-CODE => "parent_exon_key",
           -NAME => "parent_exon_key",
           -DESCRIPTION => "The exon key to identify a projected transcript's parent transcript.",
           -VALUE => "$parent_exon_key");

        print "Transcript attribute 'parent_stable_id' with value $parent_stable_id will be stored.\n";
        print "Transcript attribute 'parent_exon_key' with value $parent_exon_key will be stored.\n";

        #$transcript->add_Attributes(@trans_attributes,$parent_sid_att,$exon_key_att);
        $transcript->add_Attributes($parent_sid_att,$exon_key_att);
      }
      my $id = $outga->store($alt_gene);
      print "-gene->store returned $id\n";
    }
    else{
      print "Failed to transform ".$ref_gene->stable_id."\n";
    }
  }
}

print "Completed\n";


sub get_transcript_exon_key {
  my $transcript = shift;
  my $string = $transcript->slice->seq_region_name.":".$transcript->biotype.":".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->seq_region_strand.":";

  my $exons = sort_by_start_end_pos($transcript->get_all_Exons);
  foreach my $exon (@{$exons}) {
    $string .= ":".$exon->seq_region_start.":".$exon->seq_region_end;
  } 

  return $string;
}

sub sort_by_start_end_pos {
  my ($unsorted) = @_;

  my @sorted = sort { if ($a->seq_region_start < $b->seq_region_start) {
        return -1;
    } elsif ($a->seq_region_start == $b->seq_region_start) {
      if ($a->seq_region_end < $b->seq_region_end) {
        return-1;
      } elsif ($a->seq_region_end == $b->seq_region_end) {
        return 0;
      } elsif ($a->seq_region_end > $b->seq_region_end) {
        return 1;
      }
        return 0;
    } elsif ($a->seq_region_start > $b->seq_region_start) {
        return 1;
    }
  } @$unsorted;

  return \@sorted;
}
