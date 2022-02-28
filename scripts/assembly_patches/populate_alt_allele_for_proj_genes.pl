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

=head2
  This script assigns alt alleles for projected genes
  and assumes that the vega alt alleles are already inserted in the output database
=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

$| = 1;

my $core_host = '';
my $core_user = '';
my $core_pass = '';
my $core_port = '';
my $core_dbname = '';

my $old_host = '';
my $old_user = '';
my $old_port = '';
my $old_dbname = '';

my $existing_alt_allele = '';
my $output_file = '';
my $stable_id_prefix; # ENST or ENSMUST
my $verbose;

&GetOptions(  'dbhost:s'   => \$core_host,
              'dbuser:s'   => \$core_user,
              'dbpass:s'   => \$core_pass,
              'dbport:n'   => \$core_port,
              'dbname:s' => \$core_dbname,
              'olddbhost:s'   => \$old_host,
              'olddbuser:s'   => \$old_user,
              'olddbport:n'   => \$old_port,
              'olddbname:s' => \$old_dbname,
              'stable_id_prefix:s'    => \$stable_id_prefix,
              'existing_alt_allele:s' => \$existing_alt_allele,
              'output_file:s'    => \$output_file,
              'verbose!'    => \$verbose);
# ouptut
my $fh;
if ($output_file && $output_file ne "stdout") {
  open WRITE,">$output_file" or die "couldn't open file ".$output_file." $!";
  $fh = \*WRITE;
} else {
  $fh = \*STDOUT;
}



if (!defined $stable_id_prefix) {
  throw("Stable ID prefix not defined");
} else {
  chomp $stable_id_prefix;
  print $fh "Stable ID prefix to match is '$stable_id_prefix'\n";
}





#get core db adaptor
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $core_host,
                                                  -user   => $core_user,
                                                  -pass   => $core_pass,
                                                  -port   => $core_port,
                                                  -dbname => $core_dbname );


my $old_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $old_host,
                                                  -user   => $old_user,
                                                  -port   => $old_port,
                                                  -dbname => $old_dbname );

#get adaptors
my $ga = $db->get_GeneAdaptor();
my $ta = $db->get_TranscriptAdaptor();
my $sa = $db->get_SliceAdaptor();
my $old_ga = $old_db->get_GeneAdaptor();
my $old_ta = $old_db->get_TranscriptAdaptor();
# alt allele group adaptor
my $aaga = $db->get_AltAlleleGroupAdaptor();



# # #
# First fetch all existing alt allele groups from the output database
# There are probably some here from havana already (previous script in ensembl/misc-scripts/alt_alleles.pl)
# # #
my $alt_allele_groups = $aaga->fetch_all();

# Just for curiousity:
# Work out what the max alt_allele_id and alt_allele_group_id is
my $max_alt_allele_group_id = 0;
my $num_genes = 0;
# what we really want is this: 
#   the hash where key = gene dbid of gene on primary assembly 
#   and value = alt_alle_group_id
my %ref_genes_to_allele_group_id;
my %alt_allele_gene_ids;

foreach my $alt_allele_group ( @{$alt_allele_groups} ) {
  if ($alt_allele_group->dbID > $max_alt_allele_group_id) {
    $max_alt_allele_group_id = $alt_allele_group->dbID;
  }
  foreach my $member ( @{$alt_allele_group->get_all_members} ) {
    $num_genes++;
    $alt_allele_gene_ids{$member->[0]} = 1;
    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($member->[0]);
    print "Got gene ".$gene->stable_id." on slice ".$gene->slice->seq_region_name."\n" if $verbose;
    my $flag_list = join (',',keys %{$member->[1]});
    if ($gene->slice->is_reference) {
      # gene is on the primary assembly unit
      $ref_genes_to_allele_group_id{$gene->dbID} = $alt_allele_group->dbID;
    }
    print $fh "Alt_allele_group ".$alt_allele_group->dbID." includes gene with dbID ".$gene->dbID." on slice ".$gene->slice->seq_region_name." ($flag_list )\n";
  }
}
print $fh "\nMax max_alt_allele_group_id = ".$max_alt_allele_group_id." and includes a total of $num_genes genes\n";





# # #
# Now we want to know whether we can add genes to existing alt_allele_groups
#   or whether we need to make new groups.
# # #

# Since this script is only trying to add alt_alleles for Projected Genes, 
#   we will always have a gene on a reference slice in the alt_allele_groups we make

# Get the projected genes 
my $aa = $db->get_AnalysisAdaptor();
my @analyses = @{$aa->fetch_all()};
my @projected_logic_names;

#the projected logic names (same as the core but with proj_ at the start)
#NB: there are some logic_names that start with proj_ at transcript level
foreach my $analysis (@analyses){
  if($analysis->logic_name() =~ m/^proj_/){
    #print $fh $analysis->logic_name()."\n";
    push @projected_logic_names, $analysis->logic_name();
  }
}

# fetch projected genes by logic name
my @projected_genes;
foreach my $logic_name (sort @projected_logic_names){
  push @projected_genes, @{$ga->fetch_all_by_logic_name($logic_name)};
  print $fh $logic_name." ".scalar(@projected_genes)."\n";
}

#work out the relationships between projected gene and parent gene
foreach my $proj_gene (@projected_genes){

  my $patch_slice = $sa->fetch_by_region(undef,$proj_gene->seq_region_name());;

  #get reference slice
  my $ref_slice = get_ref_slice($patch_slice);
  if (!$ref_slice) {
    throw("Cannot find reference slice for patch slice ".$patch_slice->name());
  }

  #create reference transcripts hash
  my @ref_genes = @{$ga->fetch_all_by_Slice($ref_slice)};
  my @ref_transcripts;
  my %ref_trans_hash = ();
  foreach my $ref_gene (@ref_genes){ # make sure get transcripts outside the patch
    my @ref_trans = @{$ref_gene->get_all_Transcripts};
    foreach my $t (@ref_trans) {
      $ref_trans_hash{get_transcript_exon_key($t)} = $t->stable_id;
      push(@ref_transcripts,$t);
    }
  }

  my $orig_transcript_id = "";
  my %unique_parent_gene;

  # we want to check that all transcripts in the projected gene point back to ONE parent
  # (this has not always been the case in the past)...
  # do this more extended check:
  foreach my $proj_transcript ( @{$proj_gene->get_all_Transcripts()} ) {
  	
  	my @parent_exon_key_attribs = @{$proj_transcript->get_all_Attributes('parent_exon_key')};
    my $parent_exon_key = "";

    if (scalar(@parent_exon_key_attribs) > 1) {
      warning("Projected transcript ".$proj_transcript->stable_id." has more than 1 parent_exon_key attribute.");
    } elsif (scalar(@parent_exon_key_attribs) > 0) {
      $parent_exon_key = $parent_exon_key_attribs[0]->value();
    }

    if ($parent_exon_key) {
      if (exists($ref_trans_hash{$parent_exon_key})) {

        # fetch parent
        my $parent_g = $old_ga->fetch_by_transcript_stable_id($ref_trans_hash{$parent_exon_key});
        
        # store parent
        push(@{$unique_parent_gene{$parent_g->stable_id}},$proj_transcript);

      } else {
        warning("Projected transcript ".$proj_transcript->stable_id()." does not have a parent.");
      }
    } else {
      throw("Projected transcript ".$proj_transcript->stable_id()." does not have a parent_exon_key attribute. It should have been added during the projection stage!");
    }
  } # foreach proj_transcript

  # check only one parent
  my $parent_gene;
  my @parent_stable_ids = keys %unique_parent_gene;
  if (scalar( @parent_stable_ids ) != 1) {
    # we could swicth this to a warning and use the 'best guess' parent
    warning("Expected only one parent. Got the following parents for projected gene ".$proj_gene->stable_id.":\n  @parent_stable_ids");
    # try to choose the correct one - best guess parent
    my $has_most_projected = 0;
    foreach my $stable_id (keys %unique_parent_gene) {
      if (scalar @{$unique_parent_gene{$stable_id}} > $has_most_projected ) {
        $parent_gene = $old_ga->fetch_by_stable_id( $stable_id );
        $has_most_projected = scalar @{$unique_parent_gene{$stable_id}};
      }
    }
    warning("Set parent to be:".$parent_gene->stable_id);
  } else {
    # we have only one parent - hurrah!
    $parent_gene = $old_ga->fetch_by_stable_id( $parent_stable_ids[0] );   
  }

  if (!defined $parent_gene) {
    throw("Unable to find $orig_transcript_id in $old_dbname");
  }
  print $fh "Fetch parent $orig_transcript_id on slice ".$parent_gene->slice->name." for projected ".$proj_gene->stable_id."\n";
  my $parent_gene_string = get_gene_info($parent_gene);
  my $orig_slice = $db->get_SliceAdaptor->fetch_by_name($parent_gene->feature_Slice->name);
  my @possible_orig_genes = @{$ga->fetch_all_by_Slice($orig_slice)};
  my $orig_gene;
  print $fh "Found ".(scalar(@possible_orig_genes))." possible genes for ".$proj_gene->stable_id." on slice ".$parent_gene->feature_Slice->name."\n";
  my $count_matched = 0;
  foreach my $possible (@possible_orig_genes) {
    my $possible_gene_string = get_gene_info($possible);
    print $fh "   PAR $parent_gene_string\n   POS $possible_gene_string\n";
    if ($parent_gene_string eq $possible_gene_string) {
      $orig_gene = $possible;
      $count_matched++;
    }
  } 
  if ($count_matched != 1) {
    throw("Matched ".$count_matched." parent genes for ".$unique_parent_gene{$parent_gene->stable_id}[0]->stable_id."\n");
    $orig_gene = undef;
  } else {
    print $fh "MATCHED! ".$proj_gene->stable_id."\n";
  }
 

  # # #
  # This is the bit where we decide whether or not to add the projected gene into an
  # existing alt_allele_group or to make a new alt_allele_group
  # # # 
  if (get_gene_info($orig_gene) eq get_gene_info($proj_gene)) {
    #added orig and proj gene comparison to avoid duplicated gene_ids 
    warning("orig_gene same as proj_gene\n");
  } elsif (!($alt_allele_gene_ids{$proj_gene->dbID})) {
    print $fh "ADDING ALT ALLELE... ";
    if (exists $ref_genes_to_allele_group_id{$orig_gene->dbID} ) {
      #is the ref (parent) gene in the alt_allele_group table? 
      # Yes it is, so we are adding a member (gene) to an exisitn alt_allele_group
      # fetch the alt_allele_group
      my $alt_allele_group = $aaga->fetch_by_gene_id($orig_gene->dbID);
      my $member_found = 0;
      foreach my $allele (@{$alt_allele_group->get_all_members}) {
        my ($gene_id,$type) = @$allele;
        # see if this gene already exosts in the group
        if ($gene_id == $proj_gene->dbID) {
          $member_found = 1;
        }
        foreach my $flag ( keys %{$type} ) {
          print "\nExisting alt_allele_group includes gene $gene_id with type $flag\n" if $verbose;
        }
      }
      if ($member_found) {
        print $fh "oh no not necessary because it's already there in ".$alt_allele_group->dbID." so I'm not doing to update the group\n";
      } elsif (!($alt_allele_gene_ids{$proj_gene->dbID})) {
        # add the projected gene as a new member
        my %proj_gene_flags = ('AUTOMATICALLY_ASSIGNED' => '1');
        $alt_allele_group->add_member($proj_gene->dbID , \%proj_gene_flags);
        $alt_allele_gene_ids{$proj_gene->dbID} = 1;
        # store the changes by doing an update
        print $fh "to existing alt_allele_group";
        my $alt_allele_group_id = $aaga->update($alt_allele_group);
        print $fh " with dbID ".$alt_allele_group_id." as I see that gene ".$orig_gene->stable_id." (dbID ".$orig_gene->dbID.") on primary already stored\n";
      }
    } else {
      # We are going to make a new alt_allele_group
      # that contains two genes: projected and its parent on reference
      my %proj_gene_flags = ('AUTOMATICALLY_ASSIGNED' => '1');
      my %orig_gene_flags = ('AUTOMATICALLY_ASSIGNED' => '1', 'IS_REPRESENTATIVE' => '1');

      if (!($alt_allele_gene_ids{$proj_gene->dbID}) and !($alt_allele_gene_ids{$orig_gene->dbID})) {

       my $alt_allele_group = Bio::EnsEMBL::AltAlleleGroup->new(
                              -MEMBERS => [ [$orig_gene->dbID,\%orig_gene_flags], [$proj_gene->dbID,\%proj_gene_flags] ],
                               );
        $alt_allele_gene_ids{$proj_gene->dbID} = 1;
        foreach my $allele (@{$alt_allele_group->get_all_members}) {
          my ($gene_id,$type) = @$allele;
          foreach my $flag ( keys %{$type} ) {
            print "\nNew alt_allele_group includes gene $gene_id with type $flag\n" if $verbose;
          }
        }
        print $fh "into a new alt_allele_group";
        my $alt_allele_group_id = $aaga->store($alt_allele_group);
        print $fh " with dbID ".$alt_allele_group_id."\n";

        # add the new alt allele group to the hash
        $ref_genes_to_allele_group_id{$orig_gene->dbID} = $alt_allele_group->dbID;

      }
    }
  }
}


print $fh "DONE\n\n";
close($fh);

sub get_gene_info {
  my ($gene) = @_;
  my $string = $gene->slice->seq_region_name.":".$gene->seq_region_start.":".$gene->seq_region_end.":".$gene->seq_region_strand.":";

  my $exons = sort_by_start_end_pos($gene->get_all_Exons);
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
sub get_transcript_exon_key {
  my $transcript = shift;
  my $string = $transcript->slice->seq_region_name.":".$transcript->biotype.":".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->seq_region_strand.":";

  my $exons = sort_by_start_end_pos($transcript->get_all_Exons);
  foreach my $exon (@{$exons}) {
    $string .= ":".$exon->seq_region_start.":".$exon->seq_region_end;
  } 

  return $string;
}

sub get_ref_slice {
  my $patch_slice = shift;
  my $ref_slice;

  print "patch slice is: ".$patch_slice->name."\n";

  my @excs = $patch_slice->get_all_AssemblyExceptionFeatures();

  if (@excs) {
    foreach my $exc (@excs) {
      if (@$exc[0]) {
        if (@$exc[0]->type() !~ m/REF/) {
            $ref_slice = @$exc[0]->alternate_slice();
        }
      }
    }
  }
  print "reference slice is: ".$ref_slice->name."\n";
  return $ref_slice;
}
