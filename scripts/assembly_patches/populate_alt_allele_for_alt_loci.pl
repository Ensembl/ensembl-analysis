#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head2
  This script assigns alt allele relationships between the genes on the slices corresponding to
  the alt_seq_mapping dna align features and the genes on their corresponding hit slices.
  There needs to be an exact match in terms of biotypes, exon structure and gene start coordinates or
  a protein coverage and percent identity of the canonical transcript protein on the primary assembly
  against the one on the assembly exception above a given threshold and vice-versa AND the
  starts of the gene alleles have to lie within a GENE_START_WINDOW bp window.
  
  It requires that canonical transcripts are set.
=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(align_proteins);

$| = 1;

my $min_cov = 95;
my $min_pid = 95;
my $core_host = '';
my $core_user = '';
my $core_pass = '';
my $core_port = '';
my $core_dbname = '';

my $output_file = '';
my $verbose;

my $GENE_START_WINDOW = 100; # the gene allele starts are allowed to lie within a window of up to 100 bp long either upstream or downstream

&GetOptions(  'min_cov:n'  => \$min_cov,
              'min_pid:n'  => \$min_pid,
              'dbhost:s'   => \$core_host,
              'dbuser:s'   => \$core_user,
              'dbpass:s'   => \$core_pass,
              'dbport:n'   => \$core_port,
              'dbname:s'   => \$core_dbname,
              'output_file:s' => \$output_file,
              'verbose!'    => \$verbose);
# ouptut
my $fh;
if ($output_file && $output_file ne "stdout") {
  open WRITE,">$output_file" or die "couldn't open file ".$output_file." $!";
  $fh = \*WRITE;
} else {
  $fh = \*STDOUT;
}

# get adaptors
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $core_host,
                                                  -user   => $core_user,
                                                  -pass   => $core_pass,
                                                  -port   => $core_port,
                                                  -dbname => $core_dbname );
my $sa = $db->get_SliceAdaptor();
my $da = $db->get_DnaAlignFeatureAdaptor();
my $aaga = $db->get_AltAlleleGroupAdaptor();

my %alt_gene_flags = ('AUTOMATICALLY_ASSIGNED' => '1');
my %primary_gene_flags = ('AUTOMATICALLY_ASSIGNED' => '1', 'IS_REPRESENTATIVE' => '1');

my @alt_seq_mappings = @{$da->fetch_all_by_logic_name("alt_seq_mapping")};
foreach my $alt_seq_mapping (@alt_seq_mappings) {
  print $alt_seq_mapping->seq_region_name()." ".$alt_seq_mapping->seq_region_start()." ".$alt_seq_mapping->seq_region_end()." ".$alt_seq_mapping->hseqname()." ".$alt_seq_mapping->hstart()." ".$alt_seq_mapping->hend()."\n";
  
  # get the genes on the source region (usually hap/patch)
  my $ae_slice = $sa->fetch_by_region('toplevel',$alt_seq_mapping->seq_region_name(),$alt_seq_mapping->seq_region_start(),$alt_seq_mapping->seq_region_end(),$alt_seq_mapping->seq_region_strand());
  my @ae_genes = @{$ae_slice->get_all_Genes()};
  print "num of ae genes:".scalar(@ae_genes)."\n";

  # get the genes on the target (mapped) region (usually primary assembly)
  my $pa_slice = $sa->fetch_by_region('toplevel',$alt_seq_mapping->hseqname(),$alt_seq_mapping->hstart(),$alt_seq_mapping->hend(),$alt_seq_mapping->hseq_region_strand());
  my @pa_genes = @{$pa_slice->get_all_Genes()};
  print "num of pa genes:".scalar(@pa_genes)."\n";

  foreach my $ae_gene (@ae_genes) {
  	
  	my $ae_gene_key = get_gene_key($ae_gene,$alt_seq_mapping->hstart());
  	
  	PA_GENE: foreach my $pa_gene (@pa_genes) {
  	
  	  my $pa_gene_key = get_gene_key($pa_gene,$alt_seq_mapping->seq_region_start());

  	  my ($cov1,$pid1) = (0,0);
  	  my ($cov2,$pid2) = (0,0);
  	  my $pa_t_seq;
  	  my $ae_t_seq;
  	  my $pa_can_t = $pa_gene->canonical_transcript();
  	  my $ae_can_t = $ae_gene->canonical_transcript();
  	  
  	  if ($pa_can_t->translation() and $ae_can_t->translation()) {
  	    $pa_t_seq = $pa_can_t->translate()->seq();
  	    $ae_t_seq = $ae_can_t->translate()->seq();
  	  }

  	  if ($pa_t_seq and $ae_t_seq) {
  	  	($cov1,$pid1) = align_proteins($pa_t_seq,$ae_t_seq);
  	  	($cov2,$pid2) = align_proteins($ae_t_seq,$pa_t_seq);
  	  } else {
  	  	print "Canonical transcripts not found for gene pair ".$pa_gene->dbID()." (".$pa_gene->stable_id()."), ".$ae_gene->dbID()." (".$ae_gene->stable_id().").\n";
  	  }

  	  if (($ae_gene_key eq $pa_gene_key) or
  	      ($cov1 > $min_cov and $pid1 > $min_pid and $cov2 > $min_cov and $pid2 > $min_pid and
  	       abs(($ae_gene->seq_region_start()-$alt_seq_mapping->hstart())-($pa_gene->seq_region_start()-$alt_seq_mapping->seq_region_start())) < $GENE_START_WINDOW)) {

  	  	my $pa_gene_group = $aaga->fetch_by_gene_id($pa_gene->dbID());
  	  	my $ae_gene_group = $aaga->fetch_by_gene_id($ae_gene->dbID());
  	
  	  	if ($pa_gene_group and $ae_gene_group) {
  	  	  if ($pa_gene_group->dbID() != $ae_gene_group->dbID()) {
  	  	  	warning("FOUND alt allele relationship between gene IDs ".$pa_gene->dbID()." and ".$ae_gene->dbID()." but both of them are part of different alt allele groups so the creation of a new alt allele group is not allowed because a gene is only expected to be part of ONE alt allele group.");
  	  	  } else {
  	  	  	print "FOUND already existing alt_allele relationship for: ".$ae_gene->stable_id()." ".$ae_gene_key." and ".$pa_gene->stable_id()." ".$pa_gene_key."\n";
  	  	  }
  	  	} elsif ($pa_gene_group) {
  	  	  print "FOUND new alt_allele relationship for: ".$ae_gene->stable_id()." ".$ae_gene_key." and ".$pa_gene->stable_id()." ".$pa_gene_key." Adding the former to the latter's alt allele group ".$pa_gene_group->dbID().".\n";
  	  	  $pa_gene_group->add_member($ae_gene->dbID(),\%alt_gene_flags);
  	  	  $aaga->update($pa_gene_group);
  	  	} elsif ($ae_gene_group) {

          # if a gene on the primary assembly is found, we don't want to add more genes on the primary assembly to this alt allele group
          foreach my $ae_gene_group_gene (@{$ae_gene_group->get_all_Genes()}) {
            if ($ae_gene_group_gene->seq_region_name() !~ /CHR/) {
              print "FOUND primary assembly gene in this alt allele group. Skipping alt allele group ".$ae_gene_group->dbID()."\n";
              next PA_GENE;
            } else {
              ;
            }
          }
          
          print "FOUND new alt_allele relationship for: ".$ae_gene->stable_id()." ".$ae_gene_key." and ".$pa_gene->stable_id()." ".$pa_gene_key." Adding the latter to the former's alt allele group ".$ae_gene_group->dbID().".\n";
          $ae_gene_group->add_member($pa_gene->dbID(),\%primary_gene_flags);
  	  	  $aaga->update($ae_gene_group);

  	  	} else {
          my $new_group = Bio::EnsEMBL::AltAlleleGroup->new(-MEMBERS => [ [$pa_gene->dbID(),\%primary_gene_flags ],
                                                                          [$ae_gene->dbID(),\%alt_gene_flags] ]);
          my $new_group_id = $aaga->store($new_group);
          print "FOUND new alt_allele relationship for: ".$ae_gene->stable_id()." ".$ae_gene_key." and ".$pa_gene->stable_id()." ".$pa_gene_key." Creating a new alt allele group ".$new_group_id.".\n";
  	  	}
  	  } #if
  	} #foreach pa_gene
  } #foreach ae_gene
}
print $fh "DONE\n";
close($fh);

sub get_gene_key {
  my ($gene,$region_start) = @_;

  my $string = ($gene->seq_region_start()-$region_start).":".$gene->biotype().":".$gene->length().":".$gene->seq_region_strand();
  
  my $transcripts = sort_by_start_end_pos($gene->get_all_Transcripts());
  foreach my $transcript (@{$transcripts}) {
    $string .= ":".get_transcript_key($transcript);
  }
  return $string;
}

sub get_transcript_key {
  my $transcript = shift;
  my $string = $transcript->biotype.":".$transcript->length().":".$transcript->seq_region_strand;

  my $exons = sort_by_start_end_pos($transcript->get_all_Exons);
  foreach my $exon (@{$exons}) {
    $string .= ":".$exon->length();
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

