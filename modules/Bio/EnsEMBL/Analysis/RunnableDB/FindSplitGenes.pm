=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::FindSplitGenes - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::FindSplitGenes;  

use warnings ;
use strict;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Gene;  
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils; 
#use Bio::EnsEMBL::Analysis::Tools::Utilities; 
use Bio::EnsEMBL::Analysis::Tools::Utilities qw ( shuffle ) ; 

use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator; 
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;   

use Bio::EnsEMBL::Registry; 
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::Rule;  

use vars qw(@ISA);  

use Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator; 

@ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator ) ; 


sub new {
  my ($class,@args) = @_;  
  my $self = $class->SUPER::new(@args);  
  return $self;
}



# fetches the genes out of the QUERY database 
# ( db which contains the predictions of species we're working on ) 


# The module FindSplitGenes queries the compara database with raw sql. 
#
# It uses any random string as input id. it does not run on slices etc. 
#  
#
sub fetch_input {
  my( $self ) = @_;  

  my $dbname = 'Compara' ;  
  my @potential_stable_ids_for_recovery ; 

  # INSERT INTO `analysis` VALUES (100,'2006-11-22 16:39:35','FindSplitGenes',NULL,NULL,NULL,NULL,NULL,NULL,NULL,'FindSplitGenes',NULL,NULL,NULL); 

  my $dbc = Bio::EnsEMBL::Registry->get_DBAdaptor($dbname ,'compara') ;  
  my $gdba = Bio::EnsEMBL::Registry->get_adaptor($dbname,'compara','GenomeDB');
  my $mlssa = Bio::EnsEMBL::Registry->get_adaptor($dbname,'compara','MethodLinkSpeciesSet');

  unless ($dbc){ 
    throw ("Can't get Adaptor for comara ")  ; 
  }  

  my $method_link_type = "ENSEMBL_ORTHOLOGUES";

  # this should be the same so REVOCER_SPLIT_GENES ... can be removed ...
  my $qy_species  = $$MAIN_CONFIG{QUERY_SPECIES};  
  my $targeted_species = $$FIND_SPLIT_GENES{ANALYSIS_SETS}{$self->post_logic_name}{TARGETTED}; 

  my $informant_species = $$FIND_SPLIT_GENES{ANALYSIS_SETS}{$self->post_logic_name}{INFORMANT}; 

  #my $tg2_species = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{TARGETED};

  my $informant_perc_cov = $$FIND_SPLIT_GENES{ANALYSIS_SETS}{$self->post_logic_name}{INFORMANT_PERC_COV} ; 
  my $targeted_perc_cov  = $$FIND_SPLIT_GENES{ANALYSIS_SETS}{$self->post_logic_name}{TARGETTED_PERC_COV} ;   

  my $informant_gdb = $gdba->fetch_by_name_assembly($informant_species);
  my $targeted_gdb = $gdba->fetch_by_name_assembly($targeted_species); 
  my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type,[$informant_gdb, $targeted_gdb]);

  my $no_chr_constraint = 0;
 

  my $sql = "SELECT im.stable_id FROM
  homology h, homology_member thm, member tm, homology_member ihm, member im 
  WHERE h.homology_id=thm.homology_id AND 
  thm.member_id=tm.member_id AND 
  h.homology_id=ihm.homology_id AND 
  ihm.member_id=im.member_id AND 
  h.method_link_species_set_id = ? AND 
  tm.genome_db_id = ? AND im.genome_db_id = ? AND 
  ihm.perc_cov < ? AND thm.perc_cov > ? 
  group by im.stable_id;";
  
  my $sth = $dbc->prepare($sql);
  $sth->execute($mlss->dbID, $targeted_gdb->dbID, $informant_gdb->dbID,$informant_perc_cov, $targeted_perc_cov);
  
  $sql = "SELECT h.description,im.stable_id,im.chr_name,im.chr_start,im.chr_end,im.chr_strand,
  ihm.perc_cov,tm.stable_id,tm.chr_name,tm.chr_start,tm.chr_end,tm.chr_strand,thm.perc_cov FROM 
  homology h, homology_member thm, member tm, homology_member ihm, member im 
  WHERE h.homology_id=thm.homology_id AND 
  thm.member_id=tm.member_id AND 
  h.homology_id=ihm.homology_id AND 
  ihm.member_id=im.member_id AND 
  h.method_link_species_set_id = ? AND 
  tm.genome_db_id = ? AND im.genome_db_id = ? AND 
  im.stable_id = ?";
  
  my $stable_id;
  $sth->bind_columns(\$stable_id);
  my $sth1 = $dbc->prepare($sql);
  
  my $targeted_split_genes = 0;
  my %stable_ids;
  my ($min_start, $max_end);
  
  print "#Possible split genes in $targeted_species #using $informant_species as informant\n";
  
  while ($sth->fetch) {
    $sth1->execute($mlss->dbID, $targeted_gdb->dbID, $informant_gdb->dbID, $stable_id);
  
    next if ($sth1->rows < 2);
  
    my ($description,$im_stable_id,$im_chr_name,$im_chr_start,$im_chr_end,$im_chr_strand,
        $ihm_perc_cov,$tm_stable_id,$tm_chr_name,$tm_chr_start,$tm_chr_end,$tm_chr_strand,
        $thm_perc_cov);
  
    $sth1->bind_columns(\$description,\$im_stable_id,\$im_chr_name,
                        \$im_chr_start,\$im_chr_end,\$im_chr_strand,
                        \$ihm_perc_cov,\$tm_stable_id,\$tm_chr_name,
                        \$tm_chr_start,\$tm_chr_end,\$tm_chr_strand,
                         \$thm_perc_cov);
  
    my @split_gene_data = ();
    my $split_gene_chr; 

    while ($sth1->fetch) { 

      unless (defined $split_gene_chr) {
        $split_gene_chr = $tm_chr_name;
      }
      if ($ihm_perc_cov >$targeted_perc_cov) {# || $split_gene_chr ne $tm_chr_name) 
        @split_gene_data = ();
        $split_gene_chr = undef;
        last;
      }
      if (!$no_chr_constraint && $split_gene_chr ne $tm_chr_name) {
        @split_gene_data = ();
        $split_gene_chr = undef;
        last;
      } 

      push @split_gene_data, [$description,$im_stable_id,$im_chr_name,$im_chr_start,
                              $im_chr_end,$im_chr_strand,$ihm_perc_cov,$tm_stable_id,
                              $tm_chr_name,$tm_chr_start,$tm_chr_end,$tm_chr_strand,
                              $thm_perc_cov]; 

      $min_start = $tm_chr_start unless (defined $min_start);
      $min_start = $tm_chr_start if ($min_start > $tm_chr_start);
      $max_end = $tm_chr_end unless (defined $max_end);
      $max_end = $tm_chr_end if ($max_end < $tm_chr_end);
      $stable_ids{$im_stable_id} =1;
      $stable_ids{$tm_stable_id} =1;
  } 

  push @potential_stable_ids_for_recovery, $im_stable_id ; 

  if (scalar @split_gene_data) {
    printf "#%-5s %20s %5s %9s %9s %6s %3s %20s %5s %9s %9s %6s %3s\n", qw(desc stable_id chr start end strand cov stable_id chr start end strand cov);
    foreach my $gene_piece (@split_gene_data) {
      printf " %-5s %20s %5s %9d %9d %6s %3d %20s %5s %9d %9d %6s %3d\n",@{$gene_piece};
    }
    print "#----\n";
    $targeted_split_genes++;
  }
  %stable_ids = ();
  $min_start = undef;
  $max_end = undef; 

} 

  print "#Potentially $targeted_split_genes $targeted_species split genes\n";  

  $self->split_genes(\@potential_stable_ids_for_recovery) ; 
}  


sub run {  
  my ( $self ) = @_;  

  my @identified_stable_ids = @{$self->split_genes}; 
  my $informant_species = $$FIND_SPLIT_GENES{ANALYSIS_SETS}{$self->post_logic_name}{INFORMANT};   

  my $ga = Bio::EnsEMBL::Registry->get_adaptor($informant_species,'core','gene') ;    
  $self->output(get_longest_transcripts($ga, \@identified_stable_ids)) ; 
  print " have " . scalar (@{$self->output}) . " sequences to write\n" ;  

}


sub write_output{ 
  my ($self) = @_; 

  # sequences need to be proteins so let's translate them and check for stop codon before 

  my $written_files = $self->chunk_and_write_fasta_sequences( shuffle($self->output)) ;  
  $self->upload_input_ids( $written_files ) ;  
}


sub read_and_check_config {
  (my $self) = @_ ;   
  $self->SUPER::read_and_check_config(); 
}


sub split_genes { 
   my ( $self, $arg ) = @_; 
   $self->{split_genes} = $arg if $arg ; 
   return $self->{split_genes} ;  
}


# method has to go in Utils - used by MissingOrtholouges as well ... 

sub get_longest_transcripts{
  my ($ga, $stable_ids) = @_ ;
  my @cds ;

  GENES : for my $sid (@$stable_ids ) {
    my $gene = $ga->fetch_by_stable_id ($sid) ;
    my $maxlen = 0;
    my $nexonmax = 0;
    my $longest_transcript;

    foreach my $trans (@{$gene->get_all_Transcripts}) {
      my $len = 0;
      if ($trans->translation) {
        my @trans_exons = @{$trans->get_all_translateable_Exons()};
        foreach my $exon (@trans_exons) {
          $len+= $exon->length;
        }
        if ($len > $maxlen) {
          $maxlen = $len;
          $longest_transcript = $trans;
          $nexonmax = scalar($trans->get_all_Exons)
        }
      }
    }

    if (!defined($longest_transcript)) {
      print "No longest_transcript  transcript found for " . $gene->stable_id . "\n" ;
      next GENES ;
    }
    my $transstr = $longest_transcript->translate->seq;
    $transstr =~ s/(.{1,60})/$1\n/g;
    #print ">$sid\n$transstr\n\n" ; 
    if ($transstr =~m/\*/ ) {
      print "SKIPPING " . $longest_transcript->translation->stable_id .
        " : translation contain stop codon\n$transstr\n\n";
      next GENES ;
    }
    my $cds = $longest_transcript->translate ;
    $cds->display_id( $longest_transcript->translation->stable_id)  ;
    push @cds, $cds ;
  }
  my @rand_cds ;
  # randomize    
  return shuffle(\@cds) ;
}


1;
