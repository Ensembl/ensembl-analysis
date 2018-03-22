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

Bio::EnsEMBL::Analysis::RunnableDB::FindPartialGenes - 

=head1 SYNOPSIS

my $orthologueanalysis = Bio::EnsEMBL::Analysis::RunnableDB::FindPartialGenes->new ( 
			      -analysis   => $analysis_obj,
			     );

$orthologueanalysis->fetch_input();
$orthologueanalysis->run();
$orthologueanalysis->output();
$orthologueanalysis->write_output(); 


=head1 DESCRIPTION

This object uses information from a Ensembl Compara database to find partial 
gene predictions in an organism. It's  
use warnings ;
used to fetch the input from different databases as well as writing results 
the results to the database.a


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::FindPartialGenes;  

use strict;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Gene;  

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils 
 qw(get_transcript_with_longest_CDS get_one2one_orth_for_gene_in_other_species  ) ;  
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils; 
use Bio::EnsEMBL::Analysis::Tools::Utilities; 

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
  $self->{_partials}={};  
  return $self;
}



# fetches the genes out of the QUERY database 
# ( db which contains the predictions of species we're working on ) 

sub fetch_input {
  my( $self ) = @_; 

  # fetch all protein_coding genes from QUERY_SPECIES    
  
  $self->get_initial_geneset(
                              $$MAIN_CONFIG{QUERY_SPECIES}, 
                              $$FIND_PARTIAL_GENES{DEFAULT_GENE_BIOTYPES}
                            ); 
  print "SUMMARY : " . scalar( @{$self->genes}) . " genes fetched\n" ; 
}  


sub run {  
  my ( $self ) = @_;  

  # process all genes on this slice   
 
  # comparision species are defined in config file 
  my $qy_species  = $$MAIN_CONFIG{QUERY_SPECIES};  
  my $tg1_species = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{MAIN_ORTH}; 
  my $tg2_species = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{TEST_SPEC_1}; 
  my $tg3_species = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{TEST_SPEC_2}; 

  # loop trough QUERY genes in species (dog)
  QUERY_GENES : for my $query_gene (@{$self->genes} ) {    

    #
    # this code is mainly a cut-and-paste from Steve's find_partials.pl for mono3  
    #

    my $tg1_homol = get_one2one_orth_for_gene_in_other_species($query_gene,$tg1_species);

    if ($tg1_homol) {  

      # we found a ONE-2-ONE between QUERY (dog) and MAIN_ORTH (human) 
      # now we check if there are ortholgoues between human and test1 || human and test2 

      my $tg2_homol = get_one2one_orth_for_gene_in_other_species($tg1_homol,$tg2_species);
      my $tg3_homol = get_one2one_orth_for_gene_in_other_species($tg1_homol,$tg3_species); 

      if ($tg2_homol || $tg3_homol) {
        my ($longest_qy,$len_longest_qy,$nexon_longest_qy)    = get_transcript_with_longest_CDS($query_gene);
        my ($longest_tg1,$len_longest_tg1,$nexon_longest_tg1) = get_transcript_with_longest_CDS($tg1_homol);
        my ($longest_tg2,$len_longest_tg2,$nexon_longest_tg2) = get_transcript_with_longest_CDS($tg2_homol) if ($tg2_homol);
        my ($longest_tg3,$len_longest_tg3,$nexon_longest_tg3) = get_transcript_with_longest_CDS($tg3_homol) if ($tg3_homol);

        #print "Comparing lengths " . $qy_species . " $len_longest_qy " . $tg1_species . " $len_longest_tg1" .
               " " . $tg2_species ." " . (defined($tg2_homol) ? $len_longest_tg2 : "N/A") .
               " " . $tg3_species ." " . (defined($tg3_homol) ? $len_longest_tg3 : "N/A") . "\n";

        my $ratio_tg1_tg2 = 10;
        my $ratio_tg1_tg3 = 10;
        $ratio_tg1_tg2 = $len_longest_tg2/$len_longest_tg1 if ($tg2_homol);
        $ratio_tg1_tg3 = $len_longest_tg3/$len_longest_tg1 if ($tg3_homol);

        my $ratio_tg1_qy = $len_longest_qy/$len_longest_tg1;

        # if length_of_human_cds vs length_of_mouse_cds between 0.9 and 1.1 (+-10%)
        # OR if  length_of_human_cds vs legnght_of_monodelphis cds between 0.9 and 1.1 
        # and if predicted gene in query-species is less than 75% long than 1:1 orth in human

        # these values are set in the config 
        my $ratio_tg1_tg2_low  = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{RATIO_MAIN_vs_SPEC_1_LOW};
        my $ratio_tg1_tg2_high = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{RATIO_MAIN_vs_SPEC_1_HIGH};
        my $ratio_tg1_tg3_low  = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{RATIO_MAIN_vs_SPEC_2_LOW};
        my $ratio_tg1_tg3_high = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{RATIO_MAIN_vs_SPEC_2_HIGH};
        my $ratio_tg1_qy_75    = $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{RATIO_MAIN_vs_QUERY}; 
 


        if ($ratio_tg1_tg2 > $ratio_tg1_tg2_low  && $ratio_tg1_tg2 < $ratio_tg1_tg2_high ||
            $ratio_tg1_tg3 > $ratio_tg1_tg3_low  && $ratio_tg1_tg3 < $ratio_tg1_tg3_high )  {
          if ($ratio_tg1_qy < $ratio_tg1_qy_75 ) {
            print "Significant length difference for :" . $qy_species . "\t" . $longest_qy->stable_id .
                  "\t" . $query_gene->stable_id . " (" . ($query_gene->external_name ? $query_gene->external_name : '') . ")" .
                  " with " . $tg1_species . "\t" . $longest_tg1->stable_id .
                  "\t" . $tg1_homol->stable_id . " (" . $tg1_homol->external_name . ")\n";

            my $transcript_stable_id = $longest_tg1->stable_id;   

            if ( $$FIND_PARTIAL_GENES{ANALYSIS_SETS}{$self->post_logic_name}{IGNORE_ORTHOLOGUES_ON_MT_REGIONS} ) { 
	      if ( $tg1_homol->slice->seq_region_name =~m/MT/){  
	           print "ignore MT regions is switched on \n" ; 
	           print "\nIgnoring homolog : $transcript_stable_id 'cause it's on an MT region \n" ; 
	        } else { 
                   $self->partials($transcript_stable_id, $tg1_homol, $longest_tg1, $query_gene ) ; 
	        }
	    } else {
              $self->partials($transcript_stable_id, $tg1_homol, $longest_tg1, $query_gene ) ; 
	    }
          }
        }
      }
    }
  }
}



sub partials { 
  my ($self, $tsi, $tg1_homol, $longest_tg1_trans, $query_gene  ) = @_ ; 

  if ( $tsi ) {
    
    my $already_stored = 0 ; 
    if (!exists($self->{_partials}{$tsi})) {
       $self->{_partials}{$tsi}{gene} = $tg1_homol;
       $self->{_partials}{$tsi}{matches} = [];
       $self->{_partials}{$tsi}{transcript} = [];
    } else {   
       for my $longest_tr ( @{$self->{_partials}{$tsi}{transcript}} ) {  
         if ( $longest_tr eq  $longest_tg1_trans ) { 
            $already_stored = 1 ; 
         }
       } 
    } 

    unless ( $already_stored ) {   
      push @{$self->{_partials}{$tsi}{matches}}, $query_gene;
      push @{$self->{_partials}{$tsi}{transcript}}, $longest_tg1_trans; 
    }
  } 
  return $self->{_partials}; 
} 



sub write_output{ 
  my ($self) = @_; 

  my %identified_partials = %{ $self->partials } ;  

  # info output 
  foreach my $id (keys %identified_partials) {
    my $gene = $identified_partials{$id}{gene}; 

    
    print  "#desc             stable_id   chr     start       end strand cov".
            "            stable_id   chr     start       end strand cov\n";
    foreach my $match_gene (@{$identified_partials{$id}{matches}}) {
      printf  " %-9s %18s %3s %9d %9d     %2d   0 %18s %3s %8d %8d     %2d   0\n",
            "NA",$gene->stable_id, $gene->slice->seq_region_name, $gene->start, 
            $gene->end, $gene->strand,
            $match_gene->stable_id, $match_gene->slice->seq_region_name, $match_gene->start, 
            $match_gene->end, $match_gene->strand;
    }
    print  "#----\n";
  }  

  my @sequences_to_write ; 
 
  foreach my $id (keys %identified_partials) {   
    foreach my $longest_missing_trans (@{$identified_partials{$id}{transcript}}) {   
      #print $longest_missing_trans ."\t" . $longest_missing_trans->stable_id . "\n" ;  
      push @sequences_to_write, $longest_missing_trans; 
    }
    # this loops trough the genes in query-speices (new_species) 
    #foreach my $longest_missing_trans (@{$identified_partials{$id}{matches}}) {  
    #  print "TEST" . $longest_missing_trans ."\t" . $longest_missing_trans->stable_id . "\n" ; 
    #}
  }     

  # sequences need to be proteins so let's translate them  
  # and check for stop codon before 

  my @output ;  
  TRANSLATIONS: for my $trans ( @sequences_to_write ) {   

     if ( $trans->translation ) {  
        my $translate = $trans->translate; 
        my $transstr = $translate->seq() ;   
        $transstr =~ s/(.{1,60})/$1\n/g;
        if ($transstr =~m/\*/ ) {
          print "SKIPPING " . $trans->translation->stable_id . " : translation contain stop codon\n$transstr\n\n";
          next TRANSLATIONS ;  
        } 
       $translate->display_id($trans->translation->stable_id) ; 
       push @output, $translate ;   
     } 
  }   
  my $written_files = $self->chunk_and_write_fasta_sequences( shuffle(\@output) ) ;  
  $self->upload_input_ids( $written_files ) ;  
}


sub read_and_check_config {
  (my $self) = @_ ;   
  $self->SUPER::read_and_check_config(); 
    
  # - registry file is checked by OrthologueEvaluator 
  # - compara init file and schema-versions checked by OrtholgoueAnalysis 

  # check config for FindPartialGenes
  
   my %analysis_sets = %{$$FIND_PARTIAL_GENES{ANALYSIS_SETS}} ;  

   foreach my $logic_name ( keys %analysis_sets ) {  

     unless ( exists $$EXONERATE_CONFIG_BY_LOGIC{$logic_name} ) { 
       throw("You have defined a logic_name ".
           "\"$logic_name\"\n in our OrthologueEvaluator.pm".
           " configuration file but there is no configuration for such an analysis in the\n".
           " Exonerate2Genes-config, so i don't know where to write the genes to.\n".
           " I suggest to add a configuration for $logic_name to your ".
           "Exoneate2Genes.pm config\n") ; 
     }      

     my %subsets = %{$analysis_sets{$logic_name}};
     for my $species ( qw ( MAIN_ORTH TEST_SPEC_1 TEST_SPEC_2) ) { 
       print "analysing $subsets{$species}\n" ; 
       my $sa = Bio::EnsEMBL::Registry->get_adaptor($subsets{$species},'core',"Slice" ) ;       
       unless ( $sa ) {  
          throw ( " can't get adaptor for species $subsets{$species} - maybe typo in config ?"); 
       } 
     }

     my $analysis = $self->db->get_adaptor("Analysis")->fetch_by_logic_name($logic_name) ;  
     unless ( $analysis) {  
       throw("There's no analysis $logic_name in your database - this analysis is normally added by".
        " the setup script\n" ) ; 
     } 
     
   }
}




#
##
##
## methods which need to go into Utils : 
##
##
# 
#
#sub get_transcript_with_longest_cds { 
#  my ($ga, $stable_ids) = @_ ;
#  my @cds ;
#
#  GENES : for my $sid (@$stable_ids ) {
#    my $gene = $ga->fetch_by_stable_id ($sid) ;
#    my $maxlen = 0;
#    my $nexonmax = 0;
#    my $longest_transcript;
#
#    foreach my $trans (@{$gene->get_all_Transcripts}) {
#      my $len = 0;
#      if ($trans->translation) {
#        my @trans_exons = @{$trans->get_all_translateable_Exons()};
#        foreach my $exon (@trans_exons) {
#          $len+= $exon->length;
#        }
#        if ($len > $maxlen) {
#          $maxlen = $len;
#          $longest_transcript = $trans;
#          $nexonmax = scalar($trans->get_all_Exons)
#        }
#      }
#    }
#
#    if (!defined($longest_transcript)) {
#      print "No longest_transcript  transcript found for " . $gene->stable_id . "\n" ;
#      next GENES ;
#    }
#    my $transstr = $longest_transcript->translate->seq;
#    $transstr =~ s/(.{1,60})/$1\n/g;
#    #print ">$sid\n$transstr\n\n" ; 
#    if ($transstr =~m/\*/ ) {
#      print "SKIPPING " . $longest_transcript->translation->stable_id .
#        " : translation contain stop codon\n$transstr\n\n";
#      next GENES ;
#    }
#    my $cds = $longest_transcript->translate ;
#    $cds->display_id( $longest_transcript->translation->stable_id)  ;
#    push @cds, $cds ;
#  }
#  my @rand_cds ;
#  # randomize    
#  return shuffle(\@cds) ;
#}
#
#



1;
