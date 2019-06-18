=head1 LICENSE
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::MissingOrthologues - 

=head1 SYNOPSIS

my $orthologueanalysis = Bio::EnsEMBL::Analysis::RunnableDB::MissingOrthologues->new(
			      -analysis   => $analysis_obj,
			     );

$orthologueanalysis->fetch_input();
$orthologueanalysis->run();
$orthologueanalysis->output();
$orthologueanalysis->write_output(); 


=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::OrthologueEvaluator and is 
use warnings ;
used to fetch the input from different databases as well as writing results 
the results to the database.a


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::MissingOrthologues; 
use strict;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Gene; 

use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator; 
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;   
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw (
                                                                 get_one2one_orth_for_gene_in_other_species 
                                                                 get_transcript_with_longest_CDS
                                                                 ) ; 

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw (contains_internal_stops);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (shuffle ) ; 
use Bio::EnsEMBL::Registry; 
use vars qw(@ISA); 
use Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator; 

@ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator ) ; 


sub new {
  my ($class,@args) = @_;  
  my $self = $class->SUPER::new(@args);  
  $self->read_and_check_config ; 
  $self->verbose(1) ;  

  my $trusted_one2one_set= $$FIND_MISSING_ORTHOLOGUES{ANALYSIS_SETS}{$self->post_logic_name};    
  
  print "Searching for 1:1 orthoologues between $$trusted_one2one_set[0] and $$trusted_one2one_set[1]\n" ;  

  $self->species_1($$trusted_one2one_set[0]) ; 
  $self->species_2($$trusted_one2one_set[1]) ;  

  return $self;
}


# first analysis runs module MissingOrthologues and identifies possible targets, dumps them and 
# creates input_ids , than exonerate2genes wrapper runs with logich_name in Exoner2Genes config .... 
# --> do the logic_names in both configs need to be the same ? ?? 

# fetches the genes which have homologues in different species 

sub fetch_input {
  my( $self ) = @_;  

  # get all genes on slice for 1st species in array (species 1)
  $self->get_initial_geneset($self->species_1, $$FIND_MISSING_ORTHOLOGUES{DEFAULT_GENE_BIOTYPES}) ; 
}


sub run {
  my ($self) = @_; 

  my $query_spec = $$MAIN_CONFIG{QUERY_SPECIES};  
  my (%missing_orth,%trusted_orth); 

  QUERY_GENES : for my $gene_spec1 (@{$self->genes} ) { 
     #print " processing gene " . $gene_spec1->stable_id . "\n" ; 
     my $tg1_homol = get_one2one_orth_for_gene_in_other_species($gene_spec1 ,$self->species_2) ;    
      
     if ( $tg1_homol ) {  
       $trusted_orth{$gene_spec1->stable_id}=1; 
       foreach my $homolog_to_check  ( @{ $gene_spec1->get_all_homologous_Genes()} ) { 
         my ($check_homg, $check_homology, $check_species ) = @$homolog_to_check ;  
           if  ($check_species eq $query_spec ) { 
             #print "orth. between $query_spec and $$trusted_one2one_set[0] known\n" ; 
             next QUERY_GENES  ;
           }else { 
            $missing_orth{$gene_spec1->stable_id} =  $gene_spec1 ; 
           } 
         }
      } else { 
        #print " - no 1:1 homologue found for " . $self->species_2 . "\n" ;  
      }
  }
 
  for ( keys %missing_orth) { 
    print "missing_orthologue: $_\n"; 
  }  

  print scalar( keys %missing_orth)." missing orthologues identifed beteen " . $self->species_1 .  
        " and $query_spec ( using " . $self->species_2 ." as informant )\n"; 

   # dump sequences of identified missing orthologues   
  
   my @input_id_file_names ; 

   my @longest_transcripts;  
   for my $stable_id ( keys %missing_orth ) {    
     # check for stop codon  
     my ($ltr, $m,$n) = get_transcript_with_longest_CDS($missing_orth{$stable_id}) ;  
     unless (contains_internal_stops($ltr)){ 
       push @longest_transcripts, $ltr;
     } else { 
       print "Skipping transcript because it contains internal stops\n" ;  
     } 
   }
   my @cds ; 
   for my $tr ( @longest_transcripts ) {  
     my $cds = $tr->translate;
     $cds->display_id($tr->translation->stable_id ) ; 
     push @cds, $cds ; 
   } 
   $self->output(shuffle(\@cds))   ;
}  



sub write_output{
  my ($self) = @_;
  my $written_files = $self->chunk_and_write_fasta_sequences( $self->output) ; 
  if ( $written_files ) { 
    print scalar(@$written_files) . " fasta-files written for " . $self->species_1."\n" ;   
    $self->upload_input_ids( $written_files );   
    print STDERR "input_ids uploaded\n" ;   
  } else { 
    print STDERR "NO input_ids uploaded\n" ;    
 }
}




sub check_transcript {  
      my $tr = shift ;
      my $transstr = $tr->translate->seq;
      $transstr =~ s/(.{1,60})/$1\n/g;
      #print ">$sid\n$transstr\n\n" ; 
      if ($transstr =~m/\*/ ) {
        print "SKIPPING " . $tr->translation->stable_id .
        " : translation contain stop codon\n$transstr\n\n";
        return undef ; 
      } 
  return $tr ;  
}



sub filter_homologies {  
  my ( $all_homologies , $look_for_this_species) = @_ ;    
  my @result ;
 
  HOMOLOGIES :for my $homology ( @$all_homologies ) { 
     my @all_members = @{$homology->get_all_Members} ;
     # first object is source itself so don't process this 
     shift @all_members ;

     MA: foreach my $new_member (@all_members) { 
       my $species_name_of_orthologue = $new_member->genome_db->name ;    
       if ( $species_name_of_orthologue =~m/$look_for_this_species/) {  
         push @result, $homology ; 
       }
     } 
  }  
  return \@result ; 
} 







sub read_and_check_config {
  (my $self) = @_ ;  

  $self->SUPER::read_and_check_config(); 
    
  # - registry file is checked by OrthologueEvaluator 
  # - compara init file and schema-versions checked by OrtholgoueAnalysis 


  # check config for MissingOrthologues 
   my %config_hash = %{$FIND_MISSING_ORTHOLOGUES}; 
}







# methods which need to go into Utils : 

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
