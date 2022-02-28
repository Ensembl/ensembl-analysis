=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::FindFalseIntrons - 

=head1 SYNOPSIS

  my $fs = Bio::EnsEMBL::Analysis::RunnableDB::FindFalseIntrons->new(
			      -analysis   => $analysis_obj,
#			     );

  $fs->fetch_input();
  $fs->run();
  $fs->output();
  $fs->write_output(); 


=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Analysis::Runnable::OrthologueEvaluator and is 
 used to fetch the input from different databases as well as writing results 
 the results to the database.a

 Configuration used to run this module is contained in the following files  :

 Bio::EnsEMBL::Analysis::Config::GeneBuild::ExamineGeneSets;
 Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator;

 The analysis FindFalseIntrons uses input-ids of the following format : 

    Coordsystem:Assembly:Seq_region_name:Start:End:Strand
 
            chromosome:BROADD2:12:1333:890890:1 

            chromosome:BROADD2:3:30434100:30529400:1

=head1 METHODS


=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a '_'

=cut



package Bio::EnsEMBL::Analysis::RunnableDB::FindFalseIntrons; 
use warnings ;
use strict;  

use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExamineGeneSets; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator; 
use Bio::EnsEMBL::Analysis::Config::Databases;

use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
# Not sure if this module ever existed
#use Bio::EnsEMBL::Analysis::RunnableDB::ExamineGeneSets;
use Bio::EnsEMBL::Registry; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils 
         qw(
            get_transcript_with_longest_CDS
            get_one2one_homology_for_gene_in_other_species
            get_one2one_orth_for_gene_in_other_species
           );

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::HomologyUtils 
         qw( 
            get_gene_obj_out_of_compara_homology_object
           );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
         qw(
            list_evidence 
            convert_translateable_exons_to_exon_extended_objects
           ); 

use vars qw(@ISA);  

@ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
            Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator ) ; 


sub new {
  my ($class,@args) = @_;  
  my $self = $class->SUPER::new(@args);   
  $self->read_and_check_config();  
 # $self->verbose(50000) ;  

  $self->orthologues_species( $$FIND_FALSE_INTRONS{ANALYSIS_SETS}{$self->analysis->logic_name} ) ; 
  #print join ("\n" , @{$$FIND_FALSE_INTRONS{ANALYSIS_SETS}{$self->analysis->logic_name}} ) ; 
  return $self;
}




  # inspect all introns of the QUERY_SPECIES of genes with certain 
  # biotypes ( see OrthologueEvaluator.pm config )  

sub run { 
  my ($self) = @_;  

  my @genes = @{$self->genes} ;    

  my $cnt =0 ;
  my %input_ids_to_upload ; 
  my %transcripts_marked_for_deletion; 

  my %intron_attributes; 
  GENES: for my $g ( @genes ) {          

     $cnt++;
     print "\n\nprocessing gene-nr  $cnt / ".scalar(@genes) . " : " . $g->stable_id . "\n" ; 
     print "="x80; print "\n";
 
     # we only use the translation length of the transcript with the longest 
     # CDS - otherwise statistics will be slightly different for genebuilds
     # where genes have lots of transcripts. we also only count the nr of exons 
     # of the transcript with the longest CDS  
      
     my ($t_long,$cdsl, $n_exon_max)  = get_transcript_with_longest_CDS($g); 
 
     unless ( $t_long) {  
      print " longest transcript does not exist\n" ;  
      next GENES; 
     }     

     my ( $unique_introns , $i2t ) =  get_unique_introns_of_transcripts ( $g ) ;  
     my %unique_introns = %{$unique_introns} ; ;
     my %i2t = %{$i2t};   

     # process a list of non-redundant introns of all transcripts seperately and independent from genes / transcripts  
    
      $self->{_false_introns_alt_tr} = {} ; 
      $self->{_ftr_false_introns} = {} ; 
      $self->{_coding_ex_in_alt_tr_which_aligns_false_intron} = {};  

     my %intron_is_coding_in_homologues ; 

     for my $intron ( values %unique_introns ) {  

      if ( $intron->length < 150 ) {  

         print "\nINTRON : " . $intron->length . "< 150\n========================================\n";
 
         my $aligned_homolog_trans = $self->exonerate_intron_vs_homologs( $intron, $g, $i2t);


         if ( scalar( @$aligned_homolog_trans) > 0 ) {  
           # intron aligned against homologues protein so make input_id to recover intron 

           print "INTRON_IS_CODING_IN_HOMOLOGUES  : " . join ("," , map { $_->stable_id }  @$aligned_homolog_trans ) . "\n\n" ; 
           $intron_is_coding_in_homologues{$intron} = $aligned_homolog_trans ;  

           # recover id is : 
           #  chromosome:Btau_3.1:Un:287653469:287755066:1:other_protein:Q5W135.1:203:ENSP00000354519:345781 
 
           my @exonerateVSGenwise_ids = @{$self->get_recover_ids( $g,$t_long, $aligned_homolog_trans)} ; 

           # flag the recovered intron for the recovery with RecoverFalseIntrons.pm 

           my @iid = split /\:/, $self->input_id ;  
           my $my_slice = $self->db->get_SliceAdaptor->fetch_by_region(undef,$iid[2]) ; 

           # create simple_feature ( we could also change the method how we align the intron ( i.e. use exonerate 
           # out of a runnableDB and than use the parser as well and store the whole alignemnt in the dna_align_feature table )

           my $feature = Bio::EnsEMBL::SimpleFeature->new(  
                                                       -seq_region_start =>  $intron->seq_region_start , 
                                                       -seq_region_end   =>  $intron->seq_region_end,  
                                                       -start =>  $intron->seq_region_start , 
                                                       -end   =>  $intron->seq_region_end,  
                                                       -strand=>  $intron->seq_region_strand, 
                                                       -analysis => $self->analysis() , 
                                                       -score   => 0, 
                                                       -display_label => $intron->prev_Exon->stable_id."_".$intron->next_Exon->stable_id , 
                                                       -slice => $my_slice, 
                                                      ) ;    

             # put input-ids / recovery on hold until we checked against alt transcript if we need to recover or if 
             # FI is coding in alt trans 
             $intron_attributes{$intron}{feature} = $feature ; 
             $intron_attributes{$intron}{input_ids} = $self->get_recover_ids( $g,$t_long, $aligned_homolog_trans) ;  

#             $sfa->store($feature) ;      
#             for my $exVSgen_id ( @exonerateVSGenwise_ids ) {   
#               $input_ids_to_upload{$exVSgen_id.":".$feature->dbID}=1;     
#               $intron_attributes{$intron}{input_id} = $exVSgen_id ; 
#              print "FEATURE : $intron \t $exVSgen_id:".$feature->dbID . "\n" ; 
#             }  



         } else { 
            print "no hit -  intron did not align\n" ; 
         }
         print "\n\n" ;    
      } else {  
       # // length < xxx  
       #print "intron too long \n" ; 
      }
    }  # // unique introns  

      #
      # NOW we make decision if we want to flag for recovering, or if we flag as 'alternative transcript with introns which are coding in
      # alt. transcript and in homolog' - which basically mean it's a low quality one ...... 
      #

   
              print "Summary:\n=====================\n" ;    

              my %alt_trans_ex_coding; 
              my %false_introns = %{ $self->{_false_introns_alt_tr} } ;   
              for my $k ( keys %false_introns ) {   
                 for my $tr (  @{$false_introns{$k}} ) {   
                     $alt_trans_ex_coding{$tr->stable_id}++ ; 
                     #print $tr->stable_id ."\n"; 
                 } 
              } 
              print "\n" ;  
              foreach my $tr_coding_stable_id ( keys %alt_trans_ex_coding ) {  
                print "ALT_CODING_TR : $tr_coding_stable_id\t$alt_trans_ex_coding{$tr_coding_stable_id} introns recovered\n" ; 
              } 
              print "\n" ;  
              # false transcript->stable_id -   intron intron intron 
      
              my %false_transcripts= %{ $self->{_ftr_false_introns} } ;    
        
              for my $f_tsi ( keys %false_transcripts) {    

                 print "FalseTranscript " . $f_tsi ."\t: ".scalar(@{$false_transcripts{$f_tsi}})." possible_false_introns.\n" ;   

                 my $nr_false_introns_recovered = 0 ;  

                 for my $false_intron  (  @{$false_transcripts{$f_tsi }} ) {    
                    print "\tFI coding in HOMOLOGUES: ".join ("," , map { $_->stable_id }  @{$intron_is_coding_in_homologues{$false_intron}}) . "\n";
                    print "\tFI coding in alt.Trans : ".join (" ,", map { $_->stable_id } @{${$self->{_false_introns_alt_tr}}{$false_intron}} ) ;    

                    my %fi = %{ $self->{_coding_ex_in_alt_tr_which_aligns_false_intron}} ; 
                    print "\tFI-start : " . $fi{$false_intron}{false_intron}->seq_region_start . "\t" ; 
                    print " FI-end  : " .  $fi{$false_intron}{false_intron}->seq_region_end . "\t" ;  
                    print "\toverlapping_coding_ex_in_alt_tr: ". $fi{$false_intron}{overlapping_coding_exon}->stable_id  . "\n\n" ;    
                    $nr_false_introns_recovered++;   
                 }

                 for my $false_intron  (  @{$false_transcripts{$f_tsi }} ) {      
                   my $delete_recovery_id = 0 ;  

                   if (    (scalar( @{$false_transcripts{$f_tsi }}) == $nr_false_introns_recovered)   
                        && (scalar( @{$intron_is_coding_in_homologues{$false_intron}}) == scalar(@{$self->orthologues_species}) ) )  {  
                      # all FalseIntrons are coding in alternative transcript as well as in both Orthologues  
                      print "marked_for_deletion_alt_coding_and_coding_in_homologues : $f_tsi\n" ;  
                      #$transcripts_marked_for_deletion{coding_in_alt_tr}{$f_tsi}=1;
                      #$delete_recovery_id = 1 ; 

                   } elsif ( scalar (  @{$false_transcripts{$f_tsi }} ) == $nr_false_introns_recovered ) {    
                      # all FalseIntrons are recovered in alt trans 
                      print "transcript could be marked for deletion as it's coding in alt. trans : $f_tsi\n" ;   
                      #$transcripts_marked_for_deletion{coding_in_alt_tr}{$f_tsi}=1;
                      #$delete_recovery_id = 1 ; 

                   } else { 
                      foreach my $intron_key ( keys %intron_attributes  )  {    
                        print "intron_key $intron_key\n" ; 
                        my $simple_feature = $intron_attributes{$intron_key}{feature};  
                        print "simple_feature $simple_feature\n" ; 
                        my $feature = $intron_attributes{$intron_key}{feature}; 
                        my @recover_ids = @{ $intron_attributes{$intron_key}{input_ids} } ;   
                         for my $rci ( @recover_ids ) {  
                           print " recover_id $rci\n" ; 
                         } 
                      }  
                   }
                   if (  $delete_recovery_id ) {  
                      print "deleteing recovery_ids as FalseIntron is coding:\n" ; 
                      delete $intron_attributes{$false_intron}{feature};
                      delete $intron_attributes{$false_intron}{input_ids} ;  
                   } 
                   
                 }  

              } 

              
      #my $nr_of_false_introns_in_transcript =  keys %false_introns ; 
      #my $nr_coding_ex = keys %{ $self->{_coding_ex_in_alt_tr_which_aligns_false_intron}} ;    
#
#      print "FI in trans : $nr_of_false_introns_in_transcript \n" ;  
#      print "CODING EX   : $nr_coding_ex\n" ;  
#
#      if ( $nr_of_false_introns_in_transcript == $nr_coding_ex ) {  
#        # all false introns are recovered in alt. transcript - gene can be flagged for deletion,  
#          print "Flagged for deletion : " ; 
#        # don't upload input_ids to recover 
#      } 
#

  }  # //  GENES    

  my $sfa = $self->db->get_SimpleFeatureAdaptor() ;   

 # upload SimpleFeatures for FalseIntrons with coordinates 

  

  for my $false_intron  ( keys %intron_attributes ) { 
 
    my $feature = $intron_attributes{$false_intron}{feature}; 
    if ( $feature ) { 
      $sfa->store($feature) ;       
      for my $exVSgen_id ( @{ $intron_attributes{$false_intron}{input_ids} }){  
        #print FEATURE : $feature->dbID  . "\n" ; 
        $input_ids_to_upload{$exVSgen_id.":".$feature->dbID}=1;     
        #$intron_attributes{$false_intron}{input_id} = $exVSgen_id ; 
      }   
    }
  }  

  #
  # upload stable-id's of alt. transcripts which have been flagged for removal 
  # 
  for my $set ( keys %transcripts_marked_for_deletion ) {    
    
    for my $tsi ( keys %{$transcripts_marked_for_deletion{$set}}) { 
      my $rm_transcript = $self->db->get_TranscriptAdaptor->fetch_by_stable_id($tsi) ;  

     my $feature = Bio::EnsEMBL::SimpleFeature->new( 
                                              -seq_region_start =>  $rm_transcript->seq_region_start , 
                                              -seq_region_end   =>  $rm_transcript->seq_region_end,  
                                              -start =>  $rm_transcript->seq_region_start , 
                                              -end   =>  $rm_transcript->seq_region_end,  
                                              -strand=>  $rm_transcript->seq_region_strand, 
                                              -analysis => $self->analysis() , 
                                              -score   => 0, 
                                              -display_label => $rm_transcript->stable_id."_$set", 
                                              -slice => $rm_transcript->slice,  
                                              ) ;    
    $sfa->store($feature) ;       
    print $feature->display_id() . " stored for further investigation in SimpleFeature table \n" ; 
    } 
  } 

  # NOW check analysis and upload input_ids  
  
  # THIS CAN ALL GO INTO WRITE_OUTPUT 
     
  # check if we got a post-analysis stored in the db otherwise set one up 
       
        my $analysis_obj = new Bio::EnsEMBL::Pipeline::Analysis(
                                                -logic_name => "RecoverFalseIntrons",
                                                -program=> "exonerate",
                                                -module => "RecoverFalseIntrons",
                                                -input_id_type => 'prot_slice_form',
                                               );  
        my $submit_ana = new Bio::EnsEMBL::Pipeline::Analysis(
                                              -logic_name => "Submit_RecoverFalseIntrons", 
                                              -module => "Dummy",
                                              -input_id_type => 'prot_slice_form', 
                                             );  
        my $pa = $self->pipeline_adaptor($self->db) ; 
        my $aa = $pa->get_AnalysisAdaptor; 
        $aa->store($submit_ana);
        $aa->store($analysis_obj); 

        # set up rules / conditions for post analysis     
        print "input_ids:\n" ;   
        print join("\n", keys %input_ids_to_upload)."\n\n" ;   
  
        print "reformatting input_ids as we use the new format:\n" ;   

        my @old_input_ids = keys %input_ids_to_upload; 
        my @new_input_ids = @{reformat_input_ids(\@old_input_ids)} ; 

        print "new_input_ids:\n" ;    
        for my $nid ( @new_input_ids ) { 
           print "new_input_id " . $nid . "\n" ; 
           print "new_input_id_length  " . length($nid) . "\n" ; 
        } 
        # $self->upload_input_ids( [keys %input_ids_to_upload] , "RecoverFalseIntrons" , "prot_slice_form") ;  
        $self->upload_input_ids( \@new_input_ids , "Submit_RecoverFalseIntrons" , "prot_slice_form") ;  
        my $goal_analysis = $aa->fetch_by_logic_name("RecoverFalseIntrons") ; 
        my $rule = Bio::EnsEMBL::Pipeline::Rule->new(-goalanalysis => $goal_analysis); 

        # check if rule has already been stored ...
        my $ra = $self->pipeline_adaptor($self->db)->get_RuleAdaptor(); 
        my $store_rule = 1 ; 
      RULE: foreach my $old_rule( $ra->fetch_all ) {  
        if ( $goal_analysis->logic_name eq $old_rule->goalAnalysis->logic_name ) {  
          $store_rule = 0 ;
        }
      }  
      if ( $store_rule ) {  
        $ra->store($rule) ; 
      }  

    Bio::EnsEMBL::Registry->load_all($$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE},1,1);
    my @dba = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors()}; 
    for ( @dba ) {  
        $_->disconnect_when_inactive(1);
    } 

}  





sub reformat_input_ids { 
  my ($input_id_aref ) = @_ ;  

  my @new_ids ; 
  my %tmp; 
 
  for my $id ( @{$input_id_aref} ) {   

    #my ($cs,$asm,$sname,$start,$end,$strand,$query_stable_id,$hom_stable_id,$sf_id) = split /\:/,$id ;       
    my @components = split /\:/,$id ;       
    my $sf_id = pop  @components ;
    my $hom_stable_id = pop @components ; 
    my $key = join ( ":", @components ) ;  
    push @{$tmp{$key}{$hom_stable_id}}, $sf_id ;    
  }       
   
  for my $region_key ( keys %tmp ) {   
      my %sf_id_index ; 
      print "processing region and transcript : ";  
      print $region_key . "\n" ; 
      my %hom_to_sf_id = %{$tmp{$region_key}} ;  

      for my $stable_id ( keys %hom_to_sf_id ) { 
        print $stable_id . "\t" ; 
        my @sf_ids = @{$hom_to_sf_id{ $stable_id}};
        for my $sf_id ( @sf_ids ) {  
          print "$sf_id\t" ; 
          $sf_id_index{$sf_id} = 1; 
        } 
        print "\n" ; 
      }  
      my $sf_db_id_string = ""; 
      my $counter = 0 ; 
      for my $id ( sort keys %sf_id_index ) {   
        $sf_db_id_string .= $id .",";
        $sf_id_index{$id}=$counter ;  
        print $id . " --> " . $counter . "\n" ; 
        $counter++;   
      }  
      my $hom_string ="";
      my $sf_string =""; 
      for my $stable_id ( keys %hom_to_sf_id ) {  
        $hom_string .=$stable_id."," ; 
        for my $sf_id ( sort @{$hom_to_sf_id{ $stable_id}} ) {    
            print "\t$sf_id\t" ; 
            $sf_string .= $sf_id_index{$sf_id}. ",";
        }
        print "\n"; 
        $sf_string =~s/,$//; 
        $sf_string .= "_";
      } 
      $sf_string=~s/_$//; 
      $sf_string=~s/,$//; 
      $sf_db_id_string=~s/,$//; 
      $hom_string =~s/,$//; 
      push @new_ids , $region_key . ":" . $hom_string .":" . $sf_string . ":" . $sf_db_id_string ; 
  } 
  return \@new_ids ;  
}

sub exonerate_intron_vs_homologs {    
    my ($self,$intron_to_check, $gene, $i2t, $transcript_of_intron ) = @_ ;  
     
     my %intron_aligns_in_x_homologues ; 

     my %i2t = %{$i2t};  
     my @aligned_homolog_transcripts ;    
     my @homo_genes =  @{ $self->get_homologues_genes($gene) }  ;   
     if ( scalar(@homo_genes) == 0 ) {  
        print "no homologs found\n" ; 
     } 
     HOMOLOG: for my  $homolog_gene ( @homo_genes  ) {    


         print $homolog_gene . "\t" . $homolog_gene->stable_id . "\n" ; 

         # build array of coding exons for homologues transcript  
         my @introns_in_homolog_gene ; 
         TRANS: for my $ht ( @{$homolog_gene->get_all_Transcripts } ) {   
           push @introns_in_homolog_gene,  @{$ht->get_all_Introns}; 
         }  



         TRANS: for my $homology_trans ( @{$homolog_gene->get_all_Transcripts } ) {  
         
           my $transcript_of_intron  = $i2t{$intron_to_check};  
           my $intron_id = $transcript_of_intron->stable_id."_".$intron_to_check->prev_Exon->stable_id."_".$intron_to_check->next_Exon->stable_id; 
         
           my $query_length = exonerate_sequences($intron_to_check, $homology_trans , $intron_id, $transcript_of_intron);    
                
            
           if ( $query_length > 0 ) {   
             push @aligned_homolog_transcripts, $homology_trans ; 

               #
               # check if the intron which aligned to the CDS-exon of one of the homologs is 
               # coding sequence in one of the other alternatively spliced transcripts 
               # of the gene to investigate ? or of the other homologs ? 
               #
               
               TRANSCRIPTS: for my $tr ( @{$gene->get_all_Transcripts} ) {      
                  my $tr2 = $transcript_of_intron->stable_id ;  

                  if ( $tr->stable_id =~m/$tr2/) { 
                    next TRANSCRIPTS ;   
                  }  

                  #my $query_length = exonerate_sequences( $intron_to_check, $tr, $intron_id,$transcript_of_intron ) ;   
                  
                  # compare coordinates of introns of one transcript and exons of the alternative transcript : 
                  my @exons = @{$tr->get_all_Exons} ;                    
                  for my $e ( @exons ) {  
                    if ( $e->overlaps($intron_to_check)){  
                      print "\nExon in alternative transcript overlaps FalseIntron - Intron will not be recovered.\n";  
                      print "Exon    : " . $e->stable_id ." ( " . $tr->stable_id .  ")  overlaps_intron in $intron_id\n" ; 
                      print "Positive: "  . $tr->stable_id . " vs " . $transcript_of_intron->stable_id . "\n" ;    

                      # check if the false intron is maybe coding in one of the alternative transcripts
                      $self->intron_aligns_alternative_transcript($intron_to_check, $transcript_of_intron, $e, $tr, $homology_trans ) ;  
                     
                    } 
                  }  
               }


               # OK now identify the intron -vs- coding-exon-in-homolog pair ... 
               # get all coding exons of homolog 
               my %hom_coding_alinging_exons = %{  exonerate_intron_vs_exons($intron_to_check, $homology_trans, $intron_id, $transcript_of_intron) } ; 
               print "exonerating done \n" ; 
               for my $exon_stable_id ( keys %hom_coding_alinging_exons ) {  
                  my ( $exon, $query_length ) = @{$hom_coding_alinging_exons{$exon_stable_id}} ; 

                  for my $hom_trans_intron ( @introns_in_homolog_gene ) {  
                     if ( $hom_trans_intron->overlaps($exon) ) { 
                         print "\nXXXXXXXX\nIntron  overlaps " . $exon->stable_id . " - better not recover this as it's \nXXXXXXXXX\n" ;  
                         print $hom_trans_intron->prev_Exon->stable_id . "   " . 
                               $hom_trans_intron->next_Exon->stable_id . "   " . 
                               $exon->stable_id . "\n" ;  
                               print "\nXXXXXXXXXXXXXX\n" ; 
                     } 
                  }
               } 
         } 
       } # // foreach transcript 
     } # // foreach homolog
   return (\@aligned_homolog_transcripts)  ;        
}  

# report 


sub intron_aligns_alternative_transcript { 
    my ($self, $false_intron, $false_transcript, $overlapping_exon , $alt_transcript, $homology_trans ) = @_ ;   

    # Report for Gene X : 
    # - has xx 1:1 homologues with yy alt. spliced transcripts 
    # - has xx false-introns 
    # - false-introns are recovered in tr1, tr2 
    # - are All false-introns recovered in one transcript ?  
    # - nr of transcripts with false introns $transcript-> @{introns}  
    # - ${intron} -> recovered in false Intron  
  
   
    ${$self->{_coding_ex_in_alt_tr_which_aligns_false_intron}}{$false_intron}{false_intron}= $false_intron; 
    ${$self->{_coding_ex_in_alt_tr_which_aligns_false_intron}}{$false_intron}{overlapping_coding_exon}= $overlapping_exon; 

    push @{${ $self->{_false_introns_alt_tr}}{$false_intron}}, $alt_transcript; 
    push @{${ $self->{_intron_coding_in_hom}}{$false_intron}}, $homology_trans; 
   
   # my %tmp = %{$self->_false_introns_alt_tr} ;  
   # print join ( "\n" , keys %tmp ) ; 
   # print $tmp{$false_intron} ; 

    push @{ ${$self->{_ftr_false_introns}} {$false_transcript->stable_id} }, $false_intron ; 
    
} 




sub exonerate_sequences {  
  my ($intron_to_check , $target, $intron_id, $transcript_of_intron,$version) = @_ ;   
    
    unless ( $intron_id) {  
       $intron_id = $intron_to_check->display_id ; 
    } 
    my $out_dir = "/tmp" ; 
    my $f1= $out_dir."/".$intron_id.".query"  ; # intron seq ( short seq )  
    my $f2= $out_dir."/".$target->stable_id .".target";   # gene seq/ genome seq
    my $f3= $out_dir."/".$target->stable_id .".exonerate";   
 
  
    open (FH,">$f1") || die "Cant write file $f1\n" ;    
       # this was set before in one of the subroutines  
       # may need to worry if this gets redundant / if it's unset .......     
       
       #  write INTRON to file
       
       my $seq = $intron_to_check->seq; 
       $seq=~ s/(.{1,60})/$1\n/g;
       print FH ">Intron $intron_id\n" . $seq;  
       close(FH) ;   



       #  write TRANSCRIPT to file
       open (FH,">$f2") || die "Cant write file $f2\n" ;   
           #($seq=$target->seq->seq )=~ s/(.{1,60})/$1\n/g;
           ($seq=$target->translateable_seq )=~ s/(.{1,60})/$1\n/g;
           print FH ">".$target->stable_id . "\n$seq";  
       close(FH) ;    

        my @versions = qw ( 1.0.0 ) ;  
        my @options = ( "   " , " --exhaustive " ) ; 

        my $query_length =0  ;   
        my $cmd ;
 
        for my $version ( @versions ) {   
        for my $extra_opt ( @options ) {    
           my $exonerate_version ; 

           if ( $version =~m/0.7.1/ || $version=~m/1.0.0/) { 
              $exonerate_version = "/usr/local/ensembl/bin/exonerate-$version" ;      
           }  else  { 
              $exonerate_version = "/software/ensembl/bin/exonerate-1.4.0 " ;      
           }
           my @alignment ; 


           $cmd = $exonerate_version .  " -q $f1 -t $f2 -m affine:local --bestn 1 ";   
           $cmd .= $extra_opt ; 
           $cmd.=" > $f3"; 
           system("$cmd") ;   
           my $hit = 0 ; 
           my $raw_score ; 
           open (F,"$f3") || die " can't read $f3\n" ;   
           foreach my $l (<F>){     
             $hit = 1 if ($l=~/C4 Alignment:/) ; 
             #$hit = 0 if $l=~/vulgar/;  
             
             push @alignment , $l  if $hit ;  
                 
             if ( $l=~m/Query range:/){
               my @it = split/\s+/,$l ; 
               if ( $it[5]>$it[3]) { 
                 $query_length = ($it[5] - $it[3])+1  ;  
               } else { 
                 $query_length = ($it[3] - $it[5])+1  ;  
               } 
             }
           } 
  
          close(F);    

          #  print info 

          if ( $hit )  {            
             #print "\n$f1 \n $f2 \n $f3\n" ; 
              print "hit_with $version $extra_opt \n$cmd" ; 
             print "\nIntron_exonerate_alignment_OK :  " . $transcript_of_intron->stable_id." [ ".  $intron_to_check->prev_Exon->stable_id." <intron>  "
              .$intron_to_check->next_Exon->stable_id . "]  ALIGNS " . $target->stable_id ." ]  QLEN $query_length iLEN " . $intron_to_check->length."\n\n" ;   

             print $intron_id . "aligns perfect \n" if ( $intron_to_check->length == $query_length ); 

             #if ( $query_length % 3 == 0 ){ 
               # print "  alignment_length_mod_3\n" ; 
             #}else{ 
               # print "  alignment_frameshift\n"; 
             #}  
        
             #print "\n\n".join("",@alignment) ;    
            system("rm $f1"); system("rm $f2"); system("rm $f3");   
            return $query_length ; 
          }    else {  
           # print "no_hit \n" ; 
          } 
          # jhvx 
      } # next version  
      } # next option 
          #system("rm $f1"); system("rm $f2"); system("rm $f3");   
          return $query_length ; 
}  


sub exonerate_intron_vs_exons {  
  my ($intron_to_check , $homology_trans, $intron_id, $transcript_of_intron  ) = @_ ;   

   print "INTRON_EXONERATING to identify CODING exons in homolog.....................\n";
    
    unless ( $intron_id) {  
       $intron_id = $intron_to_check->display_id ; 
    } 
    my $out_dir = "/tmp";
    my $f1= $out_dir."/".$intron_id.".query"  ; # intron seq ( short seq )  

    my @alignment ;  

    # write INTRON seq to file 
  
    open (FH,">$f1") || die "Cant write file $f1\n" ;    
    my $seq = $intron_to_check->seq; 
    $seq=~ s/(.{1,60})/$1\n/g;
    print FH ">Intron $intron_id\n" . $seq;  
    close(FH) ;  

    my %aligning_hom_exons ; 

    my @homology_exons = @{$homology_trans->get_all_translateable_Exons};

    for my $hom_exon ( @homology_exons ) {   

      my $f2= $out_dir."/".$hom_exon->stable_id .".homology_trans";   # gene seq/ genome seq
      my $f3= $out_dir."/".$hom_exon->stable_id .".exonerate";   
 
    # write HOMOLOGS CODING EXON seq to file  
    
      open (FH,">$f2") || die "Cant write file $f2\n" ;   
        ($seq=$hom_exon->seq->seq )=~ s/(.{1,60})/$1\n/g;
        print FH ">".$hom_exon->stable_id . "\n$seq";  
      close(FH) ;    
      # used for first comparision  

      my $exonerate_version ="/usr/local/ensembl/bin/exonerate-1.0.0" ; 
      my $cmd = "$exonerate_version  -q $f1 -t $f2 -m affine:local --bestn 1 ";  
      $cmd .=" > $f3";   
      system("$cmd") ;   
      
      my $hit = 0 ; 
      my $query_length =0  ;  
      my $raw_score ;  

      open (F,"$f3") || die " can't read $f3\n" ;   
        foreach my $l (<F>){     
          if ($l=~/C4 Alignment:/) {
             $hit = 1;  
             $l = "INTRON vs CDS-EXON-HOMOLOG " . $l ;  
             print "\nINTRON vs CDS-EXON-HOMOLOG\n"; 
          } 
          push @alignment , $l  if $hit ;  
                 
          if ( $l=~m/Query range:/){
            my @it = split/\s+/,$l ; 
            if ( $it[5]>$it[3]) { 
              $query_length = ($it[5] - $it[3])+1  ;  
            } else { 
                $query_length = ($it[3] - $it[5])+1  ;  
            } 
          }
        } 
      close(F);    

      if ( $hit )  {           
          print "\ntest-intron aligns vs coding exon : " . $transcript_of_intron->stable_id.
             " [ ".  $intron_to_check->prev_Exon->stable_id." <intron>  " .$intron_to_check->next_Exon->stable_id . "]  ALIGNS " 
              . $hom_exon->stable_id ." ]  QLEN $query_length iLEN " . $intron_to_check->length."\n\n" ;   

             #print "\n\n".join("",@alignment) ;    
           $aligning_hom_exons{$hom_exon->stable_id} = [ $hom_exon , $query_length] ; 
      }else {  
        print "intron does not match exon " . $hom_exon->stable_id . "\n" ;  
      } 
          system("rm $f2"); system("rm $f3");    
   } 
   system("rm $f1"); 
  return \%aligning_hom_exons ; 
}  




sub get_unique_introns_of_transcripts { 
  my ( $g ) = @_;  

      my (%unique_introns, %i2t) ; 

      for my $tr ( @{$g->get_all_Transcripts}){    

        my @tr_exons =  @{ $tr->get_all_translateable_Exons};  
       
        # not sure if this is needed for exons but it's definitely needed for introns ! 
        if ( $tr->seq_region_strand == 1) {  
           @tr_exons = sort { $a->seq_region_start <=> $b->seq_region_start } @tr_exons ; 
        } else { 
           @tr_exons = sort { $b->seq_region_start <=> $a->seq_region_start } @tr_exons ;
        } 

         my $intron_count = 0 ;  
         foreach  ( my $i=0; $i < scalar(@tr_exons) ; $i++) {    
            
           if ( $i+1 < scalar(@tr_exons) ){  
              
              $intron_count++; 
              my $intron =  new Bio::EnsEMBL::Intron( $tr_exons[$i], $tr_exons[$i+1] ) ;    
              $i2t{$intron} = $tr;  
             
    #          push @{ $i2t{$intron} } , $tr; 
              
              my $intron_hkey = $intron->slice->name."-".$intron->seq_region_start."-".
                                   $intron->seq_region_end."-".$intron->seq_region_strand;  

              unless ( exists $unique_introns{$intron_hkey}) {  
                $unique_introns{$intron_hkey}=$intron ;  
              }
         } # all (nr_exons - 1 ) introns processed 
      }  # all exons 
   } # // transcripts  
  return ( \%unique_introns, \%i2t) ;       
}


  

sub get_recover_ids {  
  my ($self,  $gene, $t_long, $aref ) = @_ ;     

  my $g_sid = $gene->stable_id ? $gene->stable_id : "NO_ID_".$gene->dbID ; 
  my %recover_ids ;  

  my ($start, $end , $chr ) = get_coords($gene) ;

  #print " coord_system of gene : " . $gene->slice->coord_system_name()."\n" ; 
  #print " coord_system of gene : " . $gene->slice->coord_system()."\n" ; 
  #print " coord_system of gene : " . $gene->slice->coord_system->version()."\n" ; 
  #print " coord_system of gene : " . $gene->slice->coord_system->name()."\n" ; 

  my $cs = $gene->slice->coord_system->name();   
  my $csv = $gene->slice->coord_system->version();   

  unless ( $csv ) {  
      warning("Error: Could not get coord_system_version for coord_system with name : ".$cs ."\n");
  }         

  for my $homolog ( @$aref ) { 
    $recover_ids{"$cs:$csv:$chr:$start:$end:1:".$t_long->stable_id.":".$homolog->stable_id}=1;
  }    

  return [keys %recover_ids] ; 
}  


sub get_homologues_genes {
  my ($self, $gene) = @_ ; 

   #my @specs_to_check = ( $self->first_species, $self->second_species ) ; 
   my @specs_to_check = @{$$FIND_FALSE_INTRONS{ANALYSIS_SETS}{$self->analysis->logic_name}}  ; 

   my @homologs; 

#   SPECIES: foreach my $spec ( @specs_to_check ) {    
#     my $homology = get_one2one_homology_for_gene_in_other_species( $gene ,$spec  ) ;    
#
#     if ( $homology ) {  
#       my $homol_gene = get_gene_obj_out_of_compara_homology_object($homology, $spec) ;    
#
#       if ( $homol_gene ) {     
#         if ( $homol_gene->biotype=~m/protein_coding/ ) { 
#           push @homologs , $homol_gene ; 
#         } 
#       }
#     } 
#   }   
#   for ( @homologs ) {  
#      print "METHOD_1 identified homologs        : " . $_->stable_id . "\t". $_ . "\n"  ; 
#   }  

   my @hg_2 ; 
   SPECIES: foreach my $other_species ( @specs_to_check ) {     
     print "getting homolog for $other_species\n"  ; 
     my $hg =  get_one2one_orth_for_gene_in_other_species($gene, $other_species) ;   
       if ( $hg ) {     
         if ( $hg->biotype=~m/protein_coding/ ) { 
           push @hg_2 , $hg; 
         } 
       } 
   } 


   return \@hg_2; 
} 







sub read_and_check_config { 
  (my $self) = @_ ; 

   #
   # Check the compara config file in Config/OrthologueEvaluator; 
   # 

   throw("Your compara-registry-file LOCATION_OF_COMPARA_REGISTRY_FILE does not exist !!".
        "\nCheck your config !")
     unless ( -e $$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE} ) ;

     print STDERR "reading compara-config file : $$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE}\n";
    # check compara configuration file and schema-versions of dbs

    Bio::EnsEMBL::Registry->load_all($$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE},1,1);
    my @dba = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors()};

    my %tmp;
    for my $a ( @dba ) {
        push @{$tmp{ $a->get_MetaContainer->get_schema_version }} , $a->dbname  ;
    }
    if ( keys %tmp > 1 ) {
       warning("You're using databases with different schemas - this can cause problems ...\n" ) ;
       for ( keys %tmp ) {
         print "schema $_:\n". join("\n" , @{$tmp{$_}} ) . "\n\n"  ;
       }
    } 
    for my $a ( @dba ) {  
        $a->disconnect_when_inactive(1); 
    }
   
} 


# work error !!!!  
#
# substr outside of string at  
#
# /lustre/work1/ensembl/jhv/project_genestructure_comparison/cvs/ensembl_42/modules/Bio/EnsEMBL/DBSQL/SequenceAdaptor.pm line 264. 
#


sub get_coords{  
 my ($gene) = @_ ;   

  my ($start, $end, $chr ) ;
 if ( ($gene->seq_region_start - 50000 ) < 0  ) {  
    $start = $gene->slice->start ; 
 } else { 
    $start = $gene->seq_region_start - 50000 ; 
 }

 if (($gene->seq_region_start + 50000 ) > $gene->slice->end) {  
     $end = $gene->slice->length ; 
 } else {   
     $end = $gene->seq_region_end + 50000 ; 
 } 
 $chr = $gene->slice->seq_region_name ; 
 return ( $start, $end, $chr ) ;                  
}




sub fetch_input {
  my( $self ) = @_;  
  print "\n\nChecking introns for the species configured in ".
        " Config/GeneBuild/OrthologEvaluator.pm : $$MAIN_CONFIG{QUERY_SPECIES}\n" ; 

  $self->cmp_species( $$MAIN_CONFIG{QUERY_SPECIES} );  
  print "getting adaptor from registry ... " . $self->cmp_species . "\n" ; 
  my $ga = Bio::EnsEMBL::Registry->get_adaptor($self->cmp_species,"core","gene");     
  unless ( $ga ) {  
    print " I can't get a gene-adaptor for " . $self->cmp_species . "\n" . 
          " Maybe you forgot to add the alias for this to your compara registry file\n" ; 
    exit(0); 
  } 

  print "==> Fetching gene from db : " . $ga->db->dbname . " \@ " . $ga->db->host . 
        " : ".  $ga->db->port ."\n" ;

  my $sa = Bio::EnsEMBL::Registry->get_adaptor($self->cmp_species,"core","slice");   
  throw("Can't get GeneAdaptor from Registry") unless $ga ; 

  my @array = split(/:/,$self->input_id) ;  
  my $slice ;   

  my @genes_for_stats;  

  # we can use any input_id to investigate all genes in genome, or we can 
  # use a slice-name like 'chromosome:NCBI34:X:100:200' to just inspect a slice  
    
  my @limit_to_biotypes = @{$$FIND_FALSE_INTRONS{DEFAULT_GENE_BIOTYPES}};  

  if(scalar(@array) < 3 || scalar(@array) > 6) {   
    # process all genes  
    print "Fetching all genes in genome - you can also supply a slice name to limit to a certain slice ...\n" ; 
    if ( @limit_to_biotypes ) { 
      print "Limiting analysis to biotypes : ";   
      for ( @limit_to_biotypes ){
        print "$_\t" ; 
        push @genes_for_stats, @{$ga->fetch_all_by_biotype($_)} ; 
      }
      print "\n"; 
    }else{ 
      @genes_for_stats = @{$ga->fetch_all};
    } 
  }else {    
      # only process a certain slice  
      my $slice = $sa->fetch_by_name ( $self->input_id) ;   
      print "==> Getting genes on slice ".$self->input_id . "\n"  ; 
      if ( @limit_to_biotypes ){ 
        print "==> Limiting analysis to biotypes :\n";   
        for my $bt ( @limit_to_biotypes){
          print "\t- $bt\n" ; 
          push @genes_for_stats, @{$slice->get_all_Genes_by_type($bt)} ; 
        }
        print "\n\n"; 
      }else{   
        print "==> using all biotypes \n" ; 
        push @genes_for_stats, @{$slice->get_all_Genes}; 
      }
    }  
   
    # this code is just for test purposes to run on small sets 
    my $short = 0 ; 
    if ( $short ) { 
      print "SHORT==1 ONLY RUNNING ON SUBSET !!!!!!!!!!!!!\n";
      @genes_for_stats = splice(@genes_for_stats,0,10);
    }  
    print "==>  have ".scalar(@genes_for_stats)." genes fetched for ".$self->cmp_species()."\n";
    $self->genes(\@genes_for_stats) ;  
} 




use vars '$AUTOLOAD';
sub AUTOLOAD {  
 my ($self,$val) = @_;
 (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
 $self->{$routine_name}=$val if $val ; 
 return $self->{$routine_name} ; 
}
sub DESTROY {} # required due to AUTOLOAD

1;
