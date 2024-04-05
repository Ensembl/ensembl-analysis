=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::RecoverFalseIntrons - 

=head1 SYNOPSIS

  my $fs = Bio::EnsEMBL::Analysis::RunnableDB::RecoverFalseIntrons->new(
			      -analysis   => $analysis_obj,
			     );

  $fs->fetch_input();
  $fs->run();
  $fs->output();
  fs->write_output(); 


=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Analysis::Runnable::OrthologueEvaluator and is 
 used to fetch the input from different databases as well as writing results 
 the results to the database.a

 The input_id_format for this analysis requires that you ran the module FindFalseIntrons beforehand. 
 This is because false introns ( coding introns in a different species ) are stored as simple 
 features in the core /reference database and their db id is used in this analysis, to check if the 
 false intron has been fully recovered or not.  

  The input_ids are quite long for this analysis so you have to patch  the input_id_analysis table : 
  
  alter table input_id_analysis change input_id input_id varchar(950) ;

  The input_ids will look like this : 

  chromosome:BROADD2:31:26912747:27017854:1:ENSCAFT00000013647:ENSMUST00000039449,ENST00000361371,ENST00000389194:0_0_0:2190659 

  chromosome:BROADD2:31:26912747:27017854:1  ENSCAFT00000013647     ENSMUST00000039449,ENST00000361371,ENST00000389194    : 0,1_0_0 :        2190659,21906560 
     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^          ^^^^^^^^^        ^^^^^^^^^^^^^^
       slice to re-run analysis on           Trans. - StableId        Id's of homologues which hace the correct struct.    indexes of       simple_feature dbID's   
                                             of 'faulty' transcript   the protein supporting features of these trans       simple_features 
                                                                      will be used to recompute the prediction             ( see below)     ( see below ) 



    the indexes of simple features provide a 'link' between the homolog transcripts which are used for the 
    re-computation and the false introns ( simple features ) which they should 'repair'. 

    i.e. ENSMUST00000039449,ENST00000361371,ENST00000389194 : 0,1_0_0 : 2190659,2190656  

    exactly means : 

       ENSMUST00000039449 has two exons which has a non-coding counterpart in the  
       questionable transcript ENSCAFT00000013647. These have index 0 and 1  in the simple-feature 
       array string. Simple-features  2190659 and 2190656 ( intron coords ) should be overlapped by the 
       new computed gene structure. 

       ENST00000361371 _0_ 2190659   - only simple feature 2190659 should be overlapped 


=head1 METHODS

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a '_'

=cut



package Bio::EnsEMBL::Analysis::RunnableDB::RecoverFalseIntrons;  
#package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateVSGenewise;  
use warnings ;
use strict; 
use List::Util qw(max); 
use Bio::EnsEMBL::Utils::Exception qw(throw warning info stack_trace_dump);
use Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes; 
#use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedExonerate; 
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript; 
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExamineGeneSets; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator; 
use Bio::EnsEMBL::Analysis::Config::Databases;
# Not sure if this module ever existed
#use Bio::EnsEMBL::Analysis::RunnableDB::ExamineGeneSets; 
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch; 
use Bio::EnsEMBL::Registry;  
use IO::String;
use Bio::SeqIO; 
use Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Analysis::Runnable::BestTargetted; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils; 
use Bio::EnsEMBL::Analysis ; 

use Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter;

# to get the registry file to load the databases : 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator;
use Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateForGenewise;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils 
         qw(
            get_transcript_with_longest_CDS
            get_one2one_orth_for_gene_in_other_species
           );
# If you really need this module you can get it from ensembl-pipeline, tag cvs/pre_runnable_delete or use Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
#use Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::HomologyUtils 
         qw( 
            get_gene_obj_out_of_compara_homology_object
           );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
         qw(
            list_evidence 
            convert_translateable_exons_to_exon_extended_objects
           ); 
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name write_seqfile);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils
  qw(create_feature_from_gapped_pieces create_feature);

use vars qw(@ISA);  

@ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild ) ; 




sub new {
  my ($class,@args) = @_;  
  my $self = $class->SUPER::new(@args);   


   


  # old input_id_format : 
  #
  # chromosome:BROADD2:3:30384105:30579389:1:ENSCAFT00000014395:ENST00000282259:2187827
  # chromosome:BROADD2:3:30384105:30579389:1:ENSCAFT00000014395:ENST00000396135:2187827
  # chromosome:BROADD2:3:30384105:30579389:1:ENSCAFT00000014395:ENST00000396137:2187809
  # chromosome:BROADD2:3:30384105:30579389:1:ENSCAFT00000014395:ENST00000396137:2187826
  # chromosome:BROADD2:3:30384105:30579389:1:ENSCAFT00000014395:ENST00000396137:2187827
  #
  # can be changed to :
  #
  # chromosome:BROADD2:3:30384105:30579389:1:ENSCAFT00000014395:ENST00000282259[1],ENST00000396135[1],ENST00000396137[1,2,3]:2187827,2187809,2187826
  #
  # would help - so this means  :
  #
  #    ENSCAFT00000014395 contains a lot of introns to recover.
  #
  #    homolog ENST00000282259 has coding intron of simple feature [1] which is 2187827
  #    homolog ENST00000396135 has coding intron of simple feature [1] which is 2187827
  #    homolog ENST00000396137 has coding intron of simple features [1,2,3]= 2187827,2187809,2187826
  #
  # the good thing about this : all possible genes are computed in one run so we don't need another run to compare the genes, or to cluster them. 
  # 
  
  my @input = split/\:/,$self->input_id ;    
  my ($csys,$asm,$name,$start,$end,$strand,$transcript_sid,$homolog_transcript_sid, $simple_feature_index, $simple_feature_string) = split/\:/,$self->input_id ;     

#  print "csys $csys\n" ; 
#  print "asm $asm \nname $name\nstart $start\nend $end\nstrand $strand\nquery$transcript_sid\nhom_string$homolog_transcript_sid\nsimple_feature_index $simple_feature_index \nsimple_feature_string $simple_feature_string\n" ; 


  if ( scalar(@input) != 10) { 
 #   throw ("wrongly-formatted input_id ; should be : ". 
 #          "CcoordSystemName:Assembly:SeqRegionName:START:END:STRAND:PROTEIN_DB:PROTEIN_id:transcript_sid:SupportingFeatureDBId".
 #          "\nexample:   chromosome:Btau_3.1:10:16466933:16595540:1:ENSBTAT00000018735:ENST00000379915:347368" ) ; 
  } 
 
  # process input_id .....
 
  # simple features are a string of format 2344,23425,3245 
  my @simple_feature_ids = split /\,/, $simple_feature_string ; 
  $self->simple_feature_db_id(\@simple_feature_ids) ;

  # simple features are a string of format ENST00000282259,ENST00000396135,ENST00000396137 
  # this string can include different species !!!!!!
  my @stable_ids_of_homolog_transcripts = split /\,/, $homolog_transcript_sid ; 
  $self->homolog_transcript_sid(\@stable_ids_of_homolog_transcripts) ; 


  # simple feature index  stored in a string like : 0_0,1,2_0 
  # gives information which of the homolog transcripts has recovered which false intron ( which is stored in simple feature table  )
   
  $self->simple_feature_index($simple_feature_index) ;   

  $self->transcript_sid ( $transcript_sid ) ;    

  my $simgw_input_id = "$csys:$asm:$name:$start:$end:$strand" ; 
  $self->simgw_input_id ( $simgw_input_id) ; 

  print "simple_feature_ids: " . join (" ", @simple_feature_ids) . "\n" ;  
  print "homolog_trans_sid : " . join (" ", @stable_ids_of_homolog_transcripts) . "\n" ; 
  print "transcript_sid    : $transcript_sid\n" ;  
  print "simgw-id          : $simgw_input_id\n" ;  

  $self->{_ex_genes} = [] ;  
  $self->{_gw_genes} = [] ;  
  $self->{_src_genes} = [] ;  
  $self->{hom_trans_cache} = {} ; 
  $self->{sf_to_hom_trans_cache} = {} ; 
  $self->{_test_exon} = {} ; 
  $self->{_homolog_exon_alignments}= {} ; 
  $self->{_exon_pairing} = {} ; 
  $self->{_register_testex_vs_homex_coverage}= {} ; 
  $self->{_register_nr_aligned_exons_vs_all_exons_in_trans} = {} ; 
  $self->{_register_f_matrix_score} = {};

  $self->{_hit_matrix}=[];

  $self->limit_input_ids(0); 
  $self->limit_genewise_options(1); 
  $self->use_artifical_input_ids(0);   
  $self->only_run_exonerate(0);  
  $self->lightweight_analysis(0);  

  if ( $self->lightweight_analysis) { 
    $self->limit_genewise_options(1); 
    print "\n\n\t\t********* WARNINNG - Running lightweight analysis and only 1 protein per species per recovery-case - this is quicker  *******\n\n" ;  
  } 

  if ( $self->limit_genewise_options) {  
    print "\n\n\t\t********* WARNING - ONLY GENEWISE OPTIONS LIMITED ! ******** \n\n\n\n" ; 
  } 

  if ( $self->only_run_exonerate ) { 
    print "\n\n\t\t******** WARNING - ONLY RUNNING EXONERATE  !!!!****** \n\n\n" ;  
    sleep(1);
  } 

  $self->debug(0);
  $self->use_dna_align_features_from_homologs_to_add_utr(0); 
  if ( $self->use_dna_align_features_from_homologs_to_add_utr ) { 
    print "\n\n\t\t******** WARNING - USING SPECIES cDNA / EST to predict UTR !!!!****** \n\n\n" ;  
    sleep(1);
  } 

  return $self;
}



sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;
  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;
  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );

    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );

  }
  else{
    &throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;

}

sub fetch_input {   
  my ( $self) = @_ ;   
 
  # 
  # fetch questionable transcript with possible false introns from QUERY_SPECIES database 
  # query_species is the species which has the gene model which we want to check // which 
  # might contain repeetive introns which should be exons  
  # 

  my $new_gene = $self->make_new_gene_from_source_transcript($self->transcript_sid);  
  $self->src_genes($new_gene) ; 

  
  my $enst_to_simple_features  = $self->build_hash_of_recoverd_introns_per_transcript() ;  
  $self->recovered_introns_of_transcript($enst_to_simple_features) ;  

  #
  # get protein supporting features for Exons of all homologes genes which have coding exon where there is  
  # an intron in the query species 
  # 
  #  
  
  my %protein_align_feature_acc ;    
  my %dna_align_feature_acc ;    

  for my $stable_id_of_homolog_transcript ( @{$self->homolog_transcript_sid} ) {   

     my ($sp,$type) = Bio::EnsEMBL::Registry->get_species_and_object_type( $stable_id_of_homolog_transcript ) ; 
     my $ta = Bio::EnsEMBL::Registry->get_adaptor($sp,"core","transcript");       
     Bio::EnsEMBL::Registry->set_disconnect_when_inactive(); 

     print " fetching homologues transcript by stable_id : $stable_id_of_homolog_transcript\n" ; 
  
     my $one2one_orth_trans = $ta->fetch_by_stable_id( $stable_id_of_homolog_transcript ) ; 
  
     unless ( $one2one_orth_trans ) {  
          throw("no transcript with stable_id $stable_id_of_homolog_transcript found\n") ; 
     }   

     $self->set_homolog_transcript($one2one_orth_trans) ; 

     # also recompute with the translation of the one2one homologue 
     my $one2one_translation = $one2one_orth_trans->translation ;  
     $protein_align_feature_acc{$one2one_translation->stable_id} = 1 ;  

     # get all supporting features for all Exons - filter out dna_align features    

     my @e_sf = map { $_->get_all_supporting_features } @{$one2one_orth_trans->get_all_Exons}  ;   

     for my $sf_array ( @e_sf) {      
        my @reduced_sf_array ; 
        if ( $self->lightweight_analysis){ 
          @reduced_sf_array = splice @$sf_array,0,3;  
        }else {  
          @reduced_sf_array = @$sf_array;
        }
        for ( @reduced_sf_array) {   
           #print "using " . $_->hseqname . "\n" ; 
          if (ref($_)=~m/Bio::EnsEMBL::DnaPepAlignFeature/){   
            $protein_align_feature_acc{$_->hseqname} = 1 ;   
            $self->sf_to_gene($_->hseqname, $one2one_orth_trans) ;    
            if ( $_->hseqname =~m/\.+/ ) {  
              (my $tmp=$_->hseqname)=~s/\..+//g;
              $self->sf_to_gene($tmp, $one2one_orth_trans) ;   
            }
          }elsif (ref($_)=~m/Bio::EnsEMBL::DnaDnaAlignFeature/){   
            $dna_align_feature_acc{$_->hseqname} = 1 ;   
          }
        } 
     }  
     # also register the ENSMUSP whatever protein as sf to compute f_matrix later  
     if ( $one2one_orth_trans->translation ) { 
        $self->sf_to_gene($one2one_orth_trans->translation->stable_id, $one2one_orth_trans) ;    
     }
   } # next homolog transcript id  
   $self->dna_align_features_of_homologs([keys %dna_align_feature_acc]) ;  

   print "\n\n" ; 
   my %sf = %{$self->{sf_to_hom_trans_cache} } ; 
   for ( keys %sf ) {  
     print "sf : " . $_ . "\n" ; 
   } 

    
   # re-run trial with old supporting features as well ... 
   # not sure if this was a good idea - if you re-run with old supporting features you can't really
   # compare/ align against the old exons of the buggy gene .....
   
#   
#   my $q_trans = $self->db->get_TranscriptAdaptor->fetch_by_stable_id($self->transcript_sid) ; 
#   my @q_sf = map { $_->get_all_supporting_features } @{$q_trans->get_all_Exons} ;  
#     
#   for my $sf_array ( @q_sf ) {  
#     for ( @$sf_array) {  
#       if (ref($_)=~m/Bio::EnsEMBL::DnaPepAlignFeature/){  
#          $protein_align_feature_acc{$_->hseqname} = 1 ;  
#          $self->sf_to_gene($_->hseqname, $q_trans) ; 
#       }
#     } 
#   }     
#
   my @target_input_ids ;   

   my @my_input_ids = qw ( Q9GM99.1 P42859.2 P42858.1 P42859-2) ;  
   if ( $self->use_artifical_input_ids ) {  
     print "\n"x30 ; 
     print "\nWarning - using artificial input_ids"x30 ; 
     print "\n"x30 ; 
     my %tmp ; 
     @tmp{@my_input_ids}=1; 
     %protein_align_feature_acc= %tmp ; 
   } 

  for ( keys %protein_align_feature_acc ) {   
    push @target_input_ids, $self->simgw_input_id.",".$_ ;  
  }   

  if ( $self->limit_input_ids ) {  
    print "\n"x30 ; 
    print "\nWarning - input-id range limited "x30 ; 
    print "\n"x30 ; 
    print "I will run on these input_ids / proteins :\n" ;  
    print "\n\n\n\n\nWARNING - INPU_IDS LIMITED ----\n\n\n\n\n" ; 
    @target_input_ids = splice(@target_input_ids,0,5) ;  
  }

  for ( @target_input_ids ){  
    print "input_id : " . $_."\n" ; 
  } 
  #
  # get protein seq from homologues transcript stable id .... 
  # 
 
   

  print "\nRunning ExonerateTranscript-runs\n====================================\n\n" ;    

  my @protein_files ; 
  my @eg ;  
  my %genomic_seq_written;    

   #
   # we will run the exonerate comptues first and we don't delete the protein files. 
   # than we will run the genewise on the protein files, than we delete protein files. 
   # that saves us time to write the files twicce .... 
   #


  PROTEIN_ID: for my $targ_xrate_input ( @target_input_ids ) {
    print "\n\ninput_id : $targ_xrate_input\n\n" ;  

     # prepare input for ExonerateTranscript run 

     my ($slice_name, $protein_id) = split /\,/, $targ_xrate_input;  

     unless ( $genomic_seq_written{$slice_name}) { 
       $self->write_genomic_sequence_to_file($slice_name) ; 
       $genomic_seq_written{$slice_name} = 1 ; 
     }  


################

    my $seqfetcher = Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher->new(
                                                                                  -db     => [("/lustre/blastdb/Ensembl/uniprot_index/")], 
                                                                                  -format  => 'fasta', 
                                                                                 );

    my $seqx ;
    eval {
      $seqx = $seqfetcher->get_Seq_by_acc($protein_id);
    };

    if($@) {
      warning("Problem fetching sequence for id [$protein_id] with $seqfetcher $@\n");
    }

    if(!defined($seqx)){
      warning("Could not find sequence for [$protein_id] with $seqfetcher "); 
    }
    my $protein_seq ;  

    unless ( defined $seqx )  {  

       print "Could not find sequence for [$protein_id] with $seqfetcher - now trying with pfetch.... \n" ;  
       my $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new() ; 
        $protein_seq = $seq_fetcher->get_Seq_by_acc($protein_id ) ;   
  
       unless ( $protein_seq ) { 
          $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new( -options => '-a' ) ;  
          $protein_seq = $seq_fetcher->get_Seq_by_acc($protein_id ) ;   
       }  

       unless ( $protein_seq ) {  
          (my $p=$protein_id)=~s/\..+//g; 

          $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new(); 
          $protein_seq = $seq_fetcher->get_Seq_by_acc($p) ;    
          unless ( $protein_seq ) { 
            $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new( -options => '-a' ) ;  
            $protein_seq = $seq_fetcher->get_Seq_by_acc($p) ;   
          }  
       } 

       unless ( $protein_seq ) {   
          print "Failed to fetch sequence $protein_id via pfetch -a / OBDA\n"; 
          next PROTEIN_ID ; 
       }  
    } else {   
        $protein_seq = $seqx ; 
    } 
     my $pepfile      = create_file_name("protein_" ,"fa" , "/tmp" );  
     write_seqfile($protein_seq, $pepfile, 'fasta') ; 
     push @protein_files, $pepfile ; 

     print STDERR "running ExonerateTranscript on Protein :  ".$protein_id.
     " and ".$self->slice->name." length ".$self->slice->length."\n";
     print "\nUsing hardcoded program vwerion exonerate 1.4.0 !\n" ;
 
      my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                                                                         -analysis => $self->analysis,
                                                                         -query => $self->slice, 
                                                                         -target_file => $self->genomic_file, 
                                                                         -query_type => 'protein',
                                                                         -query_file => $pepfile,
                                                                         -program   => "/software/ensembl/bin/exonerate-1.4.0" ,
                                                                         -options => "--model protein2genome --bestn 1 ".
                                                                                     "--maxintron 700000 " , 
                                                                        );
    
     $r->run() ;    
     my @output = @{$r->output()} ;   
     # these are actually all transcripts not genes .... 
   
     for ( @output ) {      
       map { $_->slice($self->slice()) } @{$_->get_all_Exons} ; 
       my $g= Bio::EnsEMBL::Gene->new();  
       $g->biotype("exonerate_recomp") ;  
       $_->biotype("exonerate_recomp") ;  
       $g->add_Transcript($_);  
       $self->ex_genes($g); 
       push @eg, $g ; 
     }    
  } 
  print scalar(@eg) . " exonerate genes made\n" ;

  for ( @eg ) {  
    print_whole_gene($_) if $self->debug ;
  } 


  print "\n\nRecoverFalseIntrons: Now running BlastMiniGenewise:\n" ;  
  print "================================================\n\n" ; 

  # either the analysis is stored or we need to create a new one :   
  # needs to be re-written    
  #

 my @gw_options ; 

 if ( $self->limit_genewise_options ) {    
   print "WARNING !!!!!!!!! Genewise options limited !!! \n"x10 ; 
   push @gw_options, { '-options' => '-init endbias -splice_gtag -quiet ' } ; 
   push @gw_options, { '-options' => '-splice_gtag -quiet ' }   ;  
 } else {   

   push @gw_options, { '-options' => '-quiet' } ; 
   push @gw_options, { '-options' => '-quiet -splice_gtag ' }   ; 
   push @gw_options, { '-options' => '-quiet -splice_gtag -init endbias ' } ; 
   push @gw_options, { '-options' => '-quiet -nosplice_gtag ' } ; 

 } 

  my @gw_genes ;  

  if ( $self->only_run_exonerate ) {   
     print "WARNING - only running exonerate, not genwises .\n" ; 
     for my $g ( @eg) {  
        $self->gw_genes($g);
     } 
  } else { 
  
    my $gs = Bio::SeqIO->new(-file   => "<".$self->genomic_file , -format => "FASTA" );     
    my $genomic_seq = $gs->next_seq(); 
  
  
    # hmm i need to know which strand we're on .... need have the correct strand information in the input_id .... 
    # or does this disturb iexonerate ?  
  
  
    my $sf_strand = ${$self->get_aligned_intron_coord()}[0]->strand(); 
   
  
      for my $gw_run_options ( @gw_options ) {    
          for my $protein_file ( @protein_files )    { 
             my $ps = Bio::SeqIO->new(-file   => "<$protein_file", -format => "FASTA" );     
             my $prot_seq = $ps->next_seq() ;     
    
             my ($sr,$asm,,$chr,$start,$end,$strand,$fintron_trans,$hom_trans,$simple_feat_id) = split/:/, $self->input_id ; 
    
             my $rangefeat = new Bio::EnsEMBL::FeaturePair
             (
             -start => 1, 
             -end   => ( ($self->slice->end+1) - $self->slice->start ) , 
             -strand=> $sf_strand, 
             -slice => $self->slice
             );
             $rangefeat->strand(-1) if($strand == -1); 
    
             my $features = [$rangefeat];
    
             my $gw = Bio::EnsEMBL::Analysis::Runnable::MiniGenewise->new(  
                                                                      -query=> $self->slice, 
                                                                      -protein => $prot_seq, 
                                                                      -analysis=> $self->analysis,
                                                                      -features => $features, 
                                                                      -genewise_options => $gw_run_options,  
                                                                     );  
              if ( $sf_strand == 1 )  {  
                 $gw->reversed_sequence(0); 
              } else { 
                 $gw->reversed_sequence(1); 
              } 
              $gw->run() ; 
              my @transcripts = @{ $gw->output } ;  
    
              # Due to miniseq/fullseq possibility in Genewise, the transcripts produced 
              # only hold FeaturePairs rather than DnaPep align features - we need to 
              # convert them .... this is normally done by BlastMiniGenewise.
              # genewise only makes transcripts !  
     
              for my $t ( @transcripts ) {      
                   my $g = Bio::EnsEMBL::Gene->new(
                                 -start => $t->start , 
                                 -end   => $t->end , 
                                 -strand => $t->strand , 
                                 -slice => $self->slice, 
                                 ) ; 
    
                    $g->add_Transcript($t) ;     
                    $t->analysis($self->analysis) ; 
                    $self->gw_genes($g);
          }
        } 
      }   
    
  }   

  #
  # remove protein - and genomic files ......
  # 
  for my $pf ( @protein_files ) {  
     my $cmd = "rm $pf" ; 
     system("$cmd"); 
  } 
  
 # my $cmd = "rm " . $self->genomic_file();  
 # system("$cmd");  
}





sub print_pos {  
  my ($obj,$string , $tsf)  = @_ ;    
  #print stack_trace_dump(). "\n"; 
  print "$string " . $obj->seq_region_start . "\t" ; 
  print "" . $obj->seq_region_end   . "\t" ; 
  print "" . ( $obj->seq_region_end  - $obj->seq_region_start + 1 )   . "\t" ; 
  print "" . $obj->seq_region_strand ;   
  if ( ref($obj)=~m/Gene/) {  
    print "  ".$obj->biotype." " ; 
  } 
  if ( ref($obj)=~m/Transcript/) {  
    print "  ".$obj->biotype." " ; 
  } 
  if ( $tsf ) {
    print "  " . map { $_->display_id } @$tsf ; 
  } 
  print "\n" ;  
}


sub run { 
  my ($self) = @_ ;  

  my @gw_genes = @{$self->gw_genes} ; 
  my @ex_genes = @{$self->ex_genes} ; 
  my @src_genes = @{$self->src_genes} ;  

  # setting standard biotypes for both sets  
  if ( $self->debug) { 
     print "exonerate-genes:\n=================\n";
     for my $g ( @ex_genes ) {   
         print_whole_gene($g);
         print_pos($g, "EX_GENE") ; 
     }  
     print "genewise-genes:\n=================\n";
     for my $g ( @gw_genes ) {   
       print_whole_gene($g);
       print_pos($g, "GW_GENE") ; 
     } 
     print "src-genes:\n=================\n";
     for my $g ( @src_genes ) {  
         print_whole_gene($g);
         print_pos($g, "SRC_GENE") ; 
     } 
  } 

  my %exons ; 
  for my $gw ( @gw_genes ) { 
    for my $tw ( @{$gw->get_all_Transcripts} ) {  
      my $tw_exon = scalar(@{$tw->get_all_Exons}); 
      $exons{$tw}= scalar(@{$tw->get_all_Exons}); 
    } 
  } 
  
  for my $gw ( @ex_genes ) { 
    for my $tw ( @{$gw->get_all_Transcripts} ) {  
      my $tw_exon = scalar(@{$tw->get_all_Exons}); 
      $exons{$tw}= scalar(@{$tw->get_all_Exons}); 
    } 
  } 

  my @all_genes= (@gw_genes, @ex_genes,@src_genes ) ;   

  # use BestTargetted to check which gene strucutre is best  

  my %unique_biotypes = map {$_->biotype => 1} @all_genes;  

  # genes will be clustered by biotypes so we have to handover all biotypes we got  
  
  my $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new() ;     

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BestTargetted->new
     (
      -biotypes => [keys %unique_biotypes],
      -seqfetcher => $seq_fetcher ,
      -verbose => 1  , 
      -genes => \@all_genes, 
      -analysis => $self->analysis, 
     ); 

  $runnable->run;
  my @best_targetted_genes = @{$runnable->output};  

  my @filter_out_buggy_gene ; 
  for ( @best_targetted_genes ) {  
    if ($_->biotype=~m/buggy_gene/){   
      print "BestTargetted returned our buggy_gene which we used as input. we're skipping it. skipping buggy gene\n"; 
      next ; 
    }else {  
      push @filter_out_buggy_gene, $_ ; 
    }  
  }  
  @best_targetted_genes = @filter_out_buggy_gene;  

  if ( $self->debug ) {  
     print_best_targetted_genes(\@best_targetted_genes) ;  
  }


  # 
  # i should do a kind of  get aligned intron coords by stable id of homolog thingy ......
  # 
 
  my $ma = Bio::EnsEMBL::Registry->get_adaptor("compara","compara","member");
  Bio::EnsEMBL::Registry->set_disconnect_when_inactive();

  # dynamicly choose gene-adaptor for core db 
  # ATTENTION ! Name of the species is fetched out of the compara db than the gene-adaptor is 
  # constructed out of the config file !!!   
 

  # this checks if an intron has been removed and checks how many have been removed, too - comparing 
  # to how many coding introns there are in the homolg trnscript 
 
  my @checked_genes ;  

  my %hom_trans_to_gene ; 
  my %hom_trans_is_ccds ; 
 
  for my $homolog_transcript_id ( @{$self->homolog_transcript_sid} ) {  
     
    my ($species, $object )  =  Bio::EnsEMBL::Registry->get_species_and_object_type($homolog_transcript_id ); 
    my $ta = Bio::EnsEMBL::Registry->get_adaptor($species,"core","transcript");      
    my $ga = Bio::EnsEMBL::Registry->get_adaptor($species,"core","gene");      
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
  
    my $homolog_trans = $ta->fetch_by_stable_id($homolog_transcript_id) ;  
    my $homolog_gene = $ga->fetch_by_transcript_stable_id($homolog_transcript_id) ;  

    $hom_trans_to_gene{$homolog_transcript_id} = $homolog_gene;

    # find out if the gene / transcript is a member of the CCDS set ... 
    # this needs to be changed once we start defining CCDS membership for transcripts 

    my @dbe = @{ $homolog_gene->get_all_DBEntries } ;   

    $hom_trans_is_ccds{$homolog_trans->stable_id}="NO" ;  

    for my $dbe ( @dbe ) {  
       if ( $dbe->dbname =~m/CCDS/ )  {
         $hom_trans_is_ccds{$homolog_trans->stable_id}="YES" ; 
       } else {  
         $hom_trans_is_ccds{$homolog_trans->stable_id}="no" ; 
       } 
    } 

    # compare gene-structure with structure of homologue   

    for my $new_gene (@best_targetted_genes) {  
       $self->attach_slice($new_gene) ;  
       
       # this sub changes the biotype if an intron was removed  
       my $ch_gene = $self->check_if_false_intron_has_been_removed_in_new_trans($homolog_trans, $new_gene ) ;  

       if ( $ch_gene->biotype =~m/intron_removed/) {   
         print "match " .  $ch_gene->biotype . " matches biotype intron_removed\n" if $self->debug; 
         push @checked_genes, $ch_gene ;  
       } 
    }  
  }  

  print "\n" . scalar(@checked_genes ) . " transcripts found which have the intron removed\n\n" ;  

  my @intron_removed_genes = @checked_genes ; 
   
  # now remove redudant transcripts of all genes which have the intron removed 

  my @non_redundant_genes = @{$self->remove_redundant_genes ( \@checked_genes )} ;  


   
  # 2 methods of scoring are currently implemente ..... :   
  
  #  SCORING METHOD 1 : 
  #  - calcualting the overlap between old and new gene-exons 
  #  - building ratio 
  #  - sselecting gene with highest ratio and longest translation .....
  # 
#  
#  my $very_best_gene = $self->check_exon_overlap_between_new_gene_and_src_gene(\@non_redundant_genes) ;  
#
   my @very_final_genes ;    
#
#  if ( $self->other_genes ) {  
#     push @very_final_genes, @{$self->other_genes} ; 
#  }  
# 
#
  # SCORING METHOD 2 : 
  # - calculate exon-alignment score for aligning trans against homologe exons and use the one with highest score, 
  # than compare against score of buggy transcript and select one with highest score 
  #
  # calculate f-score for old (buggy) transcript and store them
  # 
  
  my %old_f_matrix_scores;   
  my $src_gene ;  

  for my $stable_id_of_homolog_transcript ( @{$self->homolog_transcript_sid} ) {    
    $src_gene = ${$self->src_genes}[0];
    my $src_trans = ${$src_gene->get_all_Transcripts}[0]; 
    my $homolog_trans = $self->get_homolog_transcript_by_stable_id($stable_id_of_homolog_transcript) ;    

    print " \n\nCalculating_values for old gene vs $stable_id_of_homolog_transcript\n\n" if $self->debug ; 
    my $f_matrix_score_old_gene  = $self->check_test_transcript_against_homolog($src_gene, $src_trans, $homolog_trans )  ;   
    $old_f_matrix_scores{$stable_id_of_homolog_transcript} = $f_matrix_score_old_gene ;  
    print " Transcript : " .  $homolog_trans->stable_id . " scores ......  $f_matrix_score_old_gene\n" ;  
  } 

  print "\nCalculating_scores of all genes against the homologs and set biotype of best gene to 'MAX' ....\n" ;  
  my @most_similar_to_homolog_genes = @{ $self->find_gene_which_has_exon_length_distr_most_similar_to_homolog(\@non_redundant_genes) } ; 

  my @lower_genes ; 
  my @max_genes ;  
  for my $g ( @most_similar_to_homolog_genes ) {  
    if ( $g->biotype=~m/MAX/) {   
      push @max_genes, $g; 
    } else {   
       # this is where all other genes get filtered out.
       print "skipping $g\n" ;  
       push @lower_genes, $g ; 
    } 
  } 
  my @most_sim_non_redundant ;  
  if ( scalar(@max_genes) > 1 )  {  
    print "have " . scalar(@max_genes) .  " genes with the same max score.. let's see if they're redundant.\n" ; 
    @most_sim_non_redundant = @{check_if_genes_are_redundant(\@max_genes)};    

    my $gene ; 
    if ( @most_sim_non_redundant > 1 ) {  
      print "WARNING : genes have the same score but are different - i will use the gene with less exons\n" ; 
      my $less_exon = 100000000;
      for my $g ( @most_sim_non_redundant) {       
         print_whole_gene($g);
         if ( $less_exon > scalar(@{$g->get_all_Exons})  ) {  
            $less_exon = @{$g->get_all_Exons } ; 
            $gene = $g; 
         } 
      }
      @most_sim_non_redundant=();
      push @most_sim_non_redundant, $gene ;
    }  
  } else {  
    print "setting most_sim_non_redundant to max_genes\n" ;
    print "max_genes : " . scalar(@max_genes) . "\n" ;  
    print join("\n", @max_genes ) . "\n" ; 
    @most_sim_non_redundant = @max_genes ; 
  } 
 
  print "\n\n\n" ; 
  print "Comparison of all test_genes against the homologs :\n---------------------------------------------------------\n" ; 
  my @all_scores ;  

  for my $g ( @most_sim_non_redundant ) { 
    print_whole_gene($g) ;    
    print " coverage against homolog transcript for gene $g :\n\n" ;   

    my %tmp1 = %{ $self->register_testex_vs_homex_coverage($g)};  
    my %aligned_coverage = %{ $self->register_nr_aligned_exons_vs_all_exons_in_gene($g)};    
    my %f_matrix = %{ $self->register_f_matrix_score($g)};  
    my $score ;  

    for my $k  ( keys %tmp1 ) {   
      my $score_src_trans_vs_homolog = 0 ; 
      if (  $old_f_matrix_scores{$tmp1{$k}{htrans}->stable_id} ) { 
        $score_src_trans_vs_homolog = $old_f_matrix_scores{$tmp1{$k}{htrans}->stable_id}; 
      }       

      print " homolog transcript                  : " . $tmp1{$k}{htrans}->stable_id  . "\n" ;  
      print " coverage against homolog            : " . $tmp1{$k}{val} . " (nr_aligned_dog_exons / nr_of_exons_in_homolog)\n" ;  
      print " exons aligned vs all exons in gene  : " . $aligned_coverage{$k}{val} . " (nr_aligned_dog_exons / nr_of_all_dog_exons)\n" ;  
      print " f_matrix_score of new transcript    : " . $f_matrix{$k}{val} . " (f-matrix-value)\n" ;  
      print " f_matrix_score of old trans vs hom  : " .  $score_src_trans_vs_homolog . "\n\n" ;    # fix  
      $score = $f_matrix{$k}{val} ;
      if (  $score_src_trans_vs_homolog > $f_matrix{$k}{val} ) {   
        print " the f_matrix_score of the old transcript is higher than the one of the newly built one...\n";   
        print "NO_FIX_old_score_higher_".$score_src_trans_vs_homolog. "_" . $f_matrix{$k}{val} . "\n";  
        $g->biotype("NO_FIX") ;  
      }
    } 
    if ($score ){
      $g->biotype($score) ;  
      push @all_scores, $score ; 
    } 
    push @very_final_genes, $g;
  }  


  @all_scores = sort {$b <=> $a} @all_scores ;   
  my $cnt = 1 ; 
  for my $s ( @all_scores ) {   
    for my $g ( @very_final_genes ) {   
      if ($g->biotype=~m/$s/) {  
         print "score : $s --- " . $g->biotype . "\n" ; 
         $g->biotype($cnt); 
      } 
    } 
    $cnt++; 
  }  
 
  my @genes_with_utr ; 
  print "try to add utr to " . scalar(@very_final_genes) . " genes\n" ;   

  # first try to add utr from org. species....  
  
  my @dna_align_features_of_src_gene = @{$self->src_gene_dna_align_features};   
  my @dna_align_features_of_homologs = @{$self->dna_align_features_of_homologs};    

  if ( scalar(@dna_align_features_of_src_gene)){  
    print "adding utr from org species\n"; 
    for my $vg ( @very_final_genes ) {  
      push @genes_with_utr, @{$self->add_utr($vg,\@dna_align_features_of_src_gene)} ;  
    }
  } elsif(  scalar(@dna_align_features_of_homologs) && $self->use_dna_align_features_from_homologs_to_add_utr ) { 
      print "try to use cDNA's from homologs to add utr to genes\n" ; 
      print "have " .scalar( @dna_align_features_of_homologs ) . " cdnas from homologs which i could use for utr_addition \n" ; 
      if ( $self->use_dna_align_features_from_homologs_to_add_utr ) { 
        print "try to use cDNA's from homologs to add utr to genes\n" ; 
        for my $vg ( @very_final_genes ) {  
          push @genes_with_utr, @{$self->add_utr($vg,\@dna_align_features_of_homologs)} ;  
        } 
      }
  } else {  
    @genes_with_utr = @very_final_genes ; 
  }
  print "now have " . scalar(@genes_with_utr) . " genes\n" ;   

  #
  # now we can identify the most similar transcript ( transcript with high similarity-treshold vs homolog ... 
  # if there is more than one wtih the same score we might pick the one which has the highest coverage against homolog 
  #

  my $cmd = "rm " . $self->genomic_file();   
  system("$cmd");   

  my @complete_genes ;  

  for my $g ( @genes_with_utr ) {  
    my @new_transcripts = @{$g->get_all_Transcripts}; 
    if ( scalar(@new_transcripts)  > 1 ) { 
      throw("gene has more than one transcript after processing - this should not the case. check where the other transcript is added\n"); 
    }  

    $self->attach_slice($g);   
    $g = $self->attach_stuff_to_new_gene($g); 

    # attach existing, non-fixed transcripts back 
     my $orgi_gene = $self->very_original_query_gene(); 
     my @orgi_trans_diff_slice  = @{$orgi_gene->get_all_Transcripts};
     my @orgi_trans ;   

     for ( @orgi_trans_diff_slice  ) { 
         my $t = $_->transfer($self->reference_slice) ; 
         push @orgi_trans, $t ; 
     }
 
#     if ( scalar(@orgi_trans) > 1 )  {
#       print "boof ! have more than one transcript ...... need to add other transcripts .... \n" ;  
#       for my $otr ( @orgi_trans ) {   
#         print "processing transcript : " . $otr->stable_id . "\n"; 
#         if ( $otr->stable_id eq $self->src_transcript->stable_id ) { 
#            print " it's the SAME trans : ".$otr->stable_id." and " . $self->src_transcript->stable_id . " - no need to add this - skipping\n" ; 
#         }  else {  
#
# 
#            if ( $self->check_if_transcripts_are_different([@{$g->get_all_Transcripts}, $otr])) { 
#               # transcripts are all different  
#               print "Adding old transcript: " .  $otr->stable_id . " to gene \n" ;   
#               $g = $g->transfer($otr->slice);  
#               $otr->dbID(undef);
#               $otr->adaptor(undef);  
#               $g->add_Transcript($otr) ;  
#             } else {
#               print "transcripts are redundant - not addin transcript\n" ;  
#            }
#            
#         } 
#       } 
#     }  
     push @complete_genes, $g;
  } 
  for my $c ( @complete_genes ) {  
    print "GENE :  " . $c->biotype . "\t" . $c->stable_id . "\n" ;  

    for ( @{ $c->get_all_Transcripts } ) {  
       print "\tTranscript : " .  "\t" . $_."\t" ;  

       if ( $_->stable_id ) { 
         print "\t" . $_->stable_id . "\t" ;  
       }else{ 
         print " transcript_no_stable_id\t" ; 
       }  

       if (  $_->translation ) { 
         print  $_->translation. "\t".$_->translation->stable_id   ; 
       } else {  
         print " transcript_hs_no_translation" ; 
       }  
       print "\n" ; 
    } 
  }    
  #print "adding crappy stuff for verbosity\n" ; 
  #push @complete_genes, @lower_genes ; 
  $self->output(\@complete_genes) ; 
}


sub check_if_transcripts_are_different{ 
  my ($self, $aref) = @_ ;    

  # to check by hash-key you have to make sure all transcripts are on the same slice :  
  my @transcripts_diff_slice = @$aref ;  
  my @normallised_tr; 
  for ( @$aref ) {  
     my $tmp_t = $_->transfer($self->reference_slice);  
     push @normallised_tr, $tmp_t ; 
  } 


  my %tmp ; 
  for my $t (@normallised_tr  ) { 
    my @ex= @{ $t->get_all_Exons } ;
     my $string ;
     for my $e (@ex) {
       my $exon_hk = $e->hashkey ;
       $string.=$exon_hk ;
     }  
     print $string . "\n" ; 
     $tmp{$string} = 1 ;
   }
   if ( scalar(@$aref) == scalar(keys %tmp)) {   
       print "Transcripts are different\n"; 
       return 1 ; 
   }
       print "Transcripts are duplicated.\n" ; 
  return 0;  
} 

sub check_if_genes_are_redundant { 
    my ($aref) = @_ ;  

    my @max_genes=@$aref; 
    my @non_red_genes ;
    my %tmp ;
    # little check 
    for my $g ( @max_genes) { 
      if ( scalar(@{$g->get_all_Transcripts}) > 1 ){  
        throw("gene has more than one transcript - this is wrong .... : " . $g->biotype); 
      } 
    }  
    for my $g ( @max_genes ) { 
       for my $t (@{ $g->get_all_Transcripts } ) {  
         my @ex= @{ $t->get_all_Exons } ;
         my $string ;
         for my $e (@ex) {
           my $exon_hk = $e->hashkey ;
           $string.=$exon_hk ;
         } 
         $tmp{$string} = $g ;
       }
    }  

    if ( scalar(keys %tmp)==1) {  
        print "genes are redudant \n" ; 
    }else {  
        print "genes_have_same_max_score_but_are_different\n" ;   
        for ( values %tmp ) {  
          print_whole_gene($_); 
        }
    }  
    return [values %tmp] ; 
 }



sub add_utr {  
  my ( $self, $very_best_gene , $cdna_acc_ref ) = @_ ;        

  my @hseq_names = @$cdna_acc_ref ;  

#   print "now trying to add utr\n" ; 
#   my $ta = Bio::EnsEMBL::Registry->get_adaptor($$MAIN_CONFIG{QUERY_SPECIES},"core","gene");    
#   Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
#   
#   my $transcript_sid = $transcript_sid->stable_id; 
#   my $source_gene = $ta->fetch_by_transcript_stable_id($transcript_sid) ;    
#   my @source_transcripts = @{ $source_gene->get_all_Transcripts } ;   
#
#
#   for my $st ( @source_transcripts ) { 
#
#      if ($st->three_prime_utr) {  
#         print "transcript_has_3_prim_utr\t" . $source_gene->stable_id . "\n" ; 
#      } 
#      if ($st->five_prime_utr) {  
#         print "transcript_has_5_prim_utr\t" . $source_gene->stable_id . "\n" ; 
#      } 
#   } 
#   my %cdna_acc ; 
#   for my $e ( @{$source_gene->get_all_Exons} ) {  
#      my @sf = @{ $e->get_all_supporting_features}; 
#      for my $s ( @sf ) {  
#        if ( ref($s)=~m/Bio::EnsEMBL::DnaDnaAlignFeature/) {  
#          $cdna_acc{$s->hseqname}=1; 
#        } 
#      } 
#   }    
#   print scalar(keys %cdna_acc) . " cdna's found which could be added as UTR\n" ;    
#
#  if (scalar ( keys %cdna_acc ) == 0 ) { 
#      print " no cDNAs found which could be added - now trying to fetch / use cDNA from compare-species\n" ;  
#  } 

  # get cdna seq from pfetch   
  my @cdna_seqs ; 
  my $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new() ; 
  for my $acc ( @hseq_names ) {  
    my $seq = $seq_fetcher->get_Seq_by_acc($acc) ;    
    unless ( $seq ) { 
      $seq_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new( -options => '-a' ) ;  
      $seq = $seq_fetcher->get_Seq_by_acc($acc) ;    
    } 
    unless ( $seq ) {  
      print "sequence $acc could not be fetched\n" ; 
    }else { 
      push @cdna_seqs, $seq ; 
    }  
  }  

  # now we need to align the cDNA to the genome to get a proper gene structure........    
  #  align seqs with Exonerate2Genes 

  my @cdna_genes ;  
  my @cdna_file_names ; 
  for my $cdna  ( @cdna_seqs ) { 
     print " running ExonerateTranscript on cdna \n" ;  
     my $cdna_file      = create_file_name("cdna_" ,"fa" , "/tmp" );  
     write_seqfile($cdna, $cdna_file, 'fasta') ; 
      push @cdna_file_names, $cdna_file ; 
     print STDERR "running ExonerateTranscript on cDNA and slice ".$self->slice->name." length ".$self->slice->length."\n";
     print "\nUsing hardcoded program vwerion exonerate 1.4.0 !\n" ;
                                                # CDNA !!!!  
     my $options = "--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE  --score 500 ".
                   "--saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14";  

     my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                                                                         -analysis => $self->analysis,
                                                                         -query => $self->slice, 
                                                                         -target_file => $self->genomic_file, 
                                                                         -query_type => 'dna',
                                                                         -query_file => $cdna_file,
                                                                         -program   => "/usr/local/ensembl/bin/exonerate-0.8.3", 
                                                                         #-program   => "/software/ensembl/bin/exonerate-1.4.0" ,
                                                                         -options => $options , 
                                                                        ); 
    
     $r->run() ;   

     my @output = @{$r->output()} ;     

     my $percent_ident_filter = 97 ; 
     my @o2 = @output ;  
     print "have " .scalar(@output) . " featrues before filtering $percent_ident_filter\n" ;  
     # filter results 
     my $et_filter= new Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter->new(
        -percent_id => $percent_ident_filter,
        -best_in_genome => 1, 
        -coverage => 90, 
       );
     my @filtered_cdnas = @{$et_filter->filter_results(\@output)};  
     print "have " .scalar(@filtered_cdnas) . " features after filtering by percentage_ident >= $percent_ident_filter\n" ; 
     @output = @filtered_cdnas ; 

     # these are actually all transcripts not genes .... 
      
     for ( @output ) {     
       map { $_->slice($self->slice()) } @{$_->get_all_Exons} ; 
       my $g= Bio::EnsEMBL::Gene->new();  
       $g->biotype("cdna_recomp") ;  
       $_->biotype("cdna_recomp") ;  
       $g->add_Transcript($_);  
       push @cdna_genes, $g; 
     }     
  }  # next cnda 
  my $rm_string = "rm ".join (" ", @cdna_file_names ) ;  
  if ( scalar(@cdna_file_names) > 0 ) { 
    system($rm_string);  
  }
  print scalar(@cdna_genes) . " cdna_genes from cdna made in total.\n" ; 
  my @result ;     

   #
   ########## utr builder 
   #  
   
   if ( @cdna_genes ) { 
      my $utr_builder = Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder->new(
        -db        => $self->db,
        -input_id  => 'fake', 
        -analysis  => $self->analysis, 
        -dont_use_config => 1, 
       );     
      $utr_builder->gw_genes([$very_best_gene]) ;  
      $utr_builder->INPUT_GENETYPES([$very_best_gene->biotype]); 
      $utr_builder->cDNA_GENETYPE([$cdna_genes[0]->biotype]); 
      $utr_builder->MAX_EXON_LENGTH(20000); 
      $utr_builder->MAX_INTRON_LENGTH(200000);  
      $utr_builder->query($self->slice); 
      $utr_builder->VERBOSE(1);  
      $utr_builder->UTR_GENETYPE($very_best_gene->biotype); 
      $utr_builder->BLESSED_UTR_GENETYPE("blessed") ; 
      #$utr_builder->UTR_GENETYPE('utr_added');
      # set cDNA's 
      $utr_builder->cdna_genes(\@cdna_genes); 
      $utr_builder->{evidence_sets} = {
                                'est'   => [] , #\@est_biotypes,
                                'simgw' => $utr_builder->INPUT_GENETYPES, 
                                'cdna'  => $utr_builder->cDNA_GENETYPE, 
                               };
      # fetch_input does nearly all the stuff we've just done by hand .... 
      #$utr_builder->fetch_input() ;   
      $utr_builder->run() ;    
    
      @result= @{ $utr_builder->output()} ;   
      #for ( @result) {  
      #  print "UTR gene ? $_\n" ;  
      #  print_whole_gene($_); 
      #}  
  }else {  
    # no cDNA-gene made ....  
    push @result, $very_best_gene; 
  } 
  return \@result; 
}


sub attach_stuff_to_new_gene { 
  my ( $self, $very_best_gene ) = @_ ;    

   my $transcript_sid = $self->transcript_sid() ; 
   my $ta = Bio::EnsEMBL::Registry->get_adaptor($$MAIN_CONFIG{QUERY_SPECIES},"core","gene");    
   my $tla  = Bio::EnsEMBL::Registry->get_adaptor($$MAIN_CONFIG{QUERY_SPECIES},"core","translation");    
   Bio::EnsEMBL::Registry->set_disconnect_when_inactive();

   unless ( $ta ) {  
    throw(
          "could not get gene-adaptor for $$MAIN_CONFIG{QUERY_SPECIES} - check your " . 
          "registry file / your OrthologueEvaluator-config\n"
          ); 
   }
   my $source_gene = $ta->fetch_by_transcript_stable_id($transcript_sid) ;    

   my @src_dbentries =  @{$source_gene->get_all_DBEntries}  ;  

   unless($source_gene){  
    throw("could not get transcript with stable_id $transcript_sid out of database :"  .$ta->db->dbname . "\@" . $ta->db->host. "\n")  ;  
   } 

  my @source_transcripts = @{ $source_gene->get_all_Transcripts } ; 

  my $very_best_transcript = ${$very_best_gene->get_all_Transcripts}[0]; 

  my $new_gene = Bio::EnsEMBL::Gene->new(); 
  $new_gene->stable_id($source_gene->stable_id);   
  $new_gene->biotype($very_best_gene->biotype) ; 
  my $src_status = $source_gene->status ; 
  $new_gene->status($src_status) ;  
  $new_gene->confidence($src_status) ;  

  $new_gene->description($source_gene->description);  
  $new_gene->display_xref($source_gene->display_xref);  
  $new_gene->external_db($source_gene->external_db);  
  $new_gene->created_date($source_gene->created_date);  
  $new_gene->modified_date($source_gene->modified_date);  

   for my $dbe ( @src_dbentries ) {   
     $new_gene->add_DBEntry($dbe) ; 
   } 
  my @new_exons = @{ $very_best_transcript->get_all_Exons } ;   

  my @nsf ;   

  my $original_slice = $very_best_transcript->slice ;  
  $self->reference_slice($original_slice) ;  
  #my $reference_slice ; 
  for my $st ( @source_transcripts ) {  
     print " try to attach : " . $st->stable_id . "\n" ; 

    if ( $st->stable_id=~m/$transcript_sid/ ) {  

      my @old_exons = @{ $st->get_all_Exons } ;  

      for my $oe ( @old_exons ) {   
         for my $ne ( @new_exons ) {   

           if ( $oe->seq_region_end >= $ne->seq_region_start and $oe->seq_region_start <= $ne->seq_region_end )  {  
             $ne->stable_id($oe->stable_id) ; 
             my $v =  $oe->version;
             $v++; 
             $ne->version($v) ; 
             $ne->created_date($oe->created_date) ; 
             $ne->modified_date(time) ; 
           } 
         } 
      } 

      my @tdb = @{$st->get_all_DBEntries}    ;
      # now we need to edit the transcript and remove one exon ....

       # now wee need to attach all stuff to the transcript   
       # how do we deal with UTR  ?
       # how do we deal with the exon-stable-ids etc ? 
       my $new_version = $st->version ; 
       $new_version++ ; 
       my $new_transcript = new Bio::EnsEMBL::Transcript(  
                     -exons => \@new_exons , 
                     -stable_id => $st->stable_id , 
                     -version => $new_version ,  
                     -external_name =>$st->external_name ,  
                     -external_db => $st->external_db ,
                     -external_status => $st->external_status , 
                     -display_xref => $st->display_xref , 
                     -created_date => $st->{'created_date'} , 
                     -modified_date=> time, 
                     -description => $st->description , 
                     -biotype=> $st->biotype, 
                     -status =>$st->status, 
                     -analysis =>$st->analysis, 
                     ) ;     

         for my $db ( @tdb ) {  
             $new_transcript->add_DBEntry($db) ; 
         }  

         $new_transcript->translation($very_best_transcript->translation) ; 
         $new_transcript->translation->stable_id($st->translation->stable_id) ;  
         my $tl_version = $st->translation->version ; 
         $tl_version++ ; 
         $new_transcript->translation->version($tl_version) ; 
         $new_transcript->translation->created_date($st->translation->created_date) ; 
         $new_transcript->translation->modified_date(time);  


         my @tl_db = @{$st->translation->get_all_DBEntries} ; 
         for ( @tl_db ) {   
            $new_transcript->translation->add_DBEntry($_) ; 
         }
         my @all_supporting_features ;    

         for my $f ( @{ $very_best_transcript->get_all_supporting_features} ) { 
           my $nf = $f->transfer($original_slice); 
           push @all_supporting_features, $nf ; 
         } 
                 
         $new_transcript->add_supporting_features(@all_supporting_features);    

         my $nt = $new_transcript->transfer($original_slice); 
         $new_transcript = $nt ;  
      
      # if ( $reference_slice) { 
       #  my $nt = $new_transcript->transfer($reference_slice); 
       #  $new_transcript  = $nt ;  
      # }
       $new_gene->add_Transcript($new_transcript) ; 
    } else {   
         my $nt = $st->transfer($original_slice); 
         $st = $nt ;   
         print $new_gene . "\n" ;
         if ( $new_gene->get_all_Transcripts ) {       
            if ( $self->check_if_transcripts_are_different([@{$new_gene->get_all_Transcripts}, $st])) {  
            # otherwise there's confusion where to get the start / end exon from ... not sure why ...  
              print "adding : " . $st->stable_id . " to transcript \n\n" ; 
              $st->translation->adaptor($tla); 
              $new_gene->add_Transcript($st) ; 
             }
         }else { 
              print "adding : " . $st->stable_id . " to transcript \n\n" ; 
              $st->translation->adaptor($tla); 
              $new_gene->add_Transcript($st) ; 
         } 
       # transfer old ( unchagned ) source-transcript to the same slice 
       #if ( $reference_slice) { 
       #  my $nt = $st->transfer($reference_slice); 
       #  $st = $nt ;  
       #}  
       
    }   
    # first trans was added so now get reference slice  
    #my @tmp  = @{ $new_gene->get_all_Transcripts} ;    
    #$reference_slice = $tmp[0]->slice();   
  }
  return $new_gene ; 
}  




#
# compare the exon length of the new genes against the exon length of the homolog transcripts ..... 
#


sub find_gene_which_has_exon_length_distr_most_similar_to_homolog { 
  my ( $self, $genes_to_test ) = @_ ;   

  my @most_sim_to_hom_genes ; 
  my %scores ; 
  my @test_genes = @$genes_to_test ;  

  for my $g ( @test_genes ) {  
    my $test_trans = ${$g->get_all_Transcripts}[0];  

    print "\n\t\t---------------\n\nNow computing all those funny values for $g\n"  if $self->debug;

    #
    # Only compare the test_trainscript against homolog gene structures which share the same supporting features 
    #
    # i.e. : homolog ENSG1 was build from seq1, seq2, seq3
    # i.e. : homolog ENSG2 was build from seq1, seq2
    # test_transcript was build from seq3. compare test_transcript against ENSG1 
    # 
    
    my $mm_score = 0;  

    for my $sf ( @{ $test_trans->get_all_supporting_features } ) {    

       print "have supporting feature : " . $sf->hseqname . "\n" if $self->debug;  

       my @homologues_transcripts_which_are_supported_by_sf = @{$self->sf_to_gene($sf->hseqname)} ;    

       if ( scalar (@homologues_transcripts_which_are_supported_by_sf) == 0 ) {   
          (my $sfn = $sf->hseqname)=~s/\..+//; 
           @homologues_transcripts_which_are_supported_by_sf = @{$self->sf_to_gene($sfn)} ;    
           if ( scalar (@homologues_transcripts_which_are_supported_by_sf) == 0 ) {   
             throw("No homologues transcript found for supporting_feature : " . $sf->hseqname  );
           }
       }  
       print "Supporting_feature : ".$sf->hseqname . " supports these transcripts : " ; 
       for ( @homologues_transcripts_which_are_supported_by_sf ) {  
           print "  " . $_->stable_id ; 
       }
       print "\n" ; 


       for my $homolog_trans ( @homologues_transcripts_which_are_supported_by_sf) {   
         print "\nCalculating homolog_score for " . $g . " [ " . $sf->hseqname . " ]  against_homolog: " . $homolog_trans->stable_id . "\n" ; 
         my $max_score_f_matrix  = $self->check_test_transcript_against_homolog($g, $test_trans, $homolog_trans )  ;    

         if ( $max_score_f_matrix > $mm_score ) {  
           $mm_score = $max_score_f_matrix ; 
         } 
         $scores{$g}{score} = $mm_score; 
         $scores{$g}{gene}  = $g ; 
       }
    }  
    $g->biotype("lower_homolog_score");
  }   

  
  # thsi whole part needs to be re-written as it's oold ........... 
  my $max_score = 0 ; 
  for my $k ( keys %scores ) {   
    my $sg = $scores{$k}{gene} ;  

    print "transcript " . $sg->seq_region_start . "\t" . $sg->seq_region_end . " scores : " . $scores{$k}{score} ."\n" ;  

    if ( $max_score < $scores{$k}{score}) {  
       $max_score = $scores{$k}{score}; 
    }   
  }   

  for my $g ( @test_genes ) {   
    if ( exists $scores{$g}{score}){ 
      if ($scores{$g}{score} == $max_score) { 
         my $sg = $scores{$g}{gene} ;  
         print "\nmax_scoring transcript : " . $sg->seq_region_start . "\t" . $sg->seq_region_end . " score : " . $scores{$g}{score} ."\n" ;  
         $g->biotype("MAX");
      }
    }
  } 
  return \@test_genes ; 
} 














# this routine is used to get the test_exon object from a reference string while parsing exonerate results 

sub test_exon {  
   my ($self,$test_exon) = @_ ; 
 
   unless ( exists ${$self->{_test_exon}}{$test_exon} ) { 
      ${$self->{_test_exon}}{$test_exon}=$test_exon;
   } 
   return ${$self->{_test_exon}}{$test_exon}; 
} 

#
# this routine stores the homolog exon - test exon pairings which align to each other 
#

sub exon_pairs {   
   my ($self,$homolog_exon,$test_exon) = @_ ; 

   if ( $homolog_exon && $test_exon ) {  
     if ( exists (  ${$self->{_exon_pairing}}{$homolog_exon} )) {  
       my %uniq;
       @uniq{ @{${$self->{_exon_pairing}}{$homolog_exon}} }=1;
 
       unless ( $uniq{$test_exon} ) { 
         push  @{${$self->{_exon_pairing}}{$homolog_exon}}, $test_exon;
       }
     } else {  
       push  @{${$self->{_exon_pairing}}{$homolog_exon}}, $test_exon;
     } 
   } 
   return ${$self->{_exon_pairing}}{$homolog_exon} ; 
} 


sub homolog_exon_alignments{  
   my ($self,$homolog_exon, $val) = @_ ; 

   if ( $val ) {  
     ${$self->{_homolog_exon_alignments}}{$homolog_exon}++; 
   }
   return  ${$self->{_homolog_exon_alignments}}{$homolog_exon}; 
} 




sub check_test_transcript_against_homolog { 
   my ($self,$test_gene, $test_transcript, $homolog_trans ) = @_ ;  

   if ( $self->debug ) { 
     print "\n\n################ test_gene vs HOM TRANSCRIPT " . $homolog_trans->stable_id . " ######################\n";  
     print "\n\t\tTEST_TRANSCRIPT :\n" ; 
     print_whole_trans($test_transcript) ; 
     print "\n\t\tHOMOLOG_TRANS :\n" ; 
     print_whole_trans($homolog_trans) ; 
   } 
    my $max_score_f_matrix ; 

   # prepare matrix 

   my @test_exons    = @{ $test_transcript->get_all_translateable_Exons() };   
   for (my $i=0 ; $i < scalar(@test_exons) ; $i++ ) {    
     $test_exons[$i]->{_jhv_rank}=$i ; 
     $self->test_exon($test_exons[$i]);
   } 

   my @homolog_exons = @{ $homolog_trans->get_all_translateable_Exons() };   
   for (my $i=0 ; $i < scalar(@homolog_exons) ; $i++ ) {    
     $homolog_exons[$i]->{_jhv_rank}=$i ; 
   } 
   my $max_t_rank = scalar(@test_exons) ; 
   my $max_h_rank = scalar(@homolog_exons) ; 

   for ( my $i=0; $i<$max_h_rank; $i++)  {   
     for ( my $j=0; $j<$max_t_rank; $j++ ) {    
          $self->hit_matrix($i,$j,0); 
     }
   }  

   my %exon_counterparts ; 
   my $var ;  

   # write all test exons into one big file to speed exonerate up ...

    my $test_exons_filename = "/tmp/jhv_". $$  . ".test.exon_seqs";   
    open (FH,">$test_exons_filename") || die "Cant write file $test_exons_filename\n" ; 
    for my $te ( @test_exons ) {  
      (my $seq=$te->seq->seq )=~ s/(.{1,60})/$1\n/g;
      print FH ">". "test_exon_rank " . $te->{_jhv_rank} . "\t" . $te ."\t".$te->hashkey . " " .  $te->length . "\n$seq";
    } 
    close(FH) ;  


   for my $he ( @homolog_exons) {   
      $self->exonerate_he_exon_against_test_exons($he, $test_exons_filename ) ;    
   }
   system("rm $test_exons_filename");
 
   my @hit_exons = @{$self->hit_matrix}; 

   #my $length_A = $#hit_exons ;
   #my $length_B = $#{$hit_exons[0]} ;  

   my $length_A = scalar(@homolog_exons) ; 
   my $length_B = scalar(@test_exons) ; 

   #print "length_A : $length_A " . scalar(@homolog_exons) . "\n" ; 

  if ( $self->debug ) {  
    print "\n\n" ; 
    for ( my $i=0;$i < $length_A; $i++) { 
       printf "HIT_MATRIX :\t"; 
       for (my $j = 0; $j < $length_B; $j++ ) {  
          printf "%6s", $hit_exons[$i][$j] ; 
       } 
         print "\t" . $homolog_exons[$i]->stable_id . "  $i\n" ; 
     }
     print "\n\n\n";  
  }
  # compute F-Matrix 
  
  my @f_matrix ;   

  for ( my $i=0;$i < $length_A; $i++) { 
     $f_matrix[$i][0] =  $hit_exons[$i][0]; 
  }  

  for (my $j = 0; $j < $length_B; $j++ ) {  
     $f_matrix[0][$j] =  $hit_exons[0][$j]; 
  }    


  for ( my $i=1; $i < $length_A; $i++) { 
       #print "i = $i\t" ; 
    for ( my $j= 1; $j <$length_B; $j++ ) {    
       #print "j : $j\n" ;
       my $choice_1 = $f_matrix[$i-1][$j-1] + $hit_exons[$i][$j];
       my $choice_2 = $f_matrix[$i-1][$j]   + $hit_exons[$i][$j];
       my $choice_3 = $f_matrix[$i][$j-1]   + $hit_exons[$i][$j];  
       my $max = max($choice_1,$choice_2,$choice_3);  
       #print "i = $i  j=$j  $choice_1,$choice_2,$choice_3 -- > max $max\n" ; 

       $f_matrix[$i][$j]=max($choice_1,$choice_2,$choice_3) ;   
    }
  }   

  if ( $self->debug ) { 
    for ( my $i=0; $i < $length_A; $i++) { 
      for ( my $j=0; $j < $length_B; $j++ ) {   
         printf "%5s", $f_matrix[$i][$j] ; 
      } 
        print  "\t" .$homolog_exons[$i]->stable_id . "\n" ;  
    }  
  }
  # now trace-back and see which test_gene has the most exons which align with the highest score ..... 
  
  $max_score_f_matrix = $f_matrix[$length_A-1][$length_B-1] ; 
  if ( $self->debug ) { 
    print "\n\nmax_score for aligning test_transcript to homolog is : ".$homolog_trans->stable_id." : $f_matrix[$length_A-1][$length_B-1]\n\n\n";  
    print "\n\n----- " . $homolog_trans->stable_id . "  VS  $test_transcript      ---------\n\n" ; 
  }  

  # little routine to check the pairing of the exons   OLD 


   my $exon_hits = 0 ;   
   my %exon_pairs = %{$self->{_exon_pairing}};

   HOMO_EXON: for my $he ( @homolog_exons ) {      
      unless ( $exon_pairs{$he} ) {  
         print $he->stable_id . "  DOES NOT HAVE A MATCH \n" if $self->debug ; 
         next HOMO_EXON; 
      } 

      print "\n" if ( scalar(  @{$exon_pairs{$he}}) > 1 ) ;  

      #print "\nnr_of_te_exons aligning M1: " . scalar(@{$exon_pairs{$he}}) . "\n" ; 
      #print "nr_of_te_exons aligning M2: " . $self->homolog_exon_alignments($he) . "\n\n" ; 

      if ( scalar(  @{$exon_pairs{$he}}) ==  1 ) { 
         $exon_hits++ ;  
      }elsif ( scalar(  @{$exon_pairs{$he}}) >  1 ) { 
        print "EXON HAS MORE THAN ONE HIT : " . $self->homolog_exon_alignments($he) . "\n" if $self->debug ; 
      } 

      for my $dog_exon ( @{$exon_pairs{$he}}) {  
         print $he->stable_id . "  ALIGNS  : " if ( $self->debug );
         my $diff = $he->length - $dog_exon->length ;  
         if ( $self->debug ) { 
           print " " . $dog_exon->seq_region_start . " " . $dog_exon->seq_region_end . " " .$dog_exon->length . " ". $he->length . " diff: " ;  
           print $diff . "\n" ;  
         } 
      }
      print "\n" if ( scalar(  @{$exon_pairs{$he}}) > 1 ) ; 
   }    

   print "\n\nComparison : " .  $homolog_trans->stable_id . "\t$test_gene\n" if $self->debug ; 
   #print "dog-transcripts has " . $exon_hits . " exons out of " . scalar(@test_exons) . " exons aligned : " . ( $exon_hits / scalar(@test_exons)) . "\n"  ; 
   #print "nr hom exons  : " . scalar(@homolog_exons) . "\n" ; 
   #print "nr dog exons  : " . scalar(@test_exons) . "\n" ; 
   #print "aligned exons : " . $exon_hits . "\n" ; 
   #print "hits vs homologs-coverage ratio : " . ( $exon_hits / scalar(@homolog_exons)) . "  (nr_aligned_dog_exons / nr_of_exons_in_homolog) \n";  
   #print "hits vs all-exons-in-dog-ratio  : " . ( $exon_hits / scalar(@test_exons)) . "  (nr_aligned_dog_exons / nr_of_all_dog_exons )  \n";  
   #print "exon_hits : $exon_hits    \n\n\n" ; 
   print "f_matrix max_score for aligning test_transcript to homolog is : ".$homolog_trans->stable_id." : $f_matrix[$length_A-1][$length_B-1]\n"
    if $self->debug ; 

   $self->register_testex_vs_homex_coverage($test_gene, ( $exon_hits / scalar(@homolog_exons)),$homolog_trans) ;  

   # nr of aligned dog exons vs all exons in dog - somehow a measurement how similar the structures are .... 
   $self->register_nr_aligned_exons_vs_all_exons_in_gene($test_gene, ( $exon_hits / scalar(@test_exons)),$homolog_trans) ;   

   $self->register_f_matrix_score($test_gene, $max_score_f_matrix, $homolog_trans) ; 
   
   $self->{_exon_pairing} = {};
  return $max_score_f_matrix ; 
}





 #
 # nr_aligned_dog_exons / nr_of_exons_in_homolog -measurement how much of the homolog is covered 
 #
sub register_testex_vs_homex_coverage{
  my ( $self, $test_gene, $val, $homolog_trans ) = @_ ;  
  
  if ( defined $val ) {   
   print "adding values for : $test_gene\n" if $self->debug ; 
   ${$self->{_register_testex_vs_homex_coverage}}{$test_gene}{$homolog_trans}{val} = $val ;  
   ${$self->{_register_testex_vs_homex_coverage}}{$test_gene}{$homolog_trans}{htrans} = $homolog_trans ; 
  } 
  return  ${$self->{_register_testex_vs_homex_coverage}}{$test_gene} ; 
} 



 #
 # nr of aligned dog exons vs all exons in dog - somehow a measurement how similar the structures are .... 
 #
sub register_nr_aligned_exons_vs_all_exons_in_gene{ 
  my ( $self, $test_gene, $val ,$homolog_trans) = @_ ;  
  
  if ( defined $val ) {  
    ${$self->{_register_nr_aligned_exons_vs_all_exons_in_trans}}{$test_gene}{$homolog_trans}{val} = $val ;  
    ${$self->{_register_nr_aligned_exons_vs_all_exons_in_trans}}{$test_gene}{$homolog_trans}{htrans} = $homolog_trans ; 
  } 
  return  ${$self->{_register_nr_aligned_exons_vs_all_exons_in_trans}}{$test_gene} ; 
} 


 
sub register_f_matrix_score { 
  my ($self, $test_gene, $max_score_f_matrix , $homolog_trans) = @_; 
  
  if ( defined $max_score_f_matrix) {  
    ${$self->{_register_f_matrix_score}}{$test_gene}{$homolog_trans}{val} = $max_score_f_matrix;  
    ${$self->{_register_f_matrix_score}}{$test_gene}{$homolog_trans}{htrans} = $homolog_trans ; 
  } 
  return  ${$self->{_register_f_matrix_score}}{$test_gene} ; 
} 




sub write_exon_to_file {  
   my ($exon) = @_ ;   

   my $f1= "/tmp/query_".$exon->stable_id . "_" .$$."_exon.query" ;  
   open (FH,">$f1") || die "Cant write file $f1\n" ;
   my $seq = $exon->seq->seq;
   $seq=~ s/(.{1,60})/$1\n/g;
   print FH ">homolog_exon_rank\t". $exon->{_jhv_rank} . "\t" . $exon->stable_id . "\n" . $seq;
   close(FH) ;
   return $f1 ; 
} 



sub exonerate_he_exon_against_test_exons { 
  my ($self,$he, $test_exon_filename ) = @_ ;   

  my @all_hits_for_he;  
   
  # write sequence of homolog exon to file seq to file 
  my $homolog_exon_seq_file = write_exon_to_file($he) ; 

  #print "Trying to align all test_exons against " . $he->stable_id . "\t" . $he->seq_region_start . "\t" . $he->seq_region_end . "\n" ; 

  my $result_file = "/tmp/jhv_". $he->stable_id . "_" . $$ . ".exonerate.alignment"; 

  # use option to find all exons ..... ie. match ENSMUSE00000185591  vs 
  
  my $cmd = "/software/ensembl/bin/exonerate-1.4.0 -q $homolog_exon_seq_file -t $test_exon_filename -m affine:local --bestn 1 > $result_file"; 
  system("$cmd") ; 

  # routine to parse the results out of the exonerate alignment file and write it to hit-matrix
  my $at_least_one_hit = $self->parse_exonerate_results($result_file,$he) ;   
  system("rm $result_file");

  if ( $at_least_one_hit ) { 
    #print STDERR "run produced hit : " . $he->stable_id . "\n" ; 
  } else { 
    #print STDERR "NO hit found for " . $he->stable_id . "  running with --exhaustive \n" ; 
  
  } 


  unless ( defined $at_least_one_hit ) {  
    my $cmd = "/software/ensembl/bin/exonerate-1.4.0 -q $homolog_exon_seq_file -t $test_exon_filename -m affine:local --bestn 1 --exhaustive  > $result_file"; 
    system("$cmd") ;  

    my $exhaustive_hit = $self->parse_exonerate_results($result_file,$he) ;  
    system("rm $result_file"); 

    unless ( $exhaustive_hit ) {   
      print STDERR " NO HIT even with  EXHAUSTIVE option for : " . $he->stable_id . "\n" if $self->debug;
    }
  } 
  system("rm $homolog_exon_seq_file");

}






sub parse_exonerate_results { 
  my ($self,$results_file,$he) = @_ ; 

  my $hit = 0 ;
  my $query_length =0  ;
  my $raw_score ;
  my @alignment ;  
  my $vulgar= ""; 
  my %results_hash; 

  my @result_matrix; 

  # we exonerate the sequence of the homolog exon against all coding exons of the test gene.
  # we can have multiple hits, i.e. a homolog exon can align to ONE single test-exons in different 
  # ways with different score. we only count one hit ( the best hit ) per test exon ... 

  my ($he_exon_rank, $he_stable_id, $test_exon_rank, $score,$test_exon_ref); 

  my $main_hit ; 

  open (F,"$results_file") || die " can't read $results_file\n" ;
     foreach my $l (<F>){ 
        if ($l=~/C4 Alignment:/) { 
          $hit = 1; 
          $main_hit = 1 ;
          $he_exon_rank = undef ; 
          $he_stable_id = undef ; 
          $test_exon_rank = undef ;  
          $test_exon_ref = undef ; 
          $score =undef ; 
        }
        push @alignment , $l  if $hit ;

        my @col = split/\s+/,$l;

        if ( $l=~m/Query:/ && ! defined($he_exon_rank)){
           $he_exon_rank = $col[3];  
           $he_stable_id = $col[4];  
           #print "\n\nline :$l\n"; 
           #print "QUERY : $he_exon_rank    $he_stable_id \n" ; 
        }elsif ( $l=~m/Target:/ && !defined ($test_exon_rank)){ 
           $test_exon_rank = $col[3];  
           $test_exon_ref  = $col[4];
           #print "\n\nline :$l\n"; 
           #print "TARGET : $test_exon_rank\n" ; 
        }elsif ( $l=~m/Raw score:/){ 
           $score = $col[3]; 
           #print "\n\nline :$l\n"; 
           #print "RAW SCORE : $score\n" ; 
        }   
        if ( defined $he_exon_rank && defined $test_exon_rank &&  defined $score ) {   
          #print "\ngot_a_hit : $he_exon_rank $test_exon_rank $score for:  $he_stable_id \n" ; 
          $self->hit_matrix($he_exon_rank,$test_exon_rank,$score);

          # now store the pair which aligns ( he : $te ) 
          my $test_exon_object = $self->test_exon($test_exon_ref);      
          $self->exon_pairs($he, $test_exon_object ) ;  
          $self->homolog_exon_alignments($he,1);  

          $he_exon_rank = undef ; 
          $he_stable_id = undef ; 
          $test_exon_rank = undef ; 
          $score = undef ; 
        }  
     } 
   close(F); 
    #print @alignment ; 
   unless ( $hit ) {  
        print "homolog exon : " . $he->{_jhv_rank} . " " .$he->stable_id  . " did not align at all against any test exon\n" if $self->debug; 
   }  
   return $main_hit ; 
}



sub hit_matrix { 
  my ($self, $he_exon_rank,$test_exon_rank,$score) = @_ ;  

  if ( defined $he_exon_rank && defined $test_exon_rank ) { 
     $self->{_hit_matrix}[ $he_exon_rank ][ $test_exon_rank ]= $score; 
  } 
  return $self->{_hit_matrix};
} 




sub attach_slice { 
   my ($self, $new_gene, $slice )  = @_ ;  

   my @nt = @{$new_gene->get_all_Transcripts} ;  
   throw("more than one transcript in gene structure!\n") if scalar(@nt)>1 ;     

    unless ( $slice ) { 
      $slice = $self->slice  ; 
    } 

   for my $t ( @nt ) {
     $t->slice($slice) unless $t->slice() ;  
     $t->analysis($self->analysis) ;  
     if ( $t->translation() ) { 
       $t->translation->stable_id(undef) ;  
      } 
     for my $tsf ( @{ $t->get_all_supporting_features } ) {   
        $tsf->slice($slice) unless $tsf->slice;      
        $tsf->analysis($self->analysis) unless $tsf->analysis ; 
      } 
     for my $e ( @{ $t->get_all_Exons } ) {   
       $e->slice($slice) unless $e->slice;
       for my $sf ( @{$e->get_all_supporting_features} ) { 
          $sf->analysis($self->analysis) unless $sf->analysis ; 
         $sf->slice($slice) unless $sf->slice() ; 
       } 
     } 
     if ( $t->translation() ) { 
        $t->translation->stable_id(undef) ; 
     } 
  }
  return $new_gene ; 
} 


sub print_best_targetted_genes {
  my ($aref ) = @_ ; 
  my @best_targetted_genes = @{$aref} ;  

  print "HERE's the output of BestTargetted :\n" ;  
  print "====================================\n\n"; 
   
  for my $g ( @best_targetted_genes ) {
   print_pos($g,"gene " ) ;  
       for my $t ( @{$g->get_all_Transcripts} ){ 
         my @tsf = @{$t->get_all_supporting_features} ; 
         print_pos($t,"  trans ", \@tsf ) ;    
   
         for my $e ( @{$t->get_all_Exons } ) {  
           print_pos($e,"     exon " ) ;  
         } 
         for  ( my $nr=0; $nr< scalar(@{$t->get_all_Exons }); $nr++  ) {   
           my @exons = @{ $t->get_all_Exons} ; 
           print_pos($exons[$nr],"     exon " ) ;  
           if ( $nr+1 < scalar(@exons) ) {  
             print_pos(new Bio::EnsEMBL::Intron ( $exons[$nr], $exons[$nr+1])," INTRON ") ;  
           }
         } 
         for my $i ( @{$t->get_all_Introns} ) {   
           print_pos($i, "         INTRON ");  
           print"\n";
         }
       }
      print"\n";
  }    
}



#
# make an array of non-redundant genes with the exon-hash-keys. 
# also check that all genes which will be removed from the set 
# have sf which support the same homolog transcripts 
#

sub remove_redundant_genes {   
  my ($self, $genes) = @_ ; 

  print "removing non-redundant genes\n" if $self->debug;  

  my @non_red_genes ;
  my %tmp ;

  for my $g ( @$genes) { 
    if ( scalar(@{$g->get_all_Transcripts}) > 1 ){  
      throw("gene has more than one transcript - this is wrong .... : " . $g->biotype); 
    } 
  }  

  # make an array of non-redundant genes with the exon-hash-keys
 
  for my $g ( @$genes ) {      
     #print_whole_gene($g); 
     for my $t (@{ $g->get_all_Transcripts } ) {  
       my @ex= @{ $t->get_all_Exons } ;
       my $string ;
       for my $e (@ex) {
         my $exon_hk = $e->hashkey ;
         $string.=$exon_hk ;
       } 
       # same exon_hash-keys = same gene 
       push @{$tmp{$string}}, $g ;
     }
  }  


  for my $k (keys %tmp){  

    my @g = @{$tmp{$k}} ;  

    # now check if the supporting-feature of the gene supports the same homolog transcripts 

    my %protein_align_feature_acc;  
    my %stable_ids_of_supported_transcripts;  
    my %gene_to_supported_transcripts ; 
    for my $g ( @g ) {   
      # get a list of uniq supporting features for the transcript 
      #
      #  A GENE IS ONLY REDUNDANT IF IT HAS THE SAME sf WHICH SUPPORT THE SAME TRANSCRIPTS  
      #
      #  do we really want this ?   NO !!!
      #
       my $t = ${$g->get_all_Transcripts}[0]; 
       my @e_sf = map { $_->get_all_supporting_features } @{$t->get_all_Exons}  ;  

       for my $sf_array ( @e_sf) {  
         for ( @$sf_array) {  
           if (ref($_)=~m/Bio::EnsEMBL::DnaPepAlignFeature/){   
             $protein_align_feature_acc{$_->hseqname} = 1 ;   
           }
         } 
       }   
       for my $hseqname ( keys %protein_align_feature_acc ) { 
         my @supported_transcripts = sort @{$self->sf_to_gene($hseqname)} ;
         for my $trans ( @supported_transcripts ) { 
           $stable_ids_of_supported_transcripts{$trans->stable_id} = 1 ; 
         }
       }
       my @supported_hom_stable_ids = sort keys %stable_ids_of_supported_transcripts ; 
       $gene_to_supported_transcripts{join("_",@supported_hom_stable_ids)} = $g;
    }
    for my $sid ( keys %gene_to_supported_transcripts ) {
      print "sid : $sid  ----> $gene_to_supported_transcripts{$sid} \n" if $self->debug ; 
      push @non_red_genes, $gene_to_supported_transcripts{$sid}; 
    }
      #push @non_red_genes, $g[0];
  }  
  #print "\nNON-redundant genes : " . scalar(@non_red_genes) . " left after removal of redundant strucutres ....(ALL: ".scalar(@$genes) ." )\n\n" ; 
  #for ( @non_red_genes){
  #   print_whole_gene($_);
  #}
  print "\nNON-redundant genes : " . scalar(@non_red_genes) . " left after removal of redundant strucutres ....(ALL: ".scalar(@$genes) ." )\n\n"  
   if $self->debug ; 
  return \@non_red_genes ; 
}













sub check_exon_overlap_between_new_gene_and_src_gene { 
  my ($self, $non_red_genes) = @_ ;

  my %scoring ; 
  print "check_exon_overlap_between_new_gene_and_src_gene " . scalar(@$non_red_genes) . "\n"  ;  

  # now compare the different, non-redundant  structures against the original gene .   
  
  my @src_genes = @{$self->src_genes} ;  
  my $src_gene = $src_genes[0]; 
  
  my @src_trans = @{$src_gene->get_all_Transcripts} ; 
  my @src_exons = @{ $src_trans[0]->get_all_translateable_Exons };  

  for my $ng ( @$non_red_genes ) {    

    my @tr = @{$ng->get_all_Transcripts}; 
    my @test_exons = @{ $tr[0]->get_all_translateable_Exons } ;  

    $scoring{$ng}{gene_object} = $ng ; 
    $scoring{$ng}{trans_object} = $tr[0] ; 
    $scoring{$ng}{src_exons_which_are_not_overlapped} = 0 ; 
    $scoring{$ng}{new_exons_without_overlap} = []; 
    $scoring{$ng}{overlapped_exon} = []; 
    $scoring{$ng}{nr_non_overlapping_exons} = 0 ; 



    for my $te ( @test_exons ) { 
       # looping over all new exons first ....  
       for my $se ( @src_exons ) {      
         if ( $te->seq_region_end >= $se->seq_region_start and $te->seq_region_start <= $se->seq_region_end )  {  
              #print "overlapping exons :-)\n" ; 
              #print "\nSRC_EXON  : " . $se->seq_region_start . " - " . $se->seq_region_end . "\n";
              #print "TEST_EXON : " . $te->seq_region_start . " _ " . $te->seq_region_end . "\n" ;  

              # overlap - NOW there could be some szenarios 
              # - it could be that the new exon is overlapping some shorter exons of the old prediction 
              # - it could be that one of the exons of the old gene is skipped 
              # - it could be that 
               
              push @{$scoring{$ng}{overlapped_exon}}, $se ;   
         } 
           #print "\n\n" ; 
       } 
    }   
   
    # this bit calculates how many original exons are not overlapped by the new exons.  
    my %se_exons ;  
    for my $se ( @src_exons, @test_exons  ) {        
        if (exists $se_exons{$se} ) {  
          die("this is weird - the exon keys are the same, this suggest something has gone wrong\n") ; 
        }  
        $se_exons{$se}=$se;
    } 

    my %overlapped_original_exon ;   

    print_whole_gene($ng) if $self->debug; 

    for my $se ( @src_exons ) {      
        $overlapped_original_exon{$se}{overlapping_new_exons}=[];
        for my $te ( @test_exons ) {   
         if ( $se->seq_region_end >= $te->seq_region_start and $se->seq_region_start <= $te->seq_region_end )  {   
            #print "\nmarking src_exon as OVERLAPPED       :: " .  $se_exons{$se}->seq_region_start . "\t" . $se_exons{$se}->seq_region_end."\n";
            #print "        new_exon                     :: " .  $se_exons{$te}->seq_region_start . "\t" . $se_exons{$te}->seq_region_end."\n\n";
            $overlapped_original_exon{$se}{is_overlapped}=1; 
            push @{$overlapped_original_exon{$se}{overlapping_new_exons}}, $te ; 
         }else{
            # if the exon has already been overlapped we don't change it's value  
            unless ( exists $overlapped_original_exon{$se}{is_overlapped} ) {  
              #print "marking org_exon as non-overlapped : " .  $se_exons{$se}->seq_region_start . "\t" . $se_exons{$se}->seq_region_end."\n";
              $overlapped_original_exon{$se}{is_overlapped}=0;
            } 
            if ($overlapped_original_exon{$se}{is_overlapped} == 1 ) {  
              # do nothing 
            }else{ 
              $overlapped_original_exon{$se}{is_overlapped}=0; 
            }
         } 
       }
    }   

    my $non_overlapped_exons = 0 ; 
    for my $se ( @src_exons ) {
       if ( $overlapped_original_exon{$se}{is_overlapped} == 0 ) {    
         # exon is not overlapped 
         $non_overlapped_exons++; 
       }   

      if ( scalar ( @{$overlapped_original_exon{$se}{overlapping_new_exons}}) > 1 ) {  
        print "WARN! Oh,, we might have a look at this - one old exon is overlapped by more than one new exon\n" ; 
        print "orig. exon : " . $se_exons{$se}->seq_region_start . "\t" . $se_exons{$se}->seq_region_end . "\n" ;  
        for my $te ( @{$overlapped_original_exon{$se}{overlapping_new_exons}} ) { 
           print " new. exon : " . $se_exons{$te}->seq_region_start . "\t" . $se_exons{$te}->seq_region_end . "\n" ;  
        } 
      } 
    } 


    # we should calculate as well how many original exons are overlapped by more than one new exon as this is bad, too
    # as suggests frameshifts ..... 
    #   <----------old_exon------------>
    #      <new>  <new_ex> <new_ex>         <new_ex> 

    my $new_exons_without_overlap = 0 ; 
    my %ymp ; 
    @ymp{@{$scoring{$ng}{overlapped_exon}}}=1 ; 
    for ( keys %ymp ) { 
     $new_exons_without_overlap++; 
    } 
    $scoring{$ng}{new_exons_without_overlap} = $new_exons_without_overlap ; 
    $scoring{$ng}{src_exons_which_are_not_overlapped} = $non_overlapped_exons ; 
   } 



  # now find best-scored gene(s) 

  for my $ng ( keys %scoring ) {
    
    #print "Gene : \n====================\n" ;
    #print_whole_gene($scoring{$ng}{gene_object}) ; 
    #print "nr_overlapped_exons                : " . scalar(@{$scoring{$ng}{overlapped_exon}}) . "\n" ;  
    #print "overlap_ratio                      : " . scalar(@{$scoring{$ng}{overlapped_exon}}) / scalar(@{$src_trans[0]->get_all_translateable_Exons} ) . "\n"  ;    
    $scoring{$ng}{overlap_ratio} = scalar(@{$scoring{$ng}{overlapped_exon}}) / scalar(@{$src_trans[0]->get_all_translateable_Exons} );
    #print "overlap_ratio  numbers             : " . scalar(@{$scoring{$ng}{overlapped_exon}}) . " / " . scalar(@{$src_trans[0]->get_all_translateable_Exons} ) . "\n"  ;   
    #print "translation_ratio                  : " . length($scoring{$ng}{trans_object}->translateable_seq) / length($src_trans[0]->translateable_seq ) . "\n" ; 
    #print "new_exons_without_overlap          : " . $scoring{$ng}{new_exons_without_overlap} . "\n" ; 
    # best gene seems to be the one with : 
    #  - highest nr of overlapped exons  
    #  - highest overlapped_ratio
    #  ??? nr of fixed exons ?  
  }    

  my $most_number_of_exons_fixed = -1000000; 
  my $highest_overlapping_ratio = -100000; 

  # find the gene with the maximum number of exons fixed 
  for my $ng ( keys %scoring ) {   
      if ( $scoring{$ng}{new_exons_without_overlap} >= $most_number_of_exons_fixed ) {    
        $most_number_of_exons_fixed = $scoring{$ng}{new_exons_without_overlap} ; 
      } 
      if ( $scoring{$ng}{overlap_ratio}  >= $highest_overlapping_ratio ) {     
        $highest_overlapping_ratio = $scoring{$ng}{overlap_ratio} ;
      } 
  }    


  my %best_genes ; 

  my $gene_misses_src_exon ; 
  my $nr_exons_not_overlapped = 100000 ; 

  for my $ng ( keys %scoring ) {
    if ( $scoring{$ng}{src_exons_which_are_not_overlapped} == 0  ) { 
      if ( $scoring{$ng}{overlap_ratio} == 1 ) { 
         $best_genes{$ng} = $scoring{$ng}{gene_object} ;    
         print_whole_gene($scoring{$ng}{gene_object}) ; 
         my $biotype =  "b_ratio_" . $scoring{$ng}{overlap_ratio} ; 
         $best_genes{$ng}->biotype("$biotype"); 
      } else {  
         print "skipping gene as it's overlap_ratio is not == 1. It's : $scoring{$ng}{overlap_ratio} \n" ;  
      }      
    } else {   
      if ( $nr_exons_not_overlapped  >  $scoring{$ng}{src_exons_which_are_not_overlapped} ) {  
           $nr_exons_not_overlapped =  $scoring{$ng}{src_exons_which_are_not_overlapped} ; 
           $gene_misses_src_exon = $scoring{$ng}{gene_object} ; 
      }  
      print "skipping gene - not all source exons are overlapped\n" ; 
    }  
  }    

  my @best_genes = values %best_genes ;   
  
  my $very_best_gene ;  

  if ( scalar(@best_genes) > 0 ) { 
  
    #
    # find gene with longest translation ......   
    # up to here it all worked fine but the longest translation is not always the best ..........
    # 
    
    print "\nNow selecting gene with longest translation out of " . scalar(@best_genes) . " ...\n" ;  
  
    my $max_translation_length = 0 ;  
  
    for my $bg ( @best_genes ) {   
      my $bt = ${$bg->get_all_Transcripts}[0];
         print "comparing translation length ... " .  length($bt->translateable_seq) . " vs max. length : " . $max_translation_length . "\n" ; 
      if ( length($bt->translateable_seq) > $max_translation_length ) {   
         print "translation is longer - new very best gene is : \n" ; 
         $very_best_gene = $bg ; 
         $max_translation_length = length($bt->translateable_seq); 
         print_whole_gene($very_best_gene) ;   
      }     
    }  
    if ( $very_best_gene ) { 
      print "the FINAL VERY best gene is : \n" ;  
      print_whole_gene($very_best_gene) ;   
    } else {  
      print "no very best final gene found\n" ; 
    } 
     #this is just to not loose any genes and store them for debugging 
    my @other_genes ; 
    for my $g ( @best_genes) {  
      unless ($g eq $very_best_gene ) {   
         $g->biotype("not_selected") ; 
         push @other_genes , $g ; 
      } 
    } 
    $self->other_genes(\@other_genes) ; 
    $very_best_gene->biotype("very_best_gene_old_method"); 
  }
  unless ( $very_best_gene ) {  
    print " no very best gene could be selected with method 1 as not all src_exons are overlapped\n" ; 
  }
  return $very_best_gene ; 
}








sub check_if_false_intron_has_been_removed_in_new_trans {  
   my ($self, $homolog_trans, $new_gene) = @_ ;   

   my $src_gene = $self->src_gene ;  

   my @checked_genes; 

   my @transcripts_to_check = @{$new_gene->get_all_Transcripts } ;  
   throw(" Gene has more than one transcript - this is wrong\n") if scalar(@transcripts_to_check) > 1 ;  

   my $trans_to_check = $transcripts_to_check[0]; 

   my %tmp = %{$self->recovered_introns_of_transcript};  

   my @simple_features_to_be_overlapped_by_transcript = @{ $tmp{$homolog_trans->stable_id }} ; 

   my $overlapped_simple_features_in_new_gene = 0 ; 
   my $fully_removed_simple_features_in_new_gene = 0 ; 

   my $gaps_removed = 0 ;  
   my $full_gaps_removed = 0 ;  

     SIMPLE_FEATURE: for my $sf ( @simple_features_to_be_overlapped_by_transcript ) {  
     
       if ( $self->debug) { 
         print "\n\nChecking introns of transcript ...\n" ;    
       }
       my $intron_removed;   
  
       for my $e ( @{$trans_to_check->get_all_Exons} ) {  
               if ( $self->debug) { 
                 print "checking .......\n" ;  
                 print "EXON   : " . $e->seq_region_start . "\t" .  $e->seq_region_end . "\t" . $e->seq_region_name."\t". ( $e->seq_region_end - $e->seq_region_start)."\n" ; 
                 print "GAP    : " . $sf->seq_region_start . "\t" . $sf->seq_region_end . "\t". $sf->seq_region_name."\t".( $sf->seq_region_end - $sf->seq_region_start)."\n" ;  
                 print $e->seq_region_end." >= " .$sf->seq_region_start ." && ".$e->seq_region_start . " <= ".$sf->seq_region_end . " ???\n" ; 
               }  
                 if ( $e->seq_region_end  >=   $sf->seq_region_start  &&  $e->seq_region_start  <= $sf->seq_region_end) {   
                    # we found an exon which overlaps the buyggy intron
                                                                              
                    #                        +++++++++++++++++++++++++++++++++++++++++++++++++++++
                    #                   sf_start 10                                           sf_end 20 
                    #       e_start================================  e_end  



                    


                    $overlapped_simple_features_in_new_gene++; 
                    if ( $self->debug) { 
                      print "Intron is overlapped by exon ...\n" ;    
                      print "exon: " . $e->seq_region_start . "\t" .  $e->seq_region_end . "\t" . $e->seq_region_name."\t". ( $e->seq_region_end - $e->seq_region_start)."\n" ; 
                      print "GAP : " . $sf->seq_region_start . "\t" . $sf->seq_region_end . "\t". $sf->seq_region_name."\t".( $sf->seq_region_end - $sf->seq_region_start)."\n" ;   
                    }    

                    $gaps_removed++ ; 
                    $intron_removed = 1 ; 
                    #
                    #  full overlap if : 
                    # 
                    # another check if the intron ( simple_feature ) is fully overlapped by the exon or not 
                    #
                    #                      start--EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE--end  
                    #                                  start-IIIIIIIIIIIIIIIIIIIIII-end

                    if ( $e->seq_region_start <= $sf->seq_region_start && $e->seq_region_end >=$sf->seq_region_end ) {   
                      $fully_removed_simple_features_in_new_gene ++ ; 
                      $full_gaps_removed++ ; 
                      $trans_to_check->biotype("intron_removed_".$gaps_removed."_full_".$full_gaps_removed);
                      $new_gene->biotype("intron_removed_".$gaps_removed."_full_".$full_gaps_removed); 
                    
                      if ( $self->debug )  { 
                        print "intron is fully overlapped by exon !!! \n" ; 
                        print "Exon  : " . $e->seq_region_start . 
                        "--eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee--" .
                        $e->seq_region_end . " (end)\n" ; 
                        print "  GAP :                       " .$sf->seq_region_start . "--IIIIIIIIIII--" . $sf->seq_region_end 
                       . " (end)\n" ;  
                      }
                    }  else { 
                      if ( $self->debug )  { 
                        print "intron is partially overlapped by exon !!! \n" ; 
                        print "Exon  : " . $e->seq_region_start . "--eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee--" .  
                        $e->seq_region_end . " (end)\n" ; 
                        print "  GAP :                       " .$sf->seq_region_start . "--IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII--" . 
                        $sf->seq_region_end . " (end)\n" ;  
                      }
                    }
                 } else{
                    print "Intron was NOT removed in gene ! \n" if $self->debug ; 
                    $trans_to_check->biotype("intron_removed_$gaps_removed");
                    $new_gene->biotype("intron_removed_$gaps_removed");
                 } 
       }  

      if ( $intron_removed ) {   
          if ( $self->debug ) { 
            print "intron removed in transcript :\n" ;  
            print_whole_trans($trans_to_check);  
          }
          $trans_to_check->biotype("intron_removed_".$gaps_removed."_full_".$full_gaps_removed);
          $new_gene->biotype("intron_removed_".$gaps_removed."_full_".$full_gaps_removed);
          next SIMPLE_FEATURE ; 
       }   

     }  
 
     # check if all introns are overlapped by new transcript  

     if ( scalar(@simple_features_to_be_overlapped_by_transcript) == $overlapped_simple_features_in_new_gene ) { 
       print "All introns are overlapped by gene\n" if ( $self->debug) ; 
       $trans_to_check->biotype("intron_removed_all") ; 
       $new_gene->biotype("intron_removed_all") ; 
     } 
     if ( scalar(@simple_features_to_be_overlapped_by_transcript) == $fully_removed_simple_features_in_new_gene ) { 
       print "All introns are overlapped by gene : $fully_removed_simple_features_in_new_gene\n" if ( $self->debug ) ; 
       $trans_to_check->biotype("intron_removed_all_full") ; 
       $new_gene->biotype("intron_removed_all_full") ; 
     } 
  if ( $self->debug ) {  
    print "full    gaps removed in this gene : $full_gaps_removed\n" ; 
    print "overlap gaps removed in this gene : $gaps_removed\n" ;  
  } 
  if ( $full_gaps_removed ==  0  && $gaps_removed == 0 )  {   
     print "intron has not been removed\n" if $self->debug ; 
     $new_gene->biotype("intron_NOT_removed"); 
  } 

  return $new_gene ; 
} 









sub get_aligned_intron_coord {  
  my ( $self) = @_ ;  

  my $sfa = $self->db->get_SimpleFeatureAdaptor() ;    

  my @simple_features; 
  for my $db_id ( @{ $self->simple_feature_db_id } ) { 
   
    my $sf = $sfa->fetch_by_dbID($db_id) ; 
    if ( $self->debug ) { 
      print "aligned -intron coordinates are : " . $sf->seq_region_start . "\t" . $sf->seq_region_end . "\n"  ; 
      print "aligned -intron coordinates - chromosome are : " . $sf->slice->seq_region_name . "\n";  
      print "aligned -intron -id : " .  $sf->display_id . "\n" ; 
    }  
    push @simple_features, $sf ; 
  }
  return \@simple_features ; 
}

sub print_whole_gene { 
  my ($g) = @_ ;   
  print $g . "\n" ; 
  print_pos($g,"gene :" ) ;  
  for my $t ( @{$g->get_all_Transcripts} ){ 
    print_pos($t,"  trans :" ) ;   
    my @tsf = @{$t->get_all_supporting_features} ; 
    for ( @tsf ) {  
      print "   DISPLAY : " . $_->display_id . "\n" ; 
    } 
    for my $e ( @{$t->get_all_Exons } ) { 
      print_pos($e,"     exon : " ) ;  
    } 
  } 
  print "\n\n";
}

 
sub print_whole_trans { 
  my ($t) = @_ ;    
   print_pos($t,"  trans :" ) ;    
  
    my @tsf = @{$t->get_all_supporting_features} ; 
    for ( @tsf ) {  
      print "DISPLAY : " . $_->display_id . "\n" ; 
    } 
    for my $e ( @{$t->get_all_Exons } ) { 
      print_pos($e,"     exon : " ) ;  
    } 
  print "\n\n";
}


sub ex_genes { 
   my ( $self,$aref ) = @_ ;  
   if ( $aref ) {   
     push @{$self->{_ex_genes}}, $aref ; 
   } 
   return $self->{_ex_genes} ;  
} 

sub gw_genes { 
   my ( $self,$aref ) = @_ ; 
   if ( $aref ) {   
     push @{$self->{_gw_genes}}, $aref ; 
   } 
   return $self->{_gw_genes} ;  
} 

sub src_genes { 
   my ( $self,$aref ) = @_ ; 
   if ( $aref ) {   
     push @{$self->{_src_genes}}, $aref ; 
   } 
   return $self->{_src_genes} ;  
}  

sub simple_feature_db_id { 
   my ( $self,$aref ) = @_ ; 
   $self->{_sf_dbid} = $aref if $aref ; 
   return $self->{_sf_dbid} ;  
}


sub make_new_gene_from_source_transcript{ 
  my ( $self,$transcript_sid ) = @_ ;   

   print "Reading compara-config file : $$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE}\n";
   print "defined in : ./Config/GeneBuild/Orthologue_Evaluator.pm\n" ; 
   
   Bio::EnsEMBL::Registry->load_all($$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE},1,1); 

   # get gene-adaptor for other species as well 

   my $ta = Bio::EnsEMBL::Registry->get_adaptor($$MAIN_CONFIG{QUERY_SPECIES},"core","transcript");    
   Bio::EnsEMBL::Registry->set_disconnect_when_inactive();

   my $ga = Bio::EnsEMBL::Registry->get_adaptor($$MAIN_CONFIG{QUERY_SPECIES},"core","gene");    
   Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
   

   unless ( $ta && $ga ) {  
    throw(
          "could not get gene-adaptor for $$MAIN_CONFIG{QUERY_SPECIES} - check your " . 
          "registry file / your OrthologueEvaluator-config\n"
          ); 

   } 

   # collect dna-align features for source_gene for utr-addition later  
   my $source_gene = $ga->fetch_by_transcript_stable_id($transcript_sid) ;     
   $self->very_original_query_gene($source_gene) ; 
   my %cdna_acc ; 
   for my $e ( @{$source_gene->get_all_Exons} ) {  
      my @sf = @{ $e->get_all_supporting_features}; 
      for my $s ( @sf ) {  
        if ( ref($s)=~m/Bio::EnsEMBL::DnaDnaAlignFeature/) {  
          $cdna_acc{$s->hseqname}=1; 
        } 
      } 
   }    
   print scalar(keys %cdna_acc) . " dna_align_features found for src_gene which could be added as UTR\n" ;    
   $self->src_gene_dna_align_features([keys %cdna_acc]) ; 

   my $source_transcript= $ta->fetch_by_stable_id($transcript_sid) ;     

  

   unless($source_transcript){  
    throw("could not get transcript with stable_id $transcript_sid out of database :"  .$ta->db->dbname . "\@" . $ta->db->host. "\n")  ;  
   } else { 
     $self->src_transcript($source_transcript); 
   }

  my $new_gene = Bio::EnsEMBL::Gene->new(); 
  $new_gene->add_Transcript($source_transcript) ; 
  $new_gene->biotype("buggy_gene") ;  
  return $new_gene ; 
}  




sub get_adaptor {  
  my ($self)=@_;  
  print "USING HARDCODED VALUE TEST_DB as output db\n" ; 
  return $self->get_dbadaptor('TEST_DB')->get_adaptor("Gene"); 
} 


sub write_genomic_sequence_to_file { 
   my ( $self, $slice_name  ) = @_ ; 

     my @array = split(/:/,$slice_name);

     if(@array != 6) {
       throw("Malformed slice name [$slice_name].  Format is " .
          "coord_system:version:seq_region:start:end:strand");
     } 
   
     my ($cs_name, $cs_version, $seq_region, $start, $end, $strand) = @array;
     # we want to give exonerate a bit more genomic seq 
     if($start > $end){
       my $tmp_start = $end;
       $end = $start;
       $start = $tmp_start;
     }
     print STDERR "Have coordiantes :  ".$start." ".$end." \n";
     
     # sequence_padding
     my $new_start  = $start - 5000 ; 
     my $new_end    = $end   + 5000 ;
   
     if($new_start < 1){
       $new_start = 1;
     } 

     my $sliceadp = $self->db->get_SliceAdaptor(); 

     my $slice = $sliceadp->fetch_by_region($cs_name,$seq_region, $new_start,$new_end, $strand, $cs_version); 

     $self->slice($slice) ; 
 
     if($slice->end() > $slice->seq_region_length) {
       #refetch slice, since we went past end, second call is fast
       $new_end = $slice->seq_region_length();
       $slice = $sliceadp->fetch_by_region($cs_name, $seq_region,
                                           $new_start, $new_end,
                                           $strand, $cs_version);
     }
   
      print STDERR "Have ".$slice->name." sequence to run\n";
      $self->genomic($slice);  

      # my $seq = $slice->get_repeatmasked_seq(undef,1) ;  
      # make target(genomic) and query(protein) tmp sequencefiles
      my $genfile      = create_file_name("genomic" ,"fa" , "/tmp" );    
      #write_seqfile($seq, $genfile, 'fasta') ;

      my $seq = $slice->seq(); 
      (my $str= $seq) =~ s/(.{1,60})/$1\n/g;
      $str = ">genomic_seq_header\n$str" ;  
      open(FH, ">$genfile") ;
       print FH $str;
      close(FH) ;
  
      #write_seqfile($seqio, $genfile, 'fasta') ;

      $self->genomic_file ( $genfile ) ; 
}


# building  

sub build_hash_of_recoverd_introns_per_transcript {  
  my ($self ) = @_ ;  

  my @simple_feature_index ; 
  if ($self->simple_feature_index=~m/\_/ ) { 
    @simple_feature_index = split /\_/, $self->simple_feature_index ;   
  } else { 
   push  @simple_feature_index , $self->simple_feature_index ;  
  }


  # to save some compute time we only want 2 transcripts used which revoer the same intron  
  my %recov_simple_feat_to_transcript ; 

  my @simple_feature_ids = @{$self->simple_feature_db_id};
  my @stable_ids_of_homolog_transcripts = @{$self->homolog_transcript_sid()} ;  

  my %enst_to_simple_features;

  my $sfa = $self->db->get_SimpleFeatureAdaptor() ;     

  for ( my $i =0 ; $i <@simple_feature_index ; $i++ ) {   
    # print "\n\nSTABLE_ID : " . $stable_ids_of_homolog_transcripts[$i] . ":" ;  
     my $str = $simple_feature_index[$i] ; 
     my @indexes = split /\,/, $str ; 
     for my $s_index ( @indexes ) {    
       #print "s_index : $s_index " ;  
       #print " $s_index " ;  
       #print "simple_feature_id : " . $simple_feature_ids[$s_index] . "\n"   
       #
       push @{$recov_simple_feat_to_transcript{$str}},  $stable_ids_of_homolog_transcripts[$i]; 
     
       my $sf = $sfa->fetch_by_dbID($simple_feature_ids[$s_index]) ; 
  
       push @{$enst_to_simple_features{ $stable_ids_of_homolog_transcripts[$i]}} , $sf ; 
       unless (  $simple_feature_ids[$s_index] ) {  
         throw("something is wrong with your simple_feature-indexes - the simple feature with index $s_index does not exist\n") ; 
       }  
       unless (  $sf ) { 
         throw("Simple feature with dbID :  $simple_feature_ids[$s_index] not found in simple_feature table of database : " .  
         $sfa->db->dbname . "\@" . $sfa->db->host . "\n" ) ; 
       }

     }  
  }    


  for ( keys %recov_simple_feat_to_transcript ) {   
      my %species_processed ; 
      my %tmp ; 
      my @selected_id_for_species ; 
      @tmp{@{$recov_simple_feat_to_transcript{$_}}} = 1 ;  
      
      for my $stable_id ( keys %tmp ) {   
        my $species = substr($stable_id,0,4); 
        unless (exists $species_processed{$species} ) {   
          push @selected_id_for_species, $stable_id ;  
        }
        $species_processed{$species} =1 ; 
      } 
      #$recov_simple_feat_to_transcript{$_}=[keys %tmp] ; 
      $recov_simple_feat_to_transcript{$_}=\@selected_id_for_species ;  
      print "pushing $_ " . join ( "\t", @selected_id_for_species ) . "\n"; 
  }  

  my @stable_ids_to_process ; 
  for my $id ( keys %recov_simple_feat_to_transcript ) {      
    my @tmp = @{$recov_simple_feat_to_transcript{$id}};  
    for ( @tmp ) { 
       push @stable_ids_to_process, $_ ;  
    }
  } 
  if ( $self->lightweight_analysis ) {  
    print "\n\nWarning : Using lightweight analysis and only 1 protein per species per recovery-case\n\n" ; 
    $self->homolog_transcript_sid(\@stable_ids_to_process);
  }   

 
  if ( $self->debug ) {  
    print "\n\nProcessing these ids :\n\n" ;  
    for ( keys %recov_simple_feat_to_transcript ) {   
       print "$_  " ; 
       for ( @{$recov_simple_feat_to_transcript{$_}}  ) { 
         print "$_  " ; 
       } 
       print "\n" ; 
    }
  }
  return \%enst_to_simple_features ;
}  

sub set_homolog_transcript {  
   my ( $self,  $trans ) = @_ ; 
   if ( $trans) {  
      ${$self->{hom_trans_cache}}{$trans->stable_id}= $trans ; 
   }
   return $self->{hom_trans_cache};    
} 

sub get_homolog_transcript_by_stable_id {  
   my ( $self,  $trans_stable_id ) = @_ ;  
   return ${$self->{hom_trans_cache}}{$trans_stable_id}; 
} 



sub sf_to_gene {  
   my ( $self,  $hseq_name, $orth_trans ) = @_ ;  

   
   if ( $hseq_name && $orth_trans ) {     
       # avoid to store transcripts twice 
       if ( exists ${$self->{sf_to_hom_trans_cache}}{$hseq_name} ) {    
          
          # hash-slice check if we've already stored the hseq-to-transcript relation 
          #      $xx{$hseq} = [ tr1 tr2 tr3 ]   
          my %tmp  ; 
          my @transcripts_stored =  @{${$self->{sf_to_hom_trans_cache}{$hseq_name}}} ; 
          @tmp{@transcripts_stored} = 1;    
          # now we got all transcript refs as keys 
          
          if ( exists $tmp{$orth_trans} ) {  
            # transcript already stored for this supporting feature 
          } else { 
            
           push @{ ${ $self->{sf_to_hom_trans_cache}{$hseq_name}}}, $orth_trans ; 
          } 
       } else {   
          # we got no entry for this hseq_name so let's create one 
           push @{ ${ $self->{sf_to_hom_trans_cache}{$hseq_name}}}, $orth_trans ; 
       } 
   } elsif ( $hseq_name ) { 
     # only $hseq_name is given so we return the transcripts associated to this supporting feature   
     return \@{${$self->{sf_to_hom_trans_cache}{$hseq_name}}} ; 
   } 
} 


use vars '$AUTOLOAD';
sub AUTOLOAD {
 my ($self,$val) = @_;
 (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name 
 $self->{$routine_name}=$val if defined $val ; 
 return $self->{$routine_name} ;
}
sub DESTROY {} # required due to AUTOLOAD
1; 

