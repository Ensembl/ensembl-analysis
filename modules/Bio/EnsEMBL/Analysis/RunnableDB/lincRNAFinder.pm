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

Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::lincRNAFinder;
use Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange); 
use Bio::EnsEMBL::Analysis::RunnableDB; 
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils ;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild Bio::EnsEMBL::Analysis::RunnableDB);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder
  Function  : instatiates a lincRNAFinder object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder 
  Exceptions: 
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args); 

  $self->read_and_check_config($LINCRNA_CONFIG_BY_LOGIC);    
  $self->OUTPUT_BIOTYPE($self->analysis->logic_name) if(!$self->OUTPUT_BIOTYPE); 

  return $self;
}



sub fetch_input{
  my ($self) = @_; 

  $self->query($self->fetch_sequence); 

  # get cdnas and convert them to single transcript genes  
  my $new_cdna = $self->get_genes_of_biotypes_by_db_hash_ref($self->NEW_SET_1_CDNA);   
  my @single_transcript_cdnas =  map { @{convert_to_single_transcript_gene($_)}  }  @$new_cdna ;

  # get protein_coding genes and convert them to single transcript genes  
  my $new_set_prot = $self->get_genes_of_biotypes_by_db_hash_ref($self->NEW_SET_2_PROT);
  my @single_trans_pc = map { @{convert_to_single_transcript_gene($_)}  }  @$new_set_prot;  

  # create runnable  
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder->new(
        -query => $self->query,
        -analysis => $self->analysis,
     );  

  # add hash-keys and hash-values directly to the $runnable hashref. quicker than using constructors...

  $runnable->set_1_cdna_genes(\@single_transcript_cdnas);
  $runnable->set_2_prot_genes(\@single_trans_pc);  
  $runnable->efg_simple_feature_genes($self->get_efg_simple_features); 
  
  $runnable->ignore_strand($self->CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND);
  $runnable->find_single_exon_candidates($self->FIND_SINGLE_EXON_LINCRNA_CANDIDATES);
  
  $runnable->check_cdna_overlap_with_both_K4_K36($self->CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36);
  $runnable->check_cdna_overlap_with_multiple_K36($self->CHECK_CDNA_OVERLAP_WITH_MULTI_K36);

  $runnable->maximum_translation_length_ratio($self->MAXIMUM_TRANSLATION_LENGTH_RATIO); 
  $runnable->max_translations_stored_per_gene($self->MAX_TRANSLATIONS_PER_GENE) ;
  
  $runnable->efg_clustering_with_cdna_analysis($self->create_analysis_object($self->DEBUG_LG_EFG_CLUSTERING_WITH_CDNA)) ; 
  $runnable->unclustered_efg_analysis( $self->create_analysis_object($self->DEBUG_LG_EFG_UNCLUSTERED));

  $self->runnable($runnable);  
};


sub update_efg_and_cdna_db { 
  my ($self) = @_;   

  # update efg features and change their analysis according to their clustering  
  print "Going to write simple features which have been converted to genes ...\n" ;  
  for my $rb ( @{ $self->runnable() } ) {  
  my @efg_genes  = @{ $rb->updated_efg_genes } ;    
    my $gfa = $self->get_dbadaptor($self->DEBUG_OUTPUT_DB)->get_GeneAdaptor() ; 
    for my $e ( @efg_genes ) { 
      $gfa->store($e) ; 
    }  
  } 
  print "Updated database " . $self->DEBUG_OUTPUT_DB . " with efg features 'genes'.\n" ;  

  for my $rb ( @{ $self->runnable() } ) {   
    print "Going to write cDNA genes for debugging ...\n";
    my $gfa = $self->get_dbadaptor($self->DEBUG_OUTPUT_DB)->get_GeneAdaptor() ; 
    for my $cdna (@{$rb->updated_cdnas}) {
      $gfa->store($cdna); 
    } 
  }
  print "Updated database " . $self->DEBUG_OUTPUT_DB . " with debugging cDNA 'genes'.\n" ;  
} 


sub write_output{
  my ($self) = @_; 

  if ( $self->WRITE_DEBUG_OUTPUT ) {
    $self->update_efg_and_cdna_db() ;
  }

  my $ga = $self->output_db->get_GeneAdaptor;

  logger_info("have ".@{$self->output}." genes to write"); 

  my $sucessful_count = 0 ; 

  GENE: foreach my $gene(@{$self->output}){  
  
   if ( !defined $gene->get_all_Transcripts ) {  
       throw (" gene does not have any transcripts ....\n" ) ; 
   }
 
    my @tr = @{ $gene->get_all_Transcripts };  
    my $max_ex = 0;  

    for ( @tr ) { 
       $_->status(undef);
       $_->analysis($self->analysis);  
    }    

    $gene->biotype($self->OUTPUT_BIOTYPE); 
    $gene->status(undef); 
    $gene->analysis($self->analysis);   

    eval{
      $ga->store($gene);
    }; 

    if($@){
      warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
    }else{
      $sucessful_count++;
      logger_info("STORED LINCRNA GENE ".$gene->dbID);
      print "STORED LINCRNA GENE ".$gene->dbID . " " . $gene->biotype . " " .  $self->output_db->dbname . " @ ".$self->output_db->host . " ".scalar(@tr) . " transcripts\n"  ; 
    }
  } 
     print $sucessful_count ." genes written to " . $self->output_db->dbname . " @ ".
     $self->output_db->host . "\n"  ;   

  if($sucessful_count != @{$self->output}){
    throw("Failed to write some genes");
  }
}


sub output_db{
  my ($self, $db) = @_; 

  if(!$self->{'output_db'}){
    my $db = $self->get_dbadaptor($self->OUTPUT_DB); 
    #$db->disconnect_when_inactive(1);
    $self->{output_db} = $db; 
  }
  return $self->{output_db};
}



sub get_efg_simple_features {  
  my ( $self ) = @_ ;  

  my $dba = $self->get_dbadaptor($self->EFG_FEATURE_DB); 
  # $dba->disconnect_when_inactive(0);    

  my $set_db = $dba->get_SimpleFeatureAdaptor();
  print "Fetching EFG domain data from " . $dba->dbname . "\@" . $dba->host . ".... \n";

  my @simple_features = @{ $set_db->fetch_all_by_Slice($self->query) };
  return  $self->convert_simple_features(\@simple_features); 
}




sub convert_simple_features {   
  my ( $self, $sf) = @_ ;    
  my @converted_sf;  
  my @simple_features= @$sf;  

  print scalar(@simple_features ) . " retrieved \n" ; 
  my %lg_names ;  
  my @filtered_sf ;
  for my $sf ( @simple_features ) { 
    for my $lg ( @{ $self->EFG_FEATURE_NAMES } ) { 
      if ( $sf->analysis ) {  
         if ( $sf->analysis->logic_name =~m/^$lg$/) {
           push @filtered_sf, $sf ; 
         } 
      } else {  
         throw ("sf with dbID " . $sf->dbID . " does not have a logic_name , this might be 'cause it has just been created ...! \n") ; 
      } 
    }
  } 
  @simple_features = @filtered_sf ;
  print " after filtering by logic names got " . scalar(@simple_features ) . " features left ...\n" ; 

  my $offset =  $self->EXTEND_EFG_FEATURES;
  print "Extending start/end of EFG domains by +-$offset bp.\n"; 

  for my $f ( @simple_features ) {

    my $f_start ;
    my $f_end ;

    $f_start = $f->start - $offset ;
    $f_end = $f->end+ $offset ;

    #print $f->strand . " : BEFORE S : " . $f->start . "  --->   $f_start \n" ; 
    #print $f->strand . " : BEFORE E : " . $f->end  . "  --->   $f_end  \n ";   
    #print $f->strand . " : BEFORE L : " . ($f->end - $f->start )   . "  ---> NOW :  " . ($f_end - $f_start) . " \n ";  

       my $ex = new Bio::EnsEMBL::Exon(
      -START     => $f_start,
      -END       => $f_end,
      -STRAND    => $f->strand,
      -SLICE     => $f->slice,
      -PHASE     => '0' ,
      -END_PHASE     => '0' ,
      -ANALYSIS  => $f->analysis, 
    );
    my $tran = new Bio::EnsEMBL::Transcript();
    $tran->add_Exon($ex) ; 
    $tran->analysis($f->analysis); 
    my $gene = new Bio::EnsEMBL::Gene() ;
    $gene->add_Transcript($tran) ;
    $gene->biotype("efg") ; # don't change this biotype 
    $gene->biotype("efg") ; # don't change this biotype 
    $gene->description($f->display_label);
    $gene->dbID($f->dbID) ; 
    $gene->analysis($f->analysis); 
    push @converted_sf, $gene ;
  }
  return \@converted_sf;
}


  
=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder
  Arg [2]   : hashref from config file
  Function  : call the superclass method to set all the varibles and carry
  out some sanity checking
  Returntype: N/A
  Exceptions: throws if certain key variables aren't set properly
  Example   : 

=cut

sub read_and_check_config{
  my ($self, $hash) = @_;
  $self->SUPER::read_and_check_config($hash); 

  #######
  #CHECKS
  #######
  foreach my $var(qw(NEW_SET_1_CDNA NEW_SET_2_PROT OUTPUT_DB OUTPUT_BIOTYPE EFG_FEATURE_NAMES EXTEND_EFG_FEATURES )){
    throw("RunnableDB::lincRNAFinder $var config variable is not defined") 
      unless($self->$var);
  } 

}


=head2 NEW_SET_1_CDNA 

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut


sub NEW_SET_1_CDNA {
  my ($self, $arg) = @_;
  if($arg){
    $self->{'NEW_SET_1_CDNA'} = $arg;
  }
  return $self->{'NEW_SET_1_CDNA'};
} 

sub NEW_SET_2_PROT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'NEW_SET_2_PROT'} = $arg;
  }
  return $self->{'NEW_SET_2_PROT'};
}


sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'OUTPUT_DB'} = $arg;
  }
  return $self->{'OUTPUT_DB'};
}  


sub EFG_FEATURE_DB{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'EFG_FEATURE_DB'} = $arg;
  }
  return $self->{'EFG_FEATURE_DB'};
}  


sub EFG_FEATURE_NAMES {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'EFG_FEATURE_NAMES'} = $arg;
  }
  return $self->{'EFG_FEATURE_NAMES'};
} 

sub EXTEND_EFG_FEATURES{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'EXTEND_EFG_FEATURES'} = $arg;
  }
  return $self->{EXTEND_EFG_FEATURES};
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'OUTPUT_BIOTYPE'} = $arg;
  }
  return $self->{'OUTPUT_BIOTYPE'};
} 

sub DEBUG_OUTPUT_DB{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'DEBUG_OUTPUT_DB'} = $arg;
  }
  return $self->{'DEBUG_OUTPUT_DB'};
} 

sub DEBUG_LG_EFG_CLUSTERING_WITH_CDNA{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'DEBUG_LG_EFG_CLUSTERING_WITH_CDNA'} = $arg;
  }
  return $self->{'DEBUG_LG_EFG_CLUSTERING_WITH_CDNA'};
}    

sub MAX_TRANSLATIONS_PER_GENE { 
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAX_TRANSLATIONS_PER_GENE'} = $arg;
  }
  return $self->{'MAX_TRANSLATIONS_PER_GENE'};
}

sub MAXIMUM_TRANSLATION_LENGTH_RATIO {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAXIMUM_TRANSLATION_LENGTH_RATIO'} = $arg;
  }
  return $self->{'MAXIMUM_TRANSLATION_LENGTH_RATIO'};
}

sub DEBUG_LG_EFG_UNCLUSTERED{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'DEBUG_LG_EFG_UNCLUSTERED'} = $arg;
  }
  return $self->{'DEBUG_LG_EFG_UNCLUSTERED'};
}   

sub WRITE_DEBUG_OUTPUT {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'WRITE_DEBUG_OUTPUT'} = $arg;
  }
  return $self->{'WRITE_DEBUG_OUTPUT'};
}  

sub CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND'} = $arg;
  }
  return $self->{'CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND'};
}


sub CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36 {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36'} = $arg;
  }
  return $self->{'CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36'};
}


sub CHECK_CDNA_OVERLAP_WITH_MULTI_K36 {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'CHECK_CDNA_OVERLAP_WITH_MULTI_K36'} = $arg;
  }
  return $self->{'CHECK_CDNA_OVERLAP_WITH_MULTI_K36'};
}



sub FIND_SINGLE_EXON_LINCRNA_CANDIDATES {
  my ($self, $arg) = @_;
  if (defined $arg){
    $self->{'FIND_SINGLE_EXON_LINCRNA_CANDIDATES'} = $arg;
  }
  return $self->{'FIND_SINGLE_EXON_LINCRNA_CANDIDATES'};
}


sub create_analysis_object { 
  my ($self,$logic_name ) = @_ ;   
   
  my $efg_out = $self->get_dbadaptor($self->DEBUG_OUTPUT_DB);  
  my $aa = $efg_out->get_AnalysisAdaptor() ;  
  my $analysis =  $aa->fetch_by_logic_name($logic_name) ; 
  if ( $analysis ) { 
     return $analysis ; 
  } 
  # need to create analysis object first   
  $analysis=Bio::EnsEMBL::Analysis->new( 
                                       -logic_name => $logic_name , 
                                       );  
  return $analysis ; 
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
