package Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder;

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


@ISA = qw (
           Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
           Bio::EnsEMBL::Analysis::RunnableDB
           );



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
  $self->{_all_gene_sets}={};  
  return $self;
}



sub fetch_input{
  my ($self) = @_;
  #fetch sequence/slice 
  $self->query($self->fetch_sequence); 

  #fetch genes
  $self->get_gene_sets; 


  my @efg_sf = @{$self->get_efg_simple_features()} ;  

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder->new(
        -query => $self->query,
        -analysis => $self->analysis,
     );  

  $runnable->set_1_cdna_genes($self->all_gene_sets("SET_1_CDNA"));
  $runnable->set_2_prot_genes($self->all_gene_sets("SET_2_PROT"));  
  $runnable->efg_simple_feature_genes( \@efg_sf ); 
  $runnable->efg_clustering_with_cdna_analysis($self->create_analysis_object($self->DEBUG_LG_EFG_CLUSTERING_WITH_CDNA)); 
  $runnable->unclustered_efg_analysis( $self->create_analysis_object($self->DEBUG_LG_EFG_UNCLUSTERED));  
  $runnable->maxium_translation_length_ratio($self->MAXIMUM_TRANSLATION_LENGTH_RATIO); 
  $self->runnable($runnable); 
};


sub run {
  my ($self) = @_;  

  $self->SUPER::run();
  if ( $self->DEBUG_WRITE_CLUSTERED_GENES ) {
    $self->update_efg_and_cdna_db() ;  
  }
  $self->write_output();  
}


sub update_efg_and_cdna_db { 
  my ($self) = @_;   

  # update efg features and change their analysis according to their clustering  
    print "Writing simple features which have been converted to genes ...\n" ;  
    for my $rb ( @{ $self->runnable() } ) {  
      my @efg_genes  = @{ $rb->updated_efg_genes } ;    
      my $gfa = $self->get_dbadaptor($self->DEBUG_OUTPUT_DB)->get_GeneAdaptor() ; 
      for my $e ( @efg_genes ) {  
        $gfa->store($e) ; 
      }  
    } 
    print "updated database " . $self->DEBUG_OUTPUT_DB . " - efg features, converted into genes, have been written...\n" ;  

    for my $rb ( @{ $self->runnable() } ) {   
      my $gfa = $self->get_dbadaptor($self->DEBUG_OUTPUT_DB)->get_GeneAdaptor() ; 
      for my $cdna (@{$rb->updated_cdnas}) {
        $gfa->store($cdna); 
      } 
    }
    my %sets_to_cluster = %{$self->CLUSTERING_INPUT_GENES};    
} 


sub write_output{
  my ($self) = @_; 

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
      print "STORED LINCRNA GENE ".$gene->dbID . " " . $gene->biotype . " " .  $self->output_db->dbname . " @ ".$self->output_db->host . "\n"  ; 
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
  $dba->disconnect_when_inactive(1);    

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
         if ( $sf->analysis->logic_name =~m/$lg/) {
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


  







sub get_gene_sets {
  my ($self) = @_;
  my @genes; 

   my %sets_to_cluster = %{$self->CLUSTERING_INPUT_GENES};  

    if ( keys %sets_to_cluster != 2 ) { 
        throw ("you should only have 2 sets to cluster against - you can't cluster against more sets \n" ); 
    } 
 
  # check if hash-key name is correct : 
  unless ( exists $sets_to_cluster{"SET_1_CDNA"} && exists $sets_to_cluster{"SET_2_PROT"}) {  
        throw( " configuration error - I expect to get 2 sets of genes with names \"SET_1_CDNA\" \n". 
               " and \"SET_2_PROT\" - check your config - you can't change these names!!!! \n" ) ; 
  }  
 
  foreach my $set ( keys %sets_to_cluster) { 
    my %this_set = %{$sets_to_cluster{$set}};      
    my @genes_in_set; 
    foreach my $database_db_name ( keys ( %this_set)) {   
       print "fetching $database_db_name\n" ; 
       my $set_db = $self->get_dbadaptor($database_db_name);
       #$set_db->disconnect_when_inactive(1);  
       my $slice = $self->fetch_sequence($self->input_id, $set_db);  
       my @biotypes = @{$this_set{$database_db_name}}; 
       for my $biotype  ( @biotypes ) { 
          my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
          if ( @$genes == 0 ) {  
            warning("No genes of biotype $biotype found in $set_db\n"); 
          } 
          print "Retrieved ".@$genes." of type ".$biotype."\n";
          push @genes_in_set, @$genes; 
      }  
    } 
    my @single_transcript_genes = @{$self->create_single_transcript_genes(\@genes_in_set)};
   $self->all_gene_sets($set,\@single_transcript_genes); 
  }    

  for ( keys  %{$self->all_gene_sets} ) { 
    print "$_  : "  . scalar (@{${$self->all_gene_sets}{$_}}) . " genes found\n";  
  }    
}


#
#  method to store all retrieved gene sets 
# 

sub all_gene_sets { 
  my ( $self, $set_name, $genes ) =  @_ ; 

  if ( $set_name && $genes ) {     
     ${$self->{_all_gene_sets}}{$set_name} = $genes ;  
     # return array with genes of set name 
     return ${$self->{_all_gene_sets}}{$set_name} ; 
  }elsif ( $set_name ) {   
     # return array with genes of set name 
     return ${$self->{_all_gene_sets}}{$set_name} ; 
  }  
  # return hash 
  return $self->{_all_gene_sets}; 
} 


sub create_single_transcript_genes{
  my ($self, $genes) = @_; 

  print "Creating single_transcript genes ....\n" ;
  my @single_transcript_genes; 

 GENE:foreach my $gene(@$genes){ 
    my @tr = @{$gene->get_all_Transcripts}; 
    if ( @tr == 1 ) {  
      push @single_transcript_genes, $gene ; 
    } else { 
      foreach my $transcript(@tr ) { 
        my $ng = Bio::EnsEMBL::Gene->new();  
        $ng->add_Transcript($transcript); 
        push @single_transcript_genes, $ng ; 
      }
    }
  }
  return \@single_transcript_genes; 
}

#CONFIG METHODS


=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : hashref from config file
  Function  : call the superclass method to set all the varibles and carry
  out some sanity checking
  Returntype: N/A
  Exceptions: throws if certain variables arent set properlu
  Example   : 

=cut

sub read_and_check_config{
  my ($self, $hash) = @_;
  $self->SUPER::read_and_check_config($hash); 

  #######
  #CHECKS
  #######
  foreach my $var(qw(CLUSTERING_INPUT_GENES OUTPUT_DB OUTPUT_BIOTYPE EFG_FEATURE_NAMES EXTEND_EFG_FEATURES )){
    throw("RunnableDB::lincRNAFinder $var config variable is not defined") 
      unless($self->$var);
  }
  my @keys = keys(%{$self->CLUSTERING_INPUT_GENES});
  throw("RunnableDB::lincRNAFinder CLUSTERING_INPUT_GENES has needs to contain values") if(!@keys);
  my %unique;
  foreach my $key(@keys){
    my $hashref_of_sets  = $self->CLUSTERING_INPUT_GENES->{$key};
    foreach my $set (%$hashref_of_sets ){
      if(!$unique{$set}){
        $unique{$set} = $key;
      }else{
          warning($set ." appears twice in your listing, make sure this ".
                  "isn't for the same database otherwise it will cause issue");
      }
    }
  }
  $self->OUTPUT_BIOTYPE($self->analysis->logic_name) if(!$self->OUTPUT_BIOTYPE);
}


=head2 CONFIG_ACCESSOR_METHODS

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut


sub CLUSTERING_INPUT_GENES {
  my ($self, $arg) = @_;
  if($arg){
    $self->{'CLUSTERING_INPUT_GENES'} = $arg;
  }
  return $self->{'CLUSTERING_INPUT_GENES'};
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

sub DEBUG_WRITE_CLUSTERED_GENES {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'DEBUG_WRITE_CLUSTERED_GENES'} = $arg;
  }
  return $self->{'DEBUG_WRITE_CLUSTERED_GENES'};
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
