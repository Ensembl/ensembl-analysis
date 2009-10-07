package Bio::EnsEMBL::Analysis::RunnableDB::MakeSimilarityInputIDs;

use strict;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise qw(GENEWISE_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (parse_config);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($output_logicname, $protein_count, $bmg_logicname) = rearrange
    (['OUTPUT_LOGICNAME', 'PROTEIN_COUNT', 'BMG_LOGICNAME'], @args); 

  # output logicname
  # protein count
  # bmg logicname

   my $create_analysis =  ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{creation_analysis};
   my $s_regex = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{submission_logic_name};
   my $bmg_regex = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{bmg_logic_name};
   my $def_protein_count = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{protein_count} ; 

   if ( $self->analysis->logic_name =~m/$create_analysis/i ) {     

     my $submission_logic_name = $self->analysis->logic_name;   
     $submission_logic_name =~m/$s_regex/i; 
     $submission_logic_name = $1 ; 
     $submission_logic_name="Submit_".$submission_logic_name ; 

     my $bmg  = $self->analysis->logic_name;   
     $bmg =~m/$bmg_regex/i; 
     $bmg = $1 ;  

     print "BlastMiniGenewise - logic name : $bmg\n" ; 
     print "SUBMISSION ANAL   - logic name : $submission_logic_name\n" ;  

     $self->bmg_logicname($bmg);
     $self->output_logicname($submission_logic_name);  
 
     if ($def_protein_count) { 
       $self->protein_count($def_protein_count); 
     }
  } 

 
  ### Defaults are over-ridden by parameters given in analysis table...
  my $ph = $self->parameters_hash;
  $self->protein_count($ph->{-protein_count}) if $ph->{-protein_count};
  $self->output_logicname($ph->{-output_logicname}) if $ph->{-output_logicname} ; ;
  $self->bmg_logicname($ph->{-bmg_logicname}) if $ph->{-bmg_logicname};

  ### ...which are over-ridden by constructor arguments. 
  $self->protein_count($protein_count);
  $self->output_logicname($output_logicname);
  $self->bmg_logicname($bmg_logicname);

  throw("Need an output logicname ".$self->output_logicname." and a bmg logicname ".
        $self->bmg_logicname." defined") if(!$self->output_logicname ||
                                            !$self->bmg_logicname);

  parse_config($self, $GENEWISE_CONFIG_BY_LOGIC, $self->bmg_logicname);
  return $self;
}


sub fetch_input{
  my ($self) = @_;
  print "\n\n***Fetching sequence from ".$self->db->dbname."***\n\n";
  $self->query($self->fetch_sequence($self->input_id, $self->db, $self->REPEATMASKING)); 

  
  my $dba_pip = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( 
                                                      -host => $self->db->host , 
                                                      -dbname => $self->db->dbname , 
                                                      -user => $self->db->username , 
                                                      -pass => $self->db->password , 
                                                      -port => $self->db->port ,  
                                                       );      
  $self->pipeline_adaptor($dba_pip); 
  my $output_analysis = $dba_pip->get_AnalysisAdaptor->fetch_by_logic_name($self->output_logicname); 

  if(!$output_analysis) { 
    throw("Failed to find an analysis with the logic_name ". $self->output_logicname." Cannot continue") ; 
  }
  $self->output_analysis($output_analysis); 
}

sub run{
  my ($self) = @_;
  my %kill_list =  %{$self->kill_list} if($self->USE_KILL_LIST);
  my @mask_exons;
  my @iids;
  # remove masked and killed hits as will be done in the build itself
  foreach my $type ( @{$self->BIOTYPES_TO_MASK} ) {  
    print "\nmasking Gene-type : $type\n\n" ; 
    foreach my $mask_genes ( @{ $self->gene_slice->get_all_Genes_by_type($type) } ) {
      foreach my $mask_exon ( @{ $mask_genes->get_all_Exons } ) {
        if ( $mask_exon->seqname eq $self->gene_slice->id ) {
          push @mask_exons, $mask_exon;
        }
      }
    }
  }
  # make the mask list non-redundant. Much faster when checking against features
  my @mask_regions;
  foreach my $mask_exon ( sort { $a->start <=> $b->start } @mask_exons ) {
    if ( @mask_regions and $mask_regions[-1]->{'end'} > $mask_exon->start ) {
      if ( $mask_exon->end > $mask_regions[-1]->{'end'} ) {
        $mask_regions[-1]->{'end'} = $mask_exon->end;
      }
    } else {
      push @mask_regions, { start => $mask_exon->start, end => $mask_exon->end }

    }
  }  


  my $num_seeds = 0;  
  foreach my $logicname(@{$self->PAF_LOGICNAMES}) {
    my %features;
    print "FETCHING FEATURES FOR :".$logicname."\n";
    my @features = @{$self->paf_slice->get_all_ProteinAlignFeatures($logicname)};
    print "HAVE ".@features." protein-align-features\n";
      FEATURE:foreach my $f(@features){
        next FEATURE if($self->PAF_MIN_SCORE_THRESHOLD && $f->score < $self->PAF_MIN_SCORE_THRESHOLD);
        next FEATURE if($self->PAF_UPPER_SCORE_THRESHOLD && $f->score > $self->PAF_UPPER_SCORE_THRESHOLD);
        push(@{$features{$f->hseqname}}, $f);
    }
    
    my @ids_to_ignore;
  SEQID: foreach my $sid ( keys %features ) {
      my $ex_idx = 0;
      my $count  = 0;
      
      #print STDERR "Looking at $sid\n";
    FEAT: foreach my $f ( sort { $a->start <=> $b->start } @{ $features{$sid} } ) {
        
        #printf STDERR "Feature: %d %d\n", $f->start, $f->end;
        for ( ; $ex_idx < @mask_regions ; ) {
          my $mask_exon = $mask_regions[$ex_idx];
          
          #printf STDERR " Mask exon %d %d\n", $mask_exon->{'start'}, $mask_exon->{'end'};
          if ( $mask_exon->{'start'} > $f->end ) {
            print "no exons will overlap this feature \n" ;  
            # no exons will overlap this feature
            next FEAT;
          } elsif ( $mask_exon->{'end'} >= $f->start ) {
            
            # overlap
            push @ids_to_ignore, $f->hseqname;
            print "Ignoring : " . $f->hseqname . "\n" ;  
            next SEQID;
          } else {
            $ex_idx++;
          }
        }
      }
    } 

    print "Ignoring ".@ids_to_ignore." features\n";
    foreach my $dud_id ( @ids_to_ignore, keys %kill_list ) {
      if ( exists $features{$dud_id} ) {
        delete $features{$dud_id};
      }
    }
    
    $num_seeds += scalar( keys %features );
  }
  # rule of thumb; split data so that each job constitutes one piece of
  # genomic DNA against ~20 proteins.
  #
  return () if($num_seeds == 0);

  my $num_chunks = int( $num_seeds / $self->protein_count ) 
    + 1;
  for ( my $x = 1 ; $x <= $num_chunks ; $x++ ) {
    
    #generate input id : $chr_name.1-$chr_length:$num_chunks:$x
    my $new_iid = $self->query->name . ":$num_chunks:$x";
    push @iids, $new_iid;
  }
  print "HAVE ".@iids." to write to the ref database\n"; 
  for ( @iids ) {  
      print "$_\n" ; 
  } 
  $self->output(\@iids);
}

sub write_output{
  my ($self) = @_;  

  my $sic = $self->pipeline_adaptor->get_StateInfoContainer;  

#  my $dba_pip = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( 
#                                                      -host => $self->db->host , 
#                                                      -dbname => $self->db->dbname , 
#                                                      -user => $self->db->username , 
#                                                      -pass => $self->db->password , 
#                                                      -port => $self->db->port ,  
#                                                       );    
#  my $sic = $dba_pip->get_StateInfoContainer; 

  print "output analysis is : " . $self->output_analysis."\n" ;
  foreach my $iid(@{$self->output}){ 
    print "try to store input_id : $iid\n" ; 
    eval{
      $sic->store_input_id_analysis($iid, 
                                    $self->output_analysis, 
                                    '');
    };
    throw("Failed to store ".$iid." $@") if($@);
    logger_info("Stored ".$iid);
  }
}


sub kill_list{
  my ($self, $arg) = @_;

  if($arg){
    $self->{kill_list} = $arg;
  }
  if(!$self->{kill_list}){
    my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => 'protein');
    $self->{kill_list} = $kill_list_object->get_kill_list;
  }
  return $self->{kill_list};
}


sub output_analysis{
  my ($self, $value) = @_;
  if($value && 
     $value->isa('Bio::EnsEMBL::Pipeline::Analysis')){
    $self->{output_analysis} = $value;
  }
  return $self->{'output_analysis'};
} 

sub pipeline_adaptor{
  my ($self, $value) = @_;
  if($value && 
     $value->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
    $self->{pipeline_adaptor} = $value;
  }
  return $self->{'pipeline_adaptor'};
}

sub protein_count{
  my ($self, $value) = @_;
  if($value){
    $self->{'protein_count'} = $value;
  }
  return $self->{'protein_count'};
}

sub output_logicname{
  my ($self, $value) = @_;
  if($value){
    $self->{'output_logicname'} = $value;
  }
  return $self->{'output_logicname'};
}

sub bmg_logicname{
  my ($self, $value) = @_;
  if($value){
    $self->{'bmg_logicname'} = $value;
  }
  return $self->{'bmg_logicname'};
}

sub paf_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->{paf_slice} = $slice;
  }
  if(!$self->{paf_slice}){
    my $slice = $self->fetch_sequence($self->input_id, $self->paf_source_db, $self->REPEATMASKING);
    $self->{paf_slice} = $slice;
  }
  return $self->{paf_slice};
}

sub gene_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->{gene_slice} = $slice;
  }
  if(!$self->{gene_slice}){
    my $slice = $self->fetch_sequence($self->input_id, $self->gene_source_db, $self->REPEATMASKING);
    $self->{gene_slice} = $slice;
  }
  return $self->{gene_slice};
}


sub paf_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->{paf_source_db} = $db;
  }
  if(!$self->{paf_source_db}){
    my $db = $self->get_dbadaptor($self->PAF_SOURCE_DB);
    $self->{paf_source_db} = $db;
  }
  return $self->{paf_source_db};
}

sub gene_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->{gene_source_db} = $db;
  }
  if(!$self->{gene_source_db}){
    my $db = $self->get_dbadaptor($self->GENE_SOURCE_DB);
    $self->{gene_source_db} = $db;
  }
  return $self->{gene_source_db};
}


=head2 CONFIG_ACCESSOR_METHODS

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise


sub PAF_LOGICNAMES{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_LOGICNAMES} = $arg;
  }
  return $self->{PAF_LOGICNAMES}
}

sub PAF_MIN_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_MIN_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_MIN_SCORE_THRESHOLD}
}

sub PAF_UPPER_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_UPPER_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_UPPER_SCORE_THRESHOLD}
}



sub PAF_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_SOURCE_DB} = $arg;
  }
  return $self->{PAF_SOURCE_DB}
}

sub GENE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_SOURCE_DB} = $arg;
  }
  return $self->{GENE_SOURCE_DB}
}




sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_DB} = $arg;
  }
  return $self->{OUTPUT_DB}
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_BIOTYPE} = $arg;
  }
  return $self->{OUTPUT_BIOTYPE}
}

sub GENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENEWISE_PARAMETERS} = $arg;
  }
  return $self->{GENEWISE_PARAMETERS}
}

sub MINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MINIGENEWISE_PARAMETERS}
}

sub MULTIMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MULTIMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MULTIMINIGENEWISE_PARAMETERS}
}

sub BLASTMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BLASTMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{BLASTMINIGENEWISE_PARAMETERS}
}



sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_PARAMETERS} = $arg;
  }
  return $self->{FILTER_PARAMETERS}
}



sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_OBJECT} = $arg;
  }
  return $self->{FILTER_OBJECT}
}


sub BIOTYPES_TO_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BIOTYPES_TO_MASK} = $arg;
  }
  return $self->{BIOTYPES_TO_MASK}
}


sub EXON_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXON_BASED_MASKING} = $arg;
  }
  return $self->{EXON_BASED_MASKING}
}


sub GENE_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_BASED_MASKING} = $arg;
  }
  return $self->{GENE_BASED_MASKING}
}


sub POST_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{POST_GENEWISE_MASK} = $arg;
  }
  return $self->{POST_GENEWISE_MASK}
}

sub PRE_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PRE_GENEWISE_MASK} = $arg;
  }
  return $self->{PRE_GENEWISE_MASK}
}

sub REPEATMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REPEATMASKING} = $arg;
  }
  return $self->{REPEATMASKING}
}

sub SEQFETCHER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_OBJECT} = $arg;
  }
  return $self->{SEQFETCHER_OBJECT}
}

sub SEQFETCHER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_PARAMS} = $arg;
  }
  return $self->{SEQFETCHER_PARAMS}
}

sub USE_KILL_LIST{
  my ($self, $arg) = @_;
  if($arg){
    $self->{USE_KILL_LIST} = $arg;
  }
  return $self->{USE_KILL_LIST}
}


sub LIMIT_TO_FEATURE_RANGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{LIMIT_TO_FEATURE_RANGE} = $arg;
  }
  return $self->{LIMIT_TO_FEATURE_RANGE}
}


sub FEATURE_RANGE_PADDING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FEATURE_RANGE_PADDING} = $arg;
  }
  return $self->{FEATURE_RANGE_PADDING}
}

sub WRITE_REJECTED{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{WRITE_REJECTED} = $arg;
  }
  return $self->{WRITE_REJECTED};
}

sub REJECTED_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REJECTED_BIOTYPE} = $arg;
  }
  return $self->{REJECTED_BIOTYPE};
}

sub SOFTMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SOFTMASKING} = $arg;
  }
  return $self->{SOFTMASKING}
}

sub EXONERATE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXONERATE_PARAMETERS} = $arg;
  }
  return $self->{EXONERATE_PARAMETERS}
}

sub MAKE_SIMGW_INPUT_ID_PARMAMS { 
  my ($self, $arg) = @_;
  if($arg){
    $self->{MAKE_SIMGW_INPUT_ID_PARMAMS} = $arg;
  }
  return $self->{MAKE_SIMGW_INPUT_ID_PARMAMS}
}   
1;
