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

Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator - 

=head1 SYNOPSIS

my $orthologueanalysis = Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator->new(
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

package Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator; 

use strict;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator;  
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils; 
use Bio::EnsEMBL::Registry; 
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;  
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;  
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;  
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::Rule; 
use vars qw(@ISA); 

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args); 
  $self->read_and_check_config;  
  $self->verbose(1) ; 
  return $self   ;
}


sub get_initial_geneset {   

  my ( $self, $species_alias,$biotypes) = @_ ; 

  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species_alias,'core') ;  

  unless ( $dba )  {
    throw("could not get database adaptor for species : $species_alias ") ;  
  } 

  my $sa = $dba->get_SliceAdaptor() ;

  print "Fetching genes out of ". $dba->dbname . " @ " . $dba->host ." : " . $dba->port ." : " . $self->slice_name .  "\n" ;
  my $qy_slice      = $sa->fetch_by_name($self->slice_name) ;

  my @pt_genes ;
  for my $bt ( @$biotypes ) { 
    my $genes_on_slice =  $qy_slice->get_all_Genes_by_type( $bt );
    print scalar(@$genes_on_slice) . " genes of biotype $bt found on slice ( $species_alias ) \n";  
    $self->genes( $genes_on_slice ) ; 
   }
} 


sub upload_input_ids {  
   my ( $self, $input_ids  ) = @_ ; 
 
   my $submit_analysis = "Submit_" . $self->post_logic_name;      

   print "Submit-analysis for next analysis is : $submit_analysis\n" ; 
   print "Logic-name for next analysis is      : " . $self->post_logic_name . "\n" ; 

   # otherwise the Registry overrides the adaptor - silly ... 
   #print "db before clear : " . $self->db(); 
   #print "\n" ;  
   #Bio::EnsEMBL::Registry->clear();  
   #print "db after clear : " . $self->db(); 
   #print "\n" ;  
   $self->pipeline_adaptor($self->db) ;  
 
    my  $input_id_type = "file_" . $self->post_logic_name ;  
    my $if = Bio::EnsEMBL::Pipeline::Utils::InputIDFactory->new(  
                 -db => $self->pipeline_adaptor , 
                 -logic_name => $submit_analysis , 
                 -file => 1 , 
                 -input_id_type => $input_id_type, 
               );  
     unless ($if->get_analysis($submit_analysis , $input_id_type,  1)) { 
       throw( "Cant find analysis $submit_analysis in db ") ; 
     }
     # check first if input_id is not already stored ...
      my $a = 
       $self->pipeline_adaptor->get_AnalysisAdaptor->fetch_by_logic_name($submit_analysis) ; 
       my $ia_db = $self->pipeline_adaptor->get_StateInfoContainer->list_input_id_by_Analysis($a) ;    
       my @input_ids_not_stored ;  
     
       my %tmp ; 
       @tmp{@$ia_db} = 1; 
       for my $i ( @$input_ids ) { 
         push @input_ids_not_stored, $i unless (exists $tmp{$i}) ; 
       } 
      print scalar(@$input_ids) - scalar(@input_ids_not_stored) . " input ids already stored in db " . $self->pipeline_adaptor->dbname . "\n" ;  
      $if->input_ids(\@input_ids_not_stored) ;  
      $if->store_input_ids;  
      print scalar(@input_ids_not_stored) . " input-ids uploaded into " . $self->pipeline_adaptor->dbname . "\n" ; 
}


sub chunk_and_write_fasta_sequences {
    my ( $self, $tref, $base_dir , $file_prefix, $file_suffix,$chunk_size  ) = @_ ;

    $chunk_size = 5 unless $chunk_size ; 
    print scalar(@$tref) . " sequences to write \n" ;

    return if scalar(@$tref) == 0 ; 

    unless ($base_dir) {
      my $conf = $$EXONERATE_CONFIG_BY_LOGIC{$self->post_logic_name}{QUERYSEQS};  
      $conf ? $base_dir = $conf : throw ( "There's no output-dir configured in Exonerate2Genes.pm " . 
      "for analysis " . $self->post_logic_name . " check the config\n"); 
    }
  
    $base_dir = $base_dir . "/" unless $base_dir =~m/\/$/;
    #$base_dir = $base_dir . $self->post_logic_name . "/" ;

    `mkdir $base_dir ` unless ( -e $base_dir);  

    unless ( $file_prefix ) {  
      if ($self->slice_name){
        $file_prefix = $self->slice_name ; 
      }elsif( $self->input_id ) { 
        $file_prefix = $self->input_id ; 
      }else{
        throw("Error - can't create file name because of missing slice-name or input_id\n") ; 
      }
    }
    $file_suffix = ".fa" unless $file_suffix ;

    print "writing squences to $base_dir\n" ;

    my @filenames ;

    my $wtf = 0 ;
    my @ltr = @$tref ;
    my ($seq_file , $name ) ;
    my $fcnt = 0 ;

    for ( my $i=0 ; $i<scalar(@ltr) ; $i++ ) {
      my $seq = $ltr[$i] ;
      # create new file if #$chunk_size chunks are already written
      if ($i % $chunk_size  == 0 ){

          $seq_file->close() if ($wtf == 1) ;
          $fcnt++;
          my $file_name  = $file_prefix . "_" . $fcnt . $file_suffix ;
          # these filenames are stored as input_ids in refdb  
          push @filenames, $file_name ;
          my $tmp = $base_dir . $file_name ;

          $seq_file = Bio::SeqIO->new(
                                       #-file => ">$base_dir. $file_name" ,
                                       -file => ">$tmp" ,
                                       -format => 'fasta'
                                     );

          $wtf = 1 ;
      }
          $seq_file->write_seq($seq);
    }
    $seq_file->close();
    return \@filenames ;
}


sub genes {
  my ($self, $g) = @_ ;
  push @{$self->{_genes}}, @$g if $g ;
  return $self->{_genes} ;
} 


sub species_1 {
  my ($self, $g) = @_ ;
  $self->{_species_1} = $g if $g ;
  return $self->{_species_1} ;
} 

sub species_2 {
  my ($self, $g) = @_ ;
  $self->{_species_2} = $g if $g ; 
  return $self->{_species_2} ; 
}

sub post_logic_name {  
  my ($self) = @_ ;
  my $post_logic_name = $self->input_id ; 
  $post_logic_name =~s/(^.*?)(\:.*)/$1/ ;  
  return $post_logic_name ; 
}

sub slice_name  {   
  my ($self) = @_ ;  
  my $slice_name = $self->input_id ; 
   $slice_name  =~s/(^.*?\:)(.*)/$2/   ;
 return $slice_name ; 
}

sub run {
  my ($self) = @_; 
}

sub write_output{
  my ($self,$missing_seq) = @_;
}




sub read_and_check_config {
  (my $self) = @_ ;  

  # check if config file exists 

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
}




sub pipeline_adaptor{  
    my ($self, $db )= @_ ;    
 
       if ( $db ) { 
         my $pa = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                                                            -host   => $db->host, 
                                                            -dbname => $db->dbname, 
                                                            -user   => $db->username,
                                                            -pass   => $db->password ,
                                                            -port   => $db->port,
                                                            );
         $self->{_PA}=$pa ;  
       } 
    return $self->{_PA} ; 
}



1;
