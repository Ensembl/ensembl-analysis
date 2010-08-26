
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->fetch_input();
$runnableDB->run();
$runnableDB->write_output(); #writes to DB

=head1 DESCRIPTION

Extends Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature to allow
use of compressed dna align features, useful when aligning millions of short 
Solexa reads

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaCloud;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa; 
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Analysis::Config::ExonerateSolexaCloudConfig;                 # for S3 details CHUNK seq + GENOMIC seq
use Bio::EnsEMBL::Analysis::Config::General qw(ANALYSIS_WORK_DIR);

# this module also uses Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexa 
#                       Bio::EnsEMBL::Analysis::Config::ExonerateAlignFeature       ( genomic query sequence + chunks )
#                         


use vars qw(@ISA);
@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa );


sub new {
  my ( $class, @args ) = @_;  


  my $self = $class->SUPER::new(@args);  

  $self->S3_SEQUENCE_DATA($S3_SEQUENCE_DATA);   

  $self->S3_GENOMIC_FASTA_SEQUENCE($S3_GENOMIC_FASTA_SEQUENCE);  

  $self->db->disconnect_when_inactive(1);    
  # input_id for S3 :     chunk_1453.fa::OUTPUT_DB::LANE_XYZ

  # skeletal@chunk_1453.fa::OUTPUT_DB::LANE_XYZ
  # chunk234.fa: 'file' with chunked sequence in S3 
  # LANE_XYZ   : key in CloudConfig.pm   
  # OUTPUT_DB  : key in Databases.pm  

  my $input_id = $self->input_id; 
  my ($base_batch, $rest ) = split "@",$self->input_id;  
  $self->input_id($rest) ;  
  my ($chunk,$output_db_key,$s3_sequence_data_key) = split "::",$rest; 
  $self->chunk($chunk); 
  $self->output_db_key($output_db_key); 
  $self->s3_sequence_data_key($s3_sequence_data_key);  
  $self->base_batch($base_batch);  

  # I will batch analyses togther by basebatch name, ie tissue or so to have an easer config handling later
  # set output db for write_output()     
  return $self;
}





sub get_file_from_s3 {  
  my ($self,$file_name,$bucket_name) = @_;

  # check if this file is the genomic file - we want to treat this file different  
  my $analysis_work_dir = $ANALYSIS_WORK_DIR ?  $ANALYSIS_WORK_DIR : "/tmp";   
  my $out_file = $analysis_work_dir ."/".$file_name;   

  my $is_genomic_file ;  
  if ( $self->genomic_file_name eq $file_name ) {     
    $is_genomic_file = 1; 
    if ( $out_file =~m/\.gz$/ ) {
      (my $tmp_out_file_name = $out_file) =~s/\.gz$//g; 
      if ( -e $tmp_out_file_name ) {    
          # this could cause trouble depending how the un-zipping goes
          return $tmp_out_file_name; 
      }  
    } 
  }


  my $no_secret_key ; 
  my $no_access_key ;  

  if ( !defined $self->aws_secret_key_id ||  length($self->aws_secret_key_id) == 0 ){ 
    $no_secret_key = 1; 
  }  
  if (!defined $self->aws_access_key_id ||  length($self->aws_access_key_id) == 0 ){ 
    $no_access_key= 1; 
  }   



  if ( $no_access_key == 1 && $no_secret_key == 1 ) {  
     # data seems to be public;  
     my $url = $bucket_name . ".S3.amazonaws.com/". $file_name;  
     if (! defined $is_genomic_file ) {  
        $out_file = $out_file.".".$$;
     } 
     my $command = "wget $url -O $out_file"; 
     if ( -e $out_file ) {  
       throw("output file $out_file exists already - can't download file aos other file would be overw-ritten \n");
     } 
     system($command);  
     if ( (defined $self->gzip_compression && $self->gzip_compression==1)||  $out_file =~m/\.gz$/) {  
        # file is compressed; let's un-compress it  
        my $cmd = "gunzip -fc $out_file > $out_file.part.tmp "; 
        system($cmd);    
        system("unlink $out_file");
        my $of = $out_file; 
        $out_file =~s/\.gz$//; 
        $cmd="mv $of.part.tmp $out_file" ;  
        system($cmd);  
        if (!defined $is_genomic_file ) { 
           $self->file_to_delete($out_file); 
        } 
     }   
  }else {   
    # data no public - need S3 module to fetch data  
    throw("fetching data from S3 with password not implemented now\n"); 
  } 
  return $out_file; 
} 


sub fetch_input {
  my ($self) = @_;

  # get chunk1.fa.gz from S3   - security credentials in CloudConfig if needed
  my $path_to_chunk_file = $self->get_file_from_s3($self->chunk,$self->chunk_bucket_name); 
  $self->QUERYSEQS($path_to_chunk_file);  

  # get softmasked genome sequence out of S3 - security credentials in CloudConfig if needed
  my $path_to_genomic_file = $self->get_file_from_s3($self->genomic_file_name, $self->genomic_bucket_name);  
  $self->GENOMICSEQS($path_to_genomic_file); 
  $self->SUPER::fetch_input(); 
}


sub write_output {  
  my ($self) = @_; 

  my $cmd = "unlink ".$self->file_to_delete();  
  system($cmd); 
  $self->SUPER::write_output(); 
} 


# accessor methods for CloudConfig file 

sub chunk_bucket_name {  
  my ($self) = shift; 
  return ${$self->S3_SEQUENCE_DATA}{$self->s3_sequence_data_key}{"bucket_name"}; 
}  

sub aws_access_key_id{  
  my ($self) = shift; 
  return ${$self->S3_SEQUENCE_DATA}{$self->s3_sequence_data_key}{"aws_access_key_id"}; 
}  

sub aws_secret_key_id{  
  my ($self) = shift; 
  return ${$self->S3_SEQUENCE_DATA}{$self->s3_sequence_data_key}{"aws_secret_key_id"}; 
}  

sub file_regex{  
  my ($self) = shift; 
  return ${$self->S3_SEQUENCE_DATA}{$self->s3_sequence_data_key}{"file_regex"}; 
}  

sub gzip_compression {  
  my ($self) = shift; 
  return ${$self->S3_SEQUENCE_DATA}{$self->s3_sequence_data_key}{"gzip_compression"}; 
}  


# accessor for S3_GENOMIC_FASTA_SEQUENCE
sub genomic_aws_access_key_id {  
  my ($self) = shift; 
  return ${$self->S3_GENOMIC_FASTA_SEQUENCE}{"aws_access_key_id"} ;
}

sub genomic_aws_secret_key_id{  
  my ($self) = shift; 
  return ${$self->S3_GENOMIC_FASTA_SEQUENCE}{"aws_secret_key_id"} ;
} 

sub genomic_bucket_name {  
  my ($self) = shift; 
  return ${$self->S3_GENOMIC_FASTA_SEQUENCE}{"bucket_name"} ;
} 

sub genomic_file_name {  
  my ($self) = shift; 
  return ${$self->S3_GENOMIC_FASTA_SEQUENCE}{"file_name"} ;
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
