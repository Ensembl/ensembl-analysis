
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->fetch_input();
arunnableDB->run();
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
use Bio::SeqIO;

# this module also uses Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexa 
#                       Bio::EnsEMBL::Analysis::Config::ExonerateAlignFeature       ( genomic query sequence + chunks )
#                         


use vars qw(@ISA);
@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa );


sub new {
  my ( $class, @args ) = @_;  

  push @args, ("-no_config_exception" , 1 );  

  my $self = $class->SUPER::new(@args);   
  $self->original_input_id($self->input_id); 


  my @ids = split "=",$self->input_id;   
  
  if (scalar(@ids) !=2 ){ 
     throw("Input_id format is wrong - it should be : BASE_BATCH\=chunk1.fa.gz::OUT_DB::1-100\n");
  }else {  
     # THIS runnable takes input_id format BATCH=chunk1.fa.gz::OUTPUT_DB::1-100 
     # other runnable ( ExonerateAlignFeature take   chunk1.fa.gz::OUTPUT_DB 
     # need to re-set input_id  

      my ( $chunk,$out_db,$range ) = split /::/,$ids[1]; 
      $self->input_id($chunk."::".$out_db) ;  
  }  


  $self->{_files_to_delete} = []; 
  $self->S3_SEQUENCE_DATA($S3_SEQUENCE_DATA);   

  $self->S3_GENOMIC_FASTA_SEQUENCE($S3_GENOMIC_FASTA_SEQUENCE);  
  $self->ANALYSIS_BASE_BATCH_CONFIG($ANALYSIS_BASE_BATCH_CONFIG); 

  $self->db->disconnect_when_inactive(1);     

  # input_id for S3 :    BLOOD = chunk_1453.fa :: OUTPUT_DB :: 1-1000
  # input_id for S3 :    BLOOD = chunk_1453.fa :: OUTPUT_DB :: 1001-2000

  # BLOOD      : S3 key which groups input data together, also S3_CHUNK_LOC key  
  # chunk234.fa: 'file' with chunked sequence in S3 
  # OUTPUT_DB  : key in Databases.pm  
  # range      - only take seq nr 1...1000 out of chunk file. optional.

  my $input_id = $self->original_input_id; 
  my ($base_batch, $rest ) = split "=",$self->original_input_id;   

  if ( !defined $base_batch && !defined $rest ) {  
    throw("input_id in wrong format. correct format : BLOOD = chunk_1453.fa :: OUTPUT_DB :: 1-1000");
  }  

  my ($chunk,$output_db_key,$range) = split "::",$rest;  

  if ( defined $range ) { 
     my ($nr_start_seq, $nr_end_seq) = split "-",$range; 
     $self->nr_start_seq($nr_start_seq);
     $self->nr_end_seq($nr_end_seq);
  } else {  
     warning("will run on whole chunk as no range has been supplied in input_id. To ".
             " supply a range use format :  BLOOD\@chunkXXX.fa::OUT_DB_KEY::1-1000 \n");
  } 
  $self->chunk($chunk); 
  $self->output_db_key($output_db_key);  
  $self->base_batch($base_batch);     


  if ( !defined $self->s3_sequence_data_key) {  
    throw "your batch key does not point to an s3_sequence_data_key \n";
  } 
  #print "setting bas batch to : $base_batch\n";  
  #print $self->base_batch ."\n\n";  

  # I will batch analyses togther by basebatch name, ie tissue or so to have an easer config handling later
  # set output db for write_output()     
  return $self;
}





sub get_file_from_s3 {  
  my ($self,$file_name,$bucket_name) = @_;


  print "Fetching file $file_name from $bucket_name \n"; 

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
     if ( (defined $self->gzip_compression && $self->gzip_compression==1)||  $file_name =~m/\.gz$/) {   
        print "un-compressing file\n";
        # file is compressed; let's un-compress it  
        my $cmd = "gunzip -fc $out_file > $out_file.part.tmp "; 
        system($cmd);     
        if ( -e $out_file ) { 
            system("unlink $out_file"); 
        } 
        my $of = $out_file; 
        $out_file =~s/\.gz$//; 
        $cmd="mv $of.part.tmp $out_file" ;  
        system($cmd);  
        if (!defined $is_genomic_file ) { 
           $self->files_to_delete($out_file); 
        } 
     }   
  }else {   
    # data no public - need S3 module to fetch data  
    throw("fetching data from S3 with password not implemented now\n"); 
  } 
  return $out_file; 
} 

sub extract_range_out_of_fasta {  
  my ( $self,$path_to_chunk_file  ) = @_;   

  my $inseq = Bio::SeqIO->new(-file=>"<$path_to_chunk_file", -format =>'fasta');  

  my $out_file_name = $path_to_chunk_file.".extract"; 
  
  my $outseq = Bio::SeqIO->new(-file=>">$out_file_name" , -format =>'fasta'); 
  
  my $nr_start = $self->nr_start_seq; 
  my $nr_end = $self->nr_end_seq; 

  my $seq_counter = 1; 
  my @all_seq; 
  SEQ: while ( my $seq = $inseq->next_seq ) {   
     if ( $seq_counter >= $nr_start && $seq_counter <=$nr_end ) {  
       print "extracting $seq_counter\n";  
       push @all_seq, $seq; 
     }   
     $seq_counter++;
  }
  for my $seq ( @all_seq) {  
    $outseq->write_seq($seq);
  }   
  print scalar(@all_seq) . " seqs written \n";    
  print "adding $out_file_name to list of files to deleete\n";  
  $self->files_to_delete($out_file_name); 
  return $out_file_name ; 
} 

sub fetch_input {
  my ($self) = @_;

  # get chunk1.fa.gz from S3   - security credentials in CloudConfig if needed

  my $path_to_chunk_file = $self->get_file_from_s3($self->chunk,$self->chunk_bucket_name);    

  if ( defined $self->nr_start_seq && defined $self->nr_end_seq) {  
    # range given, so we have to extract a range out of the fasta file  
    my $new_chunk_file_name = $self->extract_range_out_of_fasta($path_to_chunk_file);  
    $self->QUERYSEQS($new_chunk_file_name ) ;   
    if ( -e $path_to_chunk_file ) { 
      unlink($path_to_chunk_file);  
    }
  }else {   
   # no range given 
   $self->QUERYSEQS($path_to_chunk_file);   
  }

  # get softmasked genome sequence out of S3 - security credentials in CloudConfig if needed 

  my $path_to_genomic_file = $self->get_file_from_s3($self->genomic_file_name, $self->genomic_bucket_name);  
  $self->GENOMICSEQS($path_to_genomic_file); 
  $self->SUPER::fetch_input();  
}


sub write_output {  
  my ($self) = @_;  

  for my $file ( @{ $self->files_to_delete } ) { 
    if ( -e $file) {  
     my $cmd = "unlink $file"; 
     system($cmd);  
    } 
  } 
  $self->SUPER::write_output(); 
} 


# accessor methods for CloudConfig file 

sub chunk_bucket_name {  
  my ($self) = shift; 
  print "xxx xxx " . $self->s3_sequence_data_key."\n"; 
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

sub s3_sequence_data_key  {  
  my ($self) = shift;  
  #print "BASE BATCH IS : " . $self->base_batch ."\n\n";
  return ${$self->ANALYSIS_BASE_BATCH_CONFIG}{$self->base_batch}{S3_CHUNK_LOC}; 
} 


sub files_to_delete { 
  my ($self,$val) = @_; 

  if (defined  $val ) { 
    push @{$self->{_files_to_delete}},$val; 
  }   
  return $self->{_files_to_delete}; 
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
