
=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::S3Utils - base class for data retrieval from S3 storage 


=head1 SYNOPSIS

  Bio::EnsEMBL::Analysis::Tools::S3Utils qw(get_file_from_s3) 


  Bio::EnsEMBL::Analysis::Tools::S3Utils 

  to import all methods 

=head1 DESCRIPTION

This is a base class for Utility modules for data retreival from S3. 
The module is wrappign the s3cmd command. This command requires a 
config file which can be created with the s3cmd --configure option 

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS 

the rest of the documention details the exported static class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::S3Utils; 

use strict;
use warnings;
use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump); 


use Bio::EnsEMBL::Analysis::Config::S3Config; 
use Bio::EnsEMBL::Analysis::Config::General; 
use Digest::MD5; 
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ; 
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             get_file_from_s3
             get_file_from_s3_and_gunzip
            ) ; 



# mini routine to overwrite config vars. if direct value passed in, direct value used, 
# if not standard val in S3Conf used, if set no config ( ie $HOME/.s3cfg ) will be used 

sub s3_conf_para {   
  my ( $s3_config_file ) = @_; 
  my $s3_conf; 
  if (defined $s3_config_file ) { 
      $s3_conf = $s3_config_file; 
  }elsif ( defined $S3_CONFIG_FILE ) {  
      $s3_conf =$S3_CONFIG_FILE;
  }  
  my $use_config = "";
  if ( defined $s3_conf ) { 
    check_config_file($s3_conf);  
    $use_config = " -c $s3_conf ";
  } 
  return $use_config;
} 



sub get_file_from_s3_and_gunzip { 
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file ) = @_ ; 

  my $lf =  get_file_from_s3($s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file);  
  my $lfz = $lf."gz";
  system("mv $lf $lfz");
  my $outf = $lf.".$$.gunzip"; 
   
  gunzip $lfz=> $outf or throw(" Gunzip failed : $GunzipError");
  return $outf; 
}


sub get_file_from_s3 { 
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file ) = @_ ;  

  my $use_config = s3_conf_para($s3_config_file); 

  my $local_file = $local_dir ."/".$s3_file.".$$";
  my $cmd = "s3cmd --force $use_config get $s3_bucket/$s3_file $local_file"; 
  system($cmd);  

  if (defined $check && $check ==1 ) {     
     my $md5;
     my @lines = @{ get_md5sum($s3_bucket,$s3_file,$s3_config_file) } ;  
     for my $l ( @lines ) {  
       my @item = split /\s+/,$l; 
       if ( $item[4]=~m/$s3_file/) { 
         $md5=$item[3]; 
       }
     }  
     if ( check_file_vs_md5 ($md5,$local_file) ) {  
       return $local_file; 
     } else {  
       throw("md5sums for downloded file and local file do not match \n"); 
     } 
  }  
  return $local_file;
}



sub check_file_vs_md5{ 
  my ( $md5,$local_file) = @_  ;  

  my $ctx = Digest::MD5->new;
  open(F,$local_file); 
  $ctx->addfile(*F); 
  my $ms = $ctx->hexdigest;   
  if ( $ms eq $md5 ){  
    return 1;
  } 
  return 0 ; 
} 

sub get_md5sum {  
  my ( $s3_bucket,$s3_file,$s3_config_file ) = @_ ;   

  my $use_config = s3_conf_para($s3_config_file); 

  my @all;
  my $cmd = "s3cmd $use_config ls --list-md5 $s3_bucket/$s3_file"; 
  local *SCMD;
  open(SCMD, "$cmd 2>&1 |") or throw("couldn't open pipe s3cmd");
  while(my $line = <SCMD>){ 
    push @all, $line ; 
  } 
  return \@all; 
} 


sub check_config_file {  
  my ($file) = @_;  
  
  if (! -e $file ) {  
     throw("config file $file does not exist\n"); 
  } 
} 


=head2 id

  Arg [1]   : Bio::EnsEMBL::Feature
  Function  : Returns a string containing an appropriate label
              for the feature
  Returntype: string
  Exceptions: none
  Example   : 

=cut



sub id {
  my $feature = shift;
  my $id;

  if($feature->can('stable_id') && $feature->stable_id){
    $id = $feature->stable_id;
  }elsif($feature->can('dbID') && $feature->dbID) {
    $id = $feature->dbID;
  }else{ 
    $id = 'no-id';
  }
  if($feature->can('biotype') && $feature->biotype){
    $id .= "_".$feature->biotype;
  }
  return $id;
}



=head2 empty_Object

  Arg [1]   : Bio::EnsEMBL::Storeable or an object which inherits from it
  Arg [2]   : Boolean, whether to remove the stable id from the given object
  Function  : remove the dbID, adaptor and if appropriate the stable id
  Returntype: Bio::EnsEMBL::Storeable
  Exceptions: n/a
  Example   : empty_Object($object);

=cut



sub empty_Object{
  my ($object, $include_stable_id) = @_;
  $object->adaptor(undef);
  $object->dbID(undef);
  $object->stable_id(undef) if($object->can("stable_id") && 
                               $include_stable_id);
  return $object;
}



=head2 lies_inside_of_slice

  Arg [1]   : Bio::EnsEMBL::Feature
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : Ensures the transcript within the slice, completely on
              the lower end, it can overhang the upper end.
  Returntype: Boolean, 1 for pass, 0 for fail, i.e. lies outside of slice
              or across lower boundary
  Exceptions: none
  Example   :

=cut


sub lies_inside_of_slice{
  my ($feature, $slice) = @_;
  if($feature->start > $slice->length || 
     $feature->end < 1){
    warning(id($feature)." lies off edge if slice ".
            $slice->name);
    return 0;
  }
  if($feature->start < 1 && $feature->end > 1){
    warning(id($feature)." lies over lower boundary".
            " of slice ".$slice->name);
    return 0;
  }
  return 1;
}


1;
