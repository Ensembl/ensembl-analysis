
=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::S3Utils - base class for data retrieval from S3 storage 


=head1 SYNOPSIS

  Bio::EnsEMBL::Analysis::Tools::S3Utils qw(get_uncompressed_file_from_s3) 

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
             get_uncompressed_file_from_s3
             get_gzip_compressed_file_from_s3
            ) ; 


=head2 get_file_from_s3 

  Arg [1]   : name of S3 bucket ( String ) 
  Arg [2]   : name of file in S3 
  Arg [3]   : local directory where to store the downloaded file 
  Arg [4]   : flag if downloade file should be checked vs automatic md5sum in amazon S3 storage 
  Arg [5]   : location of S3 config file ( Set in Config/S3Config.pm  ) 
  Arg [6]   : md5sum of uncompressed file  

  Function  : Fetches a file from Amazon Simple Storage Service S3 with s3cmd. 
              Requires that the program s3cmd is properly set up ( see s3cmd --help, 
              s3cmd --configure  
              If no config file is given, it defaults to $HOME/.s3cfg 
  Returntype: Boolean, 1 for pass, 0 for fail, i.e. lies outside of slice
              or across lower boundary
  Exceptions: none
  Example   :  get_file_from_s3("s3://ensembl-cloud-cvs","test.fa","1","/home/ensembl/.s3cfg",'');

=cut

sub get_file_from_s3 { 
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file, $md5sum_uncompressed_file) = @_ ;  

  if ( $s3_bucket !~m/^s3:\/\//){ 
   throw("Your s3-bucket name is incorrect - EXAMPLE : s3://ensembl-cvs-bucket ");
  }  

  if ( $s3_file =~m/gz/ ) { 
    return get_gzip_compressed_file_from_s3($s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file, 
                                            $md5sum_uncompressed_file);
  }else {   
    # get_uncompressed_file_from_s3 does not need $md5sum_uncompressed_file as we can directly compare 
    # vs automatically stored md5sum in amazon s3 
    return get_uncompressed_file_from_s3($s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file); 
  }
}



=head2 get_gzip_compressed_file_from_s3

  Arg [1]   : name of S3 bucket ( String ) 
  Arg [2]   : name of file in S3 
  Arg [3]   : local directory where to store the downloaded file 
  Arg [4]   : flag if downloade file should be checked vs automatic md5sum in amazon S3 storage 
  Arg [5]   : location of S3 config file ( Set in Config/S3Config.pm 
  Arg [6]   : Name of file in S3 which contains the md5sums of the uncompressed files in the bucket. 

  Function  : Fetches a file from Amazon Simple Storage Service S3 with s3cmd. 
              Requires that the program s3cmd is properly set up ( see s3cmd --help, 
              s3cmd --configure  
              If no config file is given, it defaults to $HOME/.s3cfg  

              Arg[6] is needed to check if the gunzipped / uncompressed file is identical with the 
              uploaded, compresssed file. Sometimes un-compressing can fail so this assures that data is 
              consistent

  Returntype: String describing path to downloaded file 
              or across lower boundary
  Exceptions: none
  Example   :  get_gzip_compresssed_file_from_s3(
                      "s3://ensembl-cloud-cvs",
                      "test.fa", 
                      "/tmp/",
                      "1",
                      "/home/ensembl/.s3cfg",
                      'file_with_md5sums_of_uncompressed_file_in_s3_bucket'
                     );

=cut



sub get_gzip_compressed_file_from_s3 { 
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file, $md5sum_uncompressed) = @_ ;  


  # check if uncompressed file exsts 
  my $original_file_name;  
  if ( $s3_file =~m/\.gz/ ) { 
     (my $original_file_name = $s3_file)=~s/\.gz//g;    
     my $odir_local = $local_dir."/".$original_file_name;  
     if ( -e $odir_local) {  
        # get uploaded,uncompressed md5sum
        print "\n\nun-compressed  file exists!!! $odir_local\n\n";   

       if ( defined $md5sum_uncompressed ) {  
         # the stuff below is to be able to compare the checksum of the uncompressed files vs the original file
         # to be be sure the the compression/un-compression went OK 
         print "\ngetting file with md5sums of uncompressed fasta files\n"; 
         my $md5_uncompressed_uploaded = 
            get_uncompressed_file_from_s3($s3_bucket,$md5sum_uncompressed,$local_dir ,1, $s3_config_file);   

         print "\n\nfile with uncompressed md5-sums : $md5_uncompressed_uploaded\n\n";   

         # now compare md5sum of uploaded file with uncompressed md5sums with md5sum of local \nfile 
         my $md5_local = create_md5sum_for_local_file($odir_local); 
         
         print "local md5: $md5_local for $s3_file\n";   

         
         my  $md5sums_match = compare_md5_sums ( $md5_uncompressed_uploaded,$s3_file,$md5_local) ;
 
         unlink($md5_uncompressed_uploaded);
         if ( $md5sums_match == 1 ) { 
           return $odir_local; 
         }else{  
           warning("md5sum of downloaded,uncompressed file $odir_local does not match md5sum of uploaded file - Will re-download file");
           unlink($odir_local); 
         } 
       }
     } else { 
        print "uncompressed file does not exist : $odir_local\n"; 
     } 
  }  
  print "NOW fetching $s3_file\n"; 
  my $lf =  get_uncompressed_file_from_s3($s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file);    
  # file has been  downloaded correctly and automatic S3 md5sum has been checked 
  # now we have to gunzip the file 

  my $lfz = $lf."gz";
  system("mv $lf $lfz");

  my $outf = $lf.".$$.gunzip"; 
   
  gunzip $lfz=> $outf or throw(" Gunzip failed : $GunzipError");   

  print "file gunziped to $outf\n";   
  system("unlink $lfz");
  #system("mv $outf $lf"); 
  if ( defined $md5sum_uncompressed ) {   
    print "WOW have also md5sum uncompressed\n";
    # user has uploaded md5sum file for the uncompressed data. we can now check against that. 
    my $md5_local = create_md5sum_for_local_file($outf);    
    my $md5_uploaded_file = get_uncompressed_file_from_s3($s3_bucket,$md5sum_uncompressed,$local_dir ,1, $s3_config_file);

    # file names are different in uploaded md5sum file and original file  
    print "downloaded 2 files - need to compare uncompressed md5sum for file with ungzipped md5sum\n"; 
    print "md5 uploaded : $md5_uploaded_file\n";   

    my $compare_md5_sums= compare_md5_sums ( $md5_uploaded_file,$s3_file,$md5_local);

    unlink($md5_uploaded_file);  

    if ( $compare_md5_sums == 1 ) {   
      print "all OK - md5sums match\n";  
    }  
    
  }
  if ( $s3_file =~m/\.gz/){  
      ($original_file_name = $s3_file)=~s/\.gz//g; 
      my $new_name = "$local_dir/$original_file_name"; 
      if ( ! -e $new_name ) { 
         print "RENAME: $outf $new_name\n"; 
         system("mv $outf $new_name");  
      }
      return $new_name;  
  }
  return $outf; 
}


sub compare_md5_sums {  
  my ( $md5_uploaded_file,$s3_file,$md5_local)= @_;

  print "CMP md5 : md5_uploaded : $md5_uploaded_file \n"; 
  print "CMP md5 : s3  file     : $s3_file\n";  
  print "CMP md5 : md5 local    : $md5_local\n";
  
  print "file with uploaded md5-sums: $md5_uploaded_file\n";  

  open(F,"$md5_uploaded_file")  || throw ("can't open file $md5_uploaded_file\n"); 
  my $match =0; 
  FILE: while (my $line=<F>){     
     chomp($line);  
     print "MD5-file-line: $line\n";
     my ( $md5sum_uploaded,$file_name) = split /\s+/,$line;
       print "COMPARE :  $md5sum_uploaded    VS  $md5_local \n";   
     if ( $md5sum_uploaded=~m/$md5_local/ ) { 
       print "MATCH :  $md5sum_uploaded    VS  $md5_local for $file_name \n";   
       $match = 1; 
       last FILE; # md5sums match 
     }
  } 
  close(F);   
  return  $match ;
  # s3_file ' chunk0.fa.gz 
  # fn = chunk0.fa 
} 




sub get_uncompressed_file_from_s3 { 
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file ) = @_ ;  

  print "fetching file from s3 : $s3_file\n";  

  my $use_config = s3_conf_para($s3_config_file); 
  
  my $tmp_local_file = $local_dir ."/".$s3_file.".$$"; 
  my $cmd = "s3cmd --force $use_config get $s3_bucket/$s3_file $tmp_local_file";   
  system($cmd);  

  print "data downloaded and saved to $tmp_local_file\n";  

  # check integrity of downloaded file vs S3 automatically stored md5sum
  # this assures that the download went OK   

  if (defined $check && $check ==1 ) {
     print "checking md5sum of downloaded file vs. automaticly generated S3 md5sum\n";
     my $md5;
     my @lines = @{ get_s3_stored_md5sum($s3_bucket,$s3_file,$s3_config_file) } ;  

     # problem : routine returns all lines matching, 
     # ie md5.sum     md5.sum.gz      md5.sumgz   

     LINES: for my $line ( @lines ) {
       my ( $date, $time, $size, $md5string, $s3_location ) = split /\s+/,$line;   
       print "line $line\n"; 
       #2010-09-08 00:52    151635   a509fd9c8b41944d73552c36728a1b68  s3://ensembl-cloud-chunks/md5.sum 
       #print "LOC $s3_location\n"; 
       $s3_location =~s/s3:\/\///;
       my @location_string = split /\//,$s3_location; 
       my $file_name = pop @location_string ;  
       print "-$file_name-\n";  
       if ( $file_name eq $s3_file ) { 
         $md5 = $md5string ;  
         last LINES;
       } 
     }  
     if ( ! defined $md5 ) {  
        throw("md5 for file $s3_file not found\n");
     } 
     if ( check_md5_downloaded_file_vs_md5_from_s3 ($md5,$tmp_local_file) ) {   
       print "file downloaded correctly as s3-automatic md5sum match with downloaded file : $tmp_local_file\n"; 
        #   if ( ! -e $s3_file ) {  
        #    system("mv $tmp_local_file $s3_file");  
        #    print "moving $tmp_local_file $s3_file \n";
        #    exit(0); 
        return $tmp_local_file; 
     } else {  
       throw("md5sums for downloded file does not match the md5sum stored automatically in S3: $md5 - $tmp_local_file \n"); 
     } 
  } else { 
    print " not checking file \n" ;  
  }
   return $tmp_local_file;
}


sub create_md5sum_for_local_file { 
  my ($local_file) = @_  ;  

  print "Creating md5sum for $local_file ...\n"; 
  my $ctx = Digest::MD5->new;
  open(F,$local_file); 
  $ctx->addfile(*F); 
  my $ms = $ctx->hexdigest;    
  print "$local_file ==> $ms\n"; 
  return $ms; 
} 



sub check_md5_downloaded_file_vs_md5_from_s3 { 
  my ( $md5,$local_file) = @_  ;  

  my $md5_local = create_md5sum_for_local_file( $local_file); 
   print "md5 of local file ( $local_file : $md5_local ) \n";
  if ( $md5_local  eq $md5 ){   
    return 1;
  } 
  return 0 ; 
} 

sub get_s3_stored_md5sum {  
  my ( $s3_bucket,$s3_file,$s3_config_file ) = @_ ;   

  my $use_config = s3_conf_para($s3_config_file); 

  my @all;
  my $cmd = "s3cmd $use_config ls --list-md5 $s3_bucket/$s3_file";  
  print "CMD: $cmd\n";
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



=head2 s3_conf_para 

  Arg [1]   : String, path to config file ie "/home/ensembl/.s3cfg" 
  Function  : overwrite config vars for S3_CONFIG_FILE 
              If a value is directly passed into this methods, this value will be used.
              If S3_CONFIG_FILE set to undef in S3Config.pm, the default value ( ie $HOME/.s3cfg ) 
              will be used 
  Returntype: String or undef  

=cut




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




1;
