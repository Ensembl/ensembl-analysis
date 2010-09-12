
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
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $md5sum_uncompressed, $s3_config_file) = @_ ; 

  # check if uncompressed file exsts 
  my $original_file_name;  
  if ( $s3_file =~m/\.gz/ ) { 
     (my $original_file_name = $s3_file)=~s/\.gz//g;    
     my $odir_local = $local_dir."/".$original_file_name;  
     if ( -e $odir_local) {  
        # get uploaded,uncompressed md5sum
        print "un-compressed  file exists!!! $odir_local\n";  
        my $md5_uncompressed_uploaded = get_file_from_s3($s3_bucket,$md5sum_uncompressed,$local_dir ,1, $s3_config_file);   
        print "file with uploaded md5-sums : $md5_uncompressed_uploaded\n";  
        # now compare md5sum of uploaded file with uncompressed md5sums with md5sum of local \nfile 
         my $md5_local = create_md5sum_for_local_file($odir_local);
         print "local md5: $md5_local for $s3_file\n"; 
        if ( compare_md5_sums ( $md5_uncompressed_uploaded,$s3_file,$md5_local) == 1 ){   
          print "uncompressed file exits, has same md5 as uploaed md5 with uncompressed sums \n";  
          print "returning $odir_local\n"; 
          return $odir_local; 
        } 
     } else { 
        print "uncompressed file does not exist : $odir_local\n"; 
     } 
  }  
  print "NOW fetching $s3_file\n"; 
  my $lf =  get_file_from_s3($s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file);    
  # file has been  downloaded correctly and automatic S3 md5sum has been checked 
  # now we have to gunzip the file 

  my $lfz = $lf."gz";
  system("mv $lf $lfz");
  my $outf = $lf.".$$.gunzip"; 
   
  gunzip $lfz=> $outf or throw(" Gunzip failed : $GunzipError");   

  print "file gunziped to $outf\n";   

  #system("mv $outf $lf"); 

  if ( defined $md5sum_uncompressed ) {   
    print "WOW have also md5sum uncompressed\n";
    # user has uploaded md5sum file for the uncompressed data. we can now check against that. 
    my $md5_local = create_md5sum_for_local_file($outf);    
    my $md5_uploaded_file = get_file_from_s3($s3_bucket,$md5sum_uncompressed,$local_dir ,1, $s3_config_file);

    # file names are different in uploaded md5sum file and original file  
    print "downloaded 2 files - need to compare uncompressed md5sum for file with ungzipped md5sum\n"; 
    print "md5 uploaded : $md5_uploaded_file\n";  
    my $compare_md5_sums = compare_md5_sums ( $md5_uploaded_file,$s3_file,$md5_local);
    if ( defined $compare_md5_sums && $compare_md5_sums == 1 ) {   
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

  print "md5_uploaded : $md5_uploaded_file \n"; 
  print "s3  file     : $s3_file\n";  
  print "md5 local    : $md5_local\n";
  
  print "file with uploaded md5-sums: $md5_uploaded_file\n";  

  open(F,"$md5_uploaded_file")  || die "can't open file $md5_uploaded_file\n";    
  my $match ; 
  FILE: while (my $line=<F>){     
     chomp($line); 
     my ( $md5sum_uploaded,$file_name) = split /\s+/,$line;
     if ( $md5sum_uploaded=~m/$md5_local/ ) { 
       print "MATCH :  $md5sum_uploaded    VS  $md5_local \n";   
       $match = 1; 
       last FILE; # md5sums match 
     }
  } 
  close(F);  
  return $match ; 
  # s3_file ' chunk0.fa.gz 
  # fn = chunk0.fa 
} 


sub get_file_from_s3 { 
  my ( $s3_bucket,$s3_file,$local_dir ,$check, $s3_config_file ) = @_ ;  

  print "fetching file from s3 : $s3_file\n";  

  my $use_config = s3_conf_para($s3_config_file); 
  
  my $tmp_local_file = $local_dir ."/".$s3_file.".$$"; 
  # dowload file to xxx.$$
  my $cmd = "s3cmd --force $use_config get $s3_bucket/$s3_file $tmp_local_file";   

  system($cmd);  

  print "data downloaded and saved to $tmp_local_file\n";  

  # check integrity of downloaded file vs S3 automatically stored md5sum
  # this assures that the download went OK   

  if (defined $check && $check ==1 ) {
     print "checking md5sum of downloaded file vs. automatic S3 md5sum\n";
     my $md5;
     my @lines = @{ get_s3_stored_md5sum($s3_bucket,$s3_file,$s3_config_file) } ;  

     # problem : routine returns all lines matching, 
     # ie md5.sum     md5.sum.gz      md5.sumgz   

     LINES: for my $line ( @lines ) {
       my ( $date, $time, $size, $md5string, $s3_location ) = split /\s+/,$line;  
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

  my $ctx = Digest::MD5->new;
  open(F,$local_file); 
  $ctx->addfile(*F); 
  my $ms = $ctx->hexdigest;    
  return $ms; 
} 



sub check_md5_downloaded_file_vs_md5_from_s3 { 
  my ( $md5,$local_file) = @_  ;  

  my $md5_local = create_md5sum_for_local_file( $local_file);
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


\1;
