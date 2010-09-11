
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

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(get_file_from_s3) ; 


sub get_file_from_s3 { 
  my ( $s3_bucket,$s3_file,$local_file ,$s3_config_file, $check) = @_ ;  

  if (!defined $s3_config_file ) { 
    $s3_config_file = $S3_CONFIG_FILE; 
  }  
  check_config_file ($s3_config_file); 

  my $cmd = "s3cmd -c $s3_config_file get $s3_bucket/$s3_file $local_file";
  system($cmd);  

  if ( $check ==1 ) {    
     my @lines = @{ get_md5sum($s3_bucket,$s3_file,$s3_config_file) } ;  
     for ( @lines ) {  
       print $_ . "\n"; 
     } 
     my $cmd = "s3cmd -c $s3_config_file ls --list-md5sum get $s3_bucket/$s3_file $md5sum_file";
    
  } 
}


sub get_md5sum {  
  my ( $s3_bucket,$s3_file,$s3_config_file ) = @_ ;  

  my $tmp_file;  
  my $cmd = "s3cmd -c $s3_config_file ls --list-md5 $s3_bucket/$s3_file > $tmp_file "; 
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
