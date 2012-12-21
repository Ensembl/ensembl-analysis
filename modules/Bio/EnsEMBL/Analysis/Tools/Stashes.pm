# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Tools/Stashes.pm,v $
# $Revision: 1.4 $
package Bio::EnsEMBL::Analysis::Tools::Stashes; 

use warnings ;
use strict ;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Exporter ; 
use vars qw (@ISA @EXPORT); 

@ISA = qw(Exporter);
@EXPORT = qw( package_stash );

sub package_stash {
    my ($packageName) = @_;  

    my %result ;  

    local (*alias);  
    *stash = *{"${packageName}::"};   

    while (($varName, $globValue) = each %stash) {   
        # only return the config hash 
        next if $varName =~m/BEGIN/; 
        next if $varName =~m/import/; 

        *alias = $globValue;  
        $result{$varName}=$alias  if (defined($alias)) ; 
        $result{$varName}=\@alias if (*alias{ARRAY}) ; 
        $result{$varName}=\%alias if (*alias{HASH}) ;  
     }     

     if (scalar(keys %result >1) ) {    
       throw("Have more than one item exported from $packageName - you'll run into trouble\n")
     }  
     my $hash_name = (keys %result)[0];   
     return [$result{$hash_name},$hash_name] ; 
}

1;
