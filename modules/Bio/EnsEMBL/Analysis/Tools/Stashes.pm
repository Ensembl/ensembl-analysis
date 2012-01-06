=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Stashes - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Tools::Stashes; 

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
