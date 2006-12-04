
=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Utilities 

- base class which exports utility methods which don't take Bio::XX objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::Utilities qw(shuffle); 

  or 

  use Bio::EnsEMBL::Analysis::Tools:Utilities

  to get all methods

=head1 DESCRIPTION

This is a class which exports Utility methods for genebuilding and
other gene manupulation purposes. 

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::Utilities; 

use strict;
use warnings;
use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw( shuffle merge_config_details );




=head2 merge_config_details

  Arg [0]   : Array of Hashreferences
  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Function  : This func. merges the Configurations out differnt configuration-files into one Hash
  Returntype: Hashref. 
  Exceptions: throws as this method should be implemented by any child
  Example   : merge_database_configs ($DATABASES, $EXONERATE2GENES, $TRANSCRIPT_COALESCER) ; 

=cut


sub merge_config_details {
  my ($self,  @config_hashes )= @_ ;

  my %result ;

  # loop through all hrefs which are passed as input 

  foreach my $config_file ( @config_hashes ) {

    my %file = %$config_file ;

    foreach my $db_class ( keys %file ) {

      # process Exonerate2Genes.pm config (has section --> OUTDB)

      if ( exists ${$file{$db_class}}{OUTDB} ) {
        if ( defined ${$file{$db_class}}{OUTDB} && length(${$file{$db_class}}{OUTDB}{'-dbname'}) > 0  ) {
        # don't process undefiend OUT-DB's and  don't process defiened OUT-DB's which have no name
           #print "-dbname "  .${$file{$db_class}}{OUTDB}{'-dbname'}. "\n\n\n" ;

          $result{$db_class}{db} = ${$file{$db_class}}{OUTDB} ;

        }else {
         next ;
        }
      }

      # process /Conf/GeneBuild/Databases.pm 

      if (defined ( ${$file{$db_class}}{'-dbname'}) &&  length ( ${$file{$db_class}}{'-dbname'}) > 0 )  {
        # we process Databases.pm // parameteres for db-connection are ok
        $result{$db_class}{db} = \%{$file{$db_class}} ;

      } elsif (defined ( ${$file{$db_class}}{'-dbname'}) &&  length ( ${$file{$db_class}}{'-dbname'}) == 0  ) {
        next ;
      }

      # add / process data from other configs in format TranscriptCoalescer.pm 
      # and attach data to main config hash 

      for my $key (keys %{$file{$db_class}}) {
        $result{$db_class}{$key} = $file{$db_class}{$key};
      }
    }
  }
  return \%result ;
}






=head2 shuffle 

  Arg [1]   : Reference to Array
  Function  : randomizes the order of an array 
  Returntype: arrayref 
  Exceptions: none
  Example   : 

=cut

sub shuffle {
  my $tref = shift ;
  my $i = @$tref ;
  while ($i--) {
     my $j = int rand ($i+1);
     @$tref[$i,$j] = @$tref[$j,$i];
  }
  return $tref ;
}




1;
