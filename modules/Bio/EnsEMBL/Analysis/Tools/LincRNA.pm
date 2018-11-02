=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Tools::LincRNA;

=head1 SYNOPSIS

=head1 DESCRIPTION

Methods that used by more than one lincRNA related modules. 

=head1 CONTACT

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

=head1 METHODS
get_genes_of_biotypes_by_db_hash_ref
get_genes_of_biotypes

=cut


package Bio::EnsEMBL::Analysis::Tools::LincRNA;

use strict;
use warnings;
use Exporter;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             get_genes_of_biotypes_by_db_hash_ref
             get_genes_of_biotypes
            );



=head2 get_genes_of_biotypes_by_db_hash_ref

  Arg [0]   : Bio::EnsEMBL::Analysis::Tools::LincRNA
  Arg [1]   : Hashref. where the keys represent keys in the Databases-Hash,
              with a reference to an array where the values represent gene-biotypes
  Function  : Loops through the keys of the hash, builds a DB-connection to each DB specified by the key
              and fetches all genes of the biotypes specified

  Returntype: Returns an Array-reference to Bio::EnsEMBL::Gene objects which are located on the slice

=cut

sub get_genes_of_biotypes_by_db_hash_ref { 
  my ($self, $dbnames_2_biotypes) = @_;
  my @genes_to_fetch;
  foreach my $db_hash_key ( keys %$dbnames_2_biotypes )  {
    my @biotypes_to_fetch = @{$dbnames_2_biotypes->{$db_hash_key}};
    my $set_db = $self->hrdb_get_dba($self->param($db_hash_key));
    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    if($dna_dba) {
      $set_db->dnadb($dna_dba);
    }
    my $id = $self->param('iid');
    my $slice = $self->fetch_sequence($id, $set_db, undef, undef, $db_hash_key)  ;

    # implementation of fetch_all_biotypes ....  
    my $fetch_all_biotypes_flag ; 
    foreach my $biotype  ( @biotypes_to_fetch ) {  
      print "checking biotype: $biotype  \n"; 
      if ($biotype=~m/fetch_all_biotypes/ ) {    
        $fetch_all_biotypes_flag = 1 ; 
      }
    }  
    if ( $fetch_all_biotypes_flag ) {  
         print  "fetching ALL biotypes for slice out of db $db_hash_key :\n" ; 
         my $genes = $slice->get_all_Genes(undef,undef,1) ; 
         push @genes_to_fetch, @$genes;
         print scalar(@genes_to_fetch) . " genes fetched in total\n";
    } else { 
      foreach my $biotype  ( @biotypes_to_fetch ) { 
         my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
         if ( @$genes == 0 ) {
           $self->warning("No genes of biotype $biotype found in $set_db\n");
         } 
         # if ( $self->verbose ) { 
           # print  "$db_hash_key [ " . $set_db->dbname  . " ] Retrieved ".@$genes." of type ".$biotype."\n";
           # print  "DEBUG::HiveLincRNA::get_genes_of_biotypes_by_db_hash_ref: " . $db_hash_key . " Retrieved ".@$genes." of biotype ".$biotype."\n";
         # }
         push @genes_to_fetch, @$genes;
      }  
    }  
  } 
  return \@genes_to_fetch;
}




=head2 get_genes_of_biotypes 

  Arg [0]   : Bio::EnsEMBL::Analysis::Tools::LincRNA
  Arg [1]   : Database name and biotype, where the values represent gene-biotypes
  Function  : builds a DB-connection to a DB and fetches all genes of the specified biotype 
  Returntype: Returns an Array-reference to Bio::EnsEMBL::Gene objects which are located on the slice or the whole genome

=cut

sub get_genes_of_biotypes {
	my ( $self, $hbiotype, $dbName ) = @_;
	my @genes_to_fetch;
	my $set_db  = $self->hrdb_get_con($dbName);
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $set_db->dnadb($dna_dba);
  }

  if (defined($self->param('iid'))) {
 	  my $id = $self->param('iid'); 
	  my $slice   =	$self->fetch_sequence( $id, $set_db); 
		my $genes = $slice->get_all_Genes_by_type( $hbiotype, undef, 1 ); 
		if ( @$genes == 0 ) { 
			warn("No genes of biotype $hbiotype found in $set_db (it is possible) \n"); 
		} 
		print " [ "	. $set_db->dbc->dbname	. " ] Retrieved "	. @$genes . " of type "	. $hbiotype . "\n";
		push @genes_to_fetch, @$genes; 
  } else { 
  	my $sa = $set_db->get_SliceAdaptor; 
    my $count =0; 
    my $how_many = 0; 
    if ($hbiotype eq "fetch_all_biotypes") {
       foreach my $sliceC (@{$sa->fetch_all('toplevel')}) { 
         my $genes = $sliceC->get_all_Genes(undef,undef,1) ; 
         push @genes_to_fetch, @$genes;
      }
    	
    } else {
      foreach my $sliceB (@{$sa->fetch_all('toplevel')}) { 
        foreach my $g ( @{ $sliceB->get_all_Genes_by_type( $hbiotype, undef, 1 ) } ) {
      	  push @genes_to_fetch, $g;  
        }
      }
    }
  }
	return \@genes_to_fetch;
}




1;
