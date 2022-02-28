# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package ReadBaseDataFromGff; 

use warnings ;
use strict;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Exception qw( throw warning ); 
use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;

use Getopt::Long;



sub new {
  my ( $class, @args ) = @_ ;
  my $self = bless {}, $class ; 

  my ( $file, $path, $cs_name ) = rearrange ( ['FILE','PATH','CS_NAME'], @args ) ;

  $self->file($file) ;
  # setting defaults : 
  $self->path("ailMel1") ;
  $self->coord_system_name( "scaffold") ;

  # over-writing if defined as input 
  $self->path($path) if ( defined $path ) ;
  $self->coord_system_name($cs_name) if defined $cs_name ;

  return $self ;
}

sub file {
  my $self = shift ;
  $self->{file} =  shift if @_ ;
  return $self->{file} ;
}

sub path {
  my $self = shift ;
  $self->{path} =  shift if @_ ;
  return $self->{path} ;
}

sub coord_system_name {
  my $self = shift ;
  $self->{coord_system_name} =  shift if @_ ;
  return $self->{coord_system_name} ;
}


sub make_coord_sys {
  my ( $self, $coord_sys_name )  = @_ ;
  my $cs = Bio::EnsEMBL::CoordSystem->new( - NAME => $coord_sys_name, - RANK => 1 ) ;
  return $cs ; 
}  


sub get_slice { 
  my ( $self, $coord_system ) = @_ ;
  my $slice = Bio::EnsEMBL::Slice->new( - SEQ_REGION_NAME => "test_slice", - START => 1, - END => 100000, - COORD_SYSTEM => $coord_system ) ;
  return $slice ;
}

 
sub parse_gff_and_get_genes { 
  my ( $self ) = @_ ;  

  my %exon2trans ;
  my %utr2trans ;

  my $file = $self->file() ;
  my $cs = $self->make_coord_sys($self->coord_system_name()) ;
  my $slice = $self->get_slice($cs) ;
  open ( GFF, $file ) || die " cant read $file \n" ; 
  
  while ( my $line = <GFF> ) { 
    chomp($line) ;
  
    if ( $line =~ /^#/ ) {
    	next ;
    }
    
    my ( $scaffold, $analysis, $type , $g_start, $g_end , $bla, $g_strand, $g_phase , $description ) = split /\s+/, $line ; 
    my $strand = 1 ;
    if ( $g_strand =~ m/-/ ) {
      $strand  = -1 ;
    } 
  
    ## The type of file used here has a start phase but no end phase, so both are set to -1
    $g_phase = -1 ; 	
    my $end_phase = -1 ;
 
    if ($type =~m/CDS/) {
          my $newexon = new Bio::EnsEMBL::Exon (
            - start => $g_start,
            - end => $g_end,
            - strand => $strand,
            - slice  => $slice,
            - phase  => $g_phase,
            - end_phase => $end_phase,
          ) ;
  
       my @desc = split/\;/, $description ;  
       $description = $desc[0];  
       $description =~ s/Parent=// ;  
       push  @{ $exon2trans{$description}  } , $newexon ;    
    } elsif ($type =~ m/UTR/) {
          my $newexon = new Bio::EnsEMBL::Exon (
            - start => $g_start,
            - end => $g_end,
            - strand => $strand,
            - slice  => $slice,
            - phase  => $g_phase,
            - end_phase => $end_phase,
          ) ;

       my @desc = split/\;/, $description ;
       $description = $desc[0];
       $description =~ s/Parent=// ;
       push  @{ $utr2trans{$description}  } , $newexon ;
    } 
  }
  
  my @all_genes; 
  for my $tid ( sort keys %exon2trans ) {
       my @exons2 = @{ $exon2trans{$tid}} ;
       my @exons ;
       my @all_exons ;
       my $trans ;
  
       if ( $exons2[0]->strand == 1 ) {
         @exons = sort { $a->start <=> $b->start } @exons2 ;
       } else {
         @exons = sort { $b->start <=> $a->start } @exons2 ;
       }
       if ( defined $utr2trans{$tid} ) {
         push @exons, @{ $utr2trans{$tid} } ;
         if ( $exons[0]->strand == 1 ) {
           @all_exons = sort { $a->start <=> $b->start } @exons ;
         } else {
           @all_exons = sort { $b->start <=> $a->start } @exons ;
         }
         $trans = new Bio::EnsEMBL::Transcript ( - EXONS => \@all_exons ) ;
       } else {
         $trans = new Bio::EnsEMBL::Transcript ( - EXONS => \@exons ) ;
       }
       $trans->coding_region_start($exons[0]->start) ;
       $trans->coding_region_end($exons[-1]->end) ;
       $trans->slice($exons[0]->slice);  
  
       my $gene = Bio::EnsEMBL::Gene->new() ;
       $gene->add_Transcript($trans) ;
       $gene->slice($exons[0]->slice) ;
       push @all_genes, $gene ; 
  }
  return \@all_genes; 
} 	


sub get_Genes { 
  my ( $self ) = @_ ;
  return $self->parse_gff_and_get_genes() ;
} 

sub add_Translation {
  my ( $self, $gene ) = @_ ;
  my $translation = new Bio::EnsEMBL::Translation ;
  $translation->start(1) ;

  my $transcript = @{ $gene->get_all_Transcripts }[0] ;
  my @exons = @{ $transcript->get_all_Exons } ;
  my $start_exon ; my $end_exon ;

  foreach my $exon ( @exons ) {
    if ( $exon->start == $transcript->coding_region_start ) {
      $start_exon = $exon ;
    } elsif ( $exon->end == $transcript->coding_region_end ) {
      $end_exon = $exon ;
    }
  }

  $translation->start_Exon($start_exon) ;
  $translation->end_Exon($end_exon) ;
  $translation->end(length($translation->end_Exon)) ;
  $transcript->translation($translation) ;
  shift @{ $gene->get_all_Transcripts } ;
  $gene->add_Transcript($transcript) ;
  return $gene ;
}


#use vars '$AUTOLOAD';
#sub AUTOLOAD {
#  my ($self,$val) = @_;
#  ( my $routine_name = $AUTOLOAD ) =~ s/.*::// ; #trim package name
#  $self->{ $routine_name } = $val if defined $val ;
#  return $self->{ $routine_name } ;
#}


sub DESTROY {} # required due to AUTOLOAD


1; 
