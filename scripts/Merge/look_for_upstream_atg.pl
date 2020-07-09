#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
# Finds alternative ATGs for transcripts and adds them as transcripts
# attributes. By default we only look at 200 bp upstream of the UTR region.

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host;
my $user;
my $dbname;
my $port          = 3306;
my $pass;
my $upstream_dist = 200;
my @defgenetypes  = ('protein_coding');
my $coordsystem   = 'toplevel';
my @chromosomes;
my @genetypes;

my $dna_host;
my $dna_port = '3306';
my $dna_user;
my $dna_pass;
my $dna_dbname;
my $path;

GetOptions( 'dbhost|host|h:s'        => \$host,
            'dbuser|user|u:s'        => \$user,
            'dbname|db|D:s'          => \$dbname,
            'dbpass|pass|p:s'        => \$pass,
            'dbport|port|P:n'        => \$port,
            'dna_host=s'             => \$dna_host,
            'dna_port=s'             => \$dna_port,
            'dna_user=s'             => \$dna_user,
            'dna_pass=s'             => \$dna_pass,
            'dna_dbname=s'           => \$dna_dbname,
            'upstream_dist:n'        => \$upstream_dist,
            'chromosomes:s'          => \@chromosomes,
            'genetypes:s'            => \@genetypes,
            'coord_system|cs_name:s' => \$coordsystem,
            'path:s'                 => \$path );

if ( scalar(@chromosomes) ) {
  @chromosomes = split( /,/, join( ',', @chromosomes ) );
}

if ( scalar(@genetypes) ) {
  @genetypes = split( /,/, join( ',', @genetypes ) );
} else {
  @genetypes = @defgenetypes;
}

my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -pass   => $pass,
                                      -port   => $port,
                                      -path   => $path,
                                      -dbname => $dbname );

my $dnadb;
if ($dna_dbname) {
  $dnadb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dna_host,
                                        -user   => $dna_user,
                                        -pass   => $dna_pass,
                                        -port   => $dna_port,
                                        -path   => $path,
                                        -dbname => $dna_dbname );

  $db->dnadb($dnadb);
}

my $sa = $db->get_SliceAdaptor();
my $asm_exception_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();

print "Starting the checks\n";

foreach my $chromosome (@chromosomes) {

  ###
  # check for alternate sequence...
  # it may happen that we have more than one assembly exception on a
  # chromosome, and they may be of more than one 'type'
  # eg. HAP and PATCH
  ###
  my $check_slice = $sa->fetch_by_region( $coordsystem, $chromosome );
  # if $check_slice does not exist, try to fetch by name as it can be an Ensembl slice name
  if (!$check_slice) {
    $check_slice = $sa->fetch_by_name($chromosome);
  }
  my @except_feat = @{ $asm_exception_adaptor->fetch_all_by_Slice($check_slice) };
  foreach my $aef (@except_feat) {
    my $alt_slice = $aef->alternate_slice();
    warning("  Alternate slice is " . $alt_slice->name . "\n");
  }
  if ($check_slice->is_reference) {
    warning(" Reference slice $chromosome has ".(scalar(@except_feat))." assembly exceptions so please make sure you run this script for those exceptions too.\n");
  } else {
    warning(" Working with a non-reference slice $chromosome \n");
  }

  ###
  # get unique part(s) of slice
  # chr Y in human will have 2 unique parts because of the PAR
  # and for patches we also only want unique regions
  ###
  my @unique_slices = @{$sa->fetch_by_region_unique( $coordsystem, $chromosome  )};

  foreach my $slice (@unique_slices) {
    my $srname = $slice->seq_region_name();

    ###
    # Remove attribs for this seq_region
    ###
    print "Deleting existing 'ATG' transcript attributes for slice ". $slice->name ."\n";
    my $dsth = $db->dbc()->prepare("DELETE ta.* FROM seq_region sr, transcript t, transcript_attrib ta, attrib_type at WHERE sr.seq_region_id = t.seq_region_id and t.transcript_id=ta.transcript_id and at.attrib_type_id=ta.attrib_type_id AND at.code='upstream_ATG' and sr.name = '$srname'");
    $dsth->execute();


    ###
    # Now fetch the genes and look for upstream ATGs
    ###
    print "Fetching genes for " . $slice->name . "\n";
    my %genes_hash;
    my @genes;

    foreach my $genetype (@genetypes) {
      my @tmp_genes = @{ $slice->get_all_Genes_by_type($genetype) };
      print "Fetched ". (scalar(@tmp_genes)) ." $genetype genes for slice : " . $slice->name . "\n";
      push @genes, @tmp_genes;

    }

    print "Done fetching genes (fetched " . scalar(@genes) . ")\n";

    @genes = sort { $a->start <=> $b->start } @genes;

    foreach my $gene (@genes) {
      print "Starting gene dbID ".$gene->dbID."\n";

      foreach my $transcript ( @{ $gene->get_all_Transcripts } ) {
        # First as a paranoid check we want to confirm that the
        # transcript has translation
        if ( $transcript->translation ) {

          my %alt_ATG_starts;

          if ( $transcript->strand == 1 ) {
            %alt_ATG_starts =
              check_transcript( $transcript, $chromosome, $upstream_dist );

          } elsif ( $transcript->strand == -1 ) {
            %alt_ATG_starts =
              check_transcript( $transcript, $chromosome, $upstream_dist );

          } else {
            throw("Strand not recognised");
          }

          if (scalar(keys %alt_ATG_starts) > 0) {
            print "Updating transcript ",$transcript->dbID," \n";
            # We have attributes
            # Remember to first fetch all attribs we already have
            $transcript->get_all_Attributes();

            # Now make the attributes
            my $attributes = make_attributes(\%alt_ATG_starts);


            # Now add the new Attributes to the Transcript
            my $attr_adaptor = $db->get_AttributeAdaptor();
            $attr_adaptor->store_on_Transcript( $transcript->dbID,
                                                $attributes );
           # or do we do $transcript->add_Attributes(@); and then $transcript->store();
          }
        } ## end if ( $trans->translation)
      } ## end foreach my $trans ( @{ $gene...

    } ## end foreach my $gene (@genes)
  } ## end foreach my $slice
} ## end foreach my $chromosome (@chromosomes)


sub check_transcript {
  my ( $trans, $chromosome, $upstream_dist ) = @_;

  my $coding_region_start;
  my $trans_start;
  my %alt_atgs;
  my $has_utr             = 0;
  my $seq_to_check;

  if ($trans->strand == 1 ) {
    $coding_region_start = $trans->get_all_translateable_Exons->[0]->seq_region_start;
    $trans_start         = $trans->seq_region_start;
    print "Starting forward transcript with dbID : " . $trans->dbID . "\n";
    print "  this is the original forward coding region start: ", $coding_region_start, "\n";
    print "  this is the forward transcript start: ", $trans_start, "\n";
  } elsif ($trans->strand == -1) {
    $coding_region_start = $trans->get_all_translateable_Exons->[0]->seq_region_end;
    $trans_start         = $trans->seq_region_end;
    print "Starting reverse transcript with dbID " . $trans->dbID . "\n";
    print "  this is the original reverse coding region start: ", $coding_region_start, "\n";
    print "  this is the reverse transcript start: ", $trans_start, "\n";
  }


  # Check if the gene has UTR
  if ( $coding_region_start != $trans_start && defined($trans->five_prime_utr) ) {
    $has_utr = 1;
    # As the gene has UTR we only want to check within the UTR region
    my $start_exon   = $trans->translation->start_Exon;

    my $upstream_seq = $trans->five_prime_utr->seq;

    if ($trans->strand == 1 ) {
      print "  has utr. Seq on forward strand transcript is: ", $upstream_seq, "\n";
    } elsif ($trans->strand == -1) {
      print "  has utr. Seq on reverse strand transcript is: ", $upstream_seq, "\n";
    }

    if ( length($upstream_seq) < $upstream_dist ) {
      $upstream_dist = length($upstream_seq);
    }
    $seq_to_check =
      substr( $upstream_seq, length($upstream_seq) - $upstream_dist,
              $upstream_dist );
  } else {

    my $slice_to_check;
    if ($trans->strand == 1 ) {
      $slice_to_check =
      $sa->fetch_by_region( 'toplevel', $chromosome,
                            $coding_region_start - $upstream_dist,
                            $coding_region_start - 1 );
      $seq_to_check = $slice_to_check->seq;

      print "  no utr. Seq on forward strand transcript is: ", $seq_to_check, "\n";
    } elsif ($trans->strand == -1) {
      $slice_to_check =
      $sa->fetch_by_region( 'toplevel', $chromosome,
                            $coding_region_start + 1,
                            $coding_region_start + $upstream_dist );
      $seq_to_check = $slice_to_check->invert->seq();

      print "  no utr. Seq on reverse strand transcript is: ", $seq_to_check, "\n";
    }

  }

  my $i = 3;
  while ( $i <= $upstream_dist ) {
    my $startseq = substr( $seq_to_check, length($seq_to_check) - $i, 3 );
    print "  the codon is: ", $startseq, "\n";
    if ( $startseq eq "ATG" ) {
      my $alt_atg_start;
      print "  Found alternative ATG ";
      if ( $has_utr == 1 ) {
        my $cdna_alt_atg_start = ( $trans->cdna_coding_start - $i );
        # Transform the coordinates of the alternative codon into genomic coordinates
        $alt_atg_start = $cdna_alt_atg_start;
      } else {
        $alt_atg_start = "-$i";
      }
      print "at position $alt_atg_start\n";
      $alt_atgs{$alt_atg_start} = 1;
    } elsif ( $startseq eq "TAG" || $startseq eq "TAA" || $startseq eq "TGA" )
    {
      print "  Found stop\n";
      return %alt_atgs;
    }
    $i += 3;
  } ## end while ( $i <= $upstream_dist)
  return %alt_atgs;
} ## end sub check_transcript_forward



sub add_transcript_attrib {
  my ( $trans, $alt_ATG_start ) = @_;

  my $attribute =
    Bio::EnsEMBL::Attribute->new(
          -CODE        => 'upstream_ATG',
          -NAME        => 'upstream ATG',
          -DESCRIPTION => 'Alternative ATG found upstream of the defined ATG '
                        . 'used as start for the transcript',
          -VALUE => $alt_ATG_start );

  $trans->add_Attributes($attribute);

  return $trans;
}

sub make_attributes {
  my ($attrib_hash) = @_;
  my @attribs;

  foreach my $alt_ATG_start (keys %{$attrib_hash}) {
    my $attribute = Bio::EnsEMBL::Attribute->new(
          -CODE        => 'upstream_ATG',
          -NAME        => 'upstream ATG',
          -DESCRIPTION => 'Alternative ATG found upstream of the defined ATG '
                        . 'used as start for the transcript',
          -VALUE => $alt_ATG_start );

    push @attribs, $attribute;
  }
  return \@attribs;
}
