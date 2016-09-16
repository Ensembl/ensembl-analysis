#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

=head1 NAME

  copy_genes.pl

=head1 DESCRIPTION

This script takes database options and a file of gene ids or stable_ids
and copies them between two databases.  It can, if asked, split
multi-transcript genes into single genes.

The source database (from which you copy genes FROM) must contain DNA,
or else the script will throw.  If your source database does not contain
DNA, then provide details of a DNA database using the options --dnauser,
--dnahost, --dnadbname etc.

When using the in and out_config_name options it reads the equivalent
database from the Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases
file.

=head1 OPTIONS

  --all                 This option will copy all genes from
                        $sourcedbname to the output database.
  
  --merge               The source transcripts will be copied into existing genes in the output database
                        if they overlap at gene level if this parameter is true.
                        
  --skip_exon_sf        Prevent copying supporting features across exons increased speed for lost accuracy

  --split               This option will split multi-transcript genes
                        into single-transcript genes.

  --logic               This option will change the analysis of the
                        genes and transcripts being written to the output database to an
                        analysis with the specified logic_name.
                        
  --biotype             This option will change the biotype of the
                        genes and transcripts being written to the output database to the
                        biotype specified.
                        
  --source              This option will change the source of the
                        genes and transcripts being written to the output database to the
                        source specified.

  --remove_xrefs        Using this flag will remove xrefs from genes,
                        transcripts, and from translations.

  --remove_stable_ids   Using this flag will remove stable IDs from
                        genes, transcripts, translations, and from
                        exons.

  --transform_to        Transforms the genes from one coordinate system
                        to the given one.  The value must be on the
                        format "coord_system_name:version".

  --stable_id           Flag for indicating that input file contains
                        stable IDs and not gene IDs.

  --file                Read gene dbIDs out of a supplied file and only
                        copy these genes.

=head1 EXAMPLE

  perl copy_genes.pl \
    --source_config_name=CONSENSUS_DB \
    --target_config_name=UTR_DB \
    --file=gene_ids_to_copy

  or

  perl copy_genes.pl \
    --sourcehost=srchost --sourceuser=ensro \
    --sourcedbname=est_db \
    --targethost=trghost --targetuser=user --targetpass=XXX \
    --targetdbname=utr_db \
    --file=gene_ids_to_copy --stable_id


In the options, 'source' is synonymous with 'in' (--in_config_name is
the same option as --source_config_name etc.) and 'target' is synonymous
with 'out' (--targetdbname is the same as --outdbname).

=cut

use strict;
use warnings;

use Carp;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning verbose);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

my $sourcehost   = '';
my $sourceuser   = 'ensro';
my $sourcedbname = '';
my $sourcepass   = undef;
my $sourceport   = 3306;

my $outhost   = '';
my $outuser   = '';
my $outpass   = '';
my $outdbname = undef;
my $outport   = 3306;

my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnapass   = '';
my $dnadbname = undef;
my $dnaport   = 3306;

my $in_config_name;
my $out_config_name;

my $logic;
my $biotype;
my $source;
my $split = 0;
my $skip_exon_sf = 1;
my $infile;
my $all;
my $merge;
my $remove_xrefs;
my $remove_stable_ids;
my $transform_to;
my $stable_id;
my $verbose;

verbose('EXCEPTION');

GetOptions( 'inhost|sourcehost:s'                  => \$sourcehost,
            'inuser|sourceuser:s'                  => \$sourceuser,
            'indbname|sourcedbname:s'              => \$sourcedbname,
            'inport|sourceport:n'                  => \$sourceport,
            'inpass|sourcepass:s'                  => \$sourcepass,
            'in_config_name|source_config_name:s'  => \$in_config_name,
            'out_config_name|target_config_name:s' => \$out_config_name,
            'outhost|targethost:s'                 => \$outhost,
            'outuser|targetuser:s'                 => \$outuser,
            'outpass|targetpass:s'                 => \$outpass,
            'outdbname|targetdbname:s'             => \$outdbname,
            'outport|targetport:n'                 => \$outport,
            'dnahost:s'                            => \$dnahost,
            'dnauser:s'                            => \$dnauser,
            'dnadbname:s'                          => \$dnadbname,
            'dnaport:n'                            => \$dnaport,
            'logic:s'                              => \$logic,
            'biotype:s'                            => \$biotype,
            'source:s'                             => \$source,
            'split!'                               => \$split,
            'all!'                                 => \$all,
            'merge!'                               => \$merge,
            'skip_exon_sf!'                        => \$skip_exon_sf,
            'remove_xrefs!'                        => \$remove_xrefs,
            'remove_stable_ids!' => \$remove_stable_ids,
            'transform_to:s'     => \$transform_to,
            'verbose!'           => \$verbose,
            'stable_id!'         => \$stable_id,
            'file:s'             => \$infile ) ||
  throw("Error while parsing command line options");

if ($verbose) {
  verbose('WARNING');
}

my $transform_to_version;

my $attach_dna_db = 1 ;
if ($transform_to) {
  $attach_dna_db = 1;
  if ( $transform_to =~ m/:/ ) {
    ( $transform_to, $transform_to_version ) = split /:/, $transform_to;
  }
}

if ( $all && defined($infile) ) {
  throw("Specify either --all or --file, but not both");
}
elsif ( !$all && !defined($infile) ) {
  throw( "Specify either --all " .
         "(to copy all genes from $sourcedbname) " .
         "or --file (to copy a list of gene_ids)" );
}

#print "Connecting to ".$sourcedbname." at ".$sourcehost." ".$sourceuser."\n";

my $sourcedb;
my $dnadb;

if ( defined($in_config_name) ) {
  $sourcedb =
    get_db_adaptor_by_string( $in_config_name,
                              0, 0,
                              { -do_not_attach_dna_db => !$attach_dna_db
                              } );

  if ( $attach_dna_db && $dnadbname ) {
    $dnadb =
      new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                          -user   => $dnauser,
                                          -port   => $dnaport,
                                          -dbname => $dnadbname );

    $sourcedb->dnadb($dnadb);
  }
}
else {
  $sourcedb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $sourcehost,
                                        -user   => $sourceuser,
                                        -pass   => $sourcepass,
                                        -port   => $sourceport,
                                        -dbname => $sourcedbname );
  if ($attach_dna_db) {
    if ( !defined($dnadbname) ) {
      my $dna_query = q(SELECT count(1) FROM dna);
      my $dna_sth   = $sourcedb->dbc()->prepare($dna_query);
      $dna_sth->execute();

      my $dna_count = $dna_sth->fetchrow();

      if ( !$dna_count ) {
        croak( "\nYour source database does not contain DNA. " .
               "Please provide a database with genomic sequences.\n" );
      }
    }
    else {
      $dnadb =
        new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                            -user   => $dnauser,
                                            -port   => $dnaport,
                                            -dbname => $dnadbname );

      $sourcedb->dnadb($dnadb);
    }
  }
} ## end else [ if ( defined($in_config_name...))]

my $outdb;
if ($out_config_name) {
  $outdb = get_db_adaptor_by_string( $out_config_name, 0, 0,
                         { -do_not_attach_dna_db => !$attach_dna_db } );
}
else {
  $outdb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -pass   => $outpass,
                                        -port   => $outport,
                                        -dbname => $outdbname );
}

my $ga = $sourcedb->get_GeneAdaptor();

my @genes;
my @copy_genes;

if ( $remove_xrefs && $remove_stable_ids ) {
  print STDERR "Fetching genes, removing xrefs and stable ids. " .
    "(Adaptors and dbIDs also removed.)\n";
}
elsif ($remove_xrefs) {
  print STDERR "Fetching genes, removing xrefs, keeping stable ids. " .
    "(Adaptors and dbIDs also removed.)\n";
}
elsif ($remove_stable_ids) {
  print STDERR "Fetching genes, removing stable ids, keeping xrefs. " .
    "(Adaptors and dbIDs also removed.)\n";
}
else {
  print STDERR "Fetching genes, keeping xrefs and stable IDs. " .
    "(Adaptors and dbIDs removed.)\n";
}

my @gene_ids;
if ( defined($infile) ) {
  open( IN, $infile ) or die("Can not open '$infile' for reading: $!");

  while ( my $line = <IN> ) {
    chomp($line);
    push( @gene_ids, $line );
  }

  close(IN);
}
else {
  @gene_ids = @{ $ga->list_dbIDs() };
}

printf( "Got %d gene IDs\n", scalar(@gene_ids) );

my $outga = $outdb->get_GeneAdaptor();
my $analysis;
if ( defined($logic) ) {
  $analysis = $outdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
  if ( !defined($analysis) ) {
    $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $logic, );
  }
}

my $genes_processed = 0;
my $genes_stored    = 0;
my $batch_size      = 100;
my $genes_to_go     = scalar(@gene_ids);

while (@gene_ids) {
  my @gene_id_batch;
  while ( @gene_ids && scalar(@gene_id_batch) < $batch_size ) {
    my $gene_id = shift(@gene_ids);
    if ( defined($gene_id) ) {
      push( @gene_id_batch, $gene_id );
    }
  }

  my @genes;

  {
    my @gene_batch;
    if ($stable_id) {
      @gene_batch =
        @{ $ga->fetch_all_by_stable_id_list( \@gene_id_batch ) };
    }
    else {
      @gene_batch = @{ $ga->fetch_all_by_dbID_list( \@gene_id_batch ) };
    }

    if ($split) {
      foreach my $gene (@gene_batch) {
        push( @genes,
              @{convert_to_genes( $gene->get_all_Transcripts(),
                                  $gene->analysis(),
                                  $gene->biotype() ) } );
      }
    }
    else {
      @genes = @gene_batch;
    }
  }

  foreach my $gene (@genes) {
    my $old_stable_id = $gene->stable_id();

    $gene->load();    # fully_load_Gene($gene);

    if ( defined($logic) ) {
      $gene->analysis($analysis);
      foreach my $t ( @{ $gene->get_all_Transcripts() } ) {
        $t->analysis($analysis);
      }
    }
    
    if ( defined($biotype) ) {
      $gene->biotype($biotype);
      foreach my $t ( @{ $gene->get_all_Transcripts() } ) {
        $t->biotype($biotype);
      }
    }
    
    if ( defined($source) ) {
      $gene->source($source);
      foreach my $t ( @{ $gene->get_all_Transcripts() } ) {
        $t->source($source);
      }
    }    

    if ($transform_to) {
      if ($verbose) {
        printf( "Transforming '%s'\n", $old_stable_id );
      }

      my $transformed_gene =
        $gene->transform( $transform_to, $transform_to_version );

      # Only check transform if transform is successful.
      if ( defined($transformed_gene) ) {
        check_transform( $gene, $transformed_gene,
                         $transform_to_version );
      }
      $gene = $transformed_gene;
    }

    if ( defined($gene) ) {
      if (defined($merge)) {
        # check if there is any other gene on that region
        my $outdb_slice = get_feature_slice_from_db($gene,$outga->db());
        my @outdb_genes = @{$outga->fetch_all_by_Slice($outdb_slice,undef,1)};
        if (scalar(@outdb_genes) > 0) {
          # choose target gene to copy the source gene transcripts into
          my @sorted_outdb_genes = sort {$b->length() <=> $a->length()} @outdb_genes; # pick the longest
          my $target_gene = shift(@sorted_outdb_genes);

          $target_gene->load();
          $outga->remove($target_gene);

          # copy the source gene transcripts into the chosen target gene
          foreach my $source_transcript (@{$gene->get_all_Transcripts()}) {
            $source_transcript->load();

            my $new_source_transcript = $source_transcript->transfer($target_gene->slice());
            $new_source_transcript->dbID(0); # to make sure we don't delete any transcript of the output db whose id matches the id of the source db transcript
            $target_gene->add_Transcript($new_source_transcript);
          }
          empty_Gene( $target_gene, 0, 0 );
          $outga->store($target_gene,1,0,$skip_exon_sf);
        } else {
          $outga->store($gene,1,0,$skip_exon_sf);
        }

      } else {
        empty_Gene( $gene, $remove_stable_ids, $remove_xrefs );
        $outga->store($gene,1,0,$skip_exon_sf);
      }
      ++$genes_stored;
      --$genes_to_go;

      if ($verbose) {
        printf( "Stored gene '%s'.  Done %d, %d left to go..\n",
                $old_stable_id, $genes_processed + 1, $genes_to_go );
      }
    }
    else {
      printf( STDERR "Gene '%s' did not transform\n", $old_stable_id );
    }

    ++$genes_processed;

  } ## end foreach my $gene (@genes)
} ## end while (@gene_ids)

printf( "All done.  Processed %d genes, stored %d genes.\n",
        $genes_processed, $genes_stored );

sub check_transform {
  my ( $old_gene, $new_gene, $new_assembly_version ) = @_;

  my $old_transcripts = $old_gene->get_all_Transcripts();
  my $new_transcripts = $new_gene->get_all_Transcripts();

  if ( scalar(@$old_transcripts) != scalar(@$new_transcripts) ) {
    print STDERR "TRANSFORM_CHECK: old gene " .
      $old_gene->stable_id() . " has " .
      ( scalar(@$old_transcripts) ) . " transcripts but new gene has " .
      ( scalar(@$new_transcripts) ) . " transcripts\n";
  }

  # loop through and compare
  foreach my $old_transc ( @{$old_transcripts} ) {
    foreach my $new_transc ( @{$new_transcripts} ) {
      if ( $old_transc->stable_id() ne $new_transc->stable_id() ) {
        next;
      }

      # check number of exons, not sure if this helps
      if ( scalar( @{ $old_transc->get_all_Exons() } ) !=
           scalar( @{ $new_transc->get_all_Exons() } ) )
      {
        print STDERR "TRANSFORM_CHECK: old transcript " .
          $old_transc->stable_id() .
          " has " . ( scalar( @{ $old_transc->get_all_Exons() } ) ) .
          " exons but new transcript has " .
          ( scalar( @{ $new_transc->get_all_Exons() } ) ) . " exons\n";
      }

      # check translation
      if ( defined( $old_transc->translation() ) &&
           !defined($new_assembly_version) )
      {
        # We don't want to do have to deal with this if transforming
        # between assembly _VERSIONS_.
        my $new_translation = $new_transc->translate->seq();
        my $old_translation = $old_transc->translate->seq();

        if ( !defined( $new_transc->translation() ) ||
             $old_translation ne $new_translation )
        {
          print "TRANSFORM_CHECK: " .
            "old translation does not match new translation\n" .
            ">old_" . $old_transc->stable_id() .
            "\n" . $old_translation . "\n" . ">new_" .
            $new_transc->stable_id() . "\n" . $new_translation . "\n";
        }
      }
    } ## end foreach my $new_transc ( @{...})
  } ## end foreach my $old_transc ( @{...})

  return 1;
} ## end sub check_transform

sub get_feature_slice_from_db {
  my ( $feature, $db ) = @_;

  # This little helper routine returns a feature slice for a particular
  # region.  The slice will be associated with the given database.

  my $slice = $feature->feature_Slice();

  my @slices = @{
    $db->get_SliceAdaptor()->fetch_by_region_unique(
         $slice->coord_system_name(), $slice->seq_region_name(),
         $slice->start(),             $slice->end(),
         1,            $slice->coord_system()->version(),
         1 ) };

  if ( scalar(@slices) != 1 ) {
    # This will hopefully only happen if the Primary and Secondary
    # databases contain different assemblies.
    die( "!! Problem with projection for feature slice %s\n",
         $slice->name() );
  }

  return $slices[0];
}
