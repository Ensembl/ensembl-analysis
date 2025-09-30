#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene get_readthroughs_count);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(exon_overlap features_overlap overlap_length remove_short_frameshift_introns print_Transcript_and_Exons);

my ( $opt_host_secondary, $opt_port_secondary,
     $opt_user_secondary, $opt_password_secondary,
     $opt_database_secondary );
my ( $opt_host_primary,     $opt_port_primary, $opt_user_primary,
     $opt_password_primary, $opt_database_primary );
my ( $opt_host_dna,     $opt_port_dna, $opt_user_dna,
     $opt_password_dna, $opt_database_dna );
my ( $opt_host_output,     $opt_port_output, $opt_user_output,
     $opt_password_output, $opt_database_output );

my ( @opt_primary_include,  @opt_primary_exclude );
my ( @opt_secondary_include, @opt_secondary_exclude );

my $opt_primary_tag  = 'primary';
my $opt_secondary_tag = 'secondary';

my $opt_primary_xref = 1; # if 0, no primary xrefs are added; if 1, primary xrefs are added.
                          # Note that the xref attributes will always be added in any case.
my $opt_primary_gene_xref        = 'OTTG,Havana gene,ALT_GENE';
my $opt_primary_transcript_xref  = 'OTTT,Havana transcript,ALT_TRANS';
my $opt_primary_translation_xref = 'OTTP,Havana translation,MISC';

my $opt_njobs = 1;    # Default number of jobs.
my $opt_job   = 1;    # This job.
my $opt_file;   # list of genes from the primary db to be merged

my $opt_help = 0;

my %stored_genes; # to keep track of already stored genes and avoid duplicates

# make stdout unbuffered to ensure "PROCESSED + identifier" is not truncated when printfed
$| = 1;

if ( !GetOptions(
          'host_secondary:s'                    => \$opt_host_secondary,
          'port_secondary:i'                    => \$opt_port_secondary,
          'user_secondary:s'                    => \$opt_user_secondary,
          'password_secondary|pass_secondary:s'   => \$opt_password_secondary,
          'database_secondary|dbname_secondary:s' => \$opt_database_secondary,
          'host_primary:s'                     => \$opt_host_primary,
          'port_primary:i'                     => \$opt_port_primary,
          'user_primary:s'                     => \$opt_user_primary,
          'password_primary|pass_primary:s'     => \$opt_password_primary,
          'database_primary|dbname_primary:s'   => \$opt_database_primary,
          'host_dna:s'                        => \$opt_host_dna,
          'port_dna:i'                        => \$opt_port_dna,
          'user_dna:s'                        => \$opt_user_dna,
          'password_dna|pass_dna:s'           => \$opt_password_dna,
          'database_dna|dbname_dna:s'         => \$opt_database_dna,
          'host_output:s'                     => \$opt_host_output,
          'port_output:i'                     => \$opt_port_output,
          'user_output:s'                     => \$opt_user_output,
          'password_output|pass_output:s'     => \$opt_password_output,
          'database_output|dbname_output:s'   => \$opt_database_output,
          'secondary_include:s'                 => \@opt_secondary_include,
          'secondary_exclude:s'                 => \@opt_secondary_exclude,
          'primary_include:s'                  => \@opt_primary_include,
          'primary_exclude:s'                  => \@opt_primary_exclude,
          'primary_tag:s'                      => \$opt_primary_tag,
          'secondary_tag:s'                     => \$opt_secondary_tag,
          'primary_xref:i'                     => \$opt_primary_xref,
          'primary_gene_xref:s'                => \$opt_primary_gene_xref,
          'primary_transcript_xref:s'  => \$opt_primary_transcript_xref,
          'primary_translation_xref:s' => \$opt_primary_translation_xref,
          'njobs:i'                   => \$opt_njobs,
          'job:i'                     => \$opt_job,
          'file:s'                    => \$opt_file,
          'help|h|?!'                 => \$opt_help, ) ||
     $opt_help ||
     !( defined($opt_host_secondary) &&
        defined($opt_user_secondary) &&
        defined($opt_database_secondary) )
     ||
     !( defined($opt_host_primary) &&
        defined($opt_user_primary) &&
        defined($opt_database_primary) )
     ||
     !( defined($opt_host_output) &&
        defined($opt_user_output) &&
        defined($opt_database_output) ) ||
     !( $opt_njobs >= 1 && $opt_job >= 1 && $opt_job <= $opt_njobs ) )
{

  if ($opt_help) {
    pod2usage( -verbose => 2, -exitval => 0 );
  }
  else {
    pod2usage( -verbose => 0, -exitval => 'NOEXIT' );

    if ( !( defined($opt_host_secondary) &&
            defined($opt_user_secondary) &&
            defined($opt_database_secondary) ) )
    {
      die( 'Need connection parameters for Secondary database ' .
           '(host_secondary, user_secondary and database_secondary)' );
    }
    elsif ( !( defined($opt_host_primary) &&
               defined($opt_user_primary) &&
               defined($opt_database_primary) ) )
    {
      die( 'Need connection parameters for Primary database ' .
           '(host_primary, user_primary and database_primary)' );
    }
    elsif ( !( defined($opt_host_output) &&
               defined($opt_user_output) &&
               defined($opt_database_output) ) )
    {
      die( 'Need connection parameters for output database ' .
           '(host_output, user_output and database_output)' );
    }
    elsif (
       !( $opt_njobs >= 1 && $opt_job >= 1 && $opt_job <= $opt_njobs ) )
    {
      die( 'Number of jobs must be 1 or greater, ' .
           'and the current job needs to be ' .
           'between 1 and the number of jobs' );
    }
    elsif ( ( @opt_secondary_include && @opt_secondary_exclude ) ||
            ( @opt_primary_include && @opt_primary_exclude ) )
    {
      die('You may only use X_include or X_exclude, but not both');
    }
    else {
      die('Error in command line parsing');
    }
  } ## end else [ if ($opt_help) ]

} ## end if ( !GetOptions( 'host_secondary:s'...))

my $dna_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                '-no_cache' => 1,
                '-host'     => $opt_host_dna || $opt_host_secondary,
                '-port'     => $opt_port_dna || $opt_port_secondary,
                '-user'     => $opt_user_dna || $opt_user_secondary,
                '-pass' => $opt_password_dna || $opt_password_secondary,
                '-dbname' => $opt_database_dna || $opt_database_secondary,
  ) or
  die('Failed to connect to DNA database');

my $secondary_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                     '-no_cache' => 1,
                                     '-host'     => $opt_host_secondary,
                                     '-port'     => $opt_port_secondary,
                                     '-user'     => $opt_user_secondary,
                                     '-pass'   => $opt_password_secondary,
                                     '-dbname' => $opt_database_secondary,
                                     '-dnadb'  => $dna_dba, ) or
  die('Failed to connect to Secondary database');

my $primary_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new('-no_cache' => 1,
                                      '-host'     => $opt_host_primary,
                                      '-port'     => $opt_port_primary,
                                      '-user'     => $opt_user_primary,
                                      '-pass'   => $opt_password_primary,
                                      '-dbname' => $opt_database_primary,
                                      '-dnadb'  => $dna_dba, ) or
  die('Failed to connect to Primary database');

my $output_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new('-no_cache' => 1,
                                      '-host'     => $opt_host_output,
                                      '-port'     => $opt_port_output,
                                      '-user'     => $opt_user_output,
                                      '-pass'   => $opt_password_output,
                                      '-dbname' => $opt_database_output,
  ) or
  die('Failed to connect to output database');

@opt_secondary_include = split( /,/, join( ',', @opt_secondary_include ) );
@opt_secondary_exclude = split( /,/, join( ',', @opt_secondary_exclude ) );
@opt_primary_include  = split( /,/, join( ',', @opt_primary_include ) );
@opt_primary_exclude  = split( /,/, join( ',', @opt_primary_exclude ) );


if($opt_database_dna) {

  print "Optional DNA database\thost:\t".$opt_host_dna."\n".
                                   "\tport:\t".$opt_port_dna."\n".
                                   "\tuser:\t".$opt_user_dna."\n".
                                   "\tname:\t".$opt_database_dna."\n";
}


print <<DBINFO_END;
SECONDARY database\thost:\t$opt_host_secondary
                \tport:\t$opt_port_secondary
                \tuser:\t$opt_user_secondary
                \tname:\t$opt_database_secondary

PRIMARY database\thost:\t$opt_host_primary
               \tport:\t$opt_port_primary
               \tuser:\t$opt_user_primary
               \tname:\t$opt_database_primary

OUTPUT database\thost:\t$opt_host_output
               \tport:\t$opt_port_output
               \tuser:\t$opt_user_output
               \tname:\t$opt_database_output

Secondary logic name filter (include): @opt_secondary_include
Secondary logic name filter (exclude): @opt_secondary_exclude

Primary logic name filter (include): @opt_primary_include
Primary logic name filter (exclude): @opt_primary_exclude

Secondary tag:\t$opt_secondary_tag
Primary tag:\t$opt_primary_tag

Primary xref:       \t$opt_primary_xref
Primary gene xref:       \t$opt_primary_gene_xref
Primary transcript xref: \t$opt_primary_transcript_xref
Primary translation xref:\t$opt_primary_translation_xref

------------------------------------------------------------------------
DBINFO_END

my $SECONDARY_GA = $secondary_dba->get_GeneAdaptor();    # Used globally.
my $PRIMARY_GA  = $primary_dba->get_GeneAdaptor();     # Used globally.
my $OUTPUT_GA  = $output_dba->get_GeneAdaptor();     # This one too.
my $OUTPUT_AA  = $output_dba->get_AnalysisAdaptor();     # This one too.

# Create a chunk of work, i.e. a list of gene IDs.  We create uniformly
# sized chunks whose size depend on the number of available Primary
# genes and the number of jobs that are being run.  We pick the chunk
# associated with our job ID.

my @all_primary_gene_ids = ();
if ($opt_file) {
  # fetch all primary genes from the list in the file 
  open(GENEIDS,$opt_file);
  while (<GENEIDS>) {
    chomp;
    push(@all_primary_gene_ids,$_);
  }
  close(GENEIDS);
} else {
  # fetch all genes in the primary db
  @all_primary_gene_ids = sort { $a <=> $b } @{ $PRIMARY_GA->list_dbIDs() };
}

my $chunk_size  = int( scalar(@all_primary_gene_ids)/$opt_njobs );
my $chunk_first = ( $opt_job - 1 )*$chunk_size;
my $chunk_last  = $chunk_first + $chunk_size - 1;

printf( "Got %d genes, with %d jobs this means %d genes per job",
        scalar(@all_primary_gene_ids),
        $opt_njobs, $chunk_size );

# This ($chunk_spill) is just the number of jobs that will get one extra
# piece of work to do, because $chunk_size multiplied by $opt_njobs
# isn't quite equal to the number of Primary genes.

my $chunk_spill = scalar(@all_primary_gene_ids) % $opt_njobs;

if ( $chunk_spill > 0 ) {
  printf( " (the first %d jobs will get an extra gene).\n",
          $chunk_spill );

  if ( $opt_job <= $chunk_spill ) {
    $chunk_first += $opt_job - 1;
    $chunk_last  += $opt_job;
  }
  else {
    $chunk_first += $chunk_spill;
    $chunk_last  += $chunk_spill;
  }
}
else {
  print(".\n");
}

printf( "This is job %d. Will run genes indexed %d-%d.\n",
        $opt_job, $chunk_first, $chunk_last );

# Create the analyses for the merged genes and transcripts
my $merged_source = $opt_secondary_tag.'_'.$opt_primary_tag;
my $merged_gene_logic_name = $merged_source.'_gene';
my $merged_ig_gene_logic_name = $merged_source.'_ig_gene';
my $merged_transcript_logic_name = $merged_source.'_transcript';
my $primary_ig_gene_logic_name = $opt_primary_tag.'_ig_gene';
my $secondary_ig_gene_logic_name = $opt_secondary_tag.'_ig_gene';
my $secondary_lincrna_logic_name = $opt_secondary_tag.'_lincrna';
my $merged_lincrna_logic_name = $merged_source.'_lincrna';

my @logic_names = ($opt_secondary_tag,
                   $opt_primary_tag,
                   $merged_gene_logic_name,
                   $merged_ig_gene_logic_name,
                   $merged_transcript_logic_name,
                   $primary_ig_gene_logic_name,
                   $secondary_ig_gene_logic_name,
                   $secondary_lincrna_logic_name,
                   $merged_lincrna_logic_name);
                
foreach my $analysis_name (@logic_names) {
  if (!($OUTPUT_AA->fetch_by_logic_name($analysis_name))) {
    $OUTPUT_AA->store(new Bio::EnsEMBL::Analysis(-logic_name => $analysis_name));
  }
}

my %primary_genes_done;

# Process all the Primary genes in the current chunk.
foreach my $primary_gene_id (
                   @all_primary_gene_ids[ $chunk_first .. $chunk_last ] )
{
  if ( exists( $primary_genes_done{$primary_gene_id} ) ) {
    printf( "Skipping gene %d, already processed " .
              "(because of clustering).\n",
            $primary_gene_id );
    next;
  }

  # Given the gene with ID $primary_gene_id, pick out all the overlapping
  # Primary and Secondary genes (will return nothing if the gene cluster is
  # found to contain genes belonging to another LSF job).

  my ( $primary_genes, $secondary_genes ) =
    @{
    make_gene_cluster( $primary_gene_id,
                       $all_primary_gene_ids[$chunk_first] ) };

  # Filter the genes.
  my @filtered_primary_genes;
  my @filtered_secondary_genes;
  foreach my $primary_gene ( @{$primary_genes} ) {
    if ( filter_gene( $primary_gene, \@opt_primary_include,
                      \@opt_primary_exclude ) )
    {
      next;
    }
    push( @filtered_primary_genes, $primary_gene );
  }
  foreach my $secondary_gene ( @{$secondary_genes} ) {
    if ( filter_gene( $secondary_gene, \@opt_secondary_include,
                      \@opt_secondary_exclude ) )
    {
      next;
    }
    push( @filtered_secondary_genes, $secondary_gene );
  }

  # Process the filtered genes.

  # This loop needs to come before the call to process_genes() since
  # the dbIDs are changing when storing them.  It should loop over
  # @{$primary_genes}, not @filtered_primary_genes.
  foreach my $primary_gene ( @{$primary_genes} ) {
    $primary_genes_done{ $primary_gene->dbID() } = 1;
  }

  process_genes( \@filtered_primary_genes, \@filtered_secondary_genes );

  print("==\n");
} ## end foreach my $primary_gene_id ...

sub filter_gene {
  my ( $gene, $include_logic_names, $exclude_logic_names ) = @_;

  my $do_filter = 0;

  if ( @{$include_logic_names} ) {
    my $gene_logic_name = $gene->analysis()->logic_name();
    my $do_include      = 0;

    foreach my $logic_name ( @{$include_logic_names} ) {
      if ( $gene_logic_name eq $logic_name ) {
        $do_include = 1;
        last;
      }
    }

    if ( !$do_include ) {
      printf( "Not using gene %s (%d), " .
                "logic name '%s' is not included.\n",
              $gene->stable_id(), $gene->dbID(), $gene_logic_name );
    }

    $do_filter = !$do_include;
  }
  elsif ( @{$exclude_logic_names} ) {
    my $gene_logic_name = $gene->analysis()->logic_name();
    my $do_exclude      = 0;

    foreach my $logic_name ( @{$exclude_logic_names} ) {
      if ( $gene_logic_name eq $logic_name ) {
        $do_exclude = 1;
        last;
      }
    }

    if ($do_exclude) {
      printf( "Not using gene %s (%d), " .
                "logic name '%s' is excluded.\n",
              $gene->stable_id(), $gene->dbID(), $gene_logic_name );
    }

    $do_filter = $do_exclude;
  }

  return $do_filter;
} ## end sub filter_gene

sub make_gene_cluster {
  my ( $seed_gene_id, $lowest_allowed_gene_id ) = @_;

  # See if the seed gene overlaps with any other gene.  If it does,
  # repeat until we have all overlapping genes in the cluster.  If
  # the lowest Primary gene ID in the cluster is greater or equal to
  # $lowest_allowed_gene_id, then process it, otherwise skip this
  # cluster as it belongs to another LSF job.

  my $seed_gene = $PRIMARY_GA->fetch_by_dbID($seed_gene_id);

  printf( "Initiating gene cluster with Primary gene %s (%d), %s\n",
          $seed_gene->stable_id(),
          $seed_gene_id, $seed_gene->feature_Slice()->name() );

  my @cluster_queue = ($seed_gene);
  my %used_slices;

  my %primary_cluster;
  my %secondary_cluster;

  my $primary_sa  = $PRIMARY_GA->db()->get_SliceAdaptor();
  my $secondary_sa = $SECONDARY_GA->db()->get_SliceAdaptor();

  my @strands = (-1,1);

  # The cluster queue is a collection of possibly unprocessed genes in
  # this cluster of overlapping genes.  We use each gene in turn to
  # fetch overlapping gene from the Primary and the Secondary databases.
  # Any found gene, if it hasn't already been seen, is added to the
  # cluster queue.  When thu queue is empty, we are done.

  while (my $gene = shift(@cluster_queue)) {
  	
  	foreach my $current_strand (@strands) {
  	
      my $primary_slice = get_feature_slice_from_db($gene,$PRIMARY_GA->db());
      my $primary_slice_id = "primary:".$primary_slice->name()."_".$current_strand;
      if (exists($used_slices{$primary_slice_id})) {
        printf("Skipping Primary slice %s, already seen.\n",$primary_slice_id);
        next;
      }

      $used_slices{$primary_slice_id} = 1;

      # Fetch overlapping Primary genes.

      foreach my $primary_gene (
          @{ $PRIMARY_GA->fetch_all_by_Slice( $primary_slice, undef, 1 ) } )
      {
      	# all genes on the slice are fetched and we only want to look at the genes on the same strand
        if ($primary_gene->strand() == $current_strand) {
          my $add_to_result  = 1;
          my $primary_gene_id = $primary_gene->dbID();

          if ( $primary_gene->length() > 100_000_000 ) {
            printf( "Ignoring Primary gene %s (%d)," .
                    " too long (%d > 100,000,000)\n",
                    $primary_gene->stable_id(),
                    $primary_gene_id, $primary_gene->length() );
            next;
          }

          if ( !exists( $primary_cluster{$primary_gene_id} ) ) {
            if ( $primary_gene_id < $lowest_allowed_gene_id ) {
              print("Skipping gene cluster, does not belong to me.\n");
              return [ [], [] ];
            }

            push( @cluster_queue, $primary_gene );

            if ($add_to_result) {
              $primary_cluster{$primary_gene_id} = $primary_gene;

              printf( "Also fetched Primary gene %s (%d), %s\n",
                      $primary_gene->stable_id(),
                      $primary_gene_id,
                      $primary_gene->feature_Slice()->name() );
            }
          }
        } # end if primary gene strand
      } ## end foreach my $primary_gene ( @...)

      my $secondary_slice = get_feature_slice_from_db($gene,$SECONDARY_GA->db());
      my $secondary_slice_id = "secondary:".$secondary_slice->name()."_".$current_strand;
      if (exists($used_slices{$secondary_slice_id})) {
        printf("Skipping Secondary slice %s, already seen.\n",$secondary_slice_id);
        next;
      }

      $used_slices{$secondary_slice_id} = 1;

      # Fetch overlapping Secondary genes.

      foreach my $secondary_gene (
        @{ $SECONDARY_GA->fetch_all_by_Slice( $secondary_slice, undef, 1 ) } )
      {
      	# all genes on the slice are fetched and we only want to look at the genes on the same strand
      	if ($secondary_gene->strand() == $current_strand) {
      	
          my $add_to_result   = 1;
          my $secondary_gene_id = $secondary_gene->dbID();

          if ( $secondary_gene->length() > 100_000_000 ) {
            printf( "Ignoring Secondary gene %s (%d), " .
                    " too long (%d > 100,000,000)\n",
                    $secondary_gene->stable_id(),
                    $secondary_gene_id, $secondary_gene->length() );
            next;
          }

          if ( !exists( $secondary_cluster{$secondary_gene_id} ) ) {
            push( @cluster_queue, $secondary_gene );

            if ($add_to_result) {
              $secondary_cluster{$secondary_gene_id} = $secondary_gene;

              printf( "Also fetched Secondary gene %s (%d), %s\n",
                      $secondary_gene->stable_id(),
                      $secondary_gene_id,
                      $secondary_gene->feature_Slice()->name() );
            }
          }
      	} # end if secondary gene strand
      } ## end foreach my $secondary_gene ( ...)
  	} ## end foreach my $current_strand
  } ## end while ( my $gene = shift(...))

  return [ [ values(%primary_cluster) ], [ values(%secondary_cluster) ] ];
} ## end sub make_gene_cluster

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

sub process_genes {
  my ( $primary_genes, $secondary_genes ) = @_;

  # Start by fully loading each Primary gene.  Also add OTTP, OTTT and
  # OTTG xrefs to the various parts of the gene.  We also figure out
  # (and keep track of) whether the translations are coding, if the
  # genes are pseudogenes, if a gene is a single transcript gene, and
  # whether a gene is located across a region with a known assembly
  # error.

  foreach my $primary_gene ( @{$primary_genes} ) {
    $primary_gene->load();

    my $is_coding = 0;

    my $transcript_count = 0;
    foreach
      my $primary_transcript ( @{ $primary_gene->get_all_Transcripts() } )
    {
      my $primary_translation = $primary_transcript->translation();

      if ( defined($primary_translation) ) {
        $is_coding = 1;

        # Add "OTTP" xref to Primary translation.
        if ($opt_primary_xref) {
          add_primary_xref($primary_translation);
        }
        
        # store xref as an attribute
        my $attrib = Bio::EnsEMBL::Attribute->new(
                         -NAME        => 'Xref ID',
                         -CODE        => 'xref_id',
                         -VALUE       => $primary_translation->stable_id(),
                         -DESCRIPTION => 'ID of associated database reference'
        );
        $primary_translation->add_Attributes($attrib);
      }
      tag_transcript_analysis( $primary_transcript, $opt_primary_tag );
      $primary_transcript->source($opt_primary_tag);
      $primary_transcript->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($primary_transcript)));

      # Add "OTTT" xref to Primary transcript.
      if ($opt_primary_xref) {
        add_primary_xref($primary_transcript);
      }  
      
      # store xref as an attribute
      my $attrib = Bio::EnsEMBL::Attribute->new(
                       -NAME        => 'Xref ID',
                       -CODE        => 'xref_id',
                       -VALUE       => $primary_transcript->stable_id(),
                       -DESCRIPTION => 'ID of associated database reference'
      );
      $primary_transcript->add_Attributes($attrib);

      # If the transcript is part of a gene cluster then tag the gene
      unless($primary_gene->{__is_gene_cluster}) {
        $primary_gene->{__is_gene_cluster} = (
        scalar(@{ $primary_transcript->get_all_Attributes('gene_cluster') }) > 0 );
      }
      
      # If the transcript is labelled as "genome patch truncated" (attrib code='remark', value='genome patch truncated') then tag the transcript
      foreach my $primary_transcript_attrib (@{$primary_transcript->get_all_Attributes('remark')}) {
      	if ($primary_transcript_attrib->value() eq "genome patch truncated") {
      	  $primary_transcript->{__is_genome_patch_truncated} = 1;
      	}
      }
      ++$transcript_count;
    }

    $primary_gene->source($opt_primary_tag);
    $primary_gene->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($primary_gene)));

    my $has_assembly_error = (
              scalar(
                @{ $primary_gene->get_all_Attributes('NoTransRefError') }
              ) > 0 );

    $primary_gene->{__is_coding} = $is_coding;    # HACK
    $primary_gene->{__is_polymorphic_pseudogene} = $primary_gene->biotype() =~ /polymorphic_pseudogene/;

    $primary_gene->{__is_pseudogene} =            # HACK
      ( !$primary_gene->{__is_coding} &&
        $primary_gene->biotype() =~ /pseudogene/ );

    $primary_gene->{__has_ref_error} = $has_assembly_error;    # HACK
    $primary_gene->{__is_single_transcript} =
      ( $transcript_count == 1 );                             # HACK

    # Add "OTTG" xref to Primary gene.
    if ($opt_primary_xref) {
      add_primary_xref($primary_gene);
    }
    
    # store xref as an attribute
    my $attrib = Bio::EnsEMBL::Attribute->new(
                     -NAME        => 'Xref ID',
                     -CODE        => 'xref_id',
                     -VALUE       => $primary_gene->stable_id(),
                     -DESCRIPTION => 'ID of associated database reference'
    );
    $primary_gene->add_Attributes($attrib);

  } ## end foreach my $primary_gene ( @...)

  # For Secondary genes, we don't add any xrefs, but we figure out if the
  # transcripts are pseudogene transcripts or part of the CCDS dataset
  # (if made avalable).  We also see if a gene is a single transcript
  # gene, if it's coding and if it's an RNA gene or not.  All these
  # things are being used elsewhere to make decisions about merging and
  # copying.

  foreach my $secondary_gene ( @{$secondary_genes} ) {
    $secondary_gene->load();

    my $is_coding        = 0;
    my $is_pseudogene    = 0;
    my $transcript_count = 0;

    my $is_rna = (
         index( lc( $secondary_gene->analysis()->logic_name() ), 'ncrna' )
           >= 0 );

    foreach my $secondary_transcript (
                             @{ $secondary_gene->get_all_Transcripts() } )
    {
      $secondary_transcript->{__is_ccds}       = 0;    # HACK
      $secondary_transcript->{__is_pseudogene} = 0;    # HACK
      $secondary_transcript->{__is_rna} = (
         index( lc( $secondary_transcript->analysis()->logic_name() ), 'ncrna' )
           >= 0 );    # HACK

      my $secondary_translation = $secondary_transcript->translation();

      if ( defined($secondary_translation) ) {
        $secondary_transcript->{__is_ccds} = scalar(@{$secondary_transcript->get_all_Attributes('ccds_transcript')});    # HACK
        $is_coding = 1;
      }
      elsif ( $secondary_transcript->biotype() =~ /pseudogene/ ) {
        $secondary_transcript->{__is_pseudogene} = 1;    # HACK
      }

      $secondary_transcript->source($opt_secondary_tag);

      $secondary_transcript->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($secondary_transcript)));

      ++$transcript_count;
    } ## end foreach my $secondary_transcript...

    add_logic_name_suffix( $secondary_gene, $opt_secondary_tag );
    $secondary_gene->source($opt_secondary_tag);

    #my $is_rna = (
    #     index( lc( $secondary_gene->analysis()->logic_name() ), 'ncrna' )
    #       >= 0 );

    $secondary_gene->{__is_coding} = $is_coding;    # HACK
    $secondary_gene->{__is_single_transcript} =
      ( $transcript_count == 1 );                 # HACK
    $secondary_gene->{__is_rna} = $is_rna;          # HACK
    
    $secondary_gene->{__is_pseudogene} =            # HACK
      ( !$secondary_gene->{__is_coding} &&
        $secondary_gene->biotype() =~ /pseudogene/ );
        
    # Check pseudogene has one transcript  
    if ($secondary_gene->{__is_pseudogene}) {

      warning ($transcript_count." transcripts found for Secondary pseudogene: ".$secondary_gene->stable_id().
      ". Secondary pseudogenes should have a single transcript.")
      unless $transcript_count == 1;

    }
    
  } ## end foreach my $secondary_gene ( ...)

  my @transcripts_to_merge;
  my @transcripts_to_copy;
  my @transcripts_to_ignore;

  # This will keep a record of the Secondary gene stable id whenever stuff is
  # merged or copied in
  my %processed_ens_id_tracker;

SECONDARY_GENE:
  foreach my $secondary_gene ( @{$secondary_genes} ) {
    printf( "Secondary gene: %s (%d)\n",
            $secondary_gene->stable_id(), $secondary_gene->dbID() );

    my @secondary_transcripts = @{ $secondary_gene->get_all_Transcripts() };

  SECONDARY_TRANSCRIPT:
    foreach my $secondary_transcript (@secondary_transcripts) {
      printf( "\tSecondary transcript %s (%d)\n",
              $secondary_transcript->stable_id(),
              $secondary_transcript->dbID() );

    PRIMARY_GENE:
      foreach my $primary_gene ( @{$primary_genes} ) {
        printf( "\t\tPrimary gene: %s (%d)\n",
                $primary_gene->stable_id(), $primary_gene->dbID() );

        if ( $primary_gene->strand() != $secondary_gene->strand() ) {
          printf("\t\t\tNot on same strand\n");
          next PRIMARY_GENE;
        }

        my @primary_transcripts =
          @{ $primary_gene->get_all_Transcripts() };

      PRIMARY_TRANSCRIPT:
        foreach my $primary_transcript (@primary_transcripts) {
        	
          # If the primary_transcript is tagged as "genome patch truncated",
          # it should be treated as secondary and we'll try a partial merge
          my $backup_secondary_transcript = $secondary_transcript;
          my $backup_secondary_gene = $secondary_gene;
          my $backup_primary_gene = $primary_gene;
          
          if ($primary_transcript->{__is_genome_patch_truncated}) {
          	# GENOME PATCH TRUNCATED CASES
          	printf("\t\t\tSwapping Primary and Secondary transcript and trying partial merge due to 'genome patch truncated' attribute in Primary transcript %s (%d)\n",
                  $primary_transcript->stable_id(),
                  $primary_transcript->dbID() );
          	
          	$secondary_transcript = $primary_transcript;
          	$primary_transcript = $backup_secondary_transcript;
          	$secondary_gene = $primary_gene;
          	$primary_gene = $backup_secondary_gene;
          	
          	if (investigate_for_partial_merge($primary_transcript,$secondary_transcript)) {
              print("\t\t\tTo be partially merged\n");

              push( @transcripts_to_merge,
                    [ $primary_gene,  $primary_transcript,
                      $secondary_gene, $secondary_transcript ] );
            } elsif ( investigate_for_copy(
                                        $primary_gene,  $primary_transcript,
                                        $secondary_gene, $secondary_transcript
                    ) )
            {
              print( "\t\t\tSecondary transcript overlaps " .
                     "and may be copied\n" );
    
              push( @transcripts_to_copy,
                    [ $primary_gene,  $primary_transcript,
                      $secondary_gene, $secondary_transcript ] );
            }
          } else {
          	# ALL OTHER CASES
	        printf( "\t\t\tPrimary transcript %s (%d)\n",
	                $primary_transcript->stable_id(),
	                $primary_transcript->dbID() );
	
	        if ( investigate_for_merge(
	                               $primary_transcript, $secondary_transcript
	             ) )
	        {
	          print("\t\t\tTo be merged\n");
	
              push( @transcripts_to_merge,
	                [ $primary_gene,  $primary_transcript,
	                  $secondary_gene, $secondary_transcript ] );
	        }
	        elsif ( investigate_for_copy(
	                                    $primary_gene,  $primary_transcript,
	                                    $secondary_gene, $secondary_transcript
	                ) )
	        {
	          print( "\t\t\tSecondary transcript overlaps " .
	                 "and may be copied\n" );
	
	          push( @transcripts_to_copy,
	                [ $primary_gene,  $primary_transcript,
	                  $secondary_gene, $secondary_transcript ] );
	        }
          }
          # restore secondary transcript and genes
          $secondary_transcript = $backup_secondary_transcript;
          $secondary_gene = $backup_secondary_gene;
          $primary_gene = $backup_primary_gene;

        } ## end PRIMARY_TRANSCRIPT: foreach my $primary_transcript...
      } ## end PRIMARY_GENE: foreach my $primary_gene ( @...)

    } ## end SECONDARY_TRANSCRIPT: foreach my $secondary_transcript...
  } ## end SECONDARY_GENE: foreach my $secondary_gene ( ...)

  my %merged_secondary_transcripts;
  my %copied_secondary_transcripts;
  my %ignored_secondary_transcripts;

  # Merge the transcripts that have been found to be mergable.  We allow
  # merging a Secondary transcript into multiple Primary transcripts.
  foreach my $tuple (@transcripts_to_merge) {
    my ( $primary_gene,  $primary_transcript,
         $secondary_gene, $secondary_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    printf(
          "Merging %s (%d) into %s (%d) / %s (%d)\n",
          $secondary_transcript->stable_id(), $secondary_transcript->dbID(),
          $primary_transcript->stable_id(),  $primary_transcript->dbID(),
          $primary_gene->stable_id(),        $primary_gene->dbID() );

    my $merge_code =
      merge( $primary_gene, $primary_transcript, $secondary_transcript );

    $merged_secondary_transcripts{ $secondary_transcript->dbID() } = 1;

    if ( $merge_code != 1 ) {
      printf( "PROCESSED\t%d\t%s\n",
              $secondary_gene->dbID(), $secondary_gene->stable_id() );

      # Keep a record in memory of this Secondary stable id
      $processed_ens_id_tracker{$primary_gene->dbID()} = [$secondary_gene->dbID(),$secondary_gene->stable_id()];
    }
  }

  # Find all Secondary transcripts that are "siblings" to any merged
  # Secondary transcript.  Tag these for possibly copying into the
  # corresponding Primary gene.
  foreach my $tuple (@transcripts_to_merge) {
    my ( $primary_gene,  $primary_transcript,
         $secondary_gene, $secondary_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    foreach my $transcript ( @{ $secondary_gene->get_all_Transcripts() } )
    {
      my $key = $transcript->dbID();

      if ( !exists( $merged_secondary_transcripts{$key} ) ) {
        printf( "May also copy %s (%d) into %s (%d) " .
                  "due to merging of %s (%d)\n",
                $transcript->stable_id(),
                $transcript->dbID(),
                $primary_gene->stable_id(),
                $primary_gene->dbID(),
                $secondary_transcript->stable_id(),
                $secondary_transcript->dbID() );

        push( @transcripts_to_copy,
              [ $primary_gene,  $primary_transcript,
                $secondary_gene, $transcript ] );
      }
    }
  } ## end foreach my $tuple (@transcripts_to_merge)

  # Take care of the cases where an Secondary transcript has been selected
  # to be copied into multiple Primary genes.
  my %distinct_transcripts_to_copy;
  foreach my $tuple (@transcripts_to_copy) {
    my ( $primary_gene,  $primary_transcript,
         $secondary_gene, $secondary_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    # skip readthrough genes (defined as g containing any t having a 'readthrough_tra' attribute)
    if (get_readthroughs_count($primary_gene) > 0) {
      print("Not copying ".$secondary_transcript->stable_id()." into ".$primary_gene->stable_id()." because it contains readthrough transcripts\n.");
      next;
    }

    my $key = $secondary_transcript->dbID();

    if ( exists( $merged_secondary_transcripts{$key} ) ) {
      next;
    }

    if ( exists( $distinct_transcripts_to_copy{$key} ) ) {
      my $is_conflicting = 0;

      if ( $primary_gene->dbID() !=
           $distinct_transcripts_to_copy{$key}[0]->dbID() )
      {
        $is_conflicting = 1;
        printf( "Secondary transcript %s (%d) " .
                  "selected for copy into both " .
                  "%s (%d) and %s (%d)\n",
                $secondary_transcript->stable_id(),
                $secondary_transcript->dbID(),
                $primary_gene->stable_id(),
                $primary_gene->dbID(),
                $distinct_transcripts_to_copy{$key}[0]->stable_id(),
                $distinct_transcripts_to_copy{$key}[0]->dbID() );
      }

      if ( $primary_transcript->dbID() !=
           $distinct_transcripts_to_copy{$key}[1]->dbID() )
      {
      	# select by score based on exonic and non-exonic overlap 
      	my $exon_overlap = exon_overlap($secondary_transcript,$primary_transcript);
      	my $exon_non_overlap = abs($secondary_transcript->length()-$primary_transcript->length());
        my $overlap_score = $exon_overlap-$exon_non_overlap;

        if ($overlap_score > $distinct_transcripts_to_copy{$key}[4]) {
          if ($is_conflicting) {
            printf( "Will choose to copy %s (%d) into %s (%d) " .
                      "due to larger exon overlap (%dbp > %dbp)\n",
                    $secondary_transcript->stable_id(),
                    $secondary_transcript->dbID(),
                    $primary_gene->stable_id(),
                    $primary_gene->dbID(),
                    $overlap_score,
                    $distinct_transcripts_to_copy{$key}[4] );
          }

          $distinct_transcripts_to_copy{$key} = [
                                     $primary_gene,  $primary_transcript,
                                     $secondary_gene, $secondary_transcript,
                                     $overlap_score ];
        }
        elsif ($is_conflicting) {
          printf( "Will choose to copy %s (%d) into %s (%d) " .
                    "due to larger exon overlap (%dbp > %dbp)\n",
                  $secondary_transcript->stable_id(),
                  $secondary_transcript->dbID(),
                  $distinct_transcripts_to_copy{$key}[0]->stable_id(),
                  $distinct_transcripts_to_copy{$key}[0]->dbID(),
                  $distinct_transcripts_to_copy{$key}[4],
                  $overlap_score );
        }
      } ## end if ( $primary_transcript...)

    } ## end if ( exists( $distinct_transcripts_to_copy...))
    else {
      my $exon_overlap = exon_overlap($secondary_transcript,$primary_transcript);
      my $exon_non_overlap = abs($secondary_transcript->length()-$primary_transcript->length());
      my $overlap_score = $exon_overlap-$exon_non_overlap;
      $distinct_transcripts_to_copy{$key} = [
                 $primary_gene,
                 $primary_transcript,
                 $secondary_gene,
                 $secondary_transcript,
                 $overlap_score
      ];
    }

  } ## end foreach my $tuple (@transcripts_to_copy)

  # Copy the transcripts that have been selected for copying, but filter
  # out the ones with incomplete start/stop codons (unless they are
  # copied into genes with no protein coding transcripts).  Also, make
  # sure (by sorting) that we copy all CCDS transcripts first as these
  # have a small chance of turning non-coding genes into coding genes
  # (which means that if we do non-CCDS transcripts first, these may
  # loose their translations unnecessarily).
  foreach my $tuple ( sort { $b->[3]{__is_ccds} <=> $a->[3]{__is_ccds} }
                      values(%distinct_transcripts_to_copy) )
  {
    my ( $primary_gene,  $primary_transcript,
         $secondary_gene, $secondary_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $secondary_transcript->dbID();

    if ( !exists( $merged_secondary_transcripts{$key} ) &&
         !exists( $copied_secondary_transcripts{$key} ) )
    {
      if ( !defined( $primary_transcript->translation() ) ||
           !defined( $secondary_transcript->translation() ) ||
           has_complete_start_stop($secondary_transcript) )
      {
        printf( "Copying %s (%d) into %s (%d)\n",
                $secondary_transcript->stable_id(),
                $secondary_transcript->dbID(),
                $primary_gene->stable_id(),
                $primary_gene->dbID() );

        my $copy_code = copy( $primary_gene, $secondary_transcript );

        $copied_secondary_transcripts{$key} = 1;

        if ( $copy_code != 1 ) {
          printf( "PROCESSED\t%d\t%s\n",
                  $secondary_gene->dbID(), $secondary_gene->stable_id() );

          # Keep a record in memory of this stable id
          $processed_ens_id_tracker{$primary_gene->dbID()} = [$secondary_gene->dbID(),$secondary_gene->stable_id()];
        }
      }
      else {
        printf( "May ignore %s (%d) due to incomplete start/stop\n",
                $secondary_transcript->stable_id(),
                $secondary_transcript->dbID() );

        push( @transcripts_to_ignore,
              [ $primary_gene,  $primary_transcript,
                $secondary_gene, $secondary_transcript ] );
      }
    } ## end if ( !exists( $merged_secondary_transcripts...))
  } ## end foreach my $tuple ( sort { ...})

  # Ignore transcripts that should be ignored (due to incomplete
  # start/stop), but only if other transcripts from the same gene
  # have been merged or copied.  Otherwise copy them.  (It's probably
  # not necessary to sort here since all CCDS models presumably have
  # complete start/stop codons and thus wouldn't show up in this list,
  # but you never know).
  foreach my $tuple ( sort { $b->[3]{__is_ccds} <=> $a->[3]{__is_ccds} }
                      @transcripts_to_ignore )
  {
    my ( $primary_gene,  $primary_transcript,
         $secondary_gene, $secondary_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $secondary_transcript->dbID();

    if ( !exists( $merged_secondary_transcripts{$key} ) &&
         !exists( $copied_secondary_transcripts{$key} ) &&
         !exists( $ignored_secondary_transcripts{$key} ) )
    {
      my $do_ignore = 0;

      foreach
        my $transcript ( @{ $secondary_gene->get_all_Transcripts() } )
      {

        my $key2 = $transcript->dbID();

        if ( $key ne $key2 ) {
          if ( exists( $copied_secondary_transcripts{$key2} ) ||
               exists( $merged_secondary_transcripts{$key2} ) )
          {
            $do_ignore = 1;
            last;
          }
        }
      }

      my $copy_code = 0;
      if ($do_ignore) {
        printf( "Ignoring %s (%d) due to incomplete start/stop\n",
                $secondary_transcript->stable_id(),
                $secondary_transcript->dbID() );

        $ignored_secondary_transcripts{$key} = 1;
      }
      else {
        printf( "Copying %s (%d) into %s (%d) " .
                  "even though it has incomplete start/stop\n",
                $secondary_transcript->stable_id(),
                $secondary_transcript->dbID(),
                $primary_gene->stable_id(),
                $primary_gene->dbID() );

        $copy_code = copy( $primary_gene, $secondary_transcript );
      }

      if ( $copy_code != 1 ) {
        printf( "PROCESSED\t%d\t%s\n",
                $secondary_gene->dbID(), $secondary_gene->stable_id() );

        # Keep a record in memory of this stable id
        $processed_ens_id_tracker{$primary_gene->dbID()} = [$secondary_gene->dbID(),$secondary_gene->stable_id()];
      }

    } ## end if ( !exists( $merged_secondary_transcripts...))
  } ## end foreach my $tuple ( sort { ...})

  # Do a bit of houskeeping from the previous loop so that
  # %copied_secondary_transcripts correctly reflects the Secondary
  # transcripts actually copied.
  foreach my $tuple (@transcripts_to_ignore) {
    my ( $primary_gene,  $primary_transcript,
         $secondary_gene, $secondary_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $secondary_transcript->dbID();

    if ( !exists( $ignored_secondary_transcripts{$key} ) ) {
      $copied_secondary_transcripts{$key} = 1;
    }
  }

  foreach my $primary_gene ( @{$primary_genes} ) {

    my $old_dbID = $primary_gene->dbID();

    # If the biotype is bad, skip the store
    if(bad_biotype($primary_gene->biotype())) {

        print "WARNING: skipping store of Primary gene ".$primary_gene->dbID()." ".
              $primary_gene->stable_id()." due to bad biotype: ".$primary_gene->biotype()."\n";

        # If there is an Secondary gene that has been merged/copied to this Primary gene
        # put in some additional warnings that the gene will be thrown out
        if($processed_ens_id_tracker{$old_dbID}) {         
          my $ens_gene_db_id = $processed_ens_id_tracker{$old_dbID}->[0];
          my $ens_gene_stable_id = "";
          if($processed_ens_id_tracker{$old_dbID}->[1]) {
            $ens_gene_stable_id = $processed_ens_id_tracker{$old_dbID}->[1];
          } 

          print "WARNING: skipping store of Secondary gene ".$ens_gene_db_id." ".$ens_gene_stable_id.
                " due to merge/copy with Primary gene ".$primary_gene->dbID()." ".$primary_gene->stable_id().
                " with bad biotype: ".$primary_gene->biotype()."\n";          
        }

       next;
     }

    if (!(exists($stored_genes{$primary_gene->stable_id()}))) {
      empty_Gene($primary_gene);
      $OUTPUT_GA->store($primary_gene);
      printf( "STORED\t%s\told id = %d, new id = %d\n",
              $primary_gene->stable_id(),
              $old_dbID, $primary_gene->dbID() );
      $stored_genes{$primary_gene->stable_id()} = 1;
    } else {
      print("NOT STORED as it is already stored ".$primary_gene->stable_id()."\n");
    }

  }

} ## end sub process_genes


# This returns 1 if the biotype passed in is in a hash of bad biotypes
sub bad_biotype
{
  my ($biotype) = @_;

  # This would be more efficient to initialise at the top of the code,
  # but as the list will never be huge I'm putting it here as it's 
  # easier to read
  my %bad_biotypes = (
                       'artifact' => 1
                     );

  if($bad_biotypes{$biotype}) {
    return 1; 
  }

  else {
    return 0;
  }
  
} # End sub bad_biotype


sub add_primary_xref {
  my ($feature) = @_;

  my ( $external_db_name, $db_display_name, $type );

  if ( $feature->isa('Bio::EnsEMBL::Gene') ) {
    ( $external_db_name, $db_display_name, $type ) =
      split( /,/, $opt_primary_gene_xref );

    # $external_db_name = 'OTTG';
    # $db_display_name  = 'Havana gene';
    # $type             = 'ALT_GENE';
  }
  elsif ( $feature->isa('Bio::EnsEMBL::Transcript') ) {
    ( $external_db_name, $db_display_name, $type ) =
      split( /,/, $opt_primary_transcript_xref );

    # $external_db_name = 'OTTT';
    # $db_display_name  = 'Havana transcript';
    # $type             = 'ALT_TRANS';
  }
  elsif ( $feature->isa('Bio::EnsEMBL::Translation') ) {
    ( $external_db_name, $db_display_name, $type ) =
      split( /,/, $opt_primary_translation_xref );

    # $external_db_name = 'OTTP';
    # $db_display_name  = 'Havana translation';
    # $type             = 'MISC';
  }
  else {
    die("Can't add xref to unknown type of object");
  }


  unless($external_db_name && $db_display_name && $type) {

    throw("Could not assign one or all of the following: external_db_name, db_display_name or type\n".
          "If you are using the wrapper script, make sure the corresponding values are set in the\n".
          "config file. Typical values are:\n".
          "primary_gene_xref='OTTG,Havana gene,ALT_GENE'\n".
          "primary_transcript_xref='OTTT,Havana transcript,ALT_TRANS'\n".
          "primary_translation_xref='OTTP,Havana translation,MISC'\n"
         );
  }

  my $xref =
    Bio::EnsEMBL::DBEntry->new( '-dbname'          => $external_db_name,
                                '-db_display_name' => $db_display_name,
                                '-type'            => $type,
                                '-primary_id' => $feature->stable_id(),
                                '-display_id' => $feature->stable_id()
    );

  $feature->add_DBEntry($xref);
} ## end sub add_primary_xref

sub investigate_for_merge {
  my ( $primary_transcript, $secondary_transcript ) = @_;

  my @primary_introns = @{ $primary_transcript->get_all_Introns() };

  my $do_merge;

  if ( scalar(@primary_introns) == 0 ) {
    # This is a single exon transcript, compare the exons instead of
    # the (non-existing) introns.

    $do_merge = features_are_same( $primary_transcript->get_all_Exons(),
                                   $secondary_transcript->get_all_Exons()
    );

    if ( !$do_merge &&
         ( defined( $primary_transcript->translation() ) &&
           defined( $secondary_transcript->translation() ) ) )
    {
      # This is a single exon coding transcript, compare the coding
      # exons instead of the (non-existing) introns.

      $do_merge = features_are_same(
                   $primary_transcript->get_all_translateable_Exons(),
                   $secondary_transcript->get_all_translateable_Exons() );

      if ($do_merge) {
        print( "\t\t\tSpecial case: Single exon coding transcripts, " .
               "coding region is identical but UTR differs\n" );
      }
    }

    if ( !$do_merge ) {
      # A special case:  If the primary exon is exactly one stop
      # codon longer at the end, then treat them as identical.
      $do_merge = is_a_stop_codon_longer( $primary_transcript,
                                          $secondary_transcript );
      if ($do_merge) {
        print( "\t\t\t\tSpecial case: Single exon transcripts, " .
               "Secondary exon is stop codon short\n" );
      }
    }
    if ( !$do_merge ) {
      # A special case:  If the secondary exon is exactly one stop
      # codon longer at the end, then treat them as identical.
      $do_merge = is_a_stop_codon_longer( $secondary_transcript,
                                          $primary_transcript );
      if ($do_merge) {
        print( "\t\t\t\tSpecial case: Single exon transcripts, " .
               "Primary exon is stop codon short\n" );
      }
    }

  } ## end if ( scalar(@primary_introns...))
  else {
    my @secondary_introns = @{ $secondary_transcript->get_all_Introns() };

    $do_merge = features_are_same(\@primary_introns,\@secondary_introns);

    if ( !$do_merge ) {
      # A special case: If the secondary transcript contains a short frameshift
      # intron (1, 2, 4 or 5 bp as labelled by the script ensembl/misc-scripts/frameshift_transcript_attribs.pl),
      # treat it as if there was no short frameshift intron

      my @attribs = @{$secondary_transcript->get_all_Attributes('Frameshift')};
      if (scalar(@attribs) > 0) {

        my $secondary_transcript_no_frameshift = remove_short_frameshift_introns($secondary_transcript);
        @secondary_introns = @{$secondary_transcript_no_frameshift->get_all_Introns()};
      
        $do_merge = features_are_same(\@primary_introns,
                                      \@secondary_introns);
        if ($do_merge) {
          
          print("\t\t\t\tSpecial case: Secondary transcript ".$secondary_transcript_no_frameshift->stable_id()." containing short frameshift intron investigated for merge.\n");
          
          print("Primary translation/Secondary translation:\n");
          if ($primary_transcript->translation()) {
            print($primary_transcript->translation()->seq()."\n\n");
          } else {
          	print("Primary transcript has no translation\n\n");
          }
          if ($secondary_transcript->translation()) {
            print($secondary_transcript->translation()->seq()."\n\n");
          } else {
          	print("Secondary transcript has no translation\n\n");
          }
        }
      } # if scalar
    } # if !$do_merge
  }

  return $do_merge;
} ## end sub investigate_for_merge

sub investigate_for_partial_merge {
  # it returns true if a partial merge can be done
  # a partial merge can be done if the secondary transcript exons
  # are a subset of the primary transcript exons
  my ($primary_transcript,$secondary_transcript) = @_;

  return features_are_subset($primary_transcript->get_all_Exons(),
                             $secondary_transcript->get_all_Exons());
                                   
} ## end sub investigate_for_partial_merge

sub investigate_for_copy {
  my ( $primary_gene,  $primary_transcript,
       $secondary_gene, $secondary_transcript ) = @_;

  my @secondary_exons = @{ $secondary_transcript->get_all_Exons() };
  my @primary_exons  = @{ $primary_transcript->get_all_Exons() };

  my $do_copy = 0;

SECONDARY_EXON:
  foreach my $secondary_exon (@secondary_exons) {
    foreach my $primary_exon (@primary_exons) {
      if ( features_overlap( $secondary_exon, $primary_exon ) ) {
        $do_copy = 1;

        if ( scalar(@secondary_exons) == 1 &&
             $secondary_gene->{__is_single_transcript} &&
             scalar(@primary_exons) > 1 )
        {
          # We have a single exon Secondary gene that overlaps with bits
          # of a multi exon (and possibly multi transcript) Primary gene.
          #
          # The current rules are:
          #
          # * Copy if the overlapping Primary exon is identical to the
          #   Secondary exon.
          #
          # * Don't copy if the Secondary gene is RNA and the Primary gene
          #   has a different biotype.
          #
          # * Copy in any other case.

          if ( features_are_same( [$secondary_exon], [$primary_exon] ) ) {
            print( "\t\t\tSpecial case: Allowing copy of transcript " .
                   "from single exon gene into multi exon gene " .
                   "(perfect exon match)\n" );
          }
          elsif ( $secondary_gene->{__is_rna} &&
                  $secondary_gene->biotype() ne $primary_gene->biotype() )
          {
            print( "\t\t\tSpecial case: Won't copy transcript " .
                   "from single exon RNA gene into " .
                   "multi exon gene with different biotype\n" );
            $do_copy = 0;
          }
          else {
            print( "\t\t\tSpecial case: Allowing copy of transcript " .
                   "from single exon gene into multi exon gene\n" );
          }
        } ## end if ( scalar(@secondary_exons...))
        elsif ( $secondary_gene->{__is_rna} &&
                $secondary_gene->biotype() ne $primary_gene->biotype() )
        {
          # Don't allow copying of RNA gene transcripts if gene biotypes
          # do not match upi (regardless of exon counts).
          print( "\t\t\tSpecial case: Won't copy transcript " .
                 "from RNA gene into gene with different biotype\n" );
          $do_copy = 0;
        }

        if ($do_copy) {
          last SECONDARY_EXON;
        }
      } ## end if ( features_overlap(...))
    } ## end foreach my $primary_exon (@primary_exons)
  } ## end SECONDARY_EXON: foreach my $secondary_exon (@secondary_exons)

  return $do_copy;
} ## end sub investigate_for_copy

sub features_are_same {
  my ( $feature_set_a, $feature_set_b ) = @_;

  if ( scalar( @{$feature_set_a} ) == 0 ||
       ( scalar( @{$feature_set_a} ) != scalar( @{$feature_set_b} ) ) )
  {
    return 0;
  }

  for ( my $feature_index = 0;
        $feature_index < scalar( @{$feature_set_a} );
        ++$feature_index )
  {
    my $feature_a = $feature_set_a->[$feature_index];
    my $feature_b = $feature_set_b->[$feature_index];

    if (
      ( $feature_a->seq_region_start() != $feature_b->seq_region_start()
      ) ||
      ( $feature_a->seq_region_end() != $feature_b->seq_region_end() ) )
    {
      return 0;
    }
  }

  return 1;
} ## end sub features_are_same

sub features_are_subset {
# returns true if feature_set_b is a contiguous subset of feature_set_a
# first feature start and last feature end are allowed to mismatch
  my ( $feature_set_a, $feature_set_b ) = @_;

  if (scalar(@{$feature_set_a}) == 0 ||
     (scalar(@{$feature_set_a}) < scalar(@{$feature_set_b}))) {
    return 0;
  }

  my $is_subset = 1;
  my $feature_b_index = 0;
  my $contiguous_match = 0;
  
  for (my $feature_a_index = 0; $feature_a_index < scalar(@{$feature_set_a}); ++$feature_a_index) {
    my $feature_a = $feature_set_a->[$feature_a_index];
    my $feature_b = $feature_set_b->[$feature_b_index];
    my $start_match = start_features_match($feature_a,$feature_b,$feature_b_index);
    my $end_match = end_features_match($feature_a,$feature_b,$feature_b_index,scalar(@{$feature_set_b})-1);
    my $room_for_feature_b = (scalar(@{$feature_set_b})-$feature_b_index-1 <= scalar(@{$feature_set_a})-$feature_a_index-1);

    if (!$start_match or !$end_match) {
      if ($contiguous_match) {
        # if all previous features match and we stop matching having some features left to match,
        # return NO match as all the matches have to be contiguous
        return 0;
      }
      
      if (!$room_for_feature_b) {
        # no more room for feature b
        return 0;
      }
    } elsif ($start_match and $end_match) {
      # match and still room
      $contiguous_match = 1;
      $feature_b_index++;
      if ($feature_b_index > scalar(@{$feature_set_b})-1) {
      # final match
        return 1;
      }
    }
  }
  # NO match
  return 0;
} ## end sub features_are_subset

sub start_features_match () {
  # return true if feature b start lies within feature a for the first b feature
  # or feature b start matches feature a start
  
  my ($feature_a,$feature_b,$index_b) = @_;
  my $start_a = $feature_a->seq_region_start();
  my $start_b = $feature_b->seq_region_start();
  my $end_a = $feature_a->seq_region_end();

  return ( ($start_a == $start_b and $index_b > 0) or
           ($start_a <= $start_b and $start_b <= $end_a and $index_b == 0)
         );
}

sub end_features_match () {
  # return true if feature b end lies within feature a for the last b feature
  # or feature b end matches feature a end
  my ($feature_a,$feature_b,$index_b,$last_index_b) = @_;
  my $end_a = $feature_a->seq_region_end();
  my $end_b = $feature_b->seq_region_end();
  my $start_a = $feature_a->seq_region_start();

  return ( ($end_a == $end_b and $index_b < $last_index_b) or
           ($start_a <= $end_b and $end_b <= $end_a and $index_b == $last_index_b)
         );
}

sub has_complete_start_stop {
  my ($transcript) = @_;

  if ( defined( $transcript->translation() ) ) {
    my $cdna = $transcript->translateable_seq();
    my $beg  = uc( substr( $cdna, 0, 3 ) );
    my $end  = uc( substr( $cdna, -3 ) );

    # /^ATG.*T(AG|AA|GA)$/
    if ( ( $beg eq 'ATG' ) &&
       ( ( $end eq 'TAG' ) || ( $end eq 'TAA' ) || ( $end eq 'TGA' ) ) )
    {
      return 1;
    }
  }

  return 0;
}

sub is_a_stop_codon_longer {
  my ( $transcriptA, $transcriptB ) = @_;

  # Returns true (non-zero) if $transcriptA is a stop codon longer than
  # $transcriptB.

  if ( $transcriptA->length() == $transcriptB->length() + 3 &&
       defined( $transcriptA->translation() )
       && (
         ( $transcriptA->strand() != -1 &&
           ( $transcriptA->seq_region_start() ==
             $transcriptB->seq_region_start() ) )
         ||
         ( $transcriptA->strand() == -1 &&
           ( $transcriptA->seq_region_end() ==
             $transcriptB->seq_region_end() ) ) ) )
  {
    my $cdna = $transcriptA->translateable_seq();
    my $end = uc( substr( $cdna, -3 ) );

    # /T(AG|AA|GA)$/
    if ( ( $end eq 'TAG' ) || ( $end eq 'TAA' ) || ( $end eq 'TGA' ) ) {
      return 1;
    }
  }

  return 0;
} ## end sub is_a_stop_codon_longer

sub merge {
# Returns 1 if $source_transcript should be copied over
# at the end of the merge process (further processing is needed).
# Otherwise, returns 0 (further processing is not needed).

  my ( $target_gene, $target_transcript, $source_transcript ) = @_;

  printf( "Merge> Biotypes: source transcript: %s, " .
            "target transcript/gene: %s/%s\n",
          $source_transcript->biotype(),
          $target_transcript->biotype(),
          $target_gene->biotype() );

  if ( $target_gene->{__has_ref_error} ) {
    print("Merge> Target gene has assembly error\n");

    if ( $target_gene->biotype() ne $source_transcript->biotype() ) {
      printf( "Merge> Updating gene biotype from %s to %s\n",
              $target_gene->biotype(), $source_transcript->biotype() );
      $target_gene->biotype( $source_transcript->biotype() );
    }

    print( "Merge> Copying source transcript to target gene " .
           "(not merging)\n" );
    return copy( $target_gene, $source_transcript );
  }

  elsif ( $source_transcript->{__is_ccds} ) {
    if ($source_transcript->biotype() ne $target_transcript->biotype() )
    {
      printf( "Merge> Merge would demote biotype from %s to %s " .
                "(CCDS source transcript)\n",
              $source_transcript->biotype(),
              $target_transcript->biotype() );

      print( "Merge> Copying source transcript to target gene " .
             "(not merging)\n" );
      return copy( $target_gene, $source_transcript );
    }

    elsif($source_transcript->translation()->seq() ne $target_transcript->translation()->seq()) {
      printf( "Merge> Merge would alter coding sequence (CCDS source transcript)\n");
      print( "Merge> Copying source transcript to target gene " .
             "(not merging)\n" );
       return copy( $target_gene, $source_transcript );
    }

  }

  # Start by transferring the $source_transcript to the same slice as
  # the $target_gene.
  my $new_source_transcript =
    $source_transcript->transfer( $target_gene->slice() );

  # Copy all transcript related features from $source_transcript into
  # $target_transcript:
  #     supporting features
  #     intron supporting evidence
  {
    my @supporting_features = ();
    
    # only transfer the supporting features that overlap.
    # as the merge is based on intron match, there can be cases where
    # a longer Secondary transcript evidence would have been transferred
    # beyond the exon boundaries of the Primary target transcript
    foreach my $sf (@{ $new_source_transcript->get_all_supporting_features() }) {
      if (features_overlap($sf,$target_transcript)) {
        push(@supporting_features,$sf);
      }
    }

    printf( "Merge> Transferred %d supporting feature(s)\n",
            scalar(@supporting_features) );

    $target_transcript->add_supporting_features(@supporting_features);

    my @intron_support =
      @{ $new_source_transcript->get_all_IntronSupportingEvidence() };

    foreach my $intron_support (@intron_support) {
      $target_transcript->add_IntronSupportingEvidence($intron_support);
    }

    printf( "Merge> Transferred %d intron supporting evidence\n",
            scalar(@intron_support) );
  }

  # Transfer all exon related features from $source_transcript into
  # $target_transcript:
  #     supporting features
  {
    my @supporting_features;

    foreach
      my $source_exon ( @{ $new_source_transcript->get_all_Exons() } ) {
      my @exon_sf = ();
      # only transfer the supporting features that overlap.
      # as the merge is based on intron match, there can be cases where
      # a longer Secondary transcript evidence would have been transferred
      # beyond the exon boundaries of the Primary target transcript
      foreach my $sf (@{ $source_exon->get_all_supporting_features() }) {
        if (features_overlap($sf,$target_transcript)) {
          push(@exon_sf,$sf);
        }
      }
      push(@supporting_features,[@exon_sf]);
    }

    my $exon_index    = 0;
    my $feature_count = 0;
    foreach my $target_exon ( @{ $target_transcript->get_all_Exons() } )
    {
      $feature_count +=
        scalar( @{ $supporting_features[$exon_index] } );

      $target_exon->add_supporting_features(
                           @{ $supporting_features[ $exon_index++ ] } );

      add_logic_name_suffix( $target_exon, 'merged' );
    }

    printf( "Merge> Transferred %d " .
              "exon supporting evidence (%d exons)\n",
            $feature_count, $exon_index );

  }
  
  $target_gene->source($merged_source);
  $target_transcript->source($merged_source);

  $target_gene->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($target_gene)));
  $target_transcript->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($target_transcript)));

  return 0;
} ## end sub merge

sub copy {
# Returns 1 if $source_transcript should be copied over
# at the end of the merge process (further processing is needed).
# Otherwise, returns 0 (further processing is not needed).

  my ( $target_gene, $source_transcript ) = @_;

  printf( "Copy> Biotypes: source transcript: %s, target gene: %s\n",
          $source_transcript->biotype(), $target_gene->biotype() );

###############################################################################
# Case 1 - Primary gene is coding, Secondary transcript is a pseudogene
# Do not copy the Secondary transcript and the corresponding gene will be listed
# as processed and will not be copied at the end.
###############################################################################
  if ( $source_transcript->{__is_pseudogene} and
       ($target_gene->{__is_coding} or $target_gene->{__is_polymorphic_pseudogene}) )
  {
    print( "Copy> Source transcript is pseudogene, " .
           "will not copy it into a coding or polymorphic_pseudogene gene.\n" );
    print( "Copy> Deleting the Secondary annotation.\n");
    return 0;
  }

###############################################################################
# Case 2 - The Primary gene is a pseudogene. In this case it the Secondary
# transcript won't be copied and the corresponding gene will be listed
# as processed and will not be copied at the end. The only exception to
# this is when the Secondary transcript is CCDS, in which case the transcript
# will be copied and the biotype will be updated to the Secondary biotype.
# The Secondary gene will be listed as processed and will not be copied at
# the end of merge process
###############################################################################
  elsif ( $target_gene->{__is_pseudogene} ) {
  

    unless ( $source_transcript->{__is_ccds} ) {
       print( "Copy> Target gene is pseudogene, " .
              "will not copy anything into it.\n" );
       print( "Copy> Deleting the Secondary annotation.\n");
       return 0;
    }

    else {
      printf( "Copy> Updating gene biotype from %s to %s " .
              "(CCDS source transcript)\n",
              $target_gene->biotype(), $source_transcript->biotype() );

      $target_gene->biotype( $source_transcript->biotype() );
      $target_gene->{__is_coding} = 1;
    }   
  
  }

###############################################################################
# Case 3 - The Primary gene is labelled as belonging to a gene cluster. This
# tag is read from the transcripts, so at least one transcript was labelled.
# In this case the Secondary transcript will not be copied over.
# The exception to this is if the Secondary transcript is CCDS, in this case
# the Secondary transcript will be copied and the Primary gene biotype will
# updated if needed. Either way the corresponding the Secondary gene will be
# will be listed as processed and not copied at the end of the merge process.
###############################################################################
  elsif ( $target_gene->{__is_gene_cluster} ) {

    unless ( $source_transcript->{__is_ccds} ) {
      print( "Copy> Target gene is part of gene cluster, " .
           "will not copy overlapping Secondary transcripts ".
           "into it.\n" );
      print( "Copy> Deleting the Secondary annotation.\n");
      return 0;
    }

    elsif ($target_gene->biotype() ne $source_transcript->biotype()) {
      
      printf( "Copy> Updating gene biotype from %s to %s " .
              "(CCDS source transcript)\n",
              $target_gene->biotype(), $source_transcript->biotype() );

      $target_gene->biotype( $source_transcript->biotype() );
      $target_gene->{__is_coding} = 1;
    }
   
  }

###############################################################################
# Case 4 - The Primary gene has an assembly error, in this case the Secondary
# transcript will be copied in. The biotype of the Secondary transcript will
# overwrite the Primary gene biotype. This may be in the wrong place logically.
# The corresponding Secondary gene is listed as processed and will not be
# copied at the end of the merge process
###############################################################################
  elsif ( $target_gene->{__has_ref_error} ) {
    print( "Copy> Target gene has assembly error\n");

# If there is an assembly error, we would almost all the time have a pseudogene
# with only one transcript per gene, so we can update the biotype in any case.
# If it's really wrong, vega_check will complain
      printf( "Copy> Updating gene biotype from %s to %s\n",
              $target_gene->biotype(), $source_transcript->biotype() );
      $target_gene->biotype( $source_transcript->biotype() );

  }

###############################################################################
# Case 5 - The Primary gene is non coding (but not a pseudogene) and the
# Secondary transcript has a translation. In this case the translation is
# removed from the Secondary transcript and the exon phases are all set to
# -1, which means non-coding, before the transcript is copied. The only
# exception to this is when the Secondary transcript is CCDS, in this case
# the biotype of the Primary gene is overwritten with the Secondary transcript
# biotype. The corresponding Secondary gene is listed as processed and is
# not copied over at the end of the merge process
###############################################################################
  elsif ( !$target_gene->{__is_coding} &&
          defined( $source_transcript->translation() ) ) {
    unless ( $source_transcript->{__is_ccds} ) {
      print( "Copy> Removing translation from source transcript\n");

      $source_transcript->translation(undef);
      $source_transcript->dbID(undef);       # HACK
      $source_transcript->adaptor(undef);    # HACK

      print "Copy> Stripping exon phases for: ".$source_transcript->stable_id()."\n"; 
      strip_phase(\$source_transcript);

      printf( "Copy> Updating transcript biotype from %s to %s\n",
              $source_transcript->biotype(), $target_gene->biotype() );

      $source_transcript->biotype( $target_gene->biotype() );
    }

    else {
      printf( "Copy> Updating gene biotype from %s to %s " .
              "(CCDS source transcript)\n",
              $target_gene->biotype(), $source_transcript->biotype() );

      $target_gene->biotype( $source_transcript->biotype() );
      $target_gene->{__is_coding} = 1;
    }

  }

  # Transfer the $source_transcript to the same slice as
  # the $target_gene if it has not been transferred before.
  if (is_transcript_in_gene($target_gene,$source_transcript)) {
    print "Copy> Not copying because it has already been copied (or merged) or it will be merged later on ".$source_transcript->stable_id()."\n"; 
  } else {
  	my $new_source_transcript = $source_transcript->transfer( $target_gene->slice() );

    $new_source_transcript->source($opt_secondary_tag);

    $new_source_transcript->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($new_source_transcript)));
    $target_gene->add_Transcript($new_source_transcript);
    $target_gene->source($merged_source);
    $target_gene->analysis($OUTPUT_AA->fetch_by_logic_name(get_logic_name_from_biotype_source($target_gene)));
  }
  return 0;
} ## end sub copy


# Strip the phases off all the exons in a transcript. Used on Secondary
# coding transcripts that get demoted to non-coding. This will only
# be run on coding Secondary transcripts that are copied to non-coding, 
# non-pseudogene Primary genes. Yes that's confusing, just roll with it.
sub strip_phase {

  my ($transcript_to_strip) = @_;  
  my $exon_refs = $$transcript_to_strip->get_all_Exons();
  foreach my $exon (@{$exon_refs}) {
    $exon->phase(-1);
    $exon->end_phase(-1);
  }

}


sub tag_transcript_analysis {
  my ( $transcript, $suffix ) = @_;

  if ( !defined($suffix) ) { return }

  # Tag all transcript supporting features:
  foreach my $supporting_feature (
                       @{ $transcript->get_all_supporting_features() } )
  {
    add_logic_name_suffix( $supporting_feature, $suffix );
  }

  # Tag all intron supporting evidences:
  foreach my $intron_support (
                  @{ $transcript->get_all_IntronSupportingEvidence() } )
  {
    add_logic_name_suffix( $intron_support, $suffix );
  }

  # Tag all exons:
  foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
    add_logic_name_suffix( $exon, $suffix );

    # Also tag all exon supporting features:
    foreach my $supporting_feature (
                             @{ $exon->get_all_supporting_features() } )
    {
      add_logic_name_suffix( $supporting_feature, $suffix );
    }
  }

  # Finally tag the transcript itself:
  add_logic_name_suffix( $transcript, $suffix );

} ## end sub tag_transcript_analysis

sub add_logic_name_suffix {
  my ( $feature, $suffix ) = @_;

  if ( defined($feature) ) {
    my $analysis = $feature->analysis();

    if ( defined($analysis) ) {
      my $logic_name = $analysis->logic_name();

      if ( index( $logic_name, $suffix ) == -1 ) {
        my $new_analysis =
          Bio::EnsEMBL::Analysis->new(
                           '-logic_name' => $logic_name . '_' . $suffix,
                           '-module'     => $analysis->module() );
        $feature->analysis($new_analysis);
      }
    }
  }
}

sub is_transcript_in_gene {
# returns 1 if the transcript $source_transcript
# can be found in the gene $target_gene
  my ($target_gene,$source_transcript) = @_;
  
  my $transcript_key = get_transcript_exon_key($source_transcript);
  
  foreach my $transcript (@{$target_gene->get_all_Transcripts()}) {
    if (get_transcript_exon_key($transcript) eq $transcript_key) {
      return 1;
    }
  }
  return 0;
}

sub get_transcript_exon_key {
  my $transcript = shift;
  my $string = $transcript->slice()->seq_region_name().":".$transcript->biotype().":".$transcript->seq_region_start().":".$transcript->seq_region_end().":".$transcript->seq_region_strand().":".$transcript->coding_region_start()
.":".$transcript->coding_region_end()
.":".@{$transcript->get_all_translateable_Exons()}.":";

  my $exons = sort_by_start_end_pos($transcript->get_all_Exons());
  foreach my $exon (@{$exons}) {
    $string .= ":".$exon->seq_region_start().":".$exon->seq_region_end();
  }

  return $string;
}

sub sort_by_start_end_pos {
  my ($unsorted) = @_;

  my @sorted = sort { if ($a->seq_region_start < $b->seq_region_start) {
        return -1;
    } elsif ($a->seq_region_start == $b->seq_region_start) {
      if ($a->seq_region_end < $b->seq_region_end) {
        return-1;
      } elsif ($a->seq_region_end == $b->seq_region_end) {
        return 0;
      } elsif ($a->seq_region_end > $b->seq_region_end) {
        return 1;
      }
        return 0;
    } elsif ($a->seq_region_start > $b->seq_region_start) {
        return 1;
    }
  } @$unsorted;

  return \@sorted;
}

sub get_logic_name_from_biotype_source() {
# $obj must be a gene or a transcript
# $biotype must be any allowed biotype
# $source must be $merged_source or $opt_secondary_tag or $opt_primary_tag
# It returns the right logic_name for the given parameters above.
  my ($obj) = @_;

  my $biotype = $obj->biotype();
  my $source = $obj->source();

  if ($obj->isa("Bio::EnsEMBL::Gene")) {
    if ($source eq $merged_source) {
      
      if (index(lc($biotype),lc('lincRNA')) != -1) {
        return $merged_lincrna_logic_name;
      } elsif (index(lc($biotype),lc('IG_')) != -1 or
               index(lc($biotype),lc('TR_')) != -1) {
        return $merged_ig_gene_logic_name;
      } else {
        return $merged_gene_logic_name;
      }
      
    } elsif ($source eq $opt_primary_tag) {
      
      if (index(lc($biotype),lc('IG_')) != -1 or
          index(lc($biotype),lc('TR_')) != -1) {
        return $primary_ig_gene_logic_name; # yes, the gene logic name is re-used for transcript
      } else {
        return $opt_primary_tag; # yes, no lincrna logic name is used for transcript
      }
      
    } elsif ($source eq $opt_secondary_tag) {
      
      if (index(lc($biotype),lc('lincRNA')) != -1) {
        return $secondary_lincrna_logic_name;
      } elsif (index(lc($biotype),lc('IG_')) != -1 or
               index(lc($biotype),lc('TR_')) != -1) {
        return $secondary_ig_gene_logic_name; # yes, the gene logic name is re-used for transcript
      } elsif ($obj->{__is_rna}) {
        return $obj->analysis()->logic_name(); # it should already be 'ncrna'
      } else {
        return $opt_primary_tag; # yes, no lincrna logic name is used for transcript
      }
      
    }
  }
  elsif ($obj->isa("Bio::EnsEMBL::Transcript")) {
    if ($source eq $merged_source) {
        
      if (index(lc($biotype),lc('IG_')) != -1 or
          index(lc($biotype),lc('TR_')) != -1) {
        return $merged_ig_gene_logic_name; # yes, the gene logic name is re-used for transcript
      } else {
        return $merged_transcript_logic_name; # yes, no merged lincrna logic name is used for transcript
      }
      
    } elsif ($source eq $opt_primary_tag) {
      
      if (index(lc($biotype),lc('IG_')) != -1 or
          index(lc($biotype),lc('TR_')) != -1) {
        return $primary_ig_gene_logic_name; # yes, the gene logic name is re-used for transcript
      } else {
        return $opt_primary_tag; # yes, no lincrna logic name is used for transcript
      }
      
    } elsif ($source eq $opt_secondary_tag) {
      
      if (index(lc($biotype),lc('lincRNA')) != -1) {
        return $secondary_lincrna_logic_name;
      } elsif (index(lc($biotype),lc('IG_')) != -1 or
               index(lc($biotype),lc('TR_')) != -1) {
        return $secondary_ig_gene_logic_name; # yes, the gene logic name is re-used for transcript
      } elsif ($obj->{__is_rna}) {
        return 'ncrna';#$obj->analysis()->logic_name(); # it should already be 'ncrna'
      } else {
        return $opt_secondary_tag; # yes, no lincrna logic name is used for transcript
      }
      
    }
  }
  else {
    throw("Invalid object type. Valid types are Bio::EnsEMBL::Gene and Bio::EnsEMBL::Transcript.");
  }
  
}

__END__

=head1 NAME

merge.pl - A script that merges the manual annotations from Primary with
the automatic annotations from Secondary.

=head1 SYNOPSIS

merge.pl [options]

merge --help

  Options:
    --host_secondary  (required)
    --port_secondary
    --user_secondary  (required)
    --pass_secondary      or --password_secondary
    --dbname_secondary    or --database_secondary   (required)

    --host_primary   (required)
    --port_primary
    --user_primary   (required)
    --pass_primary       or --password_primary
    --dbname_primary     or --database_primary    (required)

    --host_dna
    --port_dna
    --user_dna
    --pass_dna       or --password_dna
    --dbname_dna     or --database_dna

    --host_output   (required)
    --port_output
    --user_output   (required)
    --pass_output       or --password_output
    --dbname_output     or --database_output    (required)

    --secondary_include   or --secondary_exclude
    --primary_include    or --primary_exclude

    --primary_tag
    --secondary_tag

    --primary_xref
    --primary_gene_xref
    --primary_transcript_xref
    --primary_translation_xref

    --jobs
    --chunk

  For longer explanation of options:
    --help or -h or -?

=head1 OPTIONS

=head2 Input database connection options

=over

=item B<--host_secondary>

The server host name for the Secondary database (required).

=item B<--port_secondary>

The port on the server to connect to.

=item B<--user_secondary>

The username to connect as (required).

=item B<--pass_secondary> or B<--password_secondary>

The password to authenticate with (optional, default is an empty
string).

=item B<--dbname_secondary> or B<--database_secondary>

The name of the Secondary database (required).

=item B<--host_primary>

The server host name for the Primary database (required).

=item B<--port_primary>

The port on the server to connect to.

=item B<--user_primary>

The username to connect as (required).

=item B<--pass_primary> or B<--password_primary>

The password to authenticate with (optional, default is an empty
string).

=item B<--dbname_primary> or B<--database_primary>

The name of the Primary database (required).

=back

=head2 Connection options for using a separate DNA database

These options are all optional and their default values are taken from
the values specified for the connection options for the Secondary input
database.

=over

=item B<--host_dna>

The server host name for the DNA database (optional).

=item B<--port_dna>

The port on the server to connect to (optional).

=item B<--user_dna>

The username to connect as (optional).

=item B<--pass_dna> or B<--password_dna>

The password to authenticate with (optional).

=item B<--dbname_dna> or B<--database_dna>

The name of the DNA database (optional).

=back

=head2 Output database connection options

=over

=item B<--host_output>

The server host name for the output database (required).

=item B<--port_output>

The port on the server to connect to.

=item B<--user_output>

The username to connect as (required).

=item B<--pass_output> or B<--password_output>

The password to authenticate with (optional, default is an empty
string).

=item B<--dbname_output> or B<--database_output>

The name of the output database (required).

=back

=head2 Filter options

=over

=item B<--secondary_include> or B<--secondary_exclude>

A comma-separated list of analysis logic names of genes that will be
included or excluded from the Secondary gene set.  You may use one of
these options, but not both.

=item B<--primary_include> or B<--primary_exclude>

A comma-separated list of analysis logic names of genes that will be
included or excluded from the Primary gene set.  You may use one of these
options, but not both.

=back

=head2 Options related to tagging

=over

=item B<--secondary_tag>

=item B<--primary_tag>

These options allows you to set the tags used for tagging analysis
logic names for Primary and Secondary genes, transcripts, exons and
translations, as well as the logic names assuciated with the various
types of supporting featores.

A lugic name will be suffixed by C<_> followed by the value of the tag,
so that the original logic name C<otter>, for example, would be extended
to C<otter_primary> (if ths was a lgic name that was attached to Primary
annotation).

The tag values are also used to set the C<source> of the genes and
transcrpts written to the output database.  Merged genes and transcripts
gets a combined C<source> value, e.g. C<secondary_primary>.

These options are optional.  The default values for these options are
C<primary> for B<--primary_tag> and C<secondary> for B<--secondary_tag>.

=back

=head2 Options related to Xrefs

=over

=item B<--primary_xref>

=item B<--primary_gene_xref>

=item B<--primary_transcript_xref>

=item B<--primary_translation_xref>

The values corresponding to these options should be a comma-separated
list of strings that makes up the C<external_db_name>,
C<db_display_name>, and C<type> of the xref database table.  This
information will be supplemented with the Primary stable IDs for genes,
transcripts and translations and added to the output gene set.

For example, the standard xrefs for Primary are specified as

  --primary_xref=1
  --primary_gene_xref='OTTG,Havana gene,ALT_GENE'
  --primary_transcript_xref='OTTT,Havana transcript,ALT_TRANS'
  --primary_translation_xref='OTTP,Havana translation,MISC'

These options are optional.  The default for these four options is to not add
any xref (--primary_xref=0). If you wish to add xrefs, set primary_xref to 1 and the other three
parameters to their corresponding values.

=back

=head2 LSF job related options

=over

=item B<--jobs>

This specifies the total number of jobs running this script in parallel
(optional, default is 1).

=item B<--chunk>

This specifies the piece of work that should be carried out by this
particular process (optional, must be between 1 and whatever value
specified by the B<--jobs> option, default is 1).

=back

=head2 Other options

=over

=item B<--help> or B<-h> or B<-?>

Displays this helpful text (optional).

=back

=cut
