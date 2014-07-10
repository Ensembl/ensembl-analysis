#!/usr/bin/env perl
# $Id: merge.pl,v 1.64 2014-01-29 14:57:41 ak4 Exp $

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

my ( $opt_host_ensembl, $opt_port_ensembl,
     $opt_user_ensembl, $opt_password_ensembl,
     $opt_database_ensembl );
my ( $opt_host_havana,     $opt_port_havana, $opt_user_havana,
     $opt_password_havana, $opt_database_havana );
my ( $opt_host_dna,     $opt_port_dna, $opt_user_dna,
     $opt_password_dna, $opt_database_dna );
my ( $opt_host_ccds,     $opt_port_ccds, $opt_user_ccds,
     $opt_password_ccds, $opt_database_ccds );
my ( $opt_host_output,     $opt_port_output, $opt_user_output,
     $opt_password_output, $opt_database_output );

my ( @opt_havana_include,  @opt_havana_exclude );
my ( @opt_ensembl_include, @opt_ensembl_exclude );

my $opt_havana_tag  = 'havana';
my $opt_ensembl_tag = 'ensembl';

my $opt_havana_gene_xref        = 'OTTG,Havana gene,ALT_GENE';
my $opt_havana_transcript_xref  = 'OTTT,Havana transcript,ALT_TRANS';
my $opt_havana_translation_xref = 'OTTP,Havana translation,MISC';

$opt_port_ensembl = $opt_port_havana = $opt_port_dna = $opt_port_ccds =
  $opt_port_output = 3306;

my $opt_njobs = 1;    # Default number of jobs.
my $opt_job   = 1;    # This job.

my $opt_help = 0;

if ( !GetOptions(
          'host_ensembl:s'                    => \$opt_host_ensembl,
          'port_ensembl:i'                    => \$opt_port_ensembl,
          'user_ensembl:s'                    => \$opt_user_ensembl,
          'password_ensembl|pass_ensembl:s'   => \$opt_password_ensembl,
          'database_ensembl|dbname_ensembl:s' => \$opt_database_ensembl,
          'host_havana:s'                     => \$opt_host_havana,
          'port_havana:i'                     => \$opt_port_havana,
          'user_havana:s'                     => \$opt_user_havana,
          'password_havana|pass_havana:s'     => \$opt_password_havana,
          'database_havana|dbname_havana:s'   => \$opt_database_havana,
          'host_dna:s'                        => \$opt_host_dna,
          'port_dna:i'                        => \$opt_port_dna,
          'user_dna:s'                        => \$opt_user_dna,
          'password_dna|pass_dna:s'           => \$opt_password_dna,
          'database_dna|dbname_dna:s'         => \$opt_database_dna,
          'host_ccds:s'                       => \$opt_host_ccds,
          'port_ccds:i'                       => \$opt_port_ccds,
          'user_ccds:s'                       => \$opt_user_ccds,
          'password_ccds|pass_ccds:s'         => \$opt_password_ccds,
          'database_ccds|dbname_ccds:s'       => \$opt_database_ccds,
          'host_output:s'                     => \$opt_host_output,
          'port_output:i'                     => \$opt_port_output,
          'user_output:s'                     => \$opt_user_output,
          'password_output|pass_output:s'     => \$opt_password_output,
          'database_output|dbname_output:s'   => \$opt_database_output,
          'ensembl_include:s'                 => \@opt_ensembl_include,
          'ensembl_exclude:s'                 => \@opt_ensembl_exclude,
          'havana_include:s'                  => \@opt_havana_include,
          'havana_exclude:s'                  => \@opt_havana_exclude,
          'havana_tag:s'                      => \$opt_havana_tag,
          'ensembl_tag:s'                     => \$opt_ensembl_tag,
          'havana_gene_xref:s'                => \$opt_havana_gene_xref,
          'havana_transcript_xref:s'  => \$opt_havana_transcript_xref,
          'havana_translation_xref:s' => \$opt_havana_translation_xref,
          'njobs:i'                   => \$opt_njobs,
          'job:i'                     => \$opt_job,
          'help|h|?!'                 => \$opt_help, ) ||
     $opt_help ||
     !( defined($opt_host_ensembl) &&
        defined($opt_user_ensembl) &&
        defined($opt_database_ensembl) )
     ||
     !( defined($opt_host_havana) &&
        defined($opt_user_havana) &&
        defined($opt_database_havana) )
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

    if ( !( defined($opt_host_ensembl) &&
            defined($opt_user_ensembl) &&
            defined($opt_database_ensembl) ) )
    {
      die( 'Need connection parameters for Ensembl database ' .
           '(host_ensembl, user_ensembl and database_ensembl)' );
    }
    elsif ( !( defined($opt_host_havana) &&
               defined($opt_user_havana) &&
               defined($opt_database_havana) ) )
    {
      die( 'Need connection parameters for Havana database ' .
           '(host_havana, user_havana and database_havana)' );
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
    elsif ( ( @opt_ensembl_include && @opt_ensembl_exclude ) ||
            ( @opt_havana_include && @opt_havana_exclude ) )
    {
      die('You may only use X_include or X_exclude, but not both');
    }
    else {
      die('Error in command line parsing');
    }
  } ## end else [ if ($opt_help) ]

} ## end if ( !GetOptions( 'host_ensembl:s'...))

my $dna_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                '-no_cache' => 1,
                '-host'     => $opt_host_dna || $opt_host_ensembl,
                '-port'     => $opt_port_dna || $opt_port_ensembl,
                '-user'     => $opt_user_dna || $opt_user_ensembl,
                '-pass' => $opt_password_dna || $opt_password_ensembl,
                '-dbname' => $opt_database_dna || $opt_database_ensembl,
  ) or
  die('Failed to connect to DNA database');

my $ensembl_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                     '-no_cache' => 1,
                                     '-host'     => $opt_host_ensembl,
                                     '-port'     => $opt_port_ensembl,
                                     '-user'     => $opt_user_ensembl,
                                     '-pass'   => $opt_password_ensembl,
                                     '-dbname' => $opt_database_ensembl,
                                     '-dnadb'  => $dna_dba, ) or
  die('Failed to connect to Ensembl database');

my $havana_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new('-no_cache' => 1,
                                      '-host'     => $opt_host_havana,
                                      '-port'     => $opt_port_havana,
                                      '-user'     => $opt_user_havana,
                                      '-pass'   => $opt_password_havana,
                                      '-dbname' => $opt_database_havana,
                                      '-dnadb'  => $dna_dba, ) or
  die('Failed to connect to Havana database');

my $output_dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new('-no_cache' => 1,
                                      '-host'     => $opt_host_output,
                                      '-port'     => $opt_port_output,
                                      '-user'     => $opt_user_output,
                                      '-pass'   => $opt_password_output,
                                      '-dbname' => $opt_database_output,
  ) or
  die('Failed to connect to output database');

my $ccds_dba;
if ( defined($opt_host_ccds) &&
     defined($opt_user_ccds)     &&
     defined($opt_database_ccds) &&
     $opt_database_ccds ne '' )
{
  $ccds_dba =
    Bio::EnsEMBL::DBSQL::DBAdaptor->new('-no_cache' => 1,
                                        '-host'     => $opt_host_ccds,
                                        '-port'     => $opt_port_ccds,
                                        '-user'     => $opt_user_ccds,
                                        '-pass'   => $opt_password_ccds,
                                        '-dbname' => $opt_database_ccds,
    ) or
    die('Failed to connect to CCDS database');
}
else {
  print("Not attaching CCDS database\n");
}

@opt_ensembl_include = split( /,/, join( ',', @opt_ensembl_include ) );
@opt_ensembl_exclude = split( /,/, join( ',', @opt_ensembl_exclude ) );
@opt_havana_include  = split( /,/, join( ',', @opt_havana_include ) );
@opt_havana_exclude  = split( /,/, join( ',', @opt_havana_exclude ) );


if($opt_database_dna) {

  print "Optional DNA database\thost:\t".$opt_host_dna."\n".
                                   "\tport:\t".$opt_port_dna."\n".
                                   "\tuser:\t".$opt_user_dna."\n".
                                   "\tname:\t".$opt_database_dna."\n";
}


print <<DBINFO_END;
ENSEMBL database\thost:\t$opt_host_ensembl
                \tport:\t$opt_port_ensembl
                \tuser:\t$opt_user_ensembl
                \tname:\t$opt_database_ensembl

HAVANA database\thost:\t$opt_host_havana
               \tport:\t$opt_port_havana
               \tuser:\t$opt_user_havana
               \tname:\t$opt_database_havana

CCDS database\thost:\t$opt_host_ccds
               \tport:\t$opt_port_ccds
               \tuser:\t$opt_user_ccds
               \tname:\t$opt_database_ccds

OUTPUT database\thost:\t$opt_host_output
               \tport:\t$opt_port_output
               \tuser:\t$opt_user_output
               \tname:\t$opt_database_output

Ensembl logic name filter (include): @opt_ensembl_include
Ensembl logic name filter (exclude): @opt_ensembl_exclude

Havana logic name filter (include): @opt_havana_include
Havana logic name filter (exclude): @opt_havana_exclude

Ensembl tag:\t$opt_ensembl_tag
Havana tag:\t$opt_havana_tag

Havana gene xref:       \t$opt_havana_gene_xref
Havana transcript xref: \t$opt_havana_transcript_xref
Havana translation xref:\t$opt_havana_translation_xref

------------------------------------------------------------------------
DBINFO_END

my $ENSEMBL_GA = $ensembl_dba->get_GeneAdaptor();    # Used globally.
my $HAVANA_GA  = $havana_dba->get_GeneAdaptor();     # Used globally.
my $OUTPUT_GA  = $output_dba->get_GeneAdaptor();     # This one too.
my $CCDS_TA;
if ( defined($ccds_dba) ) {
  $CCDS_TA = $ccds_dba->get_TranscriptAdaptor();     # Ditto.
}

# Create a chunk of work, i.e. a list of gene IDs.  We create uniformly
# sized chunks whose size depend on the number of available Havana
# genes and the number of jobs that are being run.  We pick the chunk
# associated with our job ID.

my @all_havana_gene_ids =
  sort { $a <=> $b } @{ $HAVANA_GA->list_dbIDs() };

my $chunk_size  = int( scalar(@all_havana_gene_ids)/$opt_njobs );
my $chunk_first = ( $opt_job - 1 )*$chunk_size;
my $chunk_last  = $chunk_first + $chunk_size - 1;

printf( "Got %d genes, with %d jobs this means %d genes per job",
        scalar(@all_havana_gene_ids),
        $opt_njobs, $chunk_size );

# This ($chunk_spill) is just the number of jobs that will get one extra
# piece of work to do, because $chunk_size multiplied by $opt_njobs
# isn't quite equal to the number of Havana genes.

my $chunk_spill = scalar(@all_havana_gene_ids) % $opt_njobs;

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

my %havana_genes_done;

# Process all the Havana genes in the current chunk.
foreach my $havana_gene_id (
                   @all_havana_gene_ids[ $chunk_first .. $chunk_last ] )
{
  if ( exists( $havana_genes_done{$havana_gene_id} ) ) {
    printf( "Skipping gene %d, already processed " .
              "(because of clustering).\n",
            $havana_gene_id );
    next;
  }

  # Given the gene with ID $havana_gene_id, pick out all the overlapping
  # Havana and Ensembl genes (will return nothing if the gene cluster is
  # found to contain genes belonging to another LSF job).

  my ( $havana_genes, $ensembl_genes ) =
    @{
    make_gene_cluster( $havana_gene_id,
                       $all_havana_gene_ids[$chunk_first] ) };

  # Filter the genes.
  my @filtered_havana_genes;
  my @filtered_ensembl_genes;
  foreach my $havana_gene ( @{$havana_genes} ) {
    if ( filter_gene( $havana_gene, \@opt_havana_include,
                      \@opt_havana_exclude ) )
    {
      next;
    }
    push( @filtered_havana_genes, $havana_gene );
  }
  foreach my $ensembl_gene ( @{$ensembl_genes} ) {
    if ( filter_gene( $ensembl_gene, \@opt_ensembl_include,
                      \@opt_ensembl_exclude ) )
    {
      next;
    }
    push( @filtered_ensembl_genes, $ensembl_gene );
  }

  # Process the filtered genes.

  # This loop needs to come before the call to process_genes() since
  # the dbIDs are changing when storing them.  It should loop over
  # @{$havana_genes}, not @filtered_havana_genes.
  foreach my $havana_gene ( @{$havana_genes} ) {
    $havana_genes_done{ $havana_gene->dbID() } = 1;
  }

  process_genes( \@filtered_havana_genes, \@filtered_ensembl_genes );

  print("==\n");
} ## end foreach my $havana_gene_id ...

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
  # the lowest Havana gene ID in the cluster is greater or equal to
  # $lowest_allowed_gene_id, then process it, otherwise skip this
  # cluster as it belongs to another LSF job.

  my $seed_gene = $HAVANA_GA->fetch_by_dbID($seed_gene_id);

  printf( "Initiating gene cluster with Havana gene %s (%d), %s\n",
          $seed_gene->stable_id(),
          $seed_gene_id, $seed_gene->feature_Slice()->name() );

  my @cluster_queue = ($seed_gene);
  my %used_slices;

  my %havana_cluster;
  my %ensembl_cluster;

  my $havana_sa  = $HAVANA_GA->db()->get_SliceAdaptor();
  my $ensembl_sa = $ENSEMBL_GA->db()->get_SliceAdaptor();

  # The cluster queue is a collection of possibly unprocessed genes in
  # this cluster of overlapping genes.  We use each gene in turn to
  # fetch overlapping gene from the Havana and the Ensembl databases.
  # Any found gene, if it hasn't already been seen, is added to the
  # cluster queue.  When thu queue is empty, we are done.

  while ( my $gene = shift(@cluster_queue) ) {
    my $havana_slice =
      get_feature_slice_from_db( $gene, $HAVANA_GA->db() );

    if ( exists( $used_slices{ 'havana:' . $havana_slice->name() } ) ) {
      printf( "Skipping Havana slice %s, already seen.\n",
              $havana_slice->name() );
      next;
    }

    $used_slices{ 'havana:' . $havana_slice->name() } = 1;

    # Fetch overlapping Havana genes.

    foreach my $havana_gene (
        @{ $HAVANA_GA->fetch_all_by_Slice( $havana_slice, undef, 1 ) } )
    {
      my $add_to_result  = 1;
      my $havana_gene_id = $havana_gene->dbID();

      if ( $havana_gene->length() > 100_000_000 ) {
        printf( "Ignoring Havana gene %s (%d)," .
                  " too long (%d > 100,000,000)\n",
                $havana_gene->stable_id(),
                $havana_gene_id, $havana_gene->length() );
        next;
      }

      if ( !exists( $havana_cluster{$havana_gene_id} ) ) {
        if ( $havana_gene_id < $lowest_allowed_gene_id ) {
          print("Skipping gene cluster, does not belong to me.\n");
          return [ [], [] ];
        }

        push( @cluster_queue, $havana_gene );

        if ($add_to_result) {
          $havana_cluster{$havana_gene_id} = $havana_gene;

          printf( "Also fetched Havana gene %s (%d), %s\n",
                  $havana_gene->stable_id(),
                  $havana_gene_id,
                  $havana_gene->feature_Slice()->name() );
        }
      }
    } ## end foreach my $havana_gene ( @...)

    my $ensembl_slice =
      get_feature_slice_from_db( $gene, $ENSEMBL_GA->db() );

    if ( exists( $used_slices{ 'ensembl:' . $ensembl_slice->name() } ) )
    {
      printf( "Skipping Ensembl slice %s, already seen.\n",
              $ensembl_slice->name() );
      next;
    }

    $used_slices{ 'ensembl:' . $ensembl_slice->name() } = 1;

    # Fetch overlapping Ensembl genes.

    foreach my $ensembl_gene (
      @{ $ENSEMBL_GA->fetch_all_by_Slice( $ensembl_slice, undef, 1 ) } )
    {
      my $add_to_result   = 1;
      my $ensembl_gene_id = $ensembl_gene->dbID();

      if ( $ensembl_gene->length() > 100_000_000 ) {
        printf( "Ignoring Ensembl gene %s (%d), " .
                  " too long (%d > 100,000,000)\n",
                $ensembl_gene->stable_id(),
                $ensembl_gene_id, $ensembl_gene->length() );
        next;
      }

      if ( !exists( $ensembl_cluster{$ensembl_gene_id} ) ) {
        push( @cluster_queue, $ensembl_gene );

        if ($add_to_result) {
          $ensembl_cluster{$ensembl_gene_id} = $ensembl_gene;

          printf( "Also fetched Ensembl gene %s (%d), %s\n",
                  $ensembl_gene->stable_id(),
                  $ensembl_gene_id,
                  $ensembl_gene->feature_Slice()->name() );
        }
      }
    } ## end foreach my $ensembl_gene ( ...)

  } ## end while ( my $gene = shift(...))

  return [ [ values(%havana_cluster) ], [ values(%ensembl_cluster) ] ];
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
         1,                           $slice->coord_system()->version(),
         1 ) };

  if ( scalar(@slices) != 1 ) {
    # This will hopefully only happen if the Havana and Ensembl
    # databases contain different assemblies.
    die( "!! Problem with projection for feature slice %s\n",
         $slice->name() );
  }

  return $slices[0];
}

sub process_genes {
  my ( $havana_genes, $ensembl_genes ) = @_;

  # Start by fully loading each Havana gene.  Also add OTTP, OTTT and
  # OTTG xrefs to the various parts of the gene.  We also figure out
  # (and keep track of) whether the translations are coding, if the
  # genes are pseudogenes, if a gene is a single transcript gene, and
  # whether a gene is located across a region with a known assembly
  # error.

  foreach my $havana_gene ( @{$havana_genes} ) {
    $havana_gene->load();

    my $is_coding = 0;

    my $transcript_count = 0;
    foreach
      my $havana_transcript ( @{ $havana_gene->get_all_Transcripts() } )
    {
      my $havana_translation = $havana_transcript->translation();

      if ( defined($havana_translation) ) {
        $is_coding = 1;

        # Add "OTTP" xref to Havana translation.
        add_havana_xref($havana_translation);
      }

      tag_transcript_analysis( $havana_transcript, $opt_havana_tag );
      $havana_transcript->source($opt_havana_tag);

      # Add "OTTT" xref to Havana transcript.
      add_havana_xref($havana_transcript);

      # If the transcript is part of a gene cluster then tag the gene
      unless($havana_gene->{__is_gene_cluster}) {
        $havana_gene->{__is_gene_cluster} = (
        scalar(@{ $havana_transcript->get_all_Attributes('gene_cluster') }) > 0 );
      }

      ++$transcript_count;
    }

    add_logic_name_suffix( $havana_gene, $opt_havana_tag );
    $havana_gene->source($opt_havana_tag);

    my $has_assembly_error = (
              scalar(
                @{ $havana_gene->get_all_Attributes('NoTransRefError') }
              ) > 0 );

    $havana_gene->{__is_coding} = $is_coding;    # HACK

    $havana_gene->{__is_pseudogene} =            # HACK
      ( !$havana_gene->{__is_coding} &&
        $havana_gene->biotype() =~ /pseudogene/ );

    # Check pseudogene has one transcript  
    if ($havana_gene->{__is_pseudogene}) {

      warning ($transcript_count." transcripts found for Havana pseudogene: ".$havana_gene->stable_id().
      ". Havana pseudogenes should have a single transcript.")
      unless $transcript_count == 1;

    }

    $havana_gene->{__has_ref_error} = $has_assembly_error;    # HACK
    $havana_gene->{__is_single_transcript} =
      ( $transcript_count == 1 );                             # HACK

    # Add "OTTG" xref to Havana gene.
    add_havana_xref($havana_gene);

  } ## end foreach my $havana_gene ( @...)

  # For Ensembl genes, we don't add any xrefs, but we figure out if the
  # transcripts are pseudogene transcripts or part of the CCDS dataset
  # (if made avalable).  We also see if a gene is a single transcript
  # gene, if it's coding and if it's an RNA gene or not.  All these
  # things are being used elsewhere to make decisions about merging and
  # copying.

  foreach my $ensembl_gene ( @{$ensembl_genes} ) {
    $ensembl_gene->load();

    my $is_coding        = 0;
    my $is_pseudogene    = 0;
    my $transcript_count = 0;

    my @ccds_transcripts;
    if ( defined($CCDS_TA) ) {
      my $ccds_slice =
        get_feature_slice_from_db( $ensembl_gene, $CCDS_TA->db() );
      @ccds_transcripts =
        @{ $CCDS_TA->fetch_all_by_Slice( $ccds_slice, 1 ) };
    }

    foreach my $ensembl_transcript (
                             @{ $ensembl_gene->get_all_Transcripts() } )
    {
      $ensembl_transcript->{__is_ccds}       = 0;    # HACK
      $ensembl_transcript->{__is_pseudogene} = 0;    # HACK

      my $ensembl_translation = $ensembl_transcript->translation();

      if ( defined($ensembl_translation) ) {
        if ( scalar(@ccds_transcripts) > 0 ) {
          my @translatable_exons =
            @{ $ensembl_transcript->get_all_translateable_Exons() };

          foreach my $ccds_transcript (@ccds_transcripts) {
            if ( features_are_same( \@translatable_exons,
                                    $ccds_transcript
                                      ->get_all_translateable_Exons( ) )
              )
            {
              $ensembl_transcript->{__is_ccds} = 1;    # HACK
              last;
            }
          }
        }

        $is_coding = 1;
      }
      elsif ( $ensembl_transcript->biotype() =~ /pseudogene/ ) {
        $ensembl_transcript->{__is_pseudogene} = 1;    # HACK
      }

      tag_transcript_analysis( $ensembl_transcript, $opt_ensembl_tag );
      $ensembl_transcript->source($opt_ensembl_tag);

      ++$transcript_count;
    } ## end foreach my $ensembl_transcript...

    add_logic_name_suffix( $ensembl_gene, $opt_ensembl_tag );
    $ensembl_gene->source($opt_ensembl_tag);

    my $is_rna = (
         index( lc( $ensembl_gene->analysis()->logic_name() ), 'ncrna' )
           >= 0 );

    $ensembl_gene->{__is_coding} = $is_coding;    # HACK
    $ensembl_gene->{__is_single_transcript} =
      ( $transcript_count == 1 );                 # HACK
    $ensembl_gene->{__is_rna} = $is_rna;          # HACK
  } ## end foreach my $ensembl_gene ( ...)

  my @transcripts_to_merge;
  my @transcripts_to_copy;
  my @transcripts_to_ignore;

  # This will keep a record of the Ensembl gene stable id whenever stuff is
  # merged or copied in
  my %processed_ens_id_tracker;

ENSEMBL_GENE:
  foreach my $ensembl_gene ( @{$ensembl_genes} ) {
    printf( "Ensembl gene: %s (%d)\n",
            $ensembl_gene->stable_id(), $ensembl_gene->dbID() );

    my @ensembl_transcripts = @{ $ensembl_gene->get_all_Transcripts() };

  ENSEMBL_TRANSCRIPT:
    foreach my $ensembl_transcript (@ensembl_transcripts) {
      printf( "\tEnsembl transcript %s (%d)\n",
              $ensembl_transcript->stable_id(),
              $ensembl_transcript->dbID() );

    HAVANA_GENE:
      foreach my $havana_gene ( @{$havana_genes} ) {
        printf( "\t\tHavana gene: %s (%d)\n",
                $havana_gene->stable_id(), $havana_gene->dbID() );

        if ( $havana_gene->strand() != $ensembl_gene->strand() ) {
          printf("\t\t\tNot on same strand\n");
          next HAVANA_GENE;
        }

        my @havana_transcripts =
          @{ $havana_gene->get_all_Transcripts() };

      HAVANA_TRANSCRIPT:
        foreach my $havana_transcript (@havana_transcripts) {
          printf( "\t\t\tHavana transcript %s (%d)\n",
                  $havana_transcript->stable_id(),
                  $havana_transcript->dbID() );

          if ( investigate_for_merge(
                                 $havana_transcript, $ensembl_transcript
               ) )
          {
            print("\t\t\tTo be merged\n");

            push( @transcripts_to_merge,
                  [ $havana_gene,  $havana_transcript,
                    $ensembl_gene, $ensembl_transcript ] );
          }
          elsif ( investigate_for_copy(
                                      $havana_gene,  $havana_transcript,
                                      $ensembl_gene, $ensembl_transcript
                  ) )
          {
            print( "\t\t\tEnsembl transcript overlaps " .
                   "and may be copied\n" );

            push( @transcripts_to_copy,
                  [ $havana_gene,  $havana_transcript,
                    $ensembl_gene, $ensembl_transcript ] );
          }

        } ## end HAVANA_TRANSCRIPT: foreach my $havana_transcript...
      } ## end HAVANA_GENE: foreach my $havana_gene ( @...)

    } ## end ENSEMBL_TRANSCRIPT: foreach my $ensembl_transcript...
  } ## end ENSEMBL_GENE: foreach my $ensembl_gene ( ...)

  my %merged_ensembl_transcripts;
  my %copied_ensembl_transcripts;
  my %ignored_ensembl_transcripts;

  # Merge the transcripts that have been found to be mergable.  We allow
  # merging an Ensembl transcript into multiple Havana transcripts.
  foreach my $tuple (@transcripts_to_merge) {
    my ( $havana_gene,  $havana_transcript,
         $ensembl_gene, $ensembl_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    printf(
          "Merging %s (%d) into %s (%d) / %s (%d)\n",
          $ensembl_transcript->stable_id(), $ensembl_transcript->dbID(),
          $havana_transcript->stable_id(),  $havana_transcript->dbID(),
          $havana_gene->stable_id(),        $havana_gene->dbID() );

    my $merge_code =
      merge( $havana_gene, $havana_transcript, $ensembl_transcript );

    $merged_ensembl_transcripts{ $ensembl_transcript->dbID() } = 1;

    if ( $merge_code != 1 ) {
      printf( "PROCESSED\t%d\t%s\n",
              $ensembl_gene->dbID(), $ensembl_gene->stable_id() );

      # Keep a record in memory of this Ensembl stable id
      $processed_ens_id_tracker{$havana_gene->dbID()} = [$ensembl_gene->dbID(),$ensembl_gene->stable_id()];
    }
  }

  # Find all Ensembl transcripts that are "siblings" to any merged
  # Ensembl transcript.  Tag these for possibly copying into the
  # corresponding Havana gene.
  foreach my $tuple (@transcripts_to_merge) {
    my ( $havana_gene,  $havana_transcript,
         $ensembl_gene, $ensembl_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    foreach my $transcript ( @{ $ensembl_gene->get_all_Transcripts() } )
    {
      my $key = $transcript->dbID();

      if ( !exists( $merged_ensembl_transcripts{$key} ) ) {
        printf( "May also copy %s (%d) into %s (%d) " .
                  "due to merging of %s (%d)\n",
                $transcript->stable_id(),
                $transcript->dbID(),
                $havana_gene->stable_id(),
                $havana_gene->dbID(),
                $ensembl_transcript->stable_id(),
                $ensembl_transcript->dbID() );

        push( @transcripts_to_copy,
              [ $havana_gene,  $havana_transcript,
                $ensembl_gene, $transcript ] );
      }
    }
  } ## end foreach my $tuple (@transcripts_to_merge)

  # Take care of the cases where an Ensembl transcript has been selected
  # to be copied into multiple Havana genes.
  my %distinct_transcripts_to_copy;
  foreach my $tuple (@transcripts_to_copy) {
    my ( $havana_gene,  $havana_transcript,
         $ensembl_gene, $ensembl_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $ensembl_transcript->dbID();

    if ( exists( $merged_ensembl_transcripts{$key} ) ) {
      next;
    }

    if ( exists( $distinct_transcripts_to_copy{$key} ) ) {
      my $is_conflicting = 0;

      if ( $havana_gene->dbID() !=
           $distinct_transcripts_to_copy{$key}[0]->dbID() )
      {
        $is_conflicting = 1;
        printf( "Ensembl transcript %s (%d) " .
                  "selected for copy into both " .
                  "%s (%d) and %s (%d)\n",
                $ensembl_transcript->stable_id(),
                $ensembl_transcript->dbID(),
                $havana_gene->stable_id(),
                $havana_gene->dbID(),
                $distinct_transcripts_to_copy{$key}[0]->stable_id(),
                $distinct_transcripts_to_copy{$key}[0]->dbID() );
      }

      if ( $havana_transcript->dbID() !=
           $distinct_transcripts_to_copy{$key}[1]->dbID() )
      {
        my $total_overlap =
          exon_overlap( $ensembl_transcript, $havana_transcript );

        if ( $total_overlap > $distinct_transcripts_to_copy{$key}[4] ) {
          if ($is_conflicting) {
            printf( "Will choose to copy %s (%d) into %s (%d) " .
                      "due to larger exon overlap (%dbp > %dbp)\n",
                    $ensembl_transcript->stable_id(),
                    $ensembl_transcript->dbID(),
                    $havana_gene->stable_id(),
                    $havana_gene->dbID(),
                    $total_overlap,
                    $distinct_transcripts_to_copy{$key}[4] );
          }

          $distinct_transcripts_to_copy{$key} = [
                                     $havana_gene,  $havana_transcript,
                                     $ensembl_gene, $ensembl_transcript,
                                     $total_overlap ];
        }
        elsif ($is_conflicting) {
          printf( "Will choose to copy %s (%d) into %s (%d) " .
                    "due to larger exon overlap (%dbp > %dbp)\n",
                  $ensembl_transcript->stable_id(),
                  $ensembl_transcript->dbID(),
                  $distinct_transcripts_to_copy{$key}[0]->stable_id(),
                  $distinct_transcripts_to_copy{$key}[0]->dbID(),
                  $distinct_transcripts_to_copy{$key}[4],
                  $total_overlap );
        }
      } ## end if ( $havana_transcript...)

    } ## end if ( exists( $distinct_transcripts_to_copy...))
    else {
      $distinct_transcripts_to_copy{$key} = [
                 $havana_gene,
                 $havana_transcript,
                 $ensembl_gene,
                 $ensembl_transcript,
                 exon_overlap( $ensembl_transcript, $havana_transcript )
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
    my ( $havana_gene,  $havana_transcript,
         $ensembl_gene, $ensembl_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $ensembl_transcript->dbID();

    if ( !exists( $merged_ensembl_transcripts{$key} ) &&
         !exists( $copied_ensembl_transcripts{$key} ) )
    {
      if ( !defined( $havana_transcript->translation() ) ||
           !defined( $ensembl_transcript->translation() ) ||
           has_complete_start_stop($ensembl_transcript) )
      {
        printf( "Copying %s (%d) into %s (%d)\n",
                $ensembl_transcript->stable_id(),
                $ensembl_transcript->dbID(),
                $havana_gene->stable_id(),
                $havana_gene->dbID() );

        my $copy_code = copy( $havana_gene, $ensembl_transcript );

        $copied_ensembl_transcripts{$key} = 1;

        if ( $copy_code != 1 ) {
          printf( "PROCESSED\t%d\t%s\n",
                  $ensembl_gene->dbID(), $ensembl_gene->stable_id() );

          # Keep a record in memory of this stable id
          $processed_ens_id_tracker{$havana_gene->dbID()} = [$ensembl_gene->dbID(),$ensembl_gene->stable_id()];
        }
      }
      else {
        printf( "May ignore %s (%d) due to incomplete start/stop\n",
                $ensembl_transcript->stable_id(),
                $ensembl_transcript->dbID() );

        push( @transcripts_to_ignore,
              [ $havana_gene,  $havana_transcript,
                $ensembl_gene, $ensembl_transcript ] );
      }
    } ## end if ( !exists( $merged_ensembl_transcripts...))
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
    my ( $havana_gene,  $havana_transcript,
         $ensembl_gene, $ensembl_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $ensembl_transcript->dbID();

    if ( !exists( $merged_ensembl_transcripts{$key} ) &&
         !exists( $copied_ensembl_transcripts{$key} ) &&
         !exists( $ignored_ensembl_transcripts{$key} ) )
    {
      my $do_ignore = 0;

      foreach
        my $transcript ( @{ $ensembl_gene->get_all_Transcripts() } )
      {

        my $key2 = $transcript->dbID();

        if ( $key ne $key2 ) {
          if ( exists( $copied_ensembl_transcripts{$key2} ) ||
               exists( $merged_ensembl_transcripts{$key2} ) )
          {
            $do_ignore = 1;
            last;
          }
        }
      }

      my $copy_code = 0;
      if ($do_ignore) {
        printf( "Ignoring %s (%d) due to incomplete start/stop\n",
                $ensembl_transcript->stable_id(),
                $ensembl_transcript->dbID() );

        $ignored_ensembl_transcripts{$key} = 1;
      }
      else {
        printf( "Copying %s (%d) into %s (%d) " .
                  "even though it has incomplete start/stop\n",
                $ensembl_transcript->stable_id(),
                $ensembl_transcript->dbID(),
                $havana_gene->stable_id(),
                $havana_gene->dbID() );

        $copy_code = copy( $havana_gene, $ensembl_transcript );
      }

      if ( $copy_code != 1 ) {
        printf( "PROCESSED\t%d\t%s\n",
                $ensembl_gene->dbID(), $ensembl_gene->stable_id() );

        # Keep a record in memory of this stable id
        $processed_ens_id_tracker{$havana_gene->dbID()} = [$ensembl_gene->dbID(),$ensembl_gene->stable_id()];
      }

    } ## end if ( !exists( $merged_ensembl_transcripts...))
  } ## end foreach my $tuple ( sort { ...})

  # Do a bit of houskeeping from the previous loop so that
  # %copied_ensembl_transcripts correctly reflects the Ensembl
  # transcripts actually copied.
  foreach my $tuple (@transcripts_to_ignore) {
    my ( $havana_gene,  $havana_transcript,
         $ensembl_gene, $ensembl_transcript
    ) = ( $tuple->[0], $tuple->[1], $tuple->[2], $tuple->[3] );

    my $key = $ensembl_transcript->dbID();

    if ( !exists( $ignored_ensembl_transcripts{$key} ) ) {
      $copied_ensembl_transcripts{$key} = 1;
    }
  }

  foreach my $havana_gene ( @{$havana_genes} ) {

    my $old_dbID = $havana_gene->dbID();

    # If the biotype is bad, skip the store
    if(bad_biotype($havana_gene->biotype())) {

        print "WARNING: skipping store of Havana gene ".$havana_gene->dbID()." ".
              $havana_gene->stable_id()." due to bad biotype: ".$havana_gene->biotype()."\n";

        # If there is an Ensembl gene that has been merged/copied to this Havana gene
        # put in some additional warnings that the gene will be thrown out
        if($processed_ens_id_tracker{$old_dbID}) {         
          my $ens_gene_db_id = $processed_ens_id_tracker{$old_dbID}->[0];
          my $ens_gene_stable_id = "";
          if($processed_ens_id_tracker{$old_dbID}->[1]) {
            $ens_gene_stable_id = $processed_ens_id_tracker{$old_dbID}->[1];
          } 

          print "WARNING: skipping store of Ensembl gene ".$ens_gene_db_id." ".$ens_gene_stable_id.
                " due to merge/copy with Havana gene ".$havana_gene->dbID()." ".$havana_gene->stable_id().
                " with bad biotype: ".$havana_gene->biotype()."\n";          
        }

       next;
     }

    $OUTPUT_GA->store($havana_gene);
    printf( "STORED\t%s\told id = %d, new id = %d\n",
            $havana_gene->stable_id(),
            $old_dbID, $havana_gene->dbID() );

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
                       'TEC' => 1,   
                       'artifact' => 1
                     );

  if($bad_biotypes{$biotype}) {
    return 1; 
  }

  else {
    return 0;
  }
  
} # End sub bad_biotype


sub add_havana_xref {
  my ($feature) = @_;

  my ( $external_db_name, $db_display_name, $type );

  if ( $feature->isa('Bio::EnsEMBL::Gene') ) {
    ( $external_db_name, $db_display_name, $type ) =
      split( /,/, $opt_havana_gene_xref );

    # $external_db_name = 'OTTG';
    # $db_display_name  = 'Havana gene';
    # $type             = 'ALT_GENE';
  }
  elsif ( $feature->isa('Bio::EnsEMBL::Transcript') ) {
    ( $external_db_name, $db_display_name, $type ) =
      split( /,/, $opt_havana_transcript_xref );

    # $external_db_name = 'OTTT';
    # $db_display_name  = 'Havana transcript';
    # $type             = 'ALT_TRANS';
  }
  elsif ( $feature->isa('Bio::EnsEMBL::Translation') ) {
    ( $external_db_name, $db_display_name, $type ) =
      split( /,/, $opt_havana_translation_xref );

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
          "havana_gene_xref='OTTG,Havana gene,ALT_GENE'\n".
          "havana_transcript_xref='OTTT,Havana transcript,ALT_TRANS'\n".
          "havana_translation_xref='OTTP,Havana translation,MISC'\n"
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
} ## end sub add_havana_xref

sub investigate_for_merge {
  my ( $havana_transcript, $ensembl_transcript ) = @_;

  my @havana_introns = @{ $havana_transcript->get_all_Introns() };

  my $do_merge;

  if ( scalar(@havana_introns) == 0 ) {
    # This is a single exon transcript, compare the exons instead of
    # the (non-existing) introns.

    $do_merge = features_are_same( $havana_transcript->get_all_Exons(),
                                   $ensembl_transcript->get_all_Exons()
    );

    if ( !$do_merge &&
         ( defined( $havana_transcript->translation() ) &&
           defined( $ensembl_transcript->translation() ) ) )
    {
      # This is a single exon coding transcript, compare the coding
      # exons instead of the (non-existing) introns.

      $do_merge = features_are_same(
                   $havana_transcript->get_all_translateable_Exons(),
                   $ensembl_transcript->get_all_translateable_Exons() );

      if ($do_merge) {
        print( "\t\t\tSpecial case: Single exon coding transcripts, " .
               "coding region is identical but UTR differs\n" );
      }
    }

    if ( !$do_merge ) {
      # A special case:  If the Havana exon is exactly one stop
      # codon longer at the end, then treat them as identical.
      $do_merge = is_a_stop_codon_longer( $havana_transcript,
                                          $ensembl_transcript );
      if ($do_merge) {
        print( "\t\t\t\tSpecial case: Single exon transcripts, " .
               "Ensembl exon is stop codon short\n" );
      }
    }
    if ( !$do_merge ) {
      # A special case:  If the Ensembl exon is exactly one stop
      # codon longer at the end, then treat them as identical.
      $do_merge = is_a_stop_codon_longer( $ensembl_transcript,
                                          $havana_transcript );
      if ($do_merge) {
        print( "\t\t\t\tSpecial case: Single exon transcripts, " .
               "Havana exon is stop codon short\n" );
      }
    }
  } ## end if ( scalar(@havana_introns...))
  else {
    my @ensembl_introns = @{ $ensembl_transcript->get_all_Introns() };

    $do_merge =
      features_are_same( \@havana_introns, \@ensembl_introns );
  }

  return $do_merge;
} ## end sub investigate_for_merge

sub investigate_for_copy {
  my ( $havana_gene,  $havana_transcript,
       $ensembl_gene, $ensembl_transcript ) = @_;

  my @ensembl_exons = @{ $ensembl_transcript->get_all_Exons() };
  my @havana_exons  = @{ $havana_transcript->get_all_Exons() };

  my $do_copy = 0;

ENSEMBL_EXON:
  foreach my $ensembl_exon (@ensembl_exons) {
    foreach my $havana_exon (@havana_exons) {
      if ( features_overlap( $ensembl_exon, $havana_exon ) ) {
        $do_copy = 1;

        if ( scalar(@ensembl_exons) == 1 &&
             $ensembl_gene->{__is_single_transcript} &&
             scalar(@havana_exons) > 1 )
        {
          # We have a single exon Ensembl gene that overlaps with bits
          # of a multi exon (and possibly multi transcript) Havana gene.
          #
          # The current rules are:
          #
          # * Copy if the overlapping Havana exon is identical to the
          #   Ensembl exon.
          #
          # * Don't copy if the Ensembl gene is RNA and the Havana gene
          #   has a different biotype.
          #
          # * Copy in any other case.

          if ( features_are_same( [$ensembl_exon], [$havana_exon] ) ) {
            print( "\t\t\tSpecial case: Allowing copy of transcript " .
                   "from single exon gene into multi exon gene " .
                   "(perfect exon match)\n" );
          }
          elsif ( $ensembl_gene->{__is_rna} &&
                  $ensembl_gene->biotype() ne $havana_gene->biotype() )
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
        } ## end if ( scalar(@ensembl_exons...))
        elsif ( $ensembl_gene->{__is_rna} &&
                $ensembl_gene->biotype() ne $havana_gene->biotype() )
        {
          # Don't allow copying of RNA gene transcripts if gene biotypes
          # do not match upi (regardless of exon counts).
          print( "\t\t\tSpecial case: Won't copy transcript " .
                 "from RNA gene into gene with different biotype\n" );
          $do_copy = 0;
        }

        if ($do_copy) {
          last ENSEMBL_EXON;
        }
      } ## end if ( features_overlap(...))
    } ## end foreach my $havana_exon (@havana_exons)
  } ## end ENSEMBL_EXON: foreach my $ensembl_exon (@ensembl_exons)

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

sub features_overlap {
  my ( $featureA, $featureB ) = @_;

  if (
     ( $featureA->seq_region_start() <= $featureB->seq_region_end() ) &&
     ( $featureA->seq_region_end() >= $featureB->seq_region_start() ) )
  {
    return 1;
  }

  return 0;
}

sub overlap_length {
  my ( $featureA, $featureB ) = @_;

  if ( !features_overlap( $featureA, $featureB ) ) {
    return 0;
  }

  my $min_end = $featureA->seq_region_end();
  if ( $featureB->seq_region_end() < $min_end ) {
    $min_end = $featureB->seq_region_end();
  }

  my $max_start = $featureA->seq_region_start();
  if ( $featureB->seq_region_start() > $max_start ) {
    $max_start = $featureB->seq_region_start();
  }

  return $min_end - $max_start + 1;
}

sub exon_overlap {
  my ( $trancriptA, $trancriptB ) = @_;

  my $overlap = 0;

  foreach my $exonA ( @{ $trancriptA->get_all_Exons() } ) {
    foreach my $exonB ( @{ $trancriptB->get_all_Exons() } ) {
      $overlap += overlap_length( $exonA, $exonB );
    }
  }

  return $overlap;
}

sub merge {
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
    my @supporting_features =
      @{ $new_source_transcript->get_all_supporting_features() };

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
      my $source_exon ( @{ $new_source_transcript->get_all_Exons() } )
    {
      push( @supporting_features,
            [ @{ $source_exon->get_all_supporting_features() } ] );
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

  add_logic_name_suffix( $target_transcript, 'merged' );
  add_logic_name_suffix( $target_gene,       'merged' );
  $target_gene->source( $opt_ensembl_tag . '_' . $opt_havana_tag );
  $target_transcript->source($opt_ensembl_tag . '_' . $opt_havana_tag );

  return 0;
} ## end sub merge

sub copy {
  my ( $target_gene, $source_transcript ) = @_;

  printf( "Copy> Biotypes: source transcript: %s, target gene: %s\n",
          $source_transcript->biotype(), $target_gene->biotype() );

###############################################################################
# Case 1 - Havana gene is coding, Ensembl transcript is a pseudogene
# Do not copy the Ensembl transcript. The corresponding Ensembl gene
# will not be listed as processed and therefore should be copied over
# at the end of merge process (if it wasn't processed elsewhere)
###############################################################################
  if ( $source_transcript->{__is_pseudogene} &&
       $target_gene->{__is_coding} )
  {
    print( "Copy> Source transcript is pseudogene, " .
           "will not copy it into a coding gene.\n" );
    print( "Copy> Deleting the Ensembl annotation.\n");
    return 0;
  }

###############################################################################
# Case 2 - The Havana gene is a pseudogene. In this case it the Ensembl
# transcript won't be copied and the corresponding gene will be listed
# as processed and will not be copied at the end. The only exception to
# this is when the Ensembl transcript is CCDS, in which case the transcript
# will be copied and the biotype will be updated to the Ensembl biotype.
# The Ensembl gene will be listed as processed and will not be copied at
# the end of merge process
###############################################################################
  elsif ( $target_gene->{__is_pseudogene} ) {
  

    unless ( $source_transcript->{__is_ccds} ) {
       print( "Copy> Target gene is pseudogene, " .
              "will not copy anything into it.\n" );
       print( "Copy> Deleting the Ensembl annotation.\n");
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
# Case 3 - The Havana gene is labelled as belonging to a gene cluster. This
# tag is read from the transcripts, so at least one transcript was labelled.
# In this case the Ensembl transcript will not be copied over.
# The exception to this is if the Ensembl transcript is CCDS, in this case
# the Ensembl transcript will be copied and the Havana gene biotype will
# updated if needed. Either way the corresponding the Ensembl gene will be
# will be listed as processed and not copied at the end of the merge process.
###############################################################################
  elsif ( $target_gene->{__is_gene_cluster} ) {

    unless ( $source_transcript->{__is_ccds} ) {
      print( "Copy> Target gene is part of gene cluster, " .
           "will not copy overlapping Ensembl transcripts ".
           "into it.\n" );
      print( "Copy> Deleting the Ensembl annotation.\n");
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
# Case 4 - The Havana gene has an assembly error, in this case the Ensembl
# transcript will be copied in and if the Ensembl gene has a translation and
# the Havana gene biotype doesn't match the biotype of the Ensembl transcript
# the biotype of the Ensembl transcript overwrites the Havana gene biotype.
# This may be in the wrong place logically. The corresponding Ensembl gene is
# listed as processed and will not be copied at the end of the merge process
###############################################################################
  elsif ( $target_gene->{__has_ref_error} ) {
    print( "Copy> Target gene has assembly error\n");

    if ( defined( $source_transcript->translation() ) &&
         $target_gene->biotype() ne $source_transcript->biotype() ) {
      printf( "Copy> Updating gene biotype from %s to %s\n",
              $target_gene->biotype(), $source_transcript->biotype() );
      $target_gene->biotype( $source_transcript->biotype() );
    }

  }

###############################################################################
# Case 5 - The Havana gene is non coding (but not a pseudogene) and the
# Ensembl transcript has a translation. In this case the translation is
# removed from the Ensembl transcript and the exon phases are all set to
# -1, which means non-coding, before the transcript is copied. The only
# exception to this is when the Ensembl transcript is CCDS, in this case
# the biotype of the Havana gene is overwritten with the Ensembl transcript
# biotype. The corresponding Ensembl gene is listed as processed and is
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

  # Start by transferring the $source_transcript to the same slice as
  # the $target_gene.
  my $new_source_transcript =
    $source_transcript->transfer( $target_gene->slice() );

  $target_gene->add_Transcript($new_source_transcript);

  add_logic_name_suffix( $new_source_transcript, 'copied' );
  add_logic_name_suffix( $target_gene,           'merged' );
  $target_gene->source( $opt_ensembl_tag . '_' . $opt_havana_tag );

  return 0;
} ## end sub copy


# Strip the phases off all the exons in a transcript. Used on Ensembl
# coding transcripts that get demoted to non-coding. This will only
# be run on coding Ensembl transcripts that are copied to non-coding, 
# non-pseudogene Havana genes. Yes that's confusing, just roll with it.
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

__END__

=head1 NAME

merge.pl - A script that merges the manual annotations from Havana with
the automatic annotations from Ensembl.

=head1 SYNOPSIS

merge.pl [options]

merge --help

  Options:
    --host_ensembl  (required)
    --port_ensembl
    --user_ensembl  (required)
    --pass_ensembl      or --password_ensembl
    --dbname_ensembl    or --database_ensembl   (required)

    --host_havana   (required)
    --port_havana
    --user_havana   (required)
    --pass_havana       or --password_havana
    --dbname_havana     or --database_havana    (required)

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

    --ensembl_include   or --ensembl_exclude
    --havana_include    or --havana_exclude

    --havana_tag
    --ensembl_tag

    --havana_gene_xref
    --havana_transcript_xref
    --havana_translation_xref

    --jobs
    --chunk

  For longer explanation of options:
    --help or -h or -?

=head1 OPTIONS

=head2 Input database connection options

=over

=item B<--host_ensembl>

The server host name for the Ensembl database (required).

=item B<--port_ensembl>

The port on the server to connect to (optional, default is 3306).

=item B<--user_ensembl>

The username to connect as (required).

=item B<--pass_ensembl> or B<--password_ensembl>

The password to authenticate with (optional, default is an empty
string).

=item B<--dbname_ensembl> or B<--database_ensembl>

The name of the Ensembl database (required).

=item B<--host_havana>

The server host name for the Havana database (required).

=item B<--port_havana>

The port on the server to connect to (optional, default is 3306).

=item B<--user_havana>

The username to connect as (required).

=item B<--pass_havana> or B<--password_havana>

The password to authenticate with (optional, default is an empty
string).

=item B<--dbname_havana> or B<--database_havana>

The name of the Havana database (required).

=back

=head2 Connection options for using a separate DNA database

These options are all optional and their default values are taken from
the values specified for the connection options for the Ensembl input
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

The port on the server to connect to (optional, default is 3306).

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

=item B<--ensembl_include> or B<--ensembl_exclude>

A comma-separated list of analysis logic names of genes that will be
included or excluded from the Ensembl gene set.  You may use one of
these options, but not both.

=item B<--havana_include> or B<--havana_exclude>

A comma-separated list of analysis logic names of genes that will be
included or excluded from the Havana gene set.  You may use one of these
options, but not both.

=back

=head2 Options related to tagging

=over

=item B<--ensembl_tag>

=item B<--havana_tag>

These options allows you to set the tags used for tagging analysis
logic names for Havana and Ensembl genes, transcripts, exons and
translations, as well as the logic names assuciated with the various
types of supporting featores.

A lugic name will be suffixed by C<_> followed by the value of the tag,
so that the original logic name C<otter>, for example, would be extended
to C<otter_havana> (if ths was a lgic name that was attached to Havana
annotation).

The tag values are also used to set the C<source> of the genes and
transcrpts written to the output database.  Merged genes and transcripts
gets a combined C<source> value, e.g. C<ensembl_havana>.

These options are optional.  The default values for these options are
C<havana> for B<--havana_tag> and C<ensembl> for B<--ensembl_tag>.

=back

=head2 Options related to Xrefs

=over

=item B<--havana_gene_xref>

=item B<--havana_transcript_xref>

=item B<--havana_translation_xref>

The values corresponding to these options should be a comma-separated
list of strings that makes up the C<external_db_name>,
C<db_display_name>, and C<type> of the xref database table.  This
information will be supplemented with the Havana stable IDs for genes,
transcripts and translations and added to the output gene set.

For example, the standard xrefs for Havana are specified as

  --havana_gene_xref='OTTG,Havana gene,ALT_GENE'
  --havana_transcript_xref='OTTT,Havana transcript,ALT_TRANS'
  --havana_translation_xref='OTTP,Havana translation,MISC'

These options are optional.  The default for these three options is the
default Havana xref information (as above).  It only needs to be changed
if you replace the Havana input data set with, for example, RefSeq.

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
