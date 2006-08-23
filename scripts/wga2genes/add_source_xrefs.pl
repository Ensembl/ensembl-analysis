#!/usr/local/bin/perl

# add_source_xrefs.pl
#
# Uses transcript supporting features to add
# appropriate DBEntries for genes using source
# database


use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $src_dbname,
    $src_dbhost,
    $src_dbuser,
    $src_dbport,
    $src_dbpass,
    $verbose,
);

$dbuser = 'ensro';
$src_dbuser = 'ensro';

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'srcdbname=s' => \$src_dbname,
            'srcdbuser=s' => \$src_dbuser,
            'srcdbhost=s' => \$src_dbhost,
            'srcdbport=s' => \$src_dbport,
            'srcdbpass=s' => \$src_dbpass,
            'verbose' => \$verbose,

);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


my $src_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $src_dbname,
	'-host' => $src_dbhost,
	'-user' => $src_dbuser,
	'-port' => $src_dbport,
	'-pass' => $src_dbpass
);

my $ga = $db->get_GeneAdaptor;
my $dbea = $db->get_DBEntryAdaptor;

$verbose and print STDERR "Fetching genes...\n";

my @genes = map { $ga->fetch_by_dbID($_) } @{$ga->list_dbIDs};

foreach my $g (@genes) {
  my (@t) = @{$g->get_all_Transcripts};

  my (@g_dbens);

  foreach my $t (@{$g->get_all_Transcripts}) {
    my ($sf) = @{$t->get_all_supporting_features};

    my $src_tr = $src_db->get_TranslationAdaptor->
        fetch_by_stable_id($sf->hseqname);

    my $src_t = $src_db->get_TranscriptAdaptor->
        fetch_by_translation_stable_id($src_tr->stable_id);

    my $src_g = $src_db->get_GeneAdaptor->
        fetch_by_transcript_stable_id($src_t->stable_id);

    push @g_dbens, Bio::EnsEMBL::DBEntry->new(-primary_id => $src_g->stable_id,
                                              -version => $src_g->version,
                                              -dbname => 'Ens_Hs_gene',
                                              -release => 1,
                                              -display_id => $src_g->stable_id);

    my $t_dben = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_t->stable_id,
                                            -version => $src_t->version,
                                            -dbname => 'Ens_Hs_transcript',
                                            -release => 1,
                                            -display_id => $src_t->stable_id);
    print "Writing xref ", $t_dben->primary_id, " for transcript ", $t->dbID, "\n";
    $dbea->store($t_dben, $t->dbID, 'Transcript');

    if (my $tr = $t->translation) {
      # make DBEntry for translation, pointing to source translation
      my $tr_dben = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_tr->stable_id,
                                               -version => $src_tr->version,
                                               -dbname => 'Ens_Hs_translation',
                                               -release => 1,
                                               -display_id => $src_tr->stable_id);
      print "Writing xref ", $tr_dben->primary_id, " for translation ", $tr->dbID, "\n";
      $dbea->store($tr_dben, $tr->dbID, 'Translation');
    }
  }

  foreach my $gref (@g_dbens) {
    print "Writing xref ", $gref->primary_id, " for gene ", $g->dbID, "\n";
    $dbea->store($gref, $g->dbID, 'Gene');
  }
}
