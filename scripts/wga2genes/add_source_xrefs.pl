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

my @genes;
if (@ARGV) {
  @genes = map { $ga->fetch_by_stable_id($_) } @ARGV;
} else {
  @genes = map { $ga->fetch_by_dbID($_) } @{$ga->list_dbIDs};
}

foreach my $g (@genes) {
  my (@t) = @{$g->get_all_Transcripts};

  my (%g_dbens);

  foreach my $t (@{$g->get_all_Transcripts}) {
    my @ens_pep_names;
    foreach my $sf (@{$t->get_all_supporting_features}) {
      if ($sf->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {
        push @ens_pep_names, $sf->hseqname;
      }
    }

    my %t_dbens;

    foreach my $ens_pep_name (@ens_pep_names) {    
            
      my $src_tr = $src_db->get_TranslationAdaptor->
          fetch_by_stable_id($ens_pep_name);
      
      if (defined $src_tr) {
        my $src_t = $src_db->get_TranscriptAdaptor->
            fetch_by_translation_stable_id($src_tr->stable_id);
      
        if (defined $src_t and not exists $t_dbens{$src_t->stable_id}) {
          $t_dbens{$src_t->stable_id} = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_t->stable_id,
                                                                   -version => $src_t->version,
                                                                   -dbname => 'Ens_Hs_transcript',
                                                                   -release => 1,
                                                                   -display_id => $src_t->stable_id);
          
          my $src_g = $src_db->get_GeneAdaptor->
              fetch_by_transcript_stable_id($src_t->stable_id);
      
          if (defined $src_g and not exists $g_dbens{$src_g->stable_id}) {
            $g_dbens{$src_g->stable_id} = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_g->stable_id,
                                                                     -version => $src_g->version,
                                                                     -dbname => 'Ens_Hs_gene',
                                                                     -release => 1,
                                                                     -display_id => $src_g->stable_id);
          }
        }
        
        if (my $tr = $t->translation) {
          # make DBEntry for translation, pointing to source translation
          my $tr_dben = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_tr->stable_id,
                                                   -version => $src_tr->version,
                                                   -dbname => 'Ens_Hs_translation',
                                                   -release => 1,
                                                   -display_id => $src_tr->stable_id);
          printf(STDERR "Writing xref %s for translation %s (%d)\n", 
                 $tr_dben->primary_id, $tr->stable_id, $tr->dbID) 
              if $verbose;
          $dbea->store($tr_dben, $tr->dbID, 'Translation');
        }
      }
    }
    
    foreach my $sid (keys %t_dbens) {
      printf(STDERR "Writing xref %s for transcript %s (%d)\n", 
             $t_dbens{$sid}->primary_id, $t->stable_id, $t->dbID)
          if $verbose;
      $dbea->store($t_dbens{$sid}, $g->dbID, 'Transcript');
    }
  }

  foreach my $sid (keys %g_dbens) {
    printf(STDERR "Writing xref %s for gene %s (%d)\n", 
           $g_dbens{$sid}->primary_id, $g->stable_id, $g->dbID)
        if $verbose;
    $dbea->store($g_dbens{$sid}, $g->dbID, 'Gene');
  }
}
