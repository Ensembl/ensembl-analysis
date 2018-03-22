#!/usr/bin/env perl

# add_source_xrefs.pl
#
# Uses transcript supporting features to add
# appropriate DBEntries for genes using source
# database

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $stable_id_map_file,
    @src_dbname,
    $src_dbhost,
    $src_dbuser,
    $src_dbport,
    $src_dbpass,
    $verbose,
);

$dbport = 3306;
$dbuser = 'ensro';
$src_dbuser = 'ensro';

GetOptions(
            'dbname|db|D=s' => \$dbname,
            'dbuser|user|u=s' => \$dbuser,
            'dbhost|host|h=s' => \$dbhost,
            'dbport|port|P=s' => \$dbport,
            'dbpass|pass|p=s' => \$dbpass,
            'stableidmap=s' => \$stable_id_map_file,
            'srcdbname=s@' => \@src_dbname,
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

my (@src_dbs, $stable_id_map);

if (defined $stable_id_map_file) {
  $stable_id_map = &read_map_file($stable_id_map_file);
} else {
  foreach my $src_dbname (@src_dbname) {
    my $src_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
        new(
            '-dbname' => $src_dbname,
            '-host' => $src_dbhost,
            '-user' => $src_dbuser,
            '-port' => $src_dbport,
            '-pass' => $src_dbpass
            );
    if (not defined $src_db) {
      die "Could not connect to $src_dbname\n";
    } else {
      push @src_dbs, $src_db;
    }
  }
}

my $ga = $db->get_GeneAdaptor;
my $dbea = $db->get_DBEntryAdaptor;

$verbose and print STDERR "Fetching genes...";

my @genes;
if (@ARGV) {
  @genes = map { $ga->fetch_by_stable_id($_) } @ARGV;
} else {
  @genes = map { $ga->fetch_by_dbID($_) } @{$ga->list_dbIDs};
}

$verbose and print STDERR "fetched ", scalar(@genes), " genes\n";


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

      my ($src_p_id, $src_p_ver);
      my ($src_t_id, $src_t_ver);
      my ($src_g_id, $src_g_ver);

      if (defined $stable_id_map) {
        if (exists $stable_id_map->{$ens_pep_name}) {
          ($src_p_id, $src_p_ver, $src_t_id, $src_t_ver, $src_g_id, $src_g_ver) = 
              @{$stable_id_map->{$ens_pep_name}};
        }
      } elsif (@src_dbs) {

        my ($src_p, $src_t, $src_g) = &find_source_transcript_and_gene($ens_pep_name);

        if (defined $src_p) {
          $src_p_id = $src_p->stable_id;
          $src_p_ver = $src_p->version;
        }
        if (defined $src_t) {
          $src_t_id = $src_t->stable_id;
          $src_t_ver = $src_t->version;
        }
        if (defined $src_g) {
          $src_g_id = $src_g->stable_id;
          $src_g_ver = $src_g->version;
        }
      }

      if (defined $src_t_id and not exists $t_dbens{$src_t_id}) {
        $t_dbens{$src_t_id} = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_t_id,
                                                         -version => $src_t_ver,
                                                         -dbname => 'Ens_Hs_transcript',
                                                         -display_id => $src_t_id);
      }          
      
      if (defined $src_g_id and not exists $g_dbens{$src_g_id}) {
        $g_dbens{$src_g_id} = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_g_id,
                                                         -version => $src_g_ver,
                                                         -dbname => 'Ens_Hs_gene',
                                                         -display_id => $src_g_id);
      }
        
      if (my $tr = $t->translation and defined $src_p_id) {
        # make DBEntry for translation, pointing to source translation
        my $tr_dben = Bio::EnsEMBL::DBEntry->new(-primary_id => $src_p_id,
                                                 -version => $src_p_ver,
                                                 -dbname => 'Ens_Hs_translation',
                                                 -display_id => $src_p_id);
          printf(STDERR "Writing xref %s for translation %s (%d)\n", 
                 $tr_dben->primary_id, $tr->stable_id, $tr->dbID) 
              if $verbose;
        $dbea->store($tr_dben, $tr->dbID, 'Translation');
      }

    }
    
    foreach my $sid (keys %t_dbens) {
      printf(STDERR "Writing xref %s for transcript %s (%d)\n", 
             $t_dbens{$sid}->primary_id, $t->stable_id, $t->dbID)
          if $verbose;
      $dbea->store($t_dbens{$sid}, $t->dbID, 'Transcript');
    }
  }

  foreach my $sid (keys %g_dbens) {
    printf(STDERR "Writing xref %s for gene %s (%d)\n", 
           $g_dbens{$sid}->primary_id, $g->stable_id, $g->dbID)
        if $verbose;
    $dbea->store($g_dbens{$sid}, $g->dbID, 'Gene');
  }
}



sub find_source_transcript_and_gene {
  my ($ensp_id) = shift;

  my ($src_p, $src_t, $src_g);

  foreach my $db (@src_dbs) {
    my $src_db;

    $src_p = $db->get_TranslationAdaptor->
        fetch_by_stable_id($ensp_id);
    next unless defined $src_p;    

    $src_db = $db;
    
    if (defined $src_p) {
      $src_t = $src_db->get_TranscriptAdaptor->
          fetch_by_translation_stable_id($src_p->stable_id);


      if (defined $src_t) {
        $src_g = $src_db->get_GeneAdaptor->
            fetch_by_transcript_stable_id($src_t->stable_id);
      }
    }
  }

  return ($src_p, $src_t, $src_g);
}



sub read_map_file {
  my ($f) = shift;

  my %idmap;

  open F, $f or die "Could not open map file $f\n";
  while(<F>) {
    chomp;
    my @l = split;
    my $pid = $l[0];
    
    if (not exists $idmap{$l[0]} or
        $l[1] > $idmap{$pid}->[1]) {
      $idmap{$pid} = \@l;
    }
  }
  close(F);

  return \%idmap;
}
