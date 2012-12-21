#!/usr/bin/env perl
# $Revision: 1.2 $

#
# Source: schema 20+ ensembl db with both an 
#   'old and 'new' assembly (eg both NCBI34 and NCBI35 human assemblies)
#   and a set of annotations that exist in the 'old' assembly - eg genes 
#   already built over NCBI34.
#
# Target: schema 20+ ensembl db with the 'new' assembly only (eg NCBI35)
#   and no genes.
#
# Action: Script reads the genes in the source db in both old and new assemblies
#   and compares them to see if their structure has changed. If it is identical,
#   then it writes them to the target db with the new assembly.
#
# WARNING: This script is not 'polished' - it is provided here so you
#   can have a starting point if you want to do this again (I find it easier.
#   to copy other people's scripts) but don't assume it will just work for you:
#   -- for instance, I wrote it to transfer ncrna's between assemblies, so
#   I have a check in there to weed out all 'ensembl' and 'pseudogenes'.

use warnings ;
use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host    = 'ecs2';
my $user    = 'ensro';
my $pass    = '';
my $dbname  = 'vivek_homo_sapiens_core_24_34';
my $port    = 3361;

my $c_host    = 'ecs2';
my $c_user    = 'ensro';
my $c_pass    = '';
my $c_port    = 3361;
my $c_dbname  = 'vivek_homo_sapiens_core_24_34';

my $t_host    = 'ecs4';
my $t_user    = 'ensadmin';
my $t_pass    = '****';
my $t_port    = 3350;
my $t_dbname  = 'vivek_homo_sapiens_23_35_pseudgene_check';

my $chr      = undef;
my $chrstart = undef;
my $chrend   = undef;

my $path     = 'NCBI35';

my $c_path   = 'NCBI34';

my $t_path   = 'NCBI35';

&GetOptions( 'host:s'    => \$host,
             'user:s'    => \$user,
             'pass:s'    => \$pass,
             'port:s'    => \$port,
             'dbname:s'  => \$dbname,
             'c_host:s'    => \$c_host,
             'c_user:s'    => \$c_user,
             'c_pass:s'    => \$c_pass,
             'c_port:s'    => \$c_port,
             'c_dbname:s'  => \$c_dbname,
             't_host:s'    => \$t_host,
             't_user:s'    => \$t_user,
             't_pass:s'    => \$t_pass,
             't_port:s'    => \$t_port,
             't_dbname:s'  => \$t_dbname,
             'chr:s'       => \$chr,
             'chrstart:n'  => \$chrstart,
             'chrend:n'    => \$chrend,
             'path:s'      => \$path,
             'c_path:s'    => \$c_path,
             't_path:s'    => \$t_path,
            );



my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
                                           -user => $user,
                                           -pass => $pass,
                                           -port => $port,
                                           -dbname => $dbname);

my $cdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $c_host,
                                             -user => $c_user,
                                             -pass => $c_pass,
                                             -port => $c_port,
                                             -dbname => $c_dbname);


my $tdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $t_host,
                                             -user => $t_user,
                                             -pass => $t_pass,
                                             -port => $t_port,
                                             -dbname => $t_dbname);

my $sgp = $sdb->get_SliceAdaptor;
my $aga = $sdb->get_GeneAdaptor;
my $pfa = $sdb->get_ProteinFeatureAdaptor;

my $c_sgp = $cdb->get_SliceAdaptor;
my $c_aga = $cdb->get_GeneAdaptor;

my $t_sgp = $tdb->get_SliceAdaptor;
my $t_aga = $tdb->get_GeneAdaptor;
my $t_pfa = $tdb->get_ProteinFeatureAdaptor;

my $vcontig = $sgp->fetch_by_region('chromosome', $chr,$chrstart,$chrend, 1, $path);
print STDERR "Fetched new slice for $chr, $chrstart, $chrend\n";

my $c_vcontig = $c_sgp->fetch_by_region('chromosome', $chr, undef, undef, undef, $c_path);
print STDERR "Fetched comparison slice for $chr\n";

my $t_vcontig = $t_sgp->fetch_by_region('chromosome', $chr,$chrstart,$chrend, 1, $t_path);
print STDERR "Fetched target slice\n";

my $genes = $vcontig->get_all_Genes();
print STDERR "Fetched ".scalar(@{$genes})." genes on source slice\n";

my %genehash;
foreach my $gene (@$genes) {
  if(($gene->type eq 'ensembl') || ($gene->type eq 'pseudogene')){
    next;
  }
  $genehash{$gene->stable_id} = $gene;
  print STDERR "caching new genes by stableid: ".$gene->stable_id."\n"; 
}

print STDERR "Cached ".scalar(keys %genehash)." genes on NCBI35\n";

my $c_genes = $c_vcontig->get_all_Genes();
print STDERR "Fetched comparison genes\n";

print STDERR "Comparing and writing ....\n";
my $nignored = 0;
my $ncgene = 0;
my $ndiffgene = 0;
CGENE: 
foreach my $c_gene (@$c_genes) {
  if(($c_gene->type eq 'ensembl') || ($c_gene->type eq 'pseudogene')){
    next CGENE;
  }

  $ncgene++;

  my $isdiff=0;

# Is it fully mapped?
  my @exons =  @{$c_gene->get_all_Exons};

  # Does the sequence-level 'thingy' change across exons? Contig?
#  my $firstseqname = $exons[0]->seq_region_name;
#  foreach my $exon (@exons) {
#    #print "Exon name " . $exon->seqname . " first seqname " . $firstseqname . "\n";
#    # changed ->seqname call to ->seq_region_name
#    if ($exon->seq_region_name ne $firstseqname) {
#      print "Ignoring gene " . $c_gene->stable_id . " which is on multiple sequences - not transferred\n";
#      $nignored++;
#      next CGENE;
#    }
#  }

# Is it at all mapped
#  if ($firstseqname ne $c_vcontig->name) {
#    print "Ignoring gene " . $c_gene->stable_id . " which is completely off path on $firstseqname - not transferred\n";
#    $nignored++;
#    next CGENE;
#  }


  my $gene = undef; 
  if (exists($genehash{$c_gene->stable_id})) {
    $gene = $genehash{$c_gene->stable_id};
    print STDERR "Gene ".$c_gene->stable_id." ( ".$c_gene->type.") was transferred\n";
    print STDERR " Gene has old coords " . $c_gene->stable_id . " from " . $c_gene->start . " to " . $c_gene->end."\n";

    # First check we have the same number of transcripts    
    my @transcripts = @{$gene->get_all_Transcripts};
    my @c_transcripts = @{$c_gene->get_all_Transcripts};

    if (scalar(@c_transcripts) != scalar(@transcripts)) {
      print STDERR "Gene " . $gene->stable_id . " has different numbers of transcripts\n";
      $isdiff=1;
    }

    my %tranhash;
    foreach my $tran (@transcripts) {
      # no longer have to sort exons
      # $tran->sort; 
      $tranhash{$tran->stable_id} = $tran;
    }

    if ($gene->strand != $c_gene->strand) {
      print STDERR " Note gene $gene->stable_id on different strands in two assemblies\n";
    }

    foreach my $c_tran (@c_transcripts) {
      #$c_tran->sort;

      if (exists($tranhash{$c_tran->stable_id})) {
        my $tran = $tranhash{$c_tran->stable_id};
        my @exons= @{$tran->get_all_Exons};
        my @c_exons= @{$c_tran->get_all_Exons};
            
        if (scalar(@exons) != scalar(@c_exons)) {
          print STDERR "Different numbers of exons in transcript " . $c_tran->stable_id . "\n";
          $isdiff=1;
        }

        my $nexon_to_comp = (scalar(@exons) > scalar(@c_exons)) ? scalar(@c_exons) : scalar(@exons);
        my $c_intron_len = 0;
        my $intron_len = 0;

        for (my $i=0;$i<$nexon_to_comp;$i++) {
          if ($exons[$i]->stable_id ne $c_exons[$i]->stable_id) {
            print STDERR  "Exon stable ids different for " . $exons[$i]->stable_id . " and " . 
                  $c_exons[$i]->stable_id . " in transcript " . $c_tran->stable_id . "\n";
            $isdiff=1;
          }
          if ($exons[$i]->length != $c_exons[$i]->length) {
            print STDERR "Exon lengths different for " . $exons[$i]->stable_id . " and " . 
                  $c_exons[$i]->stable_id . " in transcript " . $c_tran->stable_id . "\n";
            $isdiff=1;
          }
          if ($exons[$i]->seq->seq ne $c_exons[$i]->seq->seq) {
            print STDERR "Exon sequences different for " . $exons[$i]->stable_id . " and " . 
                  $c_exons[$i]->stable_id . " in transcript " . $c_tran->stable_id . "\n";
            $isdiff=1;
          }
          if ($i != 0) {
            if ($gene->strand == 1) {
              $intron_len   = $exons[$i]->start - $exons[$i-1]->end - 1;
            } else {
              $intron_len   = $exons[$i-1]->start - $exons[$i]->end - 1;
            }

            if ($c_gene->strand == 1) {
              $c_intron_len = $c_exons[$i]->start - $c_exons[$i-1]->end - 1;
            } else {
              $c_intron_len = $c_exons[$i-1]->start - $c_exons[$i]->end - 1;
            }
          }
        }
        my $intron_diff_len = abs($c_intron_len - $intron_len);
        if ($intron_diff_len > $c_intron_len/10) {
          print STDERR " Total intron length for transcript ".$c_tran->stable_id." changed by more than 10%%\n";
        }
      } else {
        print STDERR "Couldnt  find transcript " . $c_tran->stable_id . " to compare against\n";
        $isdiff=1;
      }
    }
  } else {
    print STDERR "Couldn't find gene " . $c_gene->stable_id . " to compare against\n";
    $isdiff=1;
  }
  if ($isdiff) {
    $ndiffgene++;
    print STDERR " Gene " . $c_gene->stable_id . " from " . $c_gene->start . " to " . $c_gene->end . " (old assembly coords) is different - not transferred\n";
    if($gene){
      print STDERR " -- it is mapped to Gene " . $gene->stable_id . " from " . $gene->start . " to " . $gene->end . " (new assembly coords)\n";
    }else{
       print STDERR "new gene just doesn't exist!\n";
    }
  } else {
    eval {
      print STDERR "About to write gene ".$gene->stable_id."\n"; 
      write_gene($t_aga,$t_vcontig,$gene,$pfa,$t_pfa);
    };
    if ($@) {
      print STDERR "Failed writing gene " . $gene->stable_id . "\n";
      print STDERR $@ . "\n";
    }
  }
}
print "N compared  = " . $ncgene . "\n";
print "N ignored   = " . $nignored . "\n";
print "N diff gene = " . $ndiffgene . "\n";
print "Done\n";

sub write_gene {
  my ($t_aga,$t_vcontig,$gene,$s_pfa,$t_pfa) = @_;
  print STDERR "Writing gene: " . $gene->stable_id . " from " . $gene->start . " to " . $gene->end . " (new assembly coords)\n";

  my %s_pfhash;
  $gene->slice($t_vcontig);
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    #$tran->sort;
    $tran->slice($t_vcontig);

    print "Transcript " . $tran->stable_id . "\n";

    # These lines force loads from the database to stop attempted lazy
    # loading during the write (which fail because they are to the wrong
    # db)

    if (defined($tran->translation)) {
      $s_pfhash{$tran->stable_id} = $s_pfa->fetch_by_translation_id($tran->translation->dbID);
    }

    my @exons= @{$tran->get_all_Exons};
    my $get = $tran->translation;
    $tran->_translation_id(undef);

    foreach my $exon (@exons) {
      $exon->stable_id;
      $exon->slice($t_vcontig);
      $exon->get_all_supporting_features; 
    }
  }

  # Transform gene to raw contig coords
  print "Gene " .$gene->start ." to " . $gene->end  . " type ".$gene->type."\n";
  #$gene->transform;
  #$gene->transfer($t_vcontig);
  # If you pass in a $c_gene, you should do a transfer here onto the target slice, but I don't think
  # this is necessary if you use the $gene retreived on the NCBI35 slice anyway...

  $t_aga->store($gene);

  foreach my $tran (@{$gene->get_all_Transcripts}) {
    if (defined($s_pfhash{$tran->stable_id})) {
      foreach my $pf (@{$s_pfhash{$tran->stable_id}}) {
        $pf->seqname($tran->translation->dbID);
        if (!$pf->score) { $pf->score(0) };
        if (!$pf->percent_id) { $pf->percent_id(0) };
        if (!$pf->p_value) { $pf->p_value(0) };
        $t_pfa->store($pf);
      }
    }
  }
}

