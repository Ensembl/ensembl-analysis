#!/usr/local/ensembl/bin/perl

# identify_bad_source_transcripts.pl
# 
# Scans the given database for protein_coding genes and
# transcripts that are unsuitable for projecting onto
# a low-coverage genome via a whole-genome alignment
#
# The resulting set can act as the starting point for a "kill list"  
# of genes and transcripts that will be ignored during the
# WGA2Genes process. 

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
	$dbname,
	$dbhost,
	$dbuser,
	$dbport,
	$dbpass
);

$dbuser = 'ensro';

&GetOptions(
	'dbname=s' => \$dbname,
	'dbuser=s' => \$dbuser,
	'dbhost=s' => \$dbhost,
	'dbport=s' => \$dbport,
	'dbpass=s' => \$dbpass
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);

my @slices;
if (@ARGV) {
  foreach my $sid (@ARGV) {
    push @slices, $db->get_SliceAdaptor->fetch_by_name($sid);
  }
} else {
  @slices = @{$db->get_SliceAdaptor->fetch_all('toplevel')};
}


my %reject;

foreach my $sl (sort { $a->length <=> $b->length } @slices) {
  my @orig_genes = @{$sl->get_all_Genes_by_type('protein_coding')};

  #
  # Reject Mitochondrial genes
  #
  if ($sl->seq_region_name eq 'MT') {
    map { $reject{$_->stable_id} =  "Mitochondrial protein" } @orig_genes;
    next;
  }

  my @genes;

  foreach my $g (@orig_genes) {
    my ($gst, $gen);

    my @trans;
    foreach my $t (@{$g->get_all_Transcripts}) {

      if (length($t->translateable_seq) % 3 != 0) {
        $reject{$t->stable_id} = "Non modulo-3 coding length";
      } elsif ($t->biotype ne 'protein_coding') {
        $reject{$t->stable_id} = "Non-coding transcript in protein_coding gene";
      } else {
        my ($cst, $cen);
      
        foreach my $e (@{$t->get_all_translateable_Exons}) {
          $cst = $e->start if not defined $cst or $cst > $e->start;
          $cen = $e->end   if not defined $cen or $cen < $e->end;
        }

        push @trans, {
          start => $cst,
          end   => $cen,
          tran => $t,
        }
      }
    }

    if (@trans) {
      foreach my $t (@trans) {
        $gst = $t->{start} if not defined $gst or $gst > $t->{start};
        $gen = $t->{end}   if not defined $gen or $gen < $t->{end};
      }
      
      push @genes, {
        start => $gst,
        end   => $gen,
        gene  => $g,
        transcripts => \@trans,
      }
    }
  }

  @genes = sort { $a->{start} <=> $b->{start} } @genes;

  my @clusters;
  foreach my $g (@genes) {
    if (not @clusters or $clusters[-1]->{end} < $g->{start} - 1) {
      push @clusters, {
        start => $g->{start},
        end   => $g->{end},
        genes => [$g],
      };
    } else {
      if ($g->{end} > $clusters[-1]->{end}) {
        $clusters[-1]->{end} = $g->{end};
      }
      push @{$clusters[-1]->{genes}}, $g;
    }
  }

  # we're interested in clusters that become separated 
  # when you remove on (apparently problemattic) gene

  foreach my $c (@clusters) {
    if (@{$c->{genes}} > 1) {      
      #printf("CLUSTER %s %d %d %d members\n", 
      #       $sl->seq_region_name, 
      #       $c->{start}, 
      #       $c->{end}, 
      #       scalar(@{$c->{genes}}));

      # indentify gene in the cluster where one transcript overlap
      # genomically with other genes, but other transcript(s) do
      # not

      foreach my $g (@{$c->{genes}}) {
        my @inside;

        foreach my $og (@{$c->{genes}}) {
          next if $og eq $g;
          
          if ($og->{start} > $g->{start} and 
              $og->{end}   < $g->{end}) {
            push @inside, $og;
          }
        }

        if (@inside) {
          # other genes fit inside this one. But is it just a dogdy
          # transcript that is causing this overlap? 
          my (@overlaps, @doesnt);
          foreach my $t (@{$g->{transcripts}}) {
            my $all_inside = 1;
            foreach my $in (@inside) {
              if ($in->{start} < $t->{start} or
                  $in->{end} >   $t->{end}) {
                $all_inside = 0;
                last;
              }
            }
            if ($all_inside) {
              push @overlaps, $t;
            } else {
              push @doesnt, $t;
            }
          }

          if (scalar(@overlaps) == 1 and scalar(@doesnt) > 1) {
            map { $reject{$_->{tran}->stable_id} = "Dodgy transcript" } @overlaps;
          }
        }
      }
    }
  }
}

my %reasons;
foreach my $id (keys %reject) {
  push @{$reasons{$reject{$id}}}, $id;
}

foreach my $reason (keys %reasons) {
  foreach my $id (@{$reasons{$reason}}) {
    print "$id\t$reason\n";
  }
}
