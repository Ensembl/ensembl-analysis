# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
# Utilities for scripts for doing tasks such as defining/filtering input

package buildchecks::ScriptUtils;

use Exporter;
use vars qw(@ISA @EXPORT);
use strict;


@ISA=qw(Exporter);

@EXPORT=qw(
           filter_to_chr_list
           get_chrlengths_v20
           get_chrlengths_v19
           sort_chr_names
           );


sub get_chrlengths_v20 {
  my $db   = shift;
  my $type = shift;
  my $coordsystem = shift;
  my $ignore_haplotypes = shift;

  my %chrhash;

  if ($coordsystem ne 'toplevel') {
    if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
      die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
    }

    my $query = "select seq_region.name, seq_region.length as mce from seq_region,coord_system where" .
                " seq_region.coord_system_id=coord_system.coord_system_id and" .
                " coord_system.version = '" . $type . "' and coord_system.name='$coordsystem'";
    if ($ignore_haplotypes) {
      $query .= " and seq_region.seq_region_id not in (select seq_region_id from assembly_exception where exc_type in ('HAP','PATCH_NOVEL','PATCH_FIX'))";
    }

    #print "query $query\n";
    my $sth = $db->dbc->prepare($query);

    $sth->execute;


    my $hashref;
    while (($hashref = $sth->fetchrow_hashref) && defined($hashref)) {
      $chrhash{$hashref->{'name'}} = $hashref->{mce};
       #print $hashref->{'name'} . " " . $hashref->{'mce'} . "\n";
    }
  } else {
    my $sa = $db->get_SliceAdaptor;

    my @slices = @{$sa->fetch_all('toplevel',undef,undef,1)};

    foreach my $slice (@slices) {
      #print $slice->seq_region_name . " " . $slice->length . " " . $slice->coord_system->version . "\n";
      if ($slice->coord_system->version eq $type) {
        $chrhash{$slice->seq_region_name} = $slice->length;
      } else {
        print STDERR "WARNING: " . $slice->seq_region_name . " " .
        $slice->length . " " . $slice->coord_system->version
        . " is NOT on requested assembly version but is toplevel, will skip it!\n";
      }
    }
  }
  return \%chrhash;
}

sub get_chrlengths_v19 {
  my $db = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my %chrhash;

  my $q = qq( SELECT chrom.name,max(chr_end) FROM assembly as ass, chromosome as chrom
              WHERE ass.type = '$type' and ass.chromosome_id = chrom.chromosome_id
              GROUP BY chrom.name
            );

  my $sth = $db->dbc->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  return \%chrhash;
}

sub sort_chr_names {
  my ($chrhash) = @_;

  my @sorted_names = reverse sort bychrnum keys %$chrhash;
  return \@sorted_names;
}

sub bychrnum {


  my @awords = split /_/, $a;
  my @bwords = split /_/, $b;

  my $anum = $awords[0];
  my $bnum = $bwords[0];


  #  if ($anum !~ /^chr/ || $bnum !~ /^chr/) {
  #    die "Chr name doesn't begin with chr for $a or $b";
  #  }

  $anum =~ s/chr//;
  $bnum =~ s/chr//;

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
      return $awords[1] cmp $bwords[1];
    }
  }
}

sub filter_to_chr_list {
  my ($chromosomes, $chrhash, $dbname) = @_;

  #filter to specified chromosome names only
  if (scalar(@$chromosomes)) {
    foreach my $chr (@$chromosomes) {
      my $found = 0;
      foreach my $chr_from_hash (keys %$chrhash) {
        if ($chr_from_hash =~ /^${chr}$/) {
          $found = 1;
          last;
        }
      }
      if (!$found) {
        print "Didn't find chromosome named $chr in database " . $dbname . "\n";
      }
    }
    HASH: foreach my $chr_from_hash (keys %$chrhash) {
      foreach my $chr (@$chromosomes) {
        if ($chr_from_hash =~ /^${chr}$/) {
          #print "Keeping $chr_from_hash\n";
          next HASH;
        }
      }
      #print "Removing $chr_from_hash\n";
      delete($chrhash->{$chr_from_hash});
    }
  }
  return 1;
}

return 1;
