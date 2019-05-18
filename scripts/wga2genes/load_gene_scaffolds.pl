#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Mapper::Unit;
use Bio::EnsEMBL::Mapper::Pair;

my $TOPLEVEL_ATTR_CODE = 'toplevel';

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $verbose,
    $seq_level,
    $direct_level,
    $asm_coord_sys_name,
    $asm_coord_sys_version,    
    $cmp_coord_sys_name,
    $test,
    %agps,
    %comp_slices,
);

$dbport = 3306;

GetOptions(
            'dbname|db|D=s' => \$dbname,
            'dbuser|user|u=s' => \$dbuser,
            'dbhost|host|h=s' => \$dbhost,
            'dbport|port|P=s' => \$dbport,
            'dbpass|pass|p=s' => \$dbpass,
            'asm_coord_sys_name=s' => \$asm_coord_sys_name,
            'asm_coord_sys_version=s' => \$asm_coord_sys_version,
            'cmp_coord_sys_name=s' => \$cmp_coord_sys_name, 
            'seqlevel' => \$seq_level,
            'directlevel' => \$direct_level,
            'test' => \$test,
            'verbose' => \$verbose,
            );

die "You must give a name for the assembled coord system\n"
    if not defined $asm_coord_sys_name;
die "You must give a name for the component coord system\n"
    if not defined $cmp_coord_sys_name;

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


my $cmp_coord_sys_version;
eval {
  my $cs = $db->get_CoordSystemAdaptor->fetch_by_name($cmp_coord_sys_name);
  
  if (defined $cs) {
    $cmp_coord_sys_version = $cs->version;
  }
};
die "Could not find entry for component coord system\n"
    if not defined $cmp_coord_sys_version;

$asm_coord_sys_version = $cmp_coord_sys_version if not defined $asm_coord_sys_version;

my $tl_type_id = &get_toplevel_attribute_type_id;
my $tl_attr = Bio::EnsEMBL::Attribute->new(-code => $TOPLEVEL_ATTR_CODE,
                                           -name => 'Top Level',
                                           -desc => 'Top Level Non-Redundant Sequence Region',
                                           -value => 1
                                           );
my $aa = $db->get_AttributeAdaptor;


$verbose and print STDERR "Reading AGP files...\n";
&read_agp_files(\%agps);

my $coord_sys = &write_coord_system($asm_coord_sys_name, $asm_coord_sys_version);


my $ass_sth = $db->dbc->prepare("INSERT into assembly VALUES(?,?,?,?,?,?,?)");

my $test_seq_reg_id_count = 1;

foreach my $gs_id (keys %agps) {
  $verbose and print STDERR "Processing $gs_id\n";

  my $agp = $agps{$gs_id};

  my $len  = $agp->{length};

  my $slice = Bio::EnsEMBL::Slice->new(-seq_region_name   => $gs_id,
                                       -start             => 1,
                                       -end               => $len,
                                       -seq_region_length => $len,
                                       -strand            => 1,
                                       -coord_system      => $coord_sys);
  my $seq_region_id;
  if ($test) {
    $seq_region_id = $test_seq_reg_id_count++;
  } else {
    eval {
      $seq_region_id = $db->get_SliceAdaptor->get_seq_region_id($slice);
    };
    if ($@) {
      # slice does not exist; store it
      $seq_region_id = $db->get_SliceAdaptor->store($slice);

      if (not defined $seq_region_id) {
        die "Serious error when trying to store slice $gs_id; aborting\n";
      }
      $aa->store_on_Slice($slice, [$tl_attr]);
    }
  }

  foreach my $comp (@{$agp->{components}}) {
    if (not exists $comp_slices{$comp->to->id}) {
      $comp_slices{$comp->to->id} = $db->get_SliceAdaptor->fetch_by_region($cmp_coord_sys_name,
                                                                     $comp->to->id);
    }

    if ($seq_level) {
      my $bit_slice = $db->get_SliceAdaptor->fetch_by_region($cmp_coord_sys_name,
                                                             $comp->to->id,
                                                             $comp->to->start,
                                                             $comp->to->end,
                                                             $comp->ori);
      my @bits = @{$bit_slice->project('seqlevel')};


      foreach my $bit (@bits) {
        my $gs_start = $bit->from_start + $comp->from->start - 1;
        my $gs_end   = $bit->from_end   + $comp->from->start - 1;
        if ($test) {
          printf("ASSEMBLY: %s %d %d %s %d %d %d\n", 
                 $slice->seq_region_name,
                 $gs_start,
                 $gs_end,
                 $bit->to_Slice->seq_region_name,
                 $bit->to_Slice->start,
                 $bit->to_Slice->end,
                 $bit->to_Slice->strand);
        } else {
          $ass_sth->execute($seq_region_id,
                            $db->get_SliceAdaptor->get_seq_region_id($bit->to_Slice),
                            $gs_start,
                            $gs_end,
                            $bit->to_Slice->start,
                            $bit->to_Slice->end,
                            $bit->to_Slice->strand);
        }
      }

    }
    if ($direct_level) {
      if ($test) {
        printf("ASSEMBLY: %s %d %d %s %d %d %d\n", 
               $slice->seq_region_name,
               $comp->from->start,
               $comp->from->end,
               $comp->to->start,
               $comp->to->end,
               $comp->ori);
      } else {
        $ass_sth->execute($seq_region_id,
                          $db->get_SliceAdaptor->get_seq_region_id($comp_slices{$comp->to->id}),
                          $comp->from->start,
                          $comp->from->end,
                          $comp->to->start,
                          $comp->to->end,
                          $comp->ori);
      }
    }
  } 
}

$ass_sth->finish;

# remove toplevel attribute from all components
# of gene scaffolds

$verbose and print STDERR "Removing toplevel from used...\n";

foreach my $sl_id (keys %comp_slices) {
  my $sl = $comp_slices{$sl_id};
  if ($test) {
    printf "Removed toplevel from slice %s sid %d\n", $sl->seq_region_name, $sl->get_seq_region_id;
  } else {
    $aa->remove_from_Slice($sl, $TOPLEVEL_ATTR_CODE);
  }
}

# finally, update meta table

$verbose and print STDERR "Updating meta with assembly.mapping ...\n";

my $mca = $db->get_MetaContainer;
if ($seq_level) {
  # need to get name and version of seq level coord system
  my $cs = $db->get_CoordSystemAdaptor->fetch_by_name('seqlevel');
  my $mapstring = sprintf("%s:%s|%s", 
                          $asm_coord_sys_name,
                          $asm_coord_sys_version,
                          $cs->name);
  if ($test) {
    printf("META: %s %s\n", "assembly.mapping", $mapstring);
  } else {
    $mca->store_key_value('assembly.mapping', $mapstring);
  }
}
if ($direct_level) {
  my $mapstring = sprintf("%s:%s#%s:%s", 
                          $asm_coord_sys_name,
                          $asm_coord_sys_version,
                          $cmp_coord_sys_name,
                          $cmp_coord_sys_version);
  if ($test) {
    printf("META: %s %s\n", "assembly.mapping", $mapstring);
  } else {
    $mca->store_key_value('assembly.mapping', $mapstring);
  }
}


##########################################################
# read_agp_files
##########################################################
sub read_agp_files {
  my ($agp_hash) = @_;

  while(<>) {
    chomp;
    /^\#/ and next;
    
    my @l = split(/\t/, $_);
    
    if ($l[4] ne 'N') {
      my $from = Bio::EnsEMBL::Mapper::Unit->new($l[0], $l[1], $l[2]);
      my $to  = Bio::EnsEMBL::Mapper::Unit->new($l[5], $l[6], $l[7]);
      my $pair = Bio::EnsEMBL::Mapper::Pair->new($from,
                                                 $to,
                                                 $l[8] eq '-' ? -1 : 1);

      push @{$agp_hash->{$l[0]}->{components}}, $pair;
      if (not exists $agp_hash->{$l[0]} or
          $agp_hash->{$l[0]}->{length} < $l[2]) {
        $agp_hash->{$l[0]}->{length} = $l[2];
      }
    }
  }
}


##########################################################
# write_coord_system
##########################################################
sub write_coord_system {
  my ($new_name, $new_version) = @_;

  # check to see that it does not already exist
  my $cs = $db->get_CoordSystemAdaptor->fetch_by_name($new_name, $new_version);

  if (not defined $cs) {
    $cs = Bio::EnsEMBL::CoordSystem->new(-name => $new_name,
                                         -version => $new_version,
                                         -rank => 1,
                                         -default => 1);

    if (not $test) {
      # if coord sys with rank 1 already exists, we need to shift
      # the others down
      my @coord_sys = @{$db->get_CoordSystemAdaptor->fetch_all};
      @coord_sys = sort { $a->rank <=> $b->rank } @coord_sys;
      
      if ($coord_sys[0]->rank == 1) {
        foreach my $cs (reverse @coord_sys) {
          my $st = $db->dbc->prepare("update coord_system set rank = ? where coord_system_id = ?");        
          $st->execute($cs->rank + 1, $cs->dbID);
          $st->finish;
        }
      }

      # can't use cached adaptor to store the new entry here because 
      # the ranks are cached. Need to construct a new adaptor
      my $csa = Bio::EnsEMBL::DBSQL::CoordSystemAdaptor->new($db);
      $csa->store($cs);
    }
  } else {
    warn "Warning: CoordSystem '$new_name' already exists in the database, so not writing\n";
  }

  return $cs;
}

##########################################################
# get_toplevel_attrbute_code
##########################################################
sub get_toplevel_attribute_type_id {
  my $st = $db->dbc->prepare("SELECT attrib_type_id from attrib_type where code = '$TOPLEVEL_ATTR_CODE'");
  $st->execute;
  my ($res) = @{$st->fetchrow_arrayref};
  $st->finish;

  return $res;
}
