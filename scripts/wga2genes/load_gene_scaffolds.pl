#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Mapper::Unit;
use Bio::EnsEMBL::Mapper::Pair;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $verbose,
    $seq_level,
    $asm_coord_sys_name,
    $asm_coord_sys_version,    
    $cmp_coord_sys_name,
    $test,
    %agps,
    %comp_slices,
);

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'asm_coord_sys_name=s' => \$asm_coord_sys_name,
            'asm_coord_sys_version=s' => \$asm_coord_sys_version,
            'cmp_coord_sys=s' => \$cmp_coord_sys_name, 
            'seqlevel' => \$seq_level,
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


my $tl_type_id = &get_toplevel_attribute_type_id;
my $tl_attr = Bio::EnsEMBL::Attribute->new(-code => 'toplevel',
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

    } else {
      if ($test) {
        printf("ASSEMBLY: %s %d %d %s %d %d %d\n", 
               $slice->seq_region_name,
               $comp->from->start,
               $comp->from->end,
               $comp->to->start,
               $comp->to->end,
               $comp->to->strand);
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

# finally, remove toplevel attribute from all components
# of gene scaffolds

$verbose and print STDERR "Removing toplevel from used...\n";

foreach my $sl_id (keys %comp_slices) {
  my $sl = $comp_slices{$sl_id};
  if ($test) {
    printf "Removed toplevel from slice %s\n", $sl->seq_region_name;
  } else {
    $aa->remove_from_Slice($sl, $tl_type_id);
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
    
    $agp_hash->{$l[0]}->{length} += $l[2] - $l[1] + 1;
    
    if ($l[3] ne 'N') {
      my $from = Bio::EnsEMBL::Mapper::Unit->new($l[0], $l[1], $l[2]);
      my $to  = Bio::EnsEMBL::Mapper::Unit->new($l[4], $l[5], $l[6]);
      my $pair = Bio::EnsEMBL::Mapper::Pair->new($from,
                                                 $to,
                                                 $l[7] eq '-' ? -1 : 1);

      push @{$agp_hash->{$l[0]}->{components}}, $pair;
    }
  }
}


##########################################################
# write_coord_system
##########################################################
sub write_coord_system {
  my ($new_name, $new_version) = @_;

  # check to see that it does not already exist
  my $csa = $db->get_CoordSystemAdaptor;
  my $cs = $csa->fetch_by_name($new_name, $new_version);

  if (not defined $cs) {
    $cs = Bio::EnsEMBL::CoordSystem->new(-name => $new_name,
                                         -version => $new_version,
                                         -rank => 1,
                                         -default => 1);
    $csa->store($cs);
  } else {
    warn "Warning: CoordSystem '$new_name' already exists in the database, so not writing\n";
  }

  return $cs;
}

##########################################################
# get_toplevel_attrbute_code
##########################################################
sub get_toplevel_attribute_type_id {
  my $st = $db->dbc->prepare("SELECT attrib_type_id from attrib_type where code = 'toplevel'");
  $st->execute;
  my ($res) = @{$st->fetchrow_arrayref};
  $st->finish;

  return $res;
}
