#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# This script compares the sections of the sequence at the beginning and
# end of the assembly exception (patch) with the reference equivalent.
# To cut down on the work, you should run this script against your database and
# a previously OK'd database. Then you can identify the new cases and only check them.

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use buildchecks::ScriptUtils;
use Getopt::Long;

# database
my $host;
my $user;
my $dbname;
my $port;
my $dnadbname;
my $dnahost;
my $dnaport;
my $dnauser;

# genomic location
my $coordsystem = 'chromosome';
my $coord_system_version = 'GRCh37';
my @seq_region_names;
my $all;

my $file;

# buffer
$| = 1;

GetOptions( 'host:s'                 => \$host,
            'user:s'                 => \$user,
            'dbname:s'               => \$dbname,
            'port:n'                 => \$port,
            'coord_system_version:s' => \$coord_system_version,
            'chromosomes:s'          => \@seq_region_names,
            'file:s'                 => \$file,
            'all'                    => \$all,
            'dnadbname:s'            => \$dnadbname,
            'dnahost:s'              => \$dnahost,
            'dnaport:s'              => \$dnaport,
            'dnauser:s'              => \$dnauser);

# sanity thing:
print STDERR "Database $dbname host $host port $port user $user\n";
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -dbname => $dbname
);

my $dnadb;
if ($dnadbname) {
  $dnadb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbname => $dnadbname,
                                       -host   => $dnahost,
                                       -port   => $dnaport,
                                       -user   => $dnauser);
  $db->dnadb($dnadb);
}

my $slice_adaptor = $db->get_SliceAdaptor();
my $asm_exception_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();

# outfile
my $ref_fh;
my $alt_fh;
if ($file && $file ne "stdout") {
  my $reffile = $file."_ref";
  my $altfile = $file."_alt";

  print STDERR "Printing to $reffile and $altfile\n";

  open REF,">$reffile" or die "couldn't open file ".$reffile." $!";
  open ALT,">$altfile" or die "couldn't open file ".$altfile." $!";
  $ref_fh = \*REF;
  $alt_fh = \*ALT;
} else {
  print STDERR "Printing to STDOUT\n";
  $ref_fh = \*STDOUT;
  $alt_fh = \*STDOUT;
}

# # #
# OK, now we begin to do stuff
# # #
my $chrhash = get_chrlengths_v20($db, $coord_system_version, $coordsystem);
if ($all) {
  if (scalar(@seq_region_names) == 0) {
    # we're ok, do not filter
  } else {
    throw("You have entered -all and -chr.");
  }
} else {
  if (scalar(@seq_region_names) == 0) {
    throw("Need either -all or chr name(s).");
  } else {
    filter_to_chr_list(\@seq_region_names,$chrhash,$db->dbc->dbname);
  }
}

# fetch slice on reference
foreach my $chr ( @{ sort_chr_names($chrhash) } ) {
  print STDERR "\nDoing chr $chr ...\n";
  my $chrstart  = 1;
  my $chrend    = $chrhash->{$chr};
  my $slicename = "$coordsystem:$coord_system_version:$chr:$chrstart:$chrend:1";

  # get reference slice
  my $ref_slice = $slice_adaptor->fetch_by_name($slicename);
  print STDERR "Got slice " . $ref_slice->name . "\n";

  # get all assembly excpetions
  # this gives us the bit of the reference chr where the exception lies
  # and it's by calling alt_slice that we get the HAP / PAR
  my @asm_except_feats =
    @{ $asm_exception_adaptor->fetch_all_by_Slice($ref_slice) };
  foreach my $aef (@asm_except_feats) {
      print STDERR "Got assembly exception feature"
        . " start "  . $aef->start
        . " end "    . $aef->end
        . " strand " . $aef->strand
        . " type "   . $aef->type
        . " slice "  . $aef->slice->name . "\n";

    # get the exception slice
    my $alt_slice = $aef->alternate_slice();
    print STDERR "Alternate slice is " . $alt_slice->name . "\n";

    # get a mini version of the reference slice
    my $mini_slicename = "$coordsystem:$coord_system_version:$chr:" . $aef->start . ":" . $aef->end . ":1";

    my $mini_ref_slice = $slice_adaptor->fetch_by_name($mini_slicename);
    print STDERR "Mini ref_slice is " . $mini_ref_slice->name . "\n";

    # ok, now walk in from both sides
    compare_seqs( $mini_ref_slice, $alt_slice, $ref_fh, $alt_fh );
  }
} # chr

sub compare_seqs {
  my ( $ref_slice, $alt_slice, $ref_fh, $alt_fh ) = @_;

  print STDERR "Ref start " . $ref_slice->start . " length  " . $ref_slice->length . " end " . $ref_slice->end . "\n";
  print STDERR "Alt start " . $alt_slice->start . " length  " . $alt_slice->length . " end " . $alt_slice->end . "\n";

  # subseq for start
  # coords relative to slice
  print STDERR "Ref " . $ref_slice->subseq( 1, 100 ) . "\n";
  print STDERR "Alt " . $alt_slice->subseq( 1, 100 ) . "\n";

  my @ref_seq_arr = split('',$ref_slice->subseq(1,100));
  my @alt_seq_arr = split('',$alt_slice->subseq(1,100));

  my $num_mismatches = 0;
  join('',map { $ref_seq_arr[$_] eq $alt_seq_arr[$_] ? $alt_seq_arr[$_] : $num_mismatches++ }
              0 .. $#alt_seq_arr);

  if ($num_mismatches >= 6) { # 95% identity instead of 100%
    print "WARNING: starts have different seqs\n";
  }

  # subseq for end
  # coords relative to slice
  my $ref_slice_end = $ref_slice->subseq(( $ref_slice->length - 99 ), $ref_slice->length);
  my $alt_slice_end = $alt_slice->subseq(( $alt_slice->length - 99 ), $alt_slice->length);

  print STDERR "Ref " . $ref_slice_end . "\n";
  print STDERR "Alt " . $alt_slice_end . "\n";

  my @ref_seq_end_arr = split('',$ref_slice_end);
  my @alt_seq_end_arr = split('',$alt_slice_end);

  $num_mismatches = 0;
  join('',map { $ref_seq_end_arr[$_] eq $alt_seq_end_arr[$_] ? $alt_seq_end_arr[$_] : $num_mismatches++ }
              0 .. $#alt_seq_end_arr);

  if ($num_mismatches >= 6) { # 95% identity instead of 100%
    print "WARNING: ends have different seqs\n";
  }

  reformat( $ref_fh, ">Ref_" . $ref_slice->name, $ref_slice->seq );
  reformat( $alt_fh, ">Alt_" . $alt_slice->name, $alt_slice->seq );
  return;
} ## end sub compare_seqs

sub reformat {
  my ( $fh, $id, $seq ) = @_;

  print $fh "$id\n";
  $seq =~ s/(.{1,60})/$1\n/g;
  print $fh $seq;
}

