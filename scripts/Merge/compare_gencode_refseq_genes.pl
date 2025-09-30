#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=head1 NAME

  compare_gencode_refseq_genes.pl 

=head1 DESCRIPTION

  compare_gencode_refseq_genes.pl reads in gencode genes and refseq genes, both stored in
  ensembl schema core databases, and looks for:
    - where gencode has a gene but refseq does not
    - where both groups have a gene but the biotype differs
  then writes a gene attribute into the ensembl database

  Future options: compare transcripts
  
=head1 OPTIONS

  -ensemblhost     host name for database eg. core
  -ensemblport     port number for database
  -ensembldbname   core database name
  -ensembluser     write-access user
  -ensemblpass     write-access password

  -refseqhost     host name for refseq database eg. otherfeatures
  -refseqport     port number for database
  -refseqdbname   otherfeatures database name
  -refsequser     read-only access user

  -dnahost        host name for dna database 
  -dnaport        port number for database 
  -dnadbname      core database name
  -dnauser        read-only access user 

  -chr_num         One chromosome / seq_region.name on which to run the script (optional)
                   if not specified, will run on the entire genome
                   recommend to leave this option out and run on entire genome

  -coord_system_version            Name of the default assembly / coord system to use

  -refseq_logicname      List of logic_names for which to fetch the genes (optional)
                         if not specified, will fetch all genes
                         recommend to use this

  -threshold      Percentage longer / shorter for fuzzy matching

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Attribute;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(is_slice_name get_biotype_groups);

# ensembl genes and dna 
my $ensemblhost;
my $ensembluser;  
my $ensembldbname; 
my $ensemblport; 
my $ensemblpass;

# refseq genes
my $refseqhost;
my $refsequser; 
my $refseqdbname;
my $refseqport;  
my $refseqpass;

# dna genes
my $dnahost;
my $dnauser; 
my $dnadbname; 
my $dnaport;  
my $dnapass;

# biotype groupings
my $productionhost; 
my $productionuser; #'ensro';
my $productionport; # 3306;
my $productiondbname; # 'ensembl_production';
my $productionpass;

# output
my $outfile = 'stdout';
my $verbose;

# assembly
my $coord_system_version; #'GRCh37';
my $coord_system_name; # 'toplevel' or 'chromosome';
my $refseq_logicname; # 'refseq_import' or 'refseq_human_import';

# region
my $chr_num;

# fuzziness (in percent. ie. 10 means 10%)
my $threshold = 10;

# to store or not to store
my $store;

#  ~~~~~~
#  No changes required below this point
#  ~~~~~~

# where do we get the ens structures from:
$| = 1;

&GetOptions(
  'ensemblhost:s'          => \$ensemblhost,
  'ensembluser:s'          => \$ensembluser,
  'ensembldbname:s'        => \$ensembldbname,
  'ensemblport:n'          => \$ensemblport,
  'ensemblpass:s'          => \$ensemblpass,
  'refseqhost:s'           => \$refseqhost,
  'refsequser:s'           => \$refsequser,
  'refseqport:s'           => \$refseqport,
  'refseqpass:s'           => \$refseqpass,
  'refseqdbname:s'         => \$refseqdbname,
  'dnahost:s'              => \$dnahost,
  'dnauser:s'              => \$dnauser,
  'dnadbname:s'            => \$dnadbname,
  'dnaport:s'              => \$dnaport,
  'dnapass:s'              => \$dnapass,
  'prodhost:s'              => \$productionhost,
  'produser:s'              => \$productionuser,
  'proddbname:s'            => \$productiondbname,
  'prodport:s'              => \$productionport,
  'prodpass:s'              => \$productionpass,
  'coord_system_version:s' => \$coord_system_version,
  'refseq_logicname:s'     => \$refseq_logicname,
  'outfile:s'              => \$outfile,
  'chr_num:s'              => \$chr_num,
  'threshold:n'            => \$threshold,
  'store'                  => \$store,
);


# a few checks
if (!defined $dnahost || !defined $dnauser || !defined $dnaport || !defined $dnadbname) {
  throw("dna database options required: -dnahost -dnauser -dnadbname -dnaport");
}
if (!defined $ensemblhost || !defined $ensembluser || !defined $ensemblport || !defined $ensembldbname) {
  throw("ensembl database options required: -ensemblhost -ensembluser -ensemblpass -ensembldbname -ensemblport");
}
if (!defined $refseqhost || !defined $refsequser || !defined $refseqport || !defined $refseqdbname) {
  throw("refseq database options required: -refseqhost -refsequser -refseqdbname -refseqport");
}
if (!defined $refseq_logicname) {
  throw("Please specify -refseq_logicname");
}
if (!defined $coord_system_version) {
  print STDERR "No -coord_system_version specified, so running across default\n";
}
if (!defined $chr_num) {
  print STDERR "Running across all toplevel $coord_system_version sequences\n";
}
if (!defined $threshold) {
  print STDERR "You did not define a theshold. Setting it to zero ie. exact matches required. (Recommended to use 10 to 20 percent.)\n";
  $threshold = 0;
}
my $factor = 1 + ($threshold/100);

# print out what we are working with
print STDERR "ENSEMBL database: name $ensembldbname host $ensemblhost port $ensemblport\n";
print STDERR "REFSEQ database: name $refseqdbname host $refseqhost port $refseqport\n";
print STDERR " REFSEQ logic: $refseq_logicname\n";
print STDERR "\n\n\n";


# connect to dbs and get adaptors
my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $dnahost,
  -user   => $dnauser,
  -port   => $dnaport,
  -dbname => $dnadbname,
  -pass   => $dnapass,

);

my $ensembldb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $ensemblhost,
  -user   => $ensembluser,
  -port   => $ensemblport,
  -pass   => $ensemblpass,
  -dbname => $ensembldbname
);
$ensembldb->dnadb($dnadb);
my $ensemblsa = $ensembldb->get_SliceAdaptor();
my $ensemblaa = $ensembldb->get_AttributeAdaptor();

my $refseqdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $refseqhost,
  -user   => $refsequser,
  -port   => $refseqport,
  -dbname => $refseqdbname,
  -pass   => $refseqpass,
);
$refseqdb->dnadb($dnadb);
my $refseqsa   = $refseqdb->get_SliceAdaptor();
my $refseqga   = $refseqdb->get_GeneAdaptor(); 

# production database
my $productiondb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#my $productiondb = new Bio::EnsEMBL::Production::DBSQL::DBAdaptor(
  -host   => $productionhost,
  -user   => $productionuser,
  -port   => $productionport,
  -dbname => $productiondbname,
);


# open outfile
my $fh;
if ($outfile && $outfile ne "stdout") {
  open FH,">$outfile" or die "couldn't open file ".$outfile." $!";
  $fh = \*FH;
} else {
  $fh = \*STDOUT;
}

# fetch biotype groupings
my $ensembl_biotype_groups = get_biotype_groups($ensembldb, 'core');

my $refseq_biotype_groups = get_biotype_groups($ensembldb, 'otherfeatures');


# # #
# OK, now we begin to do stuff
# # #

# fetch slice from ensembl database
my $ensembl_slice;
if (is_slice_name($chr_num)) {
  $ensembl_slice  = $ensemblsa->fetch_by_name($chr_num);
}
else {
  $ensembl_slice  = $ensemblsa->fetch_by_region($coord_system_name, $chr_num,undef,undef,undef,$coord_system_version);
}

throw("Could not fetch the region") unless ($ensembl_slice);
# now we fetch all ensembl genes
print STDERR "Fetching ENSEMBL genes...\n";
my @ensembl_genes = @{do_gene_fetch($ensembl_slice)};
print STDERR "Got " . scalar(@ensembl_genes) . " genes on ".$ensembl_slice->name."\n";


my $num_ensembl_matched = 0;
foreach my $ensembl_gene ( @ensembl_genes ) {
  # get refseq genes
  my $refseq_slice = $refseqsa->fetch_by_name( get_genomic_location($ensembl_gene) );
  my $refseq_genes = $refseqga->fetch_all_by_Slice( $refseq_slice, $refseq_logicname ); 
  print STDERR "Got " . scalar(@$refseq_genes) . " refseq genes on ".$refseq_slice->name." (Ensembl gene ".$ensembl_gene->stable_id.")\n";

  # OK now we actually compare them...
  my $attribute;
  my ($found, $matches) = confirm_refseq_in_ensembl($ensembl_gene,$refseq_genes);
  if ( scalar(@$refseq_genes) < 1 ){
    print STDERR "MISSING no refseq overlapping with ensembl stable_id ".$ensembl_gene->stable_id."\n";
    $attribute = make_attrib('No overlapping RefSeq annotation found');
  } elsif ($found < 1){
    print STDERR "MISSING overlapping refseq does not match with ensembl stable_id ".$ensembl_gene->stable_id."\n";
    $attribute = make_attrib('Overlapping RefSeq annotation not matched');
  } else {
    print STDERR "FOUND $found refseq overlapping with ensembl stable_id ".$ensembl_gene->stable_id."\n";
    $num_ensembl_matched ++;
    foreach my $id (keys %$matches) {
      if ($ensembl_biotype_groups->{$ensembl_gene->biotype} eq $refseq_biotype_groups->{$matches->{$id}->biotype}) {
        print STDERR "  matching biotype ensembl ".$ensembl_gene->biotype." (".$ensembl_biotype_groups->{$ensembl_gene->biotype} .
                     ") vs refseq ".$matches->{$id}->biotype." (".$refseq_biotype_groups->{$matches->{$id}->biotype}.")\n";
        $attribute = make_attrib("Overlapping RefSeq Gene ID ".$matches->{$id}->stable_id." matches and has similar biotype of ".$matches->{$id}->biotype);
      }  else {
        print STDERR "  nonmatched biotype ensembl ".$ensembl_gene->biotype." (".$ensembl_biotype_groups->{$ensembl_gene->biotype} .
                     ") vs refseq ".$matches->{$id}->biotype." (".$refseq_biotype_groups->{$matches->{$id}->biotype}.")\n";
        $attribute = make_attrib("Overlapping RefSeq Gene ID ".$matches->{$id}->stable_id." matches but different biotype of ".$matches->{$id}->biotype);
      }
    }
  }

  if ($store) {
    $ensemblaa->store_on_Gene($ensembl_gene, [$attribute]);
  }
}

print STDERR "Found $num_ensembl_matched of ".(scalar(@ensembl_genes))." ensembl genes matched by refseq\nDONE\n";

# # #
# now we only fetch the refseq gene logic names that we are interested in
# # #

sub make_attrib {
  my ($comment) = @_;

  my $attribute = Bio::EnsEMBL::Attribute->new
       (-CODE => 'refseq_compare',
        -NAME => 'refseq_compare',
        -DESCRIPTION => '',
        -VALUE => $comment);
  return $attribute;
}

sub confirm_refseq_in_ensembl {
  my ($ensembl_gene,$refseq_genes) = @_;

  my $consider_as_matched = 0;
  my $match_strand;
  my $match_span_length;
  my $match_median_length;
  my $match_start;
  my $match_end;
  my $match_name;
  my %matched_genes;

  # nothing fancy going on here
  # just trying to figure out if the ensembl and refseq genes are the 'same'

  REFSEQ: foreach my $refseq_gene (@{$refseq_genes}) {
    # check strand
    if ( $ensembl_gene->strand != $refseq_gene->strand ) {
      print STDERR "    strand_mismatch ".$ensembl_gene->stable_id." ".$ensembl_gene->strand." vs ".$refseq_gene->stable_id." ".$refseq_gene->strand."\n";
      $match_strand = 0;
      next REFSEQ;
    } else {
      $match_strand = 1;
    }


    # check length but be fuzzy... need extra options for single exon genes??
    # NOTE sometimes a gene's length is artificially long because of one dodgy transcript. Try to compensate for that my using the median or avergae or something...
    # BUT note that havana will often annotate one nice long transcript and many short ones so the median is also not always good to use
    my $ensembl_genomic_span_length = get_genomic_span_length($ensembl_gene);
    my $refseq_genomic_span_length = get_genomic_span_length($refseq_gene);
    if ( $ensembl_genomic_span_length*$factor < $refseq_genomic_span_length ) {
      print STDERR "    refseq_median_length_longer ".$ensembl_gene->stable_id." ".$ensembl_genomic_span_length." vs ".$refseq_gene->stable_id." ".$refseq_genomic_span_length."\n";
      $match_span_length = 0;
    } elsif ( $ensembl_genomic_span_length > $refseq_genomic_span_length*$factor ) {
      print STDERR "    ensembl_median_length_longer ".$ensembl_gene->stable_id." ".$ensembl_genomic_span_length." vs ".$refseq_gene->stable_id." ".$refseq_genomic_span_length."\n";
      $match_span_length = 0;
    } else {
      print STDERR "    median_lengths_similar ".$ensembl_gene->stable_id." ".$ensembl_genomic_span_length." vs ".$refseq_gene->stable_id." ".$refseq_genomic_span_length."\n";
      $match_span_length = 1;
    }


    my $ensembl_median_gene_length = get_gene_length_median($ensembl_gene);
    my $refseq_median_gene_length = get_gene_length_median($refseq_gene);
    if ( $ensembl_median_gene_length*$factor < $refseq_median_gene_length ) {
      print STDERR "    refseq_median_length_longer ".$ensembl_gene->stable_id." ".$ensembl_median_gene_length." vs ".$refseq_gene->stable_id." ".$refseq_median_gene_length."\n";
      $match_median_length = 0;
    } elsif ( $ensembl_median_gene_length > $refseq_median_gene_length*$factor ) {
      print STDERR "    ensembl_median_length_longer ".$ensembl_gene->stable_id." ".$ensembl_median_gene_length." vs ".$refseq_gene->stable_id." ".$refseq_median_gene_length."\n";
      $match_median_length = 0;
    } else {
      print STDERR "    median_lengths_similar ".$ensembl_gene->stable_id." ".$ensembl_median_gene_length." vs ".$refseq_gene->stable_id." ".$refseq_median_gene_length."\n";
      $match_median_length = 1;
    }


    my $ensembl_gene_length; 
    my $refseq_gene_length;
    if ( $match_median_length == 1) {
      $ensembl_gene_length = $ensembl_median_gene_length;
      $refseq_gene_length = $refseq_median_gene_length;
    } else {
      $ensembl_gene_length = $ensembl_genomic_span_length;
      $refseq_gene_length = $refseq_genomic_span_length;
    }

    # check coordinates but be fuzzy... need extra options for single exon genes??
    if ( $ensembl_gene->seq_region_start + ($ensembl_gene_length/$threshold) < $refseq_gene->seq_region_start ) {
      print STDERR "    ensembl_start_earlier ".$ensembl_gene->stable_id." ".$ensembl_gene->seq_region_start." vs ".$refseq_gene->stable_id." ".$refseq_gene->seq_region_start."\n";
      $match_start = 0;
    } elsif ($ensembl_gene->seq_region_start - ($ensembl_gene_length/$threshold) > $refseq_gene->seq_region_start ) {
      print STDERR "    ensembl_start_later ".$ensembl_gene->stable_id." ".$ensembl_gene->seq_region_start." vs ".$refseq_gene->stable_id." ".$refseq_gene->seq_region_start."\n";
      $match_start = 0;
    } else {
      print STDERR "    ensembl_start_similar ".$ensembl_gene->stable_id." ".$ensembl_gene->seq_region_start." vs ".$refseq_gene->stable_id." ".$refseq_gene->seq_region_start."\n";
      $match_start = 1;
    }


    if ($ensembl_gene->seq_region_end + ($ensembl_gene_length/$threshold) < $refseq_gene->seq_region_end ) {
      print STDERR "    ensembl_ends_earlier ".$ensembl_gene->stable_id." ".$ensembl_gene->seq_region_end." vs ".$refseq_gene->stable_id." ".$refseq_gene->seq_region_end."\n";
      $match_end = 0;
    } elsif ($ensembl_gene->seq_region_end - ($ensembl_gene_length/$threshold) > $refseq_gene->seq_region_end ) {
      print STDERR "    ensembl_ends_later ".$ensembl_gene->stable_id." ".$ensembl_gene->seq_region_end." vs ".$refseq_gene->stable_id." ".$refseq_gene->seq_region_end."\n";
      $match_end = 0;
    } else {  
      print STDERR "    ensembl_end_similar ".$ensembl_gene->stable_id." ".$ensembl_gene->seq_region_end." vs ".$refseq_gene->stable_id." ".$refseq_gene->seq_region_end."\n";
      $match_end = 1;
    }


    # check gene name
    my $refseq_xrefs = $refseq_gene->get_all_DBEntries('RefSeq_gene_name');

    # have seen cases where sometimes there is more than one xref and one is right and one wrong
    my @refseq_xrefs_shortlist;
    my @alternative_names = @{$ensembl_gene->slice->get_all_synonyms('RefSeq_genomic')};
    my $nc_name = '';
    foreach my $alt_name (@alternative_names) {
      if ($alt_name =~ /NC/) {
        $nc_name = $alt_name;
      }
    }
    foreach my $x (@{$refseq_xrefs}) {
      if ($x->primary_id =~ /$nc_name/) {
        unshift @refseq_xrefs_shortlist, $x;
      } else {
        push @refseq_xrefs_shortlist, $x;
      }
    }

    if ($ensembl_gene->external_name ne $refseq_xrefs_shortlist[0]->display_id) {
    #if ($ensembl_gene->external_name ne $refseq_gene->external_name) {
      print STDERR "    name_mismatch ".$ensembl_gene->stable_id." ".$ensembl_gene->external_name." vs ".$refseq_gene->stable_id." ".$refseq_xrefs_shortlist[0]->display_id."\n";  
      #print STDERR "    name_mismatch ".$ensembl_gene->stable_id." ".$ensembl_gene->external_name." vs ".$refseq_gene->stable_id." ".$refseq_gene->external_name."\n";  
      $match_name = 0;
    } else {
      print STDERR "    names_match ".$ensembl_gene->stable_id." ".$ensembl_gene->external_name." vs ".$refseq_gene->stable_id." ".$refseq_xrefs_shortlist[0]->display_id."\n";
      #print STDERR "    names_match ".$ensembl_gene->stable_id." ".$ensembl_gene->external_name." vs ".$refseq_gene->stable_id." ".$refseq_gene->external_name."\n";
      $match_name = 1;
    }


   # add up the scores... kinda
    if ( $match_strand == 1 && ($match_span_length == 1 || $match_median_length == 1) && $match_start == 1 && $match_end == 1) {
      $consider_as_matched ++;
      $matched_genes{$refseq_gene->stable_id} = $refseq_gene;
    } elsif ( $match_strand == 1 && $match_name == 1 ) {
      $consider_as_matched ++;
      $matched_genes{$refseq_gene->stable_id} = $refseq_gene;
    }
  }

  return ($consider_as_matched, \%matched_genes);
}


=head2 do_gene_fetch 

  Example    : my $genes = do_gene_fetch($slice, undef, \@biotypes, undef); 
  Description: Fetches genes on a slice having a specific biotype
  Returntype : Arrayref of Gene objects 
  Exceptions : none

=cut

sub do_gene_fetch {
  my ($slice, $logic) = @_;
  my $genes = $slice->get_all_Genes($logic);
  return $genes;
}

=head2 get_genomic_location 

  Example    : my $genomic_location = get_genomic_location($slice->name, $cluster);
  Description: Gets the genomic location of a cluster, as a slice name 
  Returntype : String 
  Exceptions : none

=cut

sub get_genomic_location {
  my ($gene) = @_;

  my @genomic_fields = split(":",$gene->slice->name);
  my $genomic_location = $genomic_fields[0].":".$genomic_fields[1].":".
                         $genomic_fields[2].":".$gene->seq_region_start.":".
                         $gene->seq_region_end.":".$genomic_fields[5];

  return $genomic_location;
}

sub get_genomic_span_length {
  my ($gene) = @_;
  my $median;

  my @transcript_lengths;
  foreach my $transcript (@{$gene->get_all_Transcripts}  ) {
    # this gives genomic span
    push @transcript_lengths, $transcript->end - $transcript->start;
    # this gives concatenated exon lengths
    #push @transcript_lengths, $transcript->end - $transcript->start;
  }

  my @sorted_transcript_lengths = sort by_number @transcript_lengths;
  return $sorted_transcript_lengths[-1];
}


sub get_gene_length_median {
  my ($gene) = @_;
  my $median;
  
  my @transcript_lengths;
  foreach my $transcript (@{$gene->get_all_Transcripts}  ) {
    # this gives genomic span
    push @transcript_lengths, $transcript->end - $transcript->start;
    # this gives concatenated exon lengths
    #push @transcript_lengths, $transcript->end - $transcript->start;
  }

  my $mid = int(scalar(@transcript_lengths)/2);
  my @sorted_transcript_lengths = sort by_number @transcript_lengths;
  if (@transcript_lengths % 2) {
      $median = $sorted_transcript_lengths[ $mid ];
  } else {
      $median = ($sorted_transcript_lengths[$mid-1] + $sorted_transcript_lengths[$mid])/2;
  } 
  return $median;
}


sub by_number {
  if ($a < $b) { -1 } elsif ($a > $b) { 1 } else { 0 };
}

