package Bio::EnsEMBL::Analysis::RunnableDB::HiveDownloadUniProtFiles;

use strict;
use warnings;
use feature 'say';

use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable');

sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my $self = shift;

  foreach my $query_hash (keys($self->param('queries'))) {
    $query_hash = ${$self->param('queries')}{$query_hash};
    say Dumper($query_hash);
    my $query_url = $self->build_query($query_hash);
    say "Downloading:\n".$query_url."\n";
    #`$query_url`;
  }

  say "Finished downloading UniProt files";
  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub build_query {
  my ($self,$query_hash) = @_;
  my $tax_group = 0;
  my $taxon_id;
  my $taxonomy;
  my $dest_dir;
  my $file_name;
  my $pe_string = "(";
  my $taxonomy_string = "";
  my $compress = "yes";
  my $fragment_string = "+AND+fragment:no";
  my $mito = "+NOT+organelle%3Amitochondrion";
  my $format = "fasta";

  my $full_query = "wget -q -O - \"http://www.uniprot.org/uniprot/?query=";
  my %pe_code = (
                  '1' => 'evidence+at+protein+level',
                  '2' => 'evidence+at+transcript+level',
                  '3' => 'inferred+from+homology',
                  '4' => 'predicted',
                  '5' => 'uncertain',
                );

  unless(${$query_hash}{'file_name'} && (${$query_hash}{'taxon_id'} || ${$query_hash}{'taxonomy'}) &&
        (${$query_hash}{'dest_dir'} || $self->param('dest_dir')) && ${$query_hash}{'pe_level'}) {
    die "Must provide the following keys:\nfile_name\ntaxon_id or taxonomy\ndest_dir\npe_level";
  }

  $file_name = ${$query_hash}{'file_name'};
  $dest_dir = ${$query_hash}{'dest_dir'};

  my @pe_array = split(',',${$query_hash}{'pe_level'});
  unless(scalar(@pe_array)) {
    die "Not PE levels found in value of pe_levels key. Format should be: '1,2'";
  }
  foreach my $pe_level (@pe_array) {
    unless($pe_level =~ /\d+/) {
     die "Could not parse a PE level from the following: ".$pe_level;
    }

    my $parsed_pe_level = $&;
    unless($parsed_pe_level >= 1 && $parsed_pe_level <= 5) {
     die "Parsed PE level is outside the normal range of 1-5: ".$parsed_pe_level;
   }
   $pe_string .= 'existence%3A%22'.$pe_code{$pe_level}.'%22+OR+';
  }

  $pe_string =~ s/\+OR\+$/\)/;

  say "PE string: ".$pe_string;

  if(${$query_hash}{'taxon_id'}) {
    $taxon_id = ${$query_hash}{'taxon_id'};
    $taxonomy_string = '+AND+taxonomy%3A+'.$taxon_id;
  } elsif(${$query_hash}{'taxonomy'}) {
    $taxonomy = ${$query_hash}{'taxonomy'};
    $taxonomy_string = '+AND+taxonomy%3A'.$taxonomy;
  }

  if(exists(${$query_hash}{'compress'}) && (${$query_hash}{'compress'} eq '0' || ${$query_hash}{'compress'} eq 'no')) {
    $compress = "no";
  }
  if(exists(${$query_hash}{'fragment'}) && (${$query_hash}{'fragment'} eq '1' || ${$query_hash}{'fragment'} eq 'yes')) {
    $fragment_string = "";
  }
  if(exists(${$query_hash}{'format'}) && (${$query_hash}{'format'} ne 'fasta' && ${$query_hash}{'format'} ne '')) {
    $format = ${$query_hash}{'format'};
  }
  if(exists(${$query_hash}{'mito'}) && (${$query_hash}{'mito'} eq '1' || ${$query_hash}{'mito'} eq 'yes')) {
    $mito = "";
  }

  $full_query .= $pe_string.$taxonomy_string.$fragment_string.$mito."&compress=".$compress."&format=".$format.
                 "\" > ".$dest_dir."/".$file_name;
  if($compress eq 'yes') {
    $full_query .= ".gz";
  }

  unless(-e $dest_dir) {
    `mkdir -p $dest_dir`;
  }
  say "FULL QUERY: ".$full_query;
  return($full_query);
}

1;
