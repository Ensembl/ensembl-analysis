#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

=pod

=head1 NAME

  anti_virus.pl

=head1 DESCRIPTION

  Script to identify genes that are built from viral proteins. It produces a summary of the
  affected genes, a list of gene identifiers to remove and a list of supporting protein evidence
  to add to the kill list. If you are checking a database that is already on the website
  you can also generate an HTML click list so you can easily check the  results.

  The script can work either on single or multiple databases:
     Multiple databases require a registry file and a compara database.
     Single database just needs standard database accessors.

  The profiles used to identify the genes are in the domain_kill_list.txt file, these
  consist of primary domain profiles that are used to find the putative viral genes and
  secondary domain profiles that are used as supporting evidence that a gene is viral but
  are not considered to be proof on their own.

  Genes are considered viral if any of the following are true:

  1. They contain a hit to a primary domain and they have at least one 'short' transcript.
     'Short' = genomic span of cds / cds length >= 2

  2. They do not have any short transcripts but all the domains in all the transcripts
     are either primary or secondary domains

  Once the viral genes have been identified and the kill list collated the script then
  finds any  other genes made from proteins on the kill list and applies the length filter
  to them to decide whether to kill them.

=head1 SYNOPSIS

  -regfile    registry file if you wish to run on multiple species

  -compara    compara database for multi db  mode

  -dbname     or use standard db acces for single species ( access is read only )

  -dbport

  -dbhost

  -kill_list  file with the viral domains in

  -click_list make the HTML click list, the database needs to be on the website for this to work!

  -dir        Directory to write output to

=head1 CONTACT

  http://lists.ensembl.org/mailman/listinfo/dev

=cut

use warnings ;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Getopt::Long qw(:config no_ignore_case);

my $compara = $ENV{'COMPARA'};
my $regfile = $ENV{'REGISTRY'};
my $file   = '' ;
my $dbname = '';
my $dbhost = '';
my $dbport = 3306;
my @dbas;
my $html;
my $URL;
my $output;
my $help;
my $dir = ".";

$| = 1;

GetOptions(
	   'regfile:s'   => \$regfile,
	   'compara:s'   => \$compara,
	   'dbname|db|D=s'    => \$dbname,
	   'dbport|port|P=s'    => \$dbport,
	   'dbhost|host|h=s'    => \$dbhost,
	   'kill_list:s' => \$file,
	   'click_list!' => \$html,
	   'dir:s'       => \$dir,
	   'help!'       => \$help,
	  );

unless (($regfile && $compara && $file) or ($dbname && $dbhost && $dbport && $file)){
   print "anti_virus.pl - Current variables are:
regfile    $regfile
compara    $compara
dbname     $dbname
dbport     $dbport
dbhost     $dbhost
kill_list  $file
click_list $html
dir        $dir
";
   exec('perldoc', $0);
   exit 1;
 }

# do registry stuff for multi database mode
if ($regfile && $compara) {
  Bio::EnsEMBL::Registry->load_all($regfile);
  @dbas = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors()};
} else {
  # load the single database
  print STDERR "Opening connection to datbase $dbname @ $dbhost : $dbport\n";

  my  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dbhost,
					       '-user'   => 'ensro',
					       '-dbname' => $dbname,
					       '-port'   => $dbport,
					      );
  push @dbas, $db;
}

die ("No databases found\n") unless (scalar(@dbas) > 0 );

# parse the kill list file first
my ($domains,$supporting_evidence) = kill_list($file);

foreach my $db (@dbas) {
  next unless $db->group eq "core";
  my @genes;
  my $gene_id;

  # adjust the species name to fit the ensembl URL
  $URL = $db->species;
  $URL =~ s/ /_/;

  # for multiple species change the output file name to the species name
  if (scalar(@dbas) == 1) {
    my @species_array = @{sql("SELECT meta_value 
         FROM meta 
         WHERE  meta_key = 'species.classification'
         ORDER BY meta_id
         LIMIT 2;",$db)};
    $URL = $species_array[1]."_".$species_array[0];
  }
  $output = $dir."/".$URL;
  print STDERR "Looking at $URL \n";
  print STDERR "Finding genes by domain\n";

  # Get a list of genes containing hits to the bad domains
  foreach my $domain (keys %$domains) {
    push @genes , @{sql("SELECT DISTINCT(transcript.gene_id)
    FROM protein_feature, transcript, translation
    WHERE hit_name = \'$domain\'
    AND protein_feature.translation_id = translation.translation_id
    AND translation.transcript_id = transcript.transcript_id", $db)};
  }

  next unless scalar(@genes) > 0;

  # make genes non redundant - the same gene can be found
  # several times using different domains
  @genes = @{make_nr(\@genes)};
  print STDERR scalar(@genes) . " genes found by domains\n";

  print STDERR "Filtering genes\n";

  # filter the gene list to get the most likely looking candidates
  # check the domains when assessing the transcripts
  my ($genelist,$killlist) = filter_genes(\@genes,$db,undef,1);

  print STDERR "Expanding geneset using supporting features\n";

  # next expand the search to include other genes built from your killed proteins
  my @extra_genes;

  foreach my $protein (keys %{$killlist}){
    push @extra_genes , @{sql("SELECT DISTINCT(gene.gene_id)
    FROM gene, transcript, transcript_supporting_feature
    LEFT JOIN protein_align_feature
    ON transcript_supporting_feature.feature_id = protein_align_feature.protein_align_feature_id
    WHERE protein_align_feature.protein_align_feature_id = transcript_supporting_feature.feature_id AND
    protein_align_feature.hit_name like \'$protein%\' AND
    transcript_supporting_feature.feature_type = \'protein_align_feature\' AND
    transcript_supporting_feature.transcript_id = transcript.transcript_id AND
    gene.gene_id = transcript.gene_id",$db)};
  }

  # make non redundant with respect to itself and the genes already found
  @extra_genes = @{make_nr(\@extra_genes,\@genes)};

  print STDERR scalar(@extra_genes) . " extra genes found by supporting features\n";

  print STDERR "Filtering genes\n";

  # filter the gene list to get the most likely looking candidates
  # dont check domains
  my ($final_genelist) = filter_genes(\@extra_genes,$db,$genelist,undef);

  print STDERR "Writing data\n\n";

  # write the output
  write_output($final_genelist,$output);

  # write out the killlist
  open (KILL   ,">$output.killlist") or die("Cannot open killlist file $file.killlist\n");
  foreach my $protein ( keys %$killlist ) {
    print KILL "$protein\n" unless $protein =~ /^ENS/;
  }
  close KILL;
}

exit;

#############################################

# query the database directly
sub sql {
  my ($query,$db) = @_;
  my @array;
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  while ( my $stable_id  = $sth->fetchrow_array ) {
    push @array,$stable_id;
  }
  return \@array;
}


# Parse the kill list
sub kill_list {
  my ($file) = @_;
  my %domains;
  my %supporting_evidence;
  open (KILLLIST,"$file") or die ("Cannot open kill list file for $file\n");
  print STDERR "Parsing kill list ";
  while (<KILLLIST>) {
    # parse the domains into a hash
    chomp;
    next if $_ =~ /^#/;
    next unless $_;

    if ($_ =~ /^2\t(\S+)\t(.+)$/) {
      $supporting_evidence{$1} = $2;
    } else {
      if ($_ =~ /^(\S+)\t(.+)$/) {
	$domains{$1} = $2;
      } else {
	print STDERR "\ncannot parse $_\n";
      }
    }
  }
  close KILLLIST;
  print STDERR "...done\n";
  return \%domains,\%supporting_evidence;
}


# Write out the HTML header
sub header {
  print HTML "<!DOCTYPE html PUBLIC >
<head>
  <title>Kill list click list</title>

<STYLE TYPE='text/css'>
<!--
TD{font-family: Verdana; font-size: 10pt;}
--->
</STYLE>
</HEAD>
<BODY>
<TABLE border='1'>
<CAPTION><b>Key:</b></CAPTION>
<TR><TD><b>1</b><TD>  Primary domain - only found in viral proteins.</TR>
<TR><TD><b>2</b><TD>  Secondary domain - associated with viral proteins not sufficient evidence on their own.</TR>
<TR><TD><b>Unknown</b> <TD> Domain is not currently in the kill list.</TR>
</TABLE><br><br>
<TABLE border='1'>
<CAPTION><EM>$URL Genes containing viral domains</EM></CAPTION>
<TR><TD>Num<TD>Gene ID<TD>Transcript ID<TD>Translation ID<TD>Exons<TD>Span / cds<TD>Supporting Evidence<TD>Domains</TR>\n";
}


# Gene tests:
# all transcripts must hit a bad domain
# at least one transcript must have a low ratio of span / cds
# also makes a non-redundant kill list of the protein features from the 'bad' genes
sub filter_genes {
  my ($genes,$db,$genelist,$check_domains) = @_;
  my $pseudogene;
  my $killlist;
  my $failed;
 GENE:  foreach my $gene_id (@$genes){
    my %bad_domains;
    my %dodgy_domains;
    my %total_domains;
    my %unknown_domains;
    my @bad_proteins;
    my $short_trans = 0;
    my $bad_trans = 0;
    my $dodgy_trans = 0;

    # so lets be strict, we want to delete only things where all the transcripts hit a bad domain and at least one of them has the small span
    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);

    if ($gene->biotype eq 'pseudogene'){
      $pseudogene++;
      next GENE;
    }

  TRANS:   foreach my $transcript (@{$gene->get_all_Transcripts}){

      $dodgy_domains{$transcript->dbID}{'count'} = 0;
      $bad_domains{$transcript->dbID}{'count'} = 0;
      $total_domains{$transcript->dbID} = 0;

      # identify the transcript supporing features
      my @sfs = @{$transcript->get_all_supporting_features};
      foreach my $sf (@sfs) {
	if ($sf->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {
	  # make the kill list non redundant
	  if ($sf->hseqname =~ /^(\S+)\.\d+/){
	    push @bad_proteins, $1;
	  } else {
	    push @bad_proteins, $sf->hseqname;
	  }
	}
      }

      # BAD TRANSCRIPTS HAVE THE SHORT SPAN
      if (( $transcript->coding_region_end - $transcript->coding_region_start ) / 
	  ( $transcript->cdna_coding_end - $transcript->cdna_coding_start ) <= 2 ) {
	$short_trans++;
      }

      if ($check_domains){
	
	# EACH TRANSCRIPT NEEDS TO HIT A BAD DOMAIN BEFORE THE GENE IS DELETED
	foreach my $pf (@{$transcript->translation->get_all_ProteinFeatures}) {
	  next unless $pf->analysis->logic_name eq "Superfamily" or
	    $pf->analysis->logic_name eq "Pfam" ;
	  $total_domains{$transcript->dbID}++;


	  # Bad domains
	  if ($domains->{$pf->hseqname}) {
	    $bad_trans++ unless $bad_domains{$transcript->dbID}{'count'} > 0;
	    $bad_domains{$transcript->dbID}{$pf->hseqname} = 1;
	    $bad_domains{$transcript->dbID}{'count'}++;
	    #	    print "BAD TRANS $bad_trans\n";
	  }
	
	  # Dodgy domains
	  if ($supporting_evidence->{$pf->hseqname}) {
	    $dodgy_trans++ unless ($dodgy_domains{$transcript->dbID} or $bad_domains{$transcript->dbID});
	    $dodgy_domains{$transcript->dbID}{$pf->hseqname} = 1;
	    $dodgy_domains{$transcript->dbID}{'count'}++;
	  }
	
	  # Unknown domains
	  unless ($supporting_evidence->{$pf->hseqname} or $domains->{$pf->hseqname}){
	    $unknown_domains{$transcript->dbID}{$pf->hseqname} = 1;
	  }
	}
      } else {
	$bad_trans++;
      }

      # call the transcript dodgy if all the domains it has are bad or dodgy ones
      if ( $total_domains{$transcript->dbID} == ($bad_domains{$transcript->dbID}{'count'} + $dodgy_domains{$transcript->dbID}{'count'}) 
	   && $total_domains{$transcript->dbID} > 1 ){
	$dodgy_trans++;
      }
    }

    # ADD UP ALL THE SCORES FOR THE ENTIRE GENE
    if ($short_trans && $bad_trans ==  scalar(@{$gene->get_all_Transcripts}) ) {
      $genelist->{$gene_id}{'domain'} = \%bad_domains;
      $genelist->{$gene_id}{'gene'} = $gene;
      $genelist->{$gene_id}{'dodgy_domains'} = \%dodgy_domains;
      $genelist->{$gene_id}{'expanded'} = 1 unless $check_domains;
      $genelist->{$gene_id}{'unknown'} = \%unknown_domains;
      foreach my $protein (@bad_proteins){
	$killlist->{$protein} = 1
      }
      next GENE;
    }

    # Dosent pass the span criteria but all of the transcripts hit domains of some sort
    # and *all* the domains are bad in some way
    if ( $dodgy_trans  ==  scalar(@{$gene->get_all_Transcripts}) ) {
      $genelist->{$gene_id}{'domain'} = \%bad_domains;
      $genelist->{$gene_id}{'gene'} = $gene;
      $genelist->{$gene_id}{'dodgy_domains'} = \%dodgy_domains;
      $genelist->{$gene_id}{'dodgy'} = 1;
      $genelist->{$gene_id}{'expanded'} = 1 unless $check_domains;
      $genelist->{$gene_id}{'unknown'} = \%unknown_domains;
      foreach my $protein (@bad_proteins){
	$killlist->{$protein} = 1
      }
      next GENE;
    }

    $failed++;
    $genelist->{$gene_id}{'domain'} = \%bad_domains;
    $genelist->{$gene_id}{'gene'} = $gene;
    $genelist->{$gene_id}{'expanded'} = 1 unless $check_domains;
    $genelist->{$gene_id}{'failed'} = 1;
    $genelist->{$gene_id}{'unknown'} = \%unknown_domains;
  }

  print STDERR "Eliminated $pseudogene pseudogenes\n"  if $pseudogene;
  print STDERR "$failed genes failed to pass the filters\n"  if $failed;

  return $genelist,$killlist;
}


# print to the various text / html fies
sub write_output {
  my ($genelist,$output) = @_;

  # open files
  open (SUMMARY,">$output.summary")  or die("Cannot open summary file $output.summary\n");
  open (HTML,">$output.html")  or die("Cannot open HTML file $output.html\n") if $html;
  open (GENES  ,">$output.gene_ids") or die("Cannot open gene list file $output.gene_ids\n");

  # header for the clicklist
  header() if $html;

  my $gene_num = 0;

  # Genes
  foreach my $gene_id (keys %$genelist) {
    my $count = 0;
    my $gene = $genelist->{$gene_id}{'gene'};
    next if $genelist->{$gene_id}{'failed'};
    $gene_num++;     
    print HTML "<TR><TD>$gene_num<TD><a href= 'http://www.ensembl.org/$URL/geneview?gene=" 
      . $gene->stable_id .";db=core'>".$gene->stable_id.' </a>' if $html;

    print SUMMARY  "Gene id :\t" . $gene->dbID;
    if ($gene->stable_id){
      print GENES $gene->stable_id. "\n";
    } else {
      print GENES $gene->dbID . "\n";
    }
    print SUMMARY  "\t". $gene->stable_id if $gene->stable_id;
    print SUMMARY "\n";
    print SUMMARY "Position :\t" . $gene->feature_Slice->name . "\n\n";

    # Transcripts
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      print SUMMARY  "  Transcript id :\t" . $transcript->dbID;
      print HTML"<TR><TD> <TD> "  if $html &&  $count > 0 ;
      $count++;
      print HTML "<TD><a href= 'http://www.ensembl.org/$URL/transview?;db=core;transcript=" . $transcript->stable_id .";db=core'>".$transcript->stable_id.' </a>'if $html;
      print SUMMARY  "\t". $transcript->stable_id if $transcript->stable_id;
      print SUMMARY "\n";

      # Translations
      if ( my $translation = $transcript->translation ) {
	print SUMMARY "  Translation id :\t" .$translation->dbID;
	print HTML "<TD><a href= 'http://www.ensembl.org/$URL/protview?transcript=" . $transcript->stable_id .";db=core'>".$translation->stable_id.' </a>'if $html ;
	print SUMMARY  "\t". $translation->stable_id if $translation->stable_id;
	print SUMMARY "\n";
      }

      # Exons
      print SUMMARY  "  Exons :\t" . scalar(@{$transcript->get_all_Exons}) . "\n" ;
      print SUMMARY  "  cdna length / genomic span : \t". 
	sprintf("%.2f", ( $transcript->coding_region_end - $transcript->coding_region_start ) /
		( $transcript->cdna_coding_end - $transcript->cdna_coding_start )) . "\n";
      print HTML "<TD>" . scalar(@{$transcript->get_all_Exons}) if $html;
      print HTML "<TD>". sprintf("%.2f", ( $transcript->coding_region_end - $transcript->coding_region_start ) /
				 ( $transcript->cdna_coding_end - $transcript->cdna_coding_start ) ). "<TD>" if $html;

      # protein evidence
      my @sfs = @{$transcript->get_all_supporting_features};
      my %done;
      my $features;
      foreach my $sf (@sfs) {
	if ($sf->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {
	  print SUMMARY "  Supporting feature :\t" . $sf->hseqname."\n";
	  if ($html){
	    print HTML '<a href= http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-newId+-e+[libs%3d{IMGTLIGM%20SWALL%20REFSEQ%20EMBL%20REFSEQP}-all:' 
	      . $sf->hseqname . ']+-vn+2>' . $sf->hseqname . '</a><br>'  unless $done{$sf->hseqname};
	  }
	  $done{$sf->hseqname} = 1;
	  $features .= $sf->hseqname;
	}
      }

      # Domains
      print HTML "<TD>"if $html;
      print SUMMARY "Domains :";

      # Bad
      my $bad_domains = $genelist->{$gene_id}{'domain'};
      foreach my $domain (keys %{$bad_domains->{$transcript->dbID}}) {
	if ($domains->{$domain}){
	  print SUMMARY "\t1:  " . $domain . "\t" . $domains->{$domain} . "\n";
	  print HTML "<B>1:</b> $domain ". $domains->{$domain} . "<br>" if $html;
	}
      }

      # Dodgy
      my $dodgy_domains = $genelist->{$gene_id}{'dodgy_domains'};
      foreach my $domain (keys %{$dodgy_domains->{$transcript->dbID}}) {
	if ($supporting_evidence->{$domain}){
	  print SUMMARY "\t2: ". $domain . "\t" . $supporting_evidence->{$domain} . "\n";
	  print HTML "<B>2:</b> $domain ". $supporting_evidence->{$domain} . "<br>" if $html;
	}
      }

      # Unknown
      my $unknown_domains = $genelist->{$gene_id}{'unknown'};
      foreach my $domain (keys %{$unknown_domains->{$transcript->dbID}}) {
	print SUMMARY "\tUnknown\t". $domain . "\tUnknown\n";
	print HTML "<B>Unknown:</b> $domain<br>" if $html;
      }

      if ( $genelist->{$gene_id}{'expanded'} ){
	print SUMMARY "\tTranscript built from evidence on the kill list\n";
	print HTML '<b>Transcript built from evidence on the kill list<b>' if $html;
      }

      print HTML "</TR>\n" if $html;
    }
    print SUMMARY "\n***************************************************************************************************************************\n\n";
  }

  # close the html
  print HTML "</FORM></BODY>" if $html;
  close SUMMARY;
  close HTML if $html;
  close GENES;
  return;
}


# Make an array non redundant w.r.t itself and / or a second array
sub make_nr {
  my ($array1,$array2) = @_;
  my %hash1;
  my %hash2;
  my @nr_array;
  foreach my $element (@$array1){
    $hash1{$element} = 1;
  }
  foreach my $element (@$array2){
    $hash2{$element} = 1;
  }
  foreach my $key (keys %hash1){
    push @nr_array, $key unless $hash2{$key};
  }
  return \@nr_array
}
