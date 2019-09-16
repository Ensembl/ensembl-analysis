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


=head1 NAME

  label_gencode_basic_transcripts

=head1 DESCRIPTION

  label_gencode_basic_transcripts dumps out a gtf annotation file for a given database

  NOTE: for new biotypes they must be entered in TWO places below:
  (i) in the %known_biotypes hash
  (ii) in the decision tree ie in one of these methods:
  returnBasicCodingAnnotation
  returnBasicNonCodingAnnotation
  returnBasicPseudogeneAnnotation



=head1 OPTIONS

  -host/dbhost     host name for database (gets put as host= in locator)
  -port/dbport     For RDBs, what port to connect to (port= in locator)
  -dbname          For RDBs, what name to connect to (dbname= in locator)
  -user/dbuser     For RDBs, what username to connect as (dbuser= in locator)
  -pass/dbpass     For RDBs, what password to use (dbpass= in locator)

  -dnahost/dnadbhost     host name for dna database (gets put as host= in locator)
  -dnaport/dnadbport     For RDBs, what port to connect to (port= in locator)
  -dnadbname             For RDBs, what name to connect to (dbname= in locator)
  -dnauser/dnadbuser     For RDBs, what username to connect as (dbuser= in locator)

  -path            Name of the assembly
  -verbose         Verbose options adds print statements

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw( min max );
use Carp;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_biotype_groups);

# this ewill help when debugging:
$| = 1;

# # #
# Set up some variables...
# # #
my $host   = '';
my $user;
my $pass   = '';
my $port   = 3306;
my $dbname = '';
my $dnahost;
my $dnauser;
my $dnaport   = 3306;
my $dnadbname;
my $coord_system_name = 'toplevel';
my $coord_system_version;
my $write; # boolean
my $verbose; # boolean
my $code = 'gencode_basic';

# use most recent
my $production_dbname = 'ensembl_production';
my $production_host;
my $production_port = 3306;
my $production_user;

# # #
# These are the options that can be used on the commandline
# # #
GetOptions(
  'host|h|dbhost:s'        => \$host,
  'user|u|dbuser:s'        => \$user,
  'dbname|db|D:s'          => \$dbname,
  'pass|dbpass|p:s'        => \$pass,
  'port|dbport|P:n'        => \$port,

  'dnahost|dnadbhost:s'    => \$dnahost,
  'dnauser|dnadbuser:s'    => \$dnauser,
  'dnadbname:s'            => \$dnadbname,
  'dnaport|dnadbport:n'    => \$dnaport,

  'prodhost|proddbhost:s'    => \$production_host,
  'produser|proddbuser:s'    => \$production_user,
  'proddbname:s'            => \$production_dbname,
  'prodport|proddbport:n'    => \$production_port,

  'path|cs_version:s'      => \$coord_system_version,
  'write'                  => \$write,
  'verbose'                => \$verbose,
);


# # #
# Couple of print outs / checks
# # #
if ($write) {
  print STDERR "We are going to WRITE attributes to the database\n";
} else {
  print STDERR "We are NOT going to write attributes to the database\n";
}



# # #
# Connect to databases
# # #
my $production_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host   => $production_host,
  -user   => $production_user,
  -port   => $production_port,
  -dbname => $production_dbname
);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -pass   => $pass,
  -dbname => $dbname
);
if ($dnadbname) {
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $dnahost,
                                                 -user => $dnauser,
                                                 -port => $dnaport,
                                                 -dbname => $dnadbname,
                                                 );
  $db->dnadb($dnadb);
}
my $sa  = $db->get_SliceAdaptor();
my $aa  = $db->get_AttributeAdaptor();



# # #
# Fetch biotype groupings
# # #
my $biotype2group = get_biotype_groups($production_db);

# UGH - horrible to hard code this
# but I want a more automated way of making sure that we get all the biotypes included in this script
# so I want the script to die if the production db has a new biotype in it

my $known_biotypes = {
                     'protein_coding'                     => 'coding',
                     'polymorphic_pseudogene'             => 'coding',
                     'IG_D_gene'                          => 'coding',
                     'IG_J_gene'                          => 'coding',
                     'IG_C_gene'                          => 'coding',
                     'IG_V_gene'                          => 'coding',
                     'IG_LV_gene'                         => 'coding',
                     'TR_C_gene'                          => 'coding',
                     'TR_J_gene'                          => 'coding',
                     'TR_V_gene'                          => 'coding',
                     'TR_D_gene'                          => 'coding',
                     'non_stop_decay'                     => 'coding_second_choice',
                     'nonsense_mediated_decay'            => 'coding_second_choice',
                     'rRNA'                               => 'noncoding',
                     'known_ncrna'                        => 'noncoding',
                     'snoRNA'                             => 'noncoding',
                     'snRNA'                              => 'noncoding',
                     'miRNA'                              => 'noncoding',
                     'antisense'                          => 'noncoding',
                     'Mt_tRNA'                            => 'noncoding',
                     'Mt_rRNA'                            => 'noncoding',
                     'sense_intronic'                     => 'noncoding',
                     'sense_overlapping'                  => 'noncoding',
                     'lincRNA'                            => 'noncoding',
                     'lncRNA'                             => 'noncoding',
                     'macro_lncRNA'                       => 'noncoding_second_choice',
                     'ribozyme'                           => 'noncoding_second_choice',
                     'scaRNA'                             => 'noncoding_second_choice',
                     'scRNA'                              => 'noncoding_second_choice',
                     'sRNA'                               => 'noncoding_second_choice',
                     'vaultRNA'                           => 'noncoding_second_choice',
                     'processed_transcript'               => 'noncoding_second_choice',
                     'misc_RNA'                           => 'noncoding_second_choice',
                     '3prime_overlapping_ncRNA'           => 'noncoding_second_choice',
                     'non_coding'                         => 'noncoding_second_choice',
                     'bidirectional_promoter_lncRNA'      => 'noncoding_second_choice',
                     'transcribed_processed_pseudogene'   => 'pseudogene_transcribed',
                     'transcribed_unitary_pseudogene'     => 'pseudogene_transcribed',
                     'transcribed_unprocessed_pseudogene' => 'pseudogene_transcribed',
                     'pseudogene'                         => 'pseudogene',
                     'processed_pseudogene'               => 'pseudogene',
                     'unprocessed_pseudogene'             => 'pseudogene',
                     'translated_processed_pseudogene'    => 'pseudogene',
                     'translated_unprocessed_pseudogene'  => 'pseudogene',
                     'unitary_pseudogene'                 => 'pseudogene',
                     'IG_pseudogene'                      => 'pseudogene',
                     'IG_C_pseudogene'                    => 'pseudogene',
                     'IG_D_pseudogene'                    => 'pseudogene',
                     'IG_J_pseudogene'                    => 'pseudogene',
                     'IG_pseudogene'                      => 'pseudogene',
                     'IG_V_pseudogene'                    => 'pseudogene',
                     'TR_J_pseudogene'                    => 'pseudogene',
                     'TR_V_pseudogene'                    => 'pseudogene',
                     'Mt_rRNA_pseudogene'                 => 'pseudogene',
                     'miRNA_pseudogene'                   => 'pseudogene',
                     'misc_RNA_pseudogene'                => 'pseudogene',
                     'rRNA_pseudogene'                    => 'pseudogene',
                     'scRNA_pseudogene'                   => 'pseudogene',
                     'snRNA_pseudogene'                   => 'pseudogene',
                     'snoRNA_pseudogene'                  => 'pseudogene',
                     'tRNA_pseudogene'                    => 'pseudogene',
                     'retained_intron'                    => 'problem',
                     'TEC'                                => 'problem',
                     'ambiguous_orf'                      => 'problem',
                     'disrupted_domain'                   => 'problem',
                     'LRG_gene'                           => 'do_not_use',
                     };





# # #
# Fetch the sequences we are interested in - all or subset
# # #
my @slices = @{ $sa->fetch_all( $coord_system_name, $coord_system_version, 1, undef ) };
print STDERR "Got ".( scalar( @slices) )." slices\n" if $verbose;



# delete old attribs
if ( $write ) {
  print STDERR " Deleting old attributes...\n" if $verbose;
  delete_old_attrib($db, $code);
}

# # #
# biotype check
# # #
foreach my $biotype (@{get_distinct_gene_biotypes($db)}) {
  if (!exists $known_biotypes->{$biotype}) {
    throw("Biotype $biotype not known to this script");
  }
}


# # #
# Now loop through each slices
# Clear up the old attributes, if any
# and then each gene on the slice
# # #
foreach my $slice ( @slices ) {
#next if $slice->seq_region_name ne '10';
  print STDERR "Doing slice ".$slice->seq_region_name."\n" if $verbose;
  my $gene_cnt = 0;
  my $transc_cnt = 0;

  # now look for new candidates
  foreach my $gene ( @{$slice->get_all_Genes} ) {
#next if $gene->stable_id ne 'ENSG00000099251';
    print STDERR "Gene ".$gene->stable_id."\n";
    # check biotype
    if (!exists $biotype2group->{$gene->biotype}) {
      throw("Gene biotype ".$gene->biotype." not known in production database");
    }
    if (!exists $known_biotypes->{$gene->biotype}) {
      throw("Gene biotype ".$gene->biotype." not known in this script");
    }

    $gene_cnt++;
    my $basic_transcripts = giveMeBasicAnnotationTranscripts($gene);

    foreach my $transcript ( @$basic_transcripts ) {
      $transc_cnt++;
      print STDERR "  Gene ".$gene->stable_id." has Basic transcript ".$transcript->stable_id."\n";
      if ($write) {
        store_attrib($aa, $transcript, $code);
        print STDERR "  writing ".$gene->stable_id." has Basic transcript ".$transcript->stable_id."\n";
      }
    }
  }
  print STDERR "Slice ".$slice->seq_region_name." has genes $gene_cnt with $transc_cnt basic transcripts\n";
}
print STDERR "DONE!\n\n";

sub delete_old_attrib {
  my ($db, $code) = @_;

  my $sql = q{
    DELETE ta
    FROM transcript_attrib ta, attrib_type att
    WHERE att.attrib_type_id = ta.attrib_type_id
    AND att.code = ? };

  my $sth = $db->dbc->prepare($sql);
  $sth->bind_param( 1, $code);
  $sth->execute();
  $sth->finish();

  return;
}

sub store_attrib {
  my ($aa, $transcript, $code) = @_;

  my ($attrib_type_id, $newcode, $name, $description ) = $aa->fetch_by_code($code);
  if (!$attrib_type_id) {
    throw("Unable to fetch attrib_type with code $code");
  }
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -NAME        => $name,
    -CODE        => $code,
    -VALUE       => 'GENCODE basic', # text preferred over boolean for BioMart
    -DESCRIPTION => $description
  );
  $aa->store_on_Transcript($transcript, [$attrib]);

  return;
}



# there are 5 methods in this file. the method that gives the basic transcripts is the first one (giveMeBasicAnnotationTranscripts) and it calls the other 4 methods.
# these are from electra

sub giveMeBasicAnnotationTranscripts{  # as a parameter I give a gene object and it returns a reference of a list with the gene's basic annotation transcript objects.

  my $gene=shift;
  my @transcripts=@{$gene->get_all_Transcripts};

  my $fullLengthTag;
  my @ncTrlengths;

  my $maxGenLength=0;
  my $maxGenLengthTr;

  my @ncProblemlengths;

  my @basicCodingAnnotationTranscripts=@{returnBasicCodingAnnotation($gene)};
  my @basicNonCodingAnnotationTranscripts=@{returnBasicNonCodingAnnotation($gene)};
  my @basicPseudogeneAnnotationTranscripts = @{returnBasicPseudogeneAnnotation($gene)};

  my @basicAnnotationTranscripts;
  if (@basicCodingAnnotationTranscripts) {
    @basicAnnotationTranscripts = @basicCodingAnnotationTranscripts;
  } elsif (@basicPseudogeneAnnotationTranscripts) {
    @basicAnnotationTranscripts = @basicPseudogeneAnnotationTranscripts;
  } else {
    @basicAnnotationTranscripts = @basicNonCodingAnnotationTranscripts;
  }

  if(scalar @basicAnnotationTranscripts ==0){   # if I dont find any (nor coding neither non-coding) I consider these biotypes which are in the UCSC the problem category!


    foreach my $transcript (@transcripts){

      if($known_biotypes->{$transcript->biotype} eq 'problem') {
#      if($transcript->biotype eq "retained_intron" or $transcript->biotype eq "TEC" or $transcript->biotype eq "ambiguous_orf" or
#         $transcript->biotype eq "disrupted_domain"){

        if(scalar @ncProblemlengths ==0){

          $ncProblemlengths[0]=$transcript;

        }else{

          my $tr=$ncProblemlengths[0];
          if($tr->length < $transcript->length){
            @ncProblemlengths=();
            push(@ncProblemlengths,$transcript);
          }elsif($tr->length == $transcript->length){

            push(@ncProblemlengths,$transcript);    # I end up having in this array all transcript object with same biggest length
          }

        }

      }
    }

    if(scalar @ncProblemlengths ==1 ){ # if I get 1 that has the biggest length

      push(@basicAnnotationTranscripts,$ncProblemlengths[0]);

    } elsif(scalar @ncProblemlengths >1 ){ #if we ended up having more than one transcript with same biggest length,I iterate through the list to get only 1,the one with the biggest genomic length
      $maxGenLength=0;

      foreach my $sameLengthTr (@ncProblemlengths){

        my $genLengthCurrentTr=getGenomicLengthOfTranscript($sameLengthTr);

        if($genLengthCurrentTr > $maxGenLength){

          $maxGenLength=$genLengthCurrentTr;
          $maxGenLengthTr=$sameLengthTr;
        }
      }

      push(@basicAnnotationTranscripts,$maxGenLengthTr);

    }


  }
  if (scalar(@basicAnnotationTranscripts) == 0) {
    throw("Gene ".$gene->stable_id." with ".(scalar(@transcripts))." transcripts has no transcripts tagged in Gencode Basic");
  }

  return \@basicAnnotationTranscripts;

}  #end of method!!!

=head2 returnBasicCodingAnnotation
 
 Arg [1]    : Bio::EnsEMBL::Gene
 Example    : $self->returnBasicCodingAnnotation($gene);
 Description: It returns a reference to a list containing the basic coding transcript objects of the given gene.
              Criteria to select the basic coding transcripts can be found at the end of this file.
 Returntype : listref of Bio::EnsEMBL::Transcript
 Exceptions : None.
 
=cut

sub returnBasicCodingAnnotation {

  my $gene=shift;
  my @transcripts=@{$gene->get_all_Transcripts};
  my @basicCodingAnnotation;

  my $fullLengthTag;
  my @cdsLongest;
  my @longest;

  my $maxGenLength=0;
  my $maxGenLengthTr;

  foreach my $transcript (@transcripts){  # pushes the protein-coding full length transcripts in the basicAnnotationTranscripts array

    $fullLengthTag=1;

    if($known_biotypes->{$transcript->biotype} eq 'coding') {

      my @attributes=@{$transcript->get_all_Attributes};
      foreach my $attribute (@attributes){
        if(($attribute->code eq "cds_end_NF" and $attribute->value==1) or ($attribute->code eq "cds_start_NF" and $attribute->value==1) or
           ($attribute->code eq "mRNA_end_NF" and $attribute->value==1) or ($attribute->code eq "mRNA_start_NF" and $attribute->value==1)){ # all these attributes signify non full-length transcripts

          $fullLengthTag=0;
        } # end of attrib if-statement
      } # end of foreach

      if($fullLengthTag==1) {
        push (@basicCodingAnnotation,$transcript);  # this transctipt is full-length protein-coding
      }

    } ############## end of pc full length
  }

  if (scalar (@basicCodingAnnotation) > 0 ) {

    return \@basicCodingAnnotation;

  } elsif (scalar (@basicCodingAnnotation) ==0 ){  # if NOT--->it pushes the protein_coding (from the broader pc category) longest CDS partial length transcript

    foreach my $transcript (@transcripts){

      if($known_biotypes->{$transcript->biotype} eq 'coding' or $known_biotypes->{$transcript->biotype} eq 'coding_second_choice') {

        if(scalar @cdsLongest ==0){

          $cdsLongest[0]=$transcript;

        } else {

          my $storedTranscript=$cdsLongest[0];
          my $cdsLengthTranscript=cdsLength($transcript);
          my $cdsLengthStoredTranscript=cdsLength($storedTranscript);

          if($cdsLengthStoredTranscript < $cdsLengthTranscript){
            @cdsLongest=();
            push(@cdsLongest,$transcript);
          }elsif($cdsLengthStoredTranscript == $cdsLengthTranscript){

            push(@cdsLongest,$transcript);    # I end up having in this array all transcript object with same biggest cds length
          }

        }
      }
    }

    if(scalar (@cdsLongest) >0){  # if we have more than 1 transcript with same biggest CDS we have to do some further filtering..

      foreach my $cdsLongestTranscript (@cdsLongest){

        if(scalar (@longest) ==0){

          $longest[0]=$cdsLongestTranscript;

        }else{

          my $storedTr=$longest[0];

          if($storedTr->length < $cdsLongestTranscript->length){
            @longest=();
            push(@longest,$cdsLongestTranscript);
          }elsif($storedTr->length == $cdsLongestTranscript->length){

            push(@longest,$cdsLongestTranscript);    # I end up having in this array all transcript object with same biggest length
          }
        }
      }

      if(scalar (@longest) ==1){

        push(@basicCodingAnnotation,$longest[0]);

      }elsif(scalar @longest >1){  # if we have more than one with the biggest length(exons+introns) we get the one with the longest genomic length (if there is more than one transcript with the longest genomic length we get the first one that we find!)

        foreach my $longestTr (@longest){

          my $genLengthCurrentTr=getGenomicLengthOfTranscript($longestTr);

          if($genLengthCurrentTr > $maxGenLength){

            $maxGenLength=$genLengthCurrentTr;
            $maxGenLengthTr=$longestTr;
          }
        }

        push(@basicCodingAnnotation,$maxGenLengthTr);

      }

    }elsif(scalar (@cdsLongest) ==1){

      push(@basicCodingAnnotation,$cdsLongest[0]);
    }

  }

  return \@basicCodingAnnotation;

} # end of method

sub returnBasicNonCodingAnnotation{  # i call it with a gene object and it returns a reference of a list with the basic non-coding transcript objects of the gene

  my $gene=shift;
  my @transcripts=@{$gene->get_all_Transcripts};
  my @basicNonCodingAnnotation;

  my $fullLengthTag;
  my @ncTrlengths;

  my $maxGenLength=0;
  my $maxGenLengthTr;

  my @ncProblemlengths;

  foreach my $transcript (@transcripts){  # pushes non-coding , well caracterized, FULL length transcripts

    $fullLengthTag=1;

    if($known_biotypes->{$transcript->biotype} eq 'noncoding') {

      my @attributes=@{$transcript->get_all_Attributes};
      foreach my $attribute (@attributes){
        if(($attribute->code eq "mRNA_end_NF" and $attribute->value==1) or ($attribute->code eq "mRNA_start_NF" and $attribute->value==1)){

          $fullLengthTag=0;
        }
      }

      if($fullLengthTag==1) {
         push (@basicNonCodingAnnotation,$transcript);
      }

    }
  } # for every transcript

  if (scalar(@basicNonCodingAnnotation) != 0) {   # if the transcripts are non-coding-well caracterized-full length I am almost finished!
    # return a limited number of transcripts
    @basicNonCodingAnnotation = getGivenCoverUntilExonsCovered($gene,80,@basicNonCodingAnnotation);

    if ($basicNonCodingAnnotation[0]) {
      return \@basicNonCodingAnnotation;
    } else {
      print "No basicNonCodingAnnotation from first non-coding group found.\n";
    }
  }

  if(scalar (@basicNonCodingAnnotation) ==0 ){ # if there are no nc-well caracterized-FULL LENGTH transcripts, I get the transcripts that are non-coding biggest length

    foreach my $transcript (@transcripts){

      if($known_biotypes->{$transcript->biotype} eq 'noncoding' or $known_biotypes->{$transcript->biotype} eq 'noncoding_second_choice') {

       if(scalar @ncTrlengths ==0){

         $ncTrlengths[0]=$transcript;

       }else{

         my $tr=$ncTrlengths[0];
         if($tr->length < $transcript->length){
           @ncTrlengths=();
           push(@ncTrlengths,$transcript);
          }elsif($tr->length == $transcript->length){

            push(@ncTrlengths,$transcript);    # I end up having in this array all transcript object with same biggest length
          }

        }
      }
    }

    if(scalar @ncTrlengths ==1 ){ # if I get 1 that has the biggest length

      push(@basicNonCodingAnnotation,$ncTrlengths[0]);

    }

     elsif(scalar @ncTrlengths >1 ){ #if we ended up having more than one transcript with same biggest length,I iterate through the list to get only 1,the one with the biggest genomic length
      foreach my $sameLengthTr (@ncTrlengths){

        my $genLengthCurrentTr=getGenomicLengthOfTranscript($sameLengthTr);

        if($genLengthCurrentTr > $maxGenLength){

          $maxGenLength=$genLengthCurrentTr;
          $maxGenLengthTr=$sameLengthTr;
        }
      }

      push(@basicNonCodingAnnotation,$maxGenLengthTr);

    }

  }
################################


  return \@basicNonCodingAnnotation; # the array may have 0, or 1 transcript! 0: if there is not such transcript with nc biotype, 1: the longest(exons) tr or the genomic longest (introns+exons) tr


} # end of method!

sub returnBasicPseudogeneAnnotation{ # i call it with a gene object and it returns a reference of a list with the basic coding transcript objects of the gene

  my $gene=shift;
  my @transcripts=@{$gene->get_all_Transcripts};
  my @basicPseudogeneAnnotation;

  my $transcribed_pseudogene = 0;
  foreach my $transcript (@transcripts){  # pushes all transcripts but expects only one per gene

    if($known_biotypes->{$transcript->biotype} eq 'pseudogene_transcribed') {
      # 07 August 2014
      # When a transcribed processed pseudogene, take all transcripts [requested by af2]
      $transcribed_pseudogene = 1;
      last;

    }  elsif($known_biotypes->{$transcript->biotype} eq 'pseudogene') {

      push (@basicPseudogeneAnnotation,$transcript);

    }
  } # for every transcript


  if ($transcribed_pseudogene == 1) {
     @basicPseudogeneAnnotation = @{filter_transcribed_pseudogene($gene->get_all_Transcripts, ['transcribed_processed_pseudogene','transcribed_unprocessed_pseudogene','transcribed_unitary_pseudogene'])};
  }

  if(scalar (@basicPseudogeneAnnotation) >1 && $transcribed_pseudogene != 1){   # if the transcripts are non-coding-well caracterized-full length I am finished!
    #throw("Pseudogene ".$gene->stable_id." should have exactly one transcript. This rule was requested by Havana to be imlemented from e76 onwards");
    warning("Pseudogene ".$gene->stable_id." should have exactly one transcript and not ".(scalar(@basicPseudogeneAnnotation))." transcripts");
  }

  return \@basicPseudogeneAnnotation;

}

sub getGivenCoverUntilExonsCovered { # I call it with a set of transcripts and it returns the longest transcripts covering the maximum number of exons until exon_coverage% of all exons (based on exonic length) are covered
  my ($gene,$exon_coverage,@basicTranscripts) = @_;

  my @uncoveredExons = @{$gene->get_all_Exons()};
  my @sortedBasicTranscripts = sort {getScoreExonsCoverAndLength($b,@uncoveredExons) cmp getScoreExonsCoverAndLength($a,@uncoveredExons)} @basicTranscripts;
  my @finalBasicTranscripts;
  my $transcript;
  while (scalar(@sortedBasicTranscripts) >= 1) {
    $transcript = shift(@sortedBasicTranscripts);
    if (getScoreExonsCoverAndLength($transcript,@uncoveredExons) > 0) {
      push @finalBasicTranscripts,$transcript;
    }
    if (transcriptExonsCoverGenePercentage($gene,@finalBasicTranscripts) >= $exon_coverage) {
      last;
    }
    @uncoveredExons = getUncoveredExons($gene,@finalBasicTranscripts);
    @sortedBasicTranscripts = sort {getScoreExonsCoverAndLength($b,@uncoveredExons) cmp getScoreExonsCoverAndLength($a,@uncoveredExons)} @sortedBasicTranscripts;
  }
  return @finalBasicTranscripts;
}

sub getUncoveredExons() {
  # returns an array of the exons in 'gene' which are not covered by any exon in 'transcripts'
  my ($gene,@transcripts) = @_;

  my @uncoveredExons;

  my @allTranscriptExons;
  foreach my $transcript (@transcripts) {
    push @allTranscriptExons,@{$transcript->get_all_Exons()};
  }

  my @uniqueExons = @{$gene->get_all_Exons()};

  my $found = 0;
  UNIQUE: foreach my $geneExon (@uniqueExons) {
  	$found = 0;
    foreach my $transcriptExon (@allTranscriptExons) {
      if ($transcriptExon->start() == $geneExon->start() and
          $transcriptExon->end() == $geneExon->end()) {
        $found = 1;
        last;
      }
    }
    if (!$found) {
      push @uncoveredExons,$geneExon;
    }
  }
  return @uncoveredExons;
}

sub getScoreExonsCoverAndLength {
  # returns the number of exons covered in 'exons' by the exons in 'transcript'
  my ($transcript,@exons) = @_;

  my $numExonsCovered = 0;
  UNIQUE: foreach my $exon (@exons) {
    foreach my $transcriptExon (@{$transcript->get_all_Exons()}) {
      if ($transcriptExon->start() == $exon->start() and
          $transcriptExon->end() == $exon->end()) {
        $numExonsCovered++;
        next UNIQUE;
      }
    }
  }

  my $length_score_percentage = 1-(1/$transcript->length());
  my $score = $numExonsCovered+$length_score_percentage;

  return $score;
}

sub transcriptExonsCoverGenePercentage() {
  # returns the percentage of exonic length in 'gene' covered by the exons in 'transcripts'
  my ($gene,@transcripts) = @_;

  my @allTranscriptExons;
  foreach my $transcript (@transcripts) {
    push @allTranscriptExons,@{$transcript->get_all_Exons()};
  }

  my @uniqueExons = @{$gene->get_all_Exons()};
  my $found = 0;
  my $total_coverage = 0;
  my $total_length = 0;
  my $max_coverage = 0; # different exons can overlap a gene exon so we need to get the maximum coverage
  UNIQUE: foreach my $geneExon (@uniqueExons) {
    foreach my $transcriptExon (@allTranscriptExons) {

      if ($transcriptExon->start() <= $geneExon->end() and
          $transcriptExon->end() >= $geneExon->start()) {

        my $coverage_start = max($transcriptExon->start(),$geneExon->start());
        my $coverage_end = min($transcriptExon->end(),$geneExon->end());
        my $coverage_length = $coverage_end-$coverage_start+1;
        if ($coverage_length > $max_coverage) {
          $max_coverage = $coverage_length;
        }
      } # end if transcriptExon

      if (($max_coverage == $geneExon->length()) or
            ($transcriptExon == $allTranscriptExons[-1])) {
        $total_coverage += $max_coverage;
        $total_length += $geneExon->length();
        $max_coverage = 0; # needs to be reset for the next unique exon within the gene
        next UNIQUE;
      } # end if max_coverage
    } # end foreach
  }
  return (($total_coverage*100)/$total_length);
}

sub getGenomicLengthOfTranscript{ # i call it with a transcript object and it returns its genomic length

  my $transcript=shift;
  my $genomicStart=$transcript->start;
  my $genomicEnd=$transcript->end;
  my $length;

  $length=$genomicEnd-$genomicStart;

  return $length;

}

sub cdsLength{


  my $transcript=shift;
  my $length;
  my $cdsStart=$transcript->cdna_coding_start;
  my $cdsEnd=$transcript->cdna_coding_end;

  $length=$cdsEnd-$cdsStart;

  return $length;
}


sub filter_transcribed_pseudogene {
  my ($transcripts, $biotypes) = @_;
  my %starts;
  my %ends;
  my %keep;

  foreach my $transcript (@$transcripts) {
    push @{$starts{$transcript->seq_region_start}}, $transcript;
    push @{$ends{$transcript->seq_region_end}}, $transcript;

    # collect the biotypes we want
    foreach my $biotype (@$biotypes) {
      if ($transcript->biotype eq $biotype) {
        $keep{$transcript->seq_region_start} = $transcript;
      }
    }
  }
  # collect only lowest start and highest end
  my $min_start = min(keys %starts);
  my $max_end = max(keys %ends);

  foreach my $t (@{$starts{$min_start}}) {
    $keep{$t->dbID} = $t;
  }
  foreach my $t (@{$ends{$max_end}}) {
    $keep{$t->dbID} = $t;
  }

  my @preferred = values %keep;
  return \@preferred;
}


sub sort_by_biotype_exons_and_length {
  my ($transcripts, $biotypes) = @_;
  my @keep_for_biotype;
  my @sorted;
  my %not_seen;
  my %info;


  TRANSCRIPT: foreach my $t (@$transcripts) {
    my $got_it = 0;
    foreach my $biotype (@$biotypes) {
print "tbio ".$t->biotype." vs mybio ".$biotype." \n";
      if ($t->biotype eq $biotype) {
        $got_it = 1;
        push @keep_for_biotype, $t;
       next TRANSCRIPT;
      }
    }
    if ($got_it != 1) {
      $not_seen{$t} =$t;
    }
  }

  foreach my $t (keys %not_seen) {
    my $nums_exons = @{$not_seen{$t}->get_all_Exons};
    my $length = $not_seen{$t}->length;
    push @{$info{$nums_exons}{$length}}, $not_seen{$t} ;
  }

  foreach my $numex (sort by_number_desc keys %info) {
    foreach my $len (sort by_number_desc keys %{$info{$numex}}) {
      push @sorted, @{$info{$numex}{$len}};
      print "pushed transcript with $numex exons and length $len \n";
    }
  }

  print "  keep ".@keep_for_biotype." + sorted ".(scalar(@sorted))." of ".(scalar(@$transcripts))." transcripts\n";
  return (\@keep_for_biotype, \@sorted);
}

sub by_number_desc {
     if ($a < $b) {
        return -1;
    } elsif ($a == $b) {
        return 0;
    } elsif ($a > $b) {
        return 1;
    }
}


sub get_distinct_gene_biotypes {
  my ($db)  = @_;
  my @biotypes;

  my $sql = q{
    SELECT distinct(biotype)
    FROM gene};

  my $sth = $db->dbc->prepare($sql);
  $sth->execute();
  while (my ($biotype) = $sth->fetchrow) {
    push @biotypes, $biotype;
  }
  $sth->finish();

  return \@biotypes;
}

=head1

1) Protein-coding transcripts

a)We have preference to these transcript biotypes :

"IG_D_gene" ,"IG_J_gene" ,"IG_C_gene", "IG_V_gene" ,"protein_coding",
"TR_C_gene","TR_J_gene","TR_V_gene" ,"TR_D_gene", "polymorphic_pseudogene"

 If we have these biotypes, in full-length transcripts we put them ALL
in the basic set.


b)If we don't get any of these,    we check these biotypes:

 "IG_D_gene" , "IG_J_gene" , "IG_C_gene" , "IG_V_gene"
,"protein_coding" , "TR_C_gene","TR_J_gene", "TR_V_gene" ,"TR_D_gene" ,
"non_stop_decay" ,"nonsense_mediated_decay" , "polymorphic_pseudogene"

to get ONE.  We want to put the longest CDS in the basic set.
if there is more than one with longest CDS we put the biggest length in the basic set.
if there is more than one with biggest length we put the longest genomic
length in the basic set.(if there is more than 1 of these too, I get the first I find.)



2) Non - coding transcripts

a)We have preference to these transcript biotypes :
 "antisense", "Mt_tRNA","Mt_rRNA" , "rRNA" ,"snoRNA" , "snRNA", "miRNA"
 If we have these biotypes, in full-length transcripts we put them ALL
in the basic set.


b)If we don't get any of these,  we check these biotypes:

"antisense", "Mt_tRNA", "Mt_rRNA","rRNA" , "snoRNA" ,"snRNA",
"processed_transcript", "lincRNA" ,
"3prime_overlapping_ncrna","non_coding" ,
"sense_intronic" ,"sense_overlapping", "miRNA" , "misc_RNA"
to get ONE.  We want to put the one with the biggest length in the basic set.
if there is more than one with biggest length we put the transcript with
the longest genomic length in the basic set.(if there is more than 1 of these too, I get the first I find.)

We want at least 1 transcript from both categories(basic coding and basic non-coding) if there is any
transcript falling in these groups.

3) If we don't get anything from both 1 and 2, last step is that I check the the UCSC the problem category which are the transcripts with these biotypes:

 "retained_intron", "TEC" , "ambiguous_orf" , "disrupted_domain"

I get 1 transcript which will have the biggest length.
If there is more than 1 with the biggest length I get the one with the biggest genomic length! (if there is more than 1 of these too, I get the first I find.)


=cut




=head1

This bit is from the UCSC website:

GENCODE Basic Set selection: The GENCODE Basic Set is intended to provide a simplified subset of the GENCODE transcript annotations that will be useful to the majority of users. The goal was to have a high-quality basic set that also covered all loci. Selection of GENCODE annotations for inclusion in the basic set was determined independently for the coding and non-coding transcripts at each gene locus.

Criteria for selection of coding transcripts (including polymorphic pseudogenes) at a given locus:
All full-length coding transcripts (except problem transcripts or transcripts that are nonsense-mediated decay) was included in the basic set.
If there were no transcripts meeting the above criteria, then the partial coding transcript with the largest CDS was included in the basic set (excluding problem transcripts).
Criteria for selection of non-coding transcripts at a given locus:
All full-length non-coding transcripts (except problem transcripts) with a well characterized biotype (see below) were included in the basic set.
If there were no transcripts meeting the above criteria, then the largest non-coding transcript was included in the basic set (excluding problem transcripts)..
It no transcripts were included by either the above criteria, the longest problem transcript is included.
Non-coding transcript categorization: Non-coding transcripts are categorized using their biotype and the following criteria:

well characterized: antisense, Mt_rRNA, Mt_tRNA, miRNA, rRNA, snRNA, snoRNA
poorly characterized: 3prime_overlapping_ncrna, lincRNA, misc_RNA, non_coding, processed_transcript, sense_intronic, sense_overlapping

=cut
