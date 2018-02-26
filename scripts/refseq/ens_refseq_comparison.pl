#!/usr/bin/env perl
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

use warnings;
use strict;
use feature 'say';

use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $primarydbname = '';
my $primaryuser   = '';
my $primaryhost   = '';
my $primaryport   = '';
my $primarypass   = '';

my $dnadbname = '';
my $dnauser   = '';
my $dnahost   = '';
my $dnaport   = '';

my $secondarydbname = '';
my $secondaryuser   = '';
my $secondaryhost   = '';
my $secondaryport   = '';

my $primary_logic_name = undef;
my $secondary_logic_name = undef;

my $primary_set_name = '';
my $secondary_set_name = '';

my $input_id;
my $write;
my $debug;

GetOptions(
            'primary_host=s'   => \$primaryhost,
            'primary_port=s'   => \$primaryport,
            'primary_user=s'   => \$primaryuser,
            'primary_dbname=s'     => \$primarydbname,
            'primary_pass=s'     => \$primarypass,

            'secondary_host=s'   => \$secondaryhost,
            'secondary_port=s'   => \$secondaryport,
            'secondary_user=s'   => \$secondaryuser,
            'secondary_dbname=s'     => \$secondarydbname,

            'dna_host=s'        => \$dnahost,
            'dna_port=s'        => \$dnaport,
            'dna_user=s'        => \$dnauser,
            'dna_dbname=s'      => \$dnadbname,

            'primary_set_name=s'  => \$primary_set_name,
            'secondary_set_name=s' => \$secondary_set_name,

            'primary_logic_name=s' => \$primary_logic_name,
            'secondary_logic_name=s' => \$secondary_logic_name,
            'iid=s' => \$input_id,
            'write!' => \$write,
            'debug!' => \$debug,
          );

say "Establishing db connections..." if ($debug);

my $primarydb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $primaryport,
  -user    => $primaryuser,
  -host    => $primaryhost,
  -pass    => $primarypass,
  -dbname  => $primarydbname);

my $secondarydb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -port    => $secondaryport,
    -user    => $secondaryuser,
    -host    => $secondaryhost,
    -dbname  => $secondarydbname);

say "Reading DNA from: ".$dnadbname." : ".$dnahost." : ".$dnaport if ($debug);
my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -port    => $dnaport,
    -user    => $dnauser,
    -host    => $dnahost,
    -dbname  => $dnadbname);

$primarydb->dnadb($dnadb);
$secondarydb->dnadb($dnadb);

say "...finished establishing db connections\n" if ($debug);

# These are a bunch of hashes for storing versions, sequences and dbibs
# I could have made them into a multidimensionsal hash, but this is more readable
# All use accessions as keys

my $attribute_code = 'enst_refseq_compare';
my $attribute_adaptor = $primarydb->get_AttributeAdaptor;
my (undef, undef, $attrib_name, $attrib_description) = $attribute_adaptor->fetch_by_code($attribute_code);
say "Comparing ".$primary_set_name." transcripts to overlapping ".$secondary_set_name." transcripts..." if ($debug);
my $slices;
if ($input_id) {
  $slices = [$primarydb->get_SliceAdaptor->fetch_by_name($input_id)];
}
else {
  $slices = $primarydb->get_SliceAdaptor->fetch_all('toplevel');
}

# Take all Ensembl transcripts and find all overlapping RefSeq transcripts in the core db
# Once these are found it loops through all of them and runs a series of tests to see if any match
# It will only run the tests if either both transcripts have a translation or neither has one
# It writes the results to the three output files below

foreach my $slice (@{$slices}) {
  my $primary_transcripts = $primarydb->get_TranscriptAdaptor->fetch_all_by_Slice($slice, 1, $primary_logic_name);
  foreach my $primary_transcript (@{$primary_transcripts}) {
    # Create a feature slice to get overlapping Ensembl transcripts in the core db
    my $transcript_slice = get_feature_slice_from_db($primary_transcript,$secondarydb);
    next unless($transcript_slice);

    my $secondary_transcripts = $secondarydb->get_TranscriptAdaptor->fetch_all_by_Slice($transcript_slice, 1, $secondary_logic_name);
    unless(scalar(@{$secondary_transcripts})) {
      say "No overlapping ".$secondary_set_name." transcripts found!" if ($debug);
      next;
   }

    my $db_id = $primary_transcript->dbID();
    my $protein_coding_only_match_count = 0;

    if ($debug) {
      say "===================================================================";
      say "Attempting to find ".$secondary_set_name." match for ".$primary_set_name." transcript: ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")";
      say "===================================================================";
    }
    my $primary_full_exon_string = make_exon_string($primary_transcript,1);
    my $primary_whole_seq = $primary_transcript->seq->seq();
    my $primary_cds_exon_string;
    my $primary_cds_seq;
    my $primary_pep_seq;
    my $primary_protein_coding;

    # Set variables if the ensembl transcript has a translation
    if($primary_transcript->translation()) {
      $primary_cds_exon_string = make_exon_string($primary_transcript);
      $primary_cds_seq = $primary_transcript->translateable_seq();
      $primary_pep_seq = $primary_transcript->translate()->seq();
      $primary_protein_coding = 1;
    } else {
      $primary_protein_coding = 0;
    }

    # The array will hold info about what tests the transcripts match on
    my @whole_structure_match;
    my @whole_structure_match_stable_ids;
    my @protein_coding_only_match_stable_ids;
    foreach my $secondary_transcript (@{$secondary_transcripts}) {

      #
      unless($secondary_transcript->strand == $primary_transcript->strand) {
        next;
      }

      if ($debug) {
        say "--------------------------------------------";
        say "Comparing to overlapping ".$secondary_set_name." transcript: ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID.")";
      }
      $secondary_transcript = $secondary_transcript->transfer($slice);

      # Make sure the transform worked
      unless($secondary_transcript) {
        say "Inner issue with transcript, post transform, skipping" if ($debug);
        next;
      }

      unless($secondary_transcript->strand == $primary_transcript->strand) {
        say "Transcripts are on different strands, so not comparing" if ($debug);
        next;
      }

      # Set variables if the RefSeq transcript has a translation
      my $secondary_transcript_cds_exon_string;
      my $secondary_transcript_cds_seq;
      my $secondary_transcript_pep_seq;
      my $secondary_transcript_protein_coding;
      if($secondary_transcript->translation()) {
        $secondary_transcript_cds_exon_string = make_exon_string($secondary_transcript);
        $secondary_transcript_cds_seq = $secondary_transcript->translateable_seq();
        $secondary_transcript_pep_seq = $secondary_transcript->translate()->seq();
        $secondary_transcript_protein_coding = 1;
      } else {
        $secondary_transcript_protein_coding = 0;
      }

      # This is just a check to determine if one transcript has a translation and the other doesn't
      # In this case I don't think it makes sense to do a comparison, so they get skipped
      if($secondary_transcript_protein_coding && !$primary_protein_coding) {
        say $secondary_set_name." transcript has translation while ".$primary_set_name." transcript does not. Skipping comparison. ".
            $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID."), ".$secondary_set_name.
            " transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID.")" if ($debug);
        next;
      } elsif (!$secondary_transcript_protein_coding && $primary_protein_coding){
        say $primary_set_name." transcript has translation while ".$secondary_set_name." transcript does not. Skipping comparison. ".
            $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID."), ".$secondary_set_name.
            " transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID.")" if ($debug);
        next;
      }

      # Set the variables for the RefSeq model
      my $secondary_transcript_full_exon_string = make_exon_string($secondary_transcript,1);
      my $secondary_transcript_whole_seq = $secondary_transcript->seq->seq;

      # The next section is just a load of simple conditionals to test the transcripts. I've made it clear as opposed
      # to compact. These conditionals fill in @whole_structure_match, which is then analysed later to display the
      # summary message for each RefSeq transcript. The first two tests are common to all transcript pairs and the
      # remaining three are for protein coding transcript pairs only

      # Check if the entire transcript seqs match
      if($secondary_transcript_whole_seq eq $primary_whole_seq) {
        say "Match on complete transcript sequence. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".
            $secondary_transcript->dbID."), ".$primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        $whole_structure_match[0] = 1;
      } else {
        say "No match on complete transcript sequence. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".
             $secondary_transcript->dbID."), ".$primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        $whole_structure_match[0] = 0;
      }

      # Check if the exon coord string matches across the entire transcripts
      if($primary_full_exon_string eq $secondary_transcript_full_exon_string) {
        say "Match on complete exon structure. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".
            $secondary_transcript->dbID."), ".$primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        $whole_structure_match[1] = 1;
      } else {
        say "No match on complete exon structure. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".
             $secondary_transcript->dbID."), ".$primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        $whole_structure_match[1] = 0;
      }

      # Check if both are protein coding (at this point both transcripts have a translation or neither do)
      if($secondary_transcript_protein_coding && $primary_protein_coding) {

        # Check if the CDS sequence matches
        if($secondary_transcript_cds_seq eq $primary_cds_seq) {
          say "Match on CDS sequence. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
              $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
          $whole_structure_match[2] = 1;
        } else {
          say "No match on CDS sequence.".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
              $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
          $whole_structure_match[2] = 0;
        }

        # Check if the translateable exon coords match
        if($secondary_transcript_cds_exon_string eq $primary_cds_exon_string) {
          say "Match on CDS exon structure. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
              $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
           $whole_structure_match[3] = 1;
        } else {
          say "No match on CDS exon structure. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID.
              "), ".$primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
           $whole_structure_match[3] = 0;
        }

        # Check if the peptide sequences match
        if($secondary_transcript_pep_seq eq $primary_pep_seq) {
          say "Match on peptide seq. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
              $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
          $whole_structure_match[4] = 1;
        } else {
          say "No match on peptide seq. ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
              $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
          $whole_structure_match[4] = 0;
        }

      } # if($secondary_transcript_protein_coding && $primary_protein_coding)

      # Since both transcripts do not have a translation, neither have one, otherwise they would have been skipped
      else {
        say "No CDS/translateable exon/peptide comparison as both transcripts are non-coding. ".$secondary_set_name." transcript ".
            $secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".$primary_set_name." transcript ".
            $primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
      }

      # Analyse the array, the first two elements are the general tests and the last three are the tests on
      # transcript pairs that have translations. If the sum of the elements is five then both transcripts
      # were protein coding and matched on all tests. If the sum is two and the models aren't coding then
      # they've matched on all possible tests. If they are protein coding and the sum of the final three
      # elements is three then they match on the coding level. Otherwise they do not have a full match.
      my $whole_count = 0;
      my $protein_coding_count = 0;
      for(my $i=0; $i < scalar(@whole_structure_match); $i++) {
        $whole_count += $whole_structure_match[$i];
        if($i >= 2) {
          $protein_coding_count += $whole_structure_match[$i];
        }
      }


      if($whole_count == 5) {
        say "Match on all tests (coding). ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
            $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        push(@whole_structure_match_stable_ids, $secondary_transcript->stable_id);
      } elsif((!$secondary_transcript_protein_coding && !$primary_protein_coding) && $whole_count == 2) {
        say "Match on all tests (non-coding). ".$secondary_set_name." transcript ".$secondary_transcript->stable_id." (".$secondary_transcript->dbID."), ".
            $primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        push(@whole_structure_match_stable_ids, $secondary_transcript->stable_id);
      } elsif(($secondary_transcript_protein_coding && $primary_protein_coding) && $protein_coding_count == 3) {
        say "Match on all CDS tests only (protein coding). ".$secondary_set_name." transcript ".$secondary_transcript->stable_id.
            " (".$secondary_transcript->dbID."), ".$primary_set_name." transcript ".$primary_transcript->stable_id." (".$primary_transcript->dbID.")" if ($debug);
        $protein_coding_only_match_count++;
        push(@protein_coding_only_match_stable_ids, $secondary_transcript->stable_id);
      }

    } # End foreach my $secondary_transcript

    # The next set of messages are summary messages to allow grepping to see how many of the models had at
    # least one match on either the whole transcript or CDS level

    my $report_string = "No whole transcript or CDS match found for ".$primary_set_name." transcript: ".$primary_transcript->stable_id." (".
                       $primary_transcript->dbID."), ".$secondary_set_name." slice ".$transcript_slice->name;
    my @attributes;
    # These say statements can be grepped out of the output file to get a per transcript result for the matches
    if(@whole_structure_match_stable_ids) {
      $report_string = join(':', sort @whole_structure_match_stable_ids, 'whole_transcript');
      push(@attributes, Bio::EnsEMBL::Attribute->new(
        -code => $attribute_code,
        -name => $attrib_name,
        -description => $attrib_description,
        -value => $report_string,
      ));
      say "\n".$report_string if ($debug);
    }
    if(@protein_coding_only_match_stable_ids) {
      $report_string = join(':', sort @protein_coding_only_match_stable_ids, 'cds_only');
      push(@attributes, Bio::EnsEMBL::Attribute->new(
        -code => $attribute_code,
        -name => $attrib_name,
        -description => $attrib_description,
        -value => $report_string,
      ));
      say "\n".$report_string if ($debug);
    }
    if ($write) {
      $attribute_adaptor->store_on_Transcript($primary_transcript, \@attributes);
    }

    say "===================================================================\n\n" if ($debug);

  } # End foreach primary_transcript
} # End foreach slice
say "...finished comparing ".$primary_set_name." transcripts to overlapping ".$secondary_set_name." transcripts" if ($debug);


=head2 get_feature_slice_from_db

 Arg [1]    : Bio::EnsEMBL::Feature
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Return a the genomic region of Arg[1] on the forward strand
 Returntype : Bio::EnsEMBL::Slice or undef
 Exceptions : Throws if more than one slice is returned

=cut

sub get_feature_slice_from_db {
  my ( $feature, $db ) = @_;

  # This little helper routine returns a feature slice for a particular
  # region.  The slice will be associated with the given database.

  my $slice = $feature->feature_Slice();

  if ($debug) {
    say "Slice start/end: ".$slice->start()." ".$slice->end();
    say "Feature start/end: ".$feature->start()." ".$feature->end();

    say "CSN: ".$slice->coord_system_name();
    say "SRN: ".$slice->seq_region_name();
    say "SRS: ".$slice->start();
    say "SRE: ".$slice->end();
    say "CSV: ".$slice->coord_system()->version();
  }

  my $slices = $db->get_SliceAdaptor()->fetch_by_region_unique(
         $slice->coord_system_name(),
         $slice->seq_region_name(),
         $slice->start(),
         $slice->end(),
         1,
         $slice->coord_system()->version(),
         1);

  if(scalar(@$slices) == 0) {
    say "No slice was returned, the gene crosses a patch boundary so skip";
    return;
  }elsif ( scalar(@$slices) > 1 ) {
    # This will hopefully only happen if the Havana and Ensembl
    # databases contain different assemblies.
    foreach my $problem_slice (@$slices) {
      say "Problem slice: ".$problem_slice->name;
    }
     die( "!! Problem with projection for feature slice %s\n",
          $slice->name() );

  }

  return $slices->[0];
}


=head2 make_exon_string

 Arg [1]    : Bio::EnsEMBL::Transcript
 Arg [2]    : Boolean false if you only want the coding regions, true for all exons
 Description: Create a unique string representing the structure of the transcript
 Returntype : String
 Exceptions : None

=cut

sub make_exon_string {
  my ($transcript,$utr) = @_;

  unless($transcript) {
    say "Issue with undefined transcript, skipping";
    return 1;
  }
  my $exons;
  if ($utr) {
    $exons = $transcript->get_all_Exons();
  } else {
    $exons = $transcript->get_all_translateable_Exons();
  }
  my $exon_string;
  foreach my $exon (@{$translateable_exons}) {
    $exon_string .= $exon->start.':'.$exon->end.':';
  }

  return $exon_string;
}
