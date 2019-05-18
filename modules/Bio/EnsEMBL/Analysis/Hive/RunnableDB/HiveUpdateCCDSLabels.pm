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

HiveUpdateCCDSLabels.pm

=head1 DESCRIPTION

This module:

- Deletes the existing ccds attributes 'ccds_transcript', the CCDS dna_align_features and the CCDS transcript xrefs on the specified chromosome or top-level sequence or all top-level sequences.
- Inserts a ccds_transcript attribute, a CCDS transcript as supporting feature and a CCDS transcript xref into the output database for each transcript whose CDS and translation match a CCDS transcript in the CCDS database.

=head1 OPTIONS

-chromosome                 If defined, run the script only on the genes on the specified chromosome or top-level sequence.
-assembly_path              Assembly path. For example: 'GRCh38'.
-ccds_dbname                CCDS database name.
-ccds_host                  CCDS database host.
-ccds_user                  CCDS database user name.
-output_dbname              Output database name.
-output_host                Output database host.
-output_user                Output database user name.
-output_pass                Output database user pass.
-dna_dbname                 DNA database name.
-dna_host                   DNA database host.
-dna_user                   DNA database user name.
-dna_pass                   DNA database user pass.
-reports_dir                Directory where the file containing the list of missing ccds stable IDs will be written.
-output_filename            File name for the file containing the list of missing ccds stable IDs. This will get an extra extension corresponding to the chromosome name.

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInsertCCDSLabels -ccds_dbname CCDSDBNAME -ccds_host CCDSHOST -ccds_user CCDSUSER -output_dbname OUTPUTDBNAME -output_host OUTPUTHOST -output_user OUTPUTUSER -output_pass OUTPUTPASS -dna_dbname DNADBNAME -dna_host DNAHOST -dna_user DNAUSER -dna_pass DNAPASS -reports_dir REPORTSDIR -output_filename OUTPUTFILENAME

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUpdateCCDSLabels;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command send_email);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);
use Bio::EnsEMBL::DBEntry;

sub param_defaults {
    return {
      chromosome => 'toplevel',
      assembly_path => undef,
      ccds_dbname => undef,
      ccds_host => undef,
      ccds_port => undef,
      ccds_user => undef,
      output_dbname => undef,
      output_host => undef,
      output_port => undef,
      output_user => undef,
      output_pass => undef,
      dna_dbname => undef,
      dna_host => undef,
      dna_port => undef,
      dna_user => undef,
      dna_pass => undef,
      reports_dir => undef,
      output_filename => undef
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  
  $self->param_required('assembly_path');
  $self->param_required('ccds_dbname');
  $self->param_required('ccds_host');
  $self->param_required('ccds_port');
  $self->param_required('ccds_user');
  $self->param_required('output_dbname');
  $self->param_required('output_host');
  $self->param_required('output_port');
  $self->param_required('output_user');
  $self->param_required('output_pass');
  $self->param_required('dna_dbname');
  $self->param_required('dna_host');
  $self->param_required('dna_port');
  $self->param_required('dna_user');
  $self->param_required('dna_pass');
  $self->param_required('reports_dir');
  $self->param_required('output_filename');

  # insert ccds_transcript attributes and CCDS transcripts as supporting features
  my @missing_ccds = $self->insert_ccds_labels($self->param('chromosome'),
                                        $self->param('assembly_path'),
                                        $self->param('ccds_dbname'),
                                        $self->param('ccds_host'),
                                        $self->param('ccds_port'),
                                        $self->param('ccds_user'),
                                        $self->param('output_dbname'),
                                        $self->param('output_host'),
                                        $self->param('output_port'),
                                        $self->param('output_user'),
                                        $self->param('output_pass'),
                                        $self->param('dna_dbname'),
                                        $self->param('dna_host'),
                                        $self->param('dna_port'),
                                        $self->param('dna_user'),
                                        $self->param('dna_pass'));

  # write the missing ccds stable IDs into a file
  my $missing_ccds_file = $self->param('reports_dir').$self->param('output_filename').".".$self->param_required('chromosome');
  open(my $fh,'>',$missing_ccds_file) or die "Could not open file '$missing_ccds_file' $!";
  print $fh join("\n",@missing_ccds);
  close $fh;

  return 1;
}

sub insert_ccds_labels {
# Deletes the existing ccds attributes 'ccds_transcript', the CCDS dna_align_features and the CCDS transcript xrefs on the specified chromosome or top-level sequence or all top-level sequences.
# Inserts a ccds_transcript attribute, a CCDS transcript as supporting feature and a CCDS transcript xref into the output database for each transcript whose CDS and translation match a CCDS transcript in the CCDS database.
  my ($self,
      $chromosome,
      $assembly_path,
      $ccds_dbname,
      $ccds_host,
      $ccds_port,
      $ccds_user,
      $output_dbname,
      $output_host,
      $output_port,
      $output_user,
      $output_pass,
      $dna_dbname,
      $dna_host,
      $dna_port,
      $dna_user,
      $dna_pass) = @_;

  my $dna_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                     '-no_cache' => 1,
                                     '-host'     => $dna_host,
                                     '-port'     => $dna_port,
                                     '-user'     => $dna_user,
                                     '-pass'     => $dna_pass,
                                     '-dbname'   => $dna_dbname
  ) or
  die('Failed to connect to the dna database');

  my $output_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                     '-no_cache' => 1,
                                     '-host'     => $output_host,
                                     '-port'     => $output_port,
                                     '-user'     => $output_user,
                                     '-pass'     => $output_pass,
                                     '-dbname'   => $output_dbname
  ) or
  die('Failed to connect to the output database');
  $output_dba->dnadb($dna_dba);
  
  my $ccds_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                     '-no_cache' => 1,
                                     '-host'     => $ccds_host,
                                     '-port'     => $ccds_port,
                                     '-user'     => $ccds_user,
                                     '-pass'     => '',
                                     '-dbname'   => $ccds_dbname
  ) or
  die('Failed to connect to the CCDS database');
  $ccds_dba->dnadb($dna_dba);

  my $output_ta = $output_dba->get_TranscriptAdaptor();
  my $output_sa = $output_dba->get_SliceAdaptor();
  my $ccds_ta = $ccds_dba->get_TranscriptAdaptor();
  my $ccds_sa = $ccds_dba->get_SliceAdaptor();
  
  # fetch output db slices
  my @slices;
  my $all_slices = $output_sa->fetch_all('toplevel',$assembly_path,1,undef);
  if ($chromosome ne 'toplevel') { # if chr was defined as parameter, choose one slice
    foreach my $sl (@{$all_slices}) {
      if ($sl->seq_region_name() eq $chromosome) {
        $slices[0] = $sl;
      }
    }
  } else { # if chr was not defined, choose all toplevel slices
    @slices = @{$all_slices};
  }

  foreach my $slice (@slices) {
  	
  	print("--------------- Processing slice: ".$slice->seq_region_name()."\n");
  	
  	# delete existing ccds attributes and features before doing anything else
  	delete_ccds_labels_from_all_transcripts_on_slice($slice);
  	
    # delete existing ccds attributes and features before doing anything else
    #foreach my $transcript (@{$slice->get_all_Transcripts()}) {
    #  $transcript->load();
    #  delete_ccds_attrib_and_feature($output_dba,$transcript);
    #  $transcript = undef;
    #}
  }
  
  # fetch ccds db slices
  my @ccds_slices;
  my $ccds_all_slices = $ccds_sa->fetch_all('toplevel',$assembly_path,1,undef);
  if ($chromosome ne 'toplevel') { # if chr was defined as parameter, choose one slice
    foreach my $sl (@{$ccds_all_slices}) {
      if ($sl->seq_region_name() eq $chromosome) {
        $ccds_slices[0] = $sl;
      }
    }
  } else { # if chr was not defined, choose all toplevel slices
    @ccds_slices = @{$ccds_all_slices};
  }

  my @missing_ccds = ();

  foreach my $slice (@ccds_slices) {
    
    print("--------------- Processing CCDS slice: ".$slice->seq_region_name()."\n");

    foreach my $ccds_transcript (@{$slice->get_all_Transcripts()}) {
      $ccds_transcript->load();
      my $ccds_found = 0;
  	
      # get output db transcripts overlapping ccds transcript
  	  my $output_slice = get_feature_slice_from_db($ccds_transcript,$output_ta->db());
      my @output_transcripts = @{$output_ta->fetch_all_by_Slice($output_slice,1,undef,undef,undef,'protein_coding')};

      # get ccds translation
      my $ccds_translation = $ccds_transcript->translation();
      my $ccds_translation_seq;
      if ($ccds_translation) {
        $ccds_translation_seq = $ccds_translation->seq();
      } else {
        $self->throw($ccds_transcript->stable_id()." does not have a translation");
      }

      # check if any overlapping transcript matches ccds and add attribute and feature
      foreach my $output_transcript (@output_transcripts) {

        my $output_translation = $output_transcript->translation();
        my $output_translation_seq;
        if ($output_translation) {
          $output_translation_seq = $output_translation->seq();
        } else {
          $self->throw($output_transcript->stable_id()." does not have a translation");
        }

        my @translateable_exons = @{$output_transcript->get_all_translateable_Exons()};

        if (features_are_same(\@translateable_exons,$ccds_transcript->get_all_translateable_Exons())) {
          if ($output_translation_seq eq $ccds_translation_seq) {
            $self->add_ccds_transcript_attrib($output_dba,$output_transcript,$ccds_transcript->stable_id());
            $self->add_ccds_supporting_feature($output_dba,$ccds_transcript,$output_transcript);
            $self->add_ccds_transcript_xref($output_dba,$ccds_transcript,$output_transcript);
            $ccds_found = 1;
          }
        }
      } # end foreach output_transcript

      if (!($ccds_found)) {
        push(@missing_ccds,$ccds_transcript->stable_id());
      }
      $ccds_transcript = undef;
    } # end foreach ccds_transcript
  } # end foreach CCDS slice
  return @missing_ccds;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub add_ccds_transcript_attrib {
# Inserts a 'ccds_transcript' attribute into the transcript 'transcript' whose value is the CCDS stable id 'ccds_stable_id'
  my ($self, $db_adaptor,$transcript,$ccds_stable_id) = @_;
  
  my $attribute_adaptor = $db_adaptor->get_AttributeAdaptor();
  my $attrib_code = 'ccds_transcript';
  my ($attrib_type_id,$newcode,$name,$description) = $attribute_adaptor->fetch_by_code($attrib_code);
  
  if (!$attrib_type_id) {
    $self->throw("Unable to fetch attrib_type with code $attrib_code");
  }

  my $ccds_attribute = Bio::EnsEMBL::Attribute->new(
                                                     -NAME        => $name,
                                                     -CODE        => $attrib_code,
                                                     -VALUE       => $ccds_stable_id,
                                                     -DESCRIPTION => $description
                                                   );

  $attribute_adaptor->store_on_Transcript($transcript,[$ccds_attribute]);
}

sub add_ccds_supporting_feature {
# Inserts a ccds transcript supporting feature (including dna align feature) associated with the transcript 'transcript' whose hit name and the rest of its parameters are based on the CCDS 'ccds_transcript'
  my ($self, $dba,$ccds_transcript,$transcript) = @_;

  my @exon_features;
  my @features;

  my @exons = @{$transcript->get_all_translateable_Exons};

  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons;
  }

  my $exon_start = 1;
  my $exon_length = 0;

  foreach my $exon (@exons) {
    my $fp = new Bio::EnsEMBL::FeaturePair();
    $fp->start   ($exon->start);
    $fp->end     ($exon->end);
    $fp->strand  ($exon->strand);
    $fp->seqname ($exon->slice->seq_region_name);
    $fp->hseqname($ccds_transcript->stable_id());
    $fp->hstart  ($exon_start);
    $exon_length = $exon->end-$exon->start;
    $fp->hend    ($exon_start+$exon_length);

    $exon_start += $exon_length+1;
    $fp->hstrand(1);

    push(@exon_features,$fp);
  }
  my $tsf = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@exon_features);
  $tsf->seqname($transcript->slice()->seq_region_name());
  $tsf->slice($transcript->slice());
  $tsf->score(100);
  $tsf->analysis($transcript->analysis());

  push (@features, $tsf);

  my $tsf_adaptor = $dba->get_TranscriptSupportingFeatureAdaptor();
  $tsf_adaptor->store($transcript->dbID(),\@features);
}

sub add_ccds_transcript_xref {
# Inserts a ccds transcript xref associated with the transcript 'transcript' whose display_label and the rest of its parameters are based on the CCDS 'ccds_transcript'
  my ($self, $db_adaptor,$ccds_transcript,$transcript) = @_;
  
  my $dbe_adaptor = $db_adaptor->get_DBEntryAdaptor();
  my ($sid_without_version,$sid_version) = split(/\./,$ccds_transcript->stable_id());
  
  my $ccds_xref = new Bio::EnsEMBL::DBEntry(
                                             -adaptor => $dbe_adaptor,
                                             -primary_id => $sid_without_version,
                                             -version => $sid_version,
                                             -dbname  => 'CCDS',
                                             -display_id => $ccds_transcript->stable_id(),
                                             -description => '',
                                             -priority => 240,
                                             -db_display_name => 'CCDS',
                                             -info_type => 'DIRECT',
                                             -type => 'MISC'
                                           );
  $ccds_xref->status('XREF');
  $dbe_adaptor->store($ccds_xref,$transcript->dbID(),'Transcript');
}

sub delete_ccds_labels_from_all_transcripts_on_slice {
# Deletes the CCDS attributes 'ccds_transcript', the CCDS dna_align_features and the CCDS transcript xrefs from the slice 'slice' in the database associated with the adaptor 'db_adaptor'
  my ($slice) = shift();
  
  my $sa = $slice->adaptor();
  my $sdbc = $sa->dbc();
  
  my $seq_region_name = $slice->seq_region_name();
  my $coord_system_name = $slice->coord_system_name();
  my $coord_system_version = $slice->coord_system->version();
  
  $sdbc->do('DELETE ta FROM transcript t,transcript_attrib ta,seq_region sr,coord_system cs 
                                 WHERE t.transcript_id=ta.transcript_id AND
                                       t.seq_region_id=sr.seq_region_id AND
                                       sr.name="'.$seq_region_name.'" AND
                                       cs.name="'.$coord_system_name.'" AND
                                       cs.version="'.$coord_system_version.'" AND
                                       cs.coord_system_id=sr.coord_system_id AND
                                       ta.attrib_type_id IN (SELECT attrib_type_id FROM attrib_type WHERE code="ccds_transcript");');

  $sdbc->do('DELETE daf,tsf FROM dna_align_feature daf,transcript_supporting_feature tsf,transcript t,seq_region sr,coord_system cs 
                                 WHERE hit_name LIKE "CCDS%" AND
                                       tsf.feature_id=daf.dna_align_feature_id AND
                                       tsf.feature_type="dna_align_feature" AND
                                       tsf.transcript_id=t.transcript_id AND
                                       t.seq_region_id=sr.seq_region_id AND
                                       sr.name="'.$seq_region_name.'" AND
                                       cs.name="'.$coord_system_name.'" AND
                                       cs.version="'.$coord_system_version.'" AND
                                       cs.coord_system_id=sr.coord_system_id;');                         
  
  $sdbc->do('DELETE x,ox FROM xref x,object_xref ox,transcript t,seq_region sr,coord_system cs,external_db ed
                                 WHERE x.xref_id=ox.xref_id AND
                                       t.transcript_id=ox.ensembl_id AND
                                       ensembl_object_type="Transcript" AND
                                       t.seq_region_id=sr.seq_region_id AND
                                       ed.external_db_id=x.external_db_id AND
                                       ed.db_name="CCDS" AND
                                       sr.name="'.$seq_region_name.'" AND
                                       cs.name="'.$coord_system_name.'" AND
                                       cs.version="'.$coord_system_version.'" AND
                                       cs.coord_system_id=sr.coord_system_id;');                         
}

sub get_feature_slice_from_db {
  my ( $feature, $db ) = @_;

  # This little helper routine returns a feature slice for a particular
  # region.  The slice will be associated with the given database.

  my $slice = $feature->feature_Slice();

  my @slices = @{
    $db->get_SliceAdaptor()->fetch_by_region_unique(
         $slice->coord_system_name(), $slice->seq_region_name(),
         $slice->start(),             $slice->end(),
         1,            $slice->coord_system()->version(),
         1 ) };

  if ( scalar(@slices) != 1 ) {
    # This will hopefully only happen if the Primary and Secondary
    # databases contain different assemblies.
    die( "!! Problem with projection for feature slice %s\n",
         $slice->name() );
  }

  return $slices[0];
}

sub features_are_same {
  my ( $feature_set_a, $feature_set_b ) = @_;

  if ( scalar( @{$feature_set_a} ) == 0 ||
       ( scalar( @{$feature_set_a} ) != scalar( @{$feature_set_b} ) ) )
  {
    return 0;
  }

  for ( my $feature_index = 0;
        $feature_index < scalar( @{$feature_set_a} );
        ++$feature_index )
  {
    my $feature_a = $feature_set_a->[$feature_index];
    my $feature_b = $feature_set_b->[$feature_index];

    if (
      ( $feature_a->seq_region_start() != $feature_b->seq_region_start()
      ) ||
      ( $feature_a->seq_region_end() != $feature_b->seq_region_end() ) )
    {
      return 0;
    }
  }

  return 1;
} ## end sub features_are_same

1;
