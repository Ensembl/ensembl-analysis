=head1 LICENSE

 Copyright [2019] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::WriteCDSCoordsToFile

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::WriteCDSCoordsToFile->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module dumps cds coord data from Ensembl dbs to file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::WriteCDSCoordsToFile;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name is_canonical_splice);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _branch_to_flow_to => 1,
    use_generic_output_type => 0,
    generic_output_type => 'cds',
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for cds coord dumping
 Returntype : None

=cut

sub fetch_input {
  my ($self) = @_;

  say "Fetching input";
  # For the combine files option we don't actually need to do anything in fetch input
  unless($self->param_required('dump_type') eq 'combine_files') {
    my $target_dba = $self->hrdb_get_dba($self->param('target_db'));

    my $dna_dba;
    if($self->param('use_genome_flatfile')) {
      unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
        $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
      }
      setup_fasta(
                   -FASTA => $self->param_required('genome_file'),
                 );
    } else {
      $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
      $target_dba->dnadb($dna_dba);
    }

    my $gene_adaptor = $target_dba->get_GeneAdaptor();
    my $gene_ids = $self->input_id();
    my $genes = [];
    foreach my $gene_id (@$gene_ids) {
      push(@$genes,$gene_adaptor->fetch_by_dbID($gene_id));
    }
    $self->param('input_genes',$genes);
    $self->hrdb_set_con($target_dba,'target_db');
  } # End unless
}


sub run {
  my ($self) = @_;

  if($self->param_required('dump_type') eq 'gene') {
    say "Dumping cds  to file from genes";
    my $genes = $self->param('input_genes');
    $self->dump_genes($genes);
  } elsif($self->param('dump_type') eq 'combine_files') {
    $self->combine_files();
  } else {
    $self->throw("Dump type not recognised. Dump type: ".$self->param('dump_type'));
  }

}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  unless($self->param('dump_type') eq 'combine_files') {
    my $output_file = $self->param('_output_file');
    my $accu_entry = "cds_file";
    $self->dataflow_output_id({$accu_entry => $output_file}, $self->param('_branch_to_flow_to'));
  }
}


sub dump_genes {
  my ($self,$genes) = @_;

  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my $output_dir  = $self->param_required('cds_dir');
  my $output_file = $self->create_filename(undef,'cds',$output_dir,1);
  my $temp_name = "$output_file";
  say "Creating output file: ".$output_file;
  my $input_type  = $self->param_required('input_type');

  my $cds_string_hash = {};
  open(OUT,">".$output_file);
  foreach my $gene (@$genes) {
    my $strand = $gene->strand;
    say "Gene: ".$gene->dbID.":".$gene->strand;
    my $transcripts = $gene->get_all_Transcripts;
    my $seq_region_name = $gene->seq_region_name;
    foreach my $transcript (@$transcripts) {
      # Note I decided the cds length is crucial for different isoforms, but it comes at the
      # cost of having to call translation. Otherwise $transcript->coding_region_start would
      # have worked here and I assume should be more lightweight (but then who knows)
      unless($transcript->translation) {
        next;
      }

      my $cds_start = $transcript->coding_region_start;
      my $cds_end = $transcript->coding_region_end;
      my $cds_length = length($transcript->translateable_seq);
      my $cds_string = $cds_start.":".$cds_end.":".$strand.":".$cds_length;
      if($cds_string_hash->{$seq_region_name}->{$cds_string}) {
        $cds_string_hash->{$seq_region_name}->{$cds_string} += 1;
      } else {
        $cds_string_hash->{$seq_region_name}->{$cds_string} = 1;
      }
    } # end foreach my $transcript (@$transcripts)
  } # foreach my $gene (@$genes)

  foreach my $seq_region (keys(%$cds_string_hash)) {
    my $out_line = $seq_region.":".$input_type.",";
    my $cds_coords = $cds_string_hash->{$seq_region};
    foreach my $cds (keys(%$cds_coords)) {
      my $count = $cds_coords->{$cds};
      $out_line .= $cds.":".$count.",";
    }
    say OUT $out_line;
  }
  close OUT;

  $self->param('_output_file',$temp_name);
}


sub combine_files {
  my ($self) = @_;

  my $full_cds_hash = {};
  my $initial_cds_files = $self->param_required('cds_file');
  foreach my $cds_file (@$initial_cds_files) {
    unless(open(IN,$cds_file)) {
      $self->throw("Could not open cds file: ".$cds_file);
    }
    while(<IN>) {
      chomp $_;
      my @line = split(',',$_);
      my $header = $line[0];
      my ($seq_region,$type) = split(':',$header);
      if($self->param('use_generic_output_type')) {
        $type = $self->param('generic_output_type');
      }

      foreach(my $i=1; $i<scalar(@line); $i++) {
        my $cds_string = $line[$i];
        my ($start,$end,$strand,$length,$count) = split(':',$cds_string);
        my $cds = $start.':'.$end.':'.$strand.":".$length;
        if($full_cds_hash->{$seq_region}->{$type}->{$cds}) {
          $full_cds_hash->{$seq_region}->{$type}->{$cds} += $count;
        } else {
          $full_cds_hash->{$seq_region}->{$type}->{$cds} = $count;
	}
      }
    }
    close IN;
  }

  open(OUT,">".$self->param_required('combined_cds_file'));
  foreach my $seq_region (keys(%$full_cds_hash)) {
    foreach my $type (keys(%{$full_cds_hash->{$seq_region}})) {
      my $line = $seq_region.":".$type.",";
      foreach my $cds (keys(%{$full_cds_hash->{$seq_region}->{$type}})) {
        my $count = $full_cds_hash->{$seq_region}->{$type}->{$cds};
        $line .= $cds.":".$count.",";
      }
      say OUT $line;
    }
  }
  close OUT;
}


sub create_filename {
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


1;
