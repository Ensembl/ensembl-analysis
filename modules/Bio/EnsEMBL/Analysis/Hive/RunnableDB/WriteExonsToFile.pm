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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::WriteExonsToFile;

use warnings;
use strict;
use feature 'say';
use File::Spec::Functions;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name is_canonical_splice);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _branch_to_flow_to => 1,
    use_generic_output_type => 0,
    generic_output_type => 'exon',
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
#  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));
#  $self->hrdb_set_con($target_dba,'target_db');
}


sub run {
  my ($self) = @_;

  my $input_dbs = $self->param_required('input_dbs');

  my $seq_region_names = {};
  my $exon_string_hash = {};
  foreach my $input_id (@$input_dbs) {
    my $target_dba = $self->hrdb_get_dba($input_id);
    my $seq_region_sql = "SELECT seq_region_id,name FROM seq_region";
    my $sth = $target_dba->dbc->prepare($seq_region_sql);
    $sth->execute();
    while(my @seq_region_array = $sth->fetchrow_array) {
      $seq_region_names->{$seq_region_array[0]} = $seq_region_array[1];
    }

    my $exon_sql = "SELECT seq_region_start,seq_region_end,seq_region_strand,seq_region_id FROM exon";
    $sth = $target_dba->dbc->prepare($exon_sql);
    $sth->execute();
    while(my @exon_array = $sth->fetchrow_array) {
      my ($seq_region_start,$seq_region_end,$seq_region_strand,$seq_region_id) = @exon_array;
      my $seq_region_name = $seq_region_names->{$seq_region_id};
      my $exon_string = $seq_region_start.":".$seq_region_end.":".$seq_region_strand;

      if($exon_string_hash->{$seq_region_name}->{$exon_string}) {
        $exon_string_hash->{$seq_region_name}->{$exon_string} += 1;
      } else {
        $exon_string_hash->{$seq_region_name}->{$exon_string} = 1;
      }
    }
  }

  my $output_dir  = $self->param_required('exons_dir');
  my $output_file = catfile($output_dir,'final.exon');

  say "Creating output file: ".$output_file;
  my $input_type  = $self->param_required('input_type');
  open(OUT,">".$output_file);
  foreach my $seq_region (keys(%$exon_string_hash)) {
    my $out_line = $seq_region.":".$input_type.",";
    my $exon_strings = $exon_string_hash->{$seq_region};
    foreach my $exon_string (keys(%$exon_strings)) {
      my $count = $exon_strings->{$exon_string};
      $out_line .= $exon_string.":".$count.",";
    }
    say OUT $out_line;
  }
  close OUT;
}


sub write_output {
  my ($self) = @_;
}

1;
