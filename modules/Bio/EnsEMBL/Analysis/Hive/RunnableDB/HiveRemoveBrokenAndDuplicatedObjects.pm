=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveDuplicatedObjects

=head1 SYNOPSIS


=head1 DESCRIPTION

Remove any models which is broken, like missing the transcript and exons or duplicated models
These are created when a job failed in the WRITE state because of memory or time limit
The 'feature_type' parameter command two helper functions, process_<feature_type> and is_ok_<feature_type>,
which need to be implemented for any 'feature_type'

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(hashkey_Object);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default options:
								feature_type => 'gene',
								check_support => 1,
								coord_system => 'toplevel',
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
	my ($self) = @_;

	return {
		%{$self->SUPER::param_defaults},
		feature_type => 'gene',
		check_support => 1,
		coord_system => 'toplevel',
	}
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch the input objects for flagging
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
	my ($self) = @_;

	my $db = $self->get_database_by_name('target_db');
	my $slice_adaptor = $db->get_SliceAdaptor;
	my $slices;
	if ($self->param_is_defined('iid')) {
		$slices = [$slice_adaptor->fetch_by_name($self->param('iid'))];
	}
	else {
		$slices = $slice_adaptor->fetch_all($self->param('coord_system'));
	}
	my @objects;
	if ($self->param('feature_type') eq 'gene') {
		foreach my $slice (@$slices) {
			push(@objects, sort {$a->start <=> $b->start || $b->end <=> $a->end} @{$slice->get_all_Genes(undef, undef, 1)});
		}
	}
	elsif ($self->param('feature_type') eq 'prediction') {
		foreach my $slice (@$slices) {
			push(@objects, sort {$a->start <=> $b->start || $b->end <=> $a->end} @{$slice->get_all_PredictionTranscripts(undef, 1)});
		}
	}
	else {
		$self->throw('Unknown feature_type "'.$self->param('feature_type').'", you may want to create a method process_'.$self->param('feature_type').' in '.__PACKAGE__);
	}
	$self->param('objects', \@objects);
}


=head2 run

 Arg [1]    : None
 Description: Flag the object which are duplicated or broken
 Returntype : None
 Exceptions : None

=cut

sub run {
	my ($self) = @_;

	my $objects = $self->param('objects');
	my %objects_to_delete;
	my $method = 'process_'.$self->param('feature_type');
  my $has_features = 'is_ok_'.$self->param('feature_type');
	for (my $i = 0; $i < @$objects; $i++) {
    if ($self->$has_features($objects->[$i])) {
      for (my $j = $i+1; $j < @$objects; $j++) {
        if ($self->$has_features($objects->[$j])) {
          if ($objects->[$i]->slice == $objects->[$j]->slice) {
            if ($objects->[$i]->end < $objects->[$j]->start) {
              last;
            }
            elsif ($objects->[$i]->start == $objects->[$j]->end) {
              my $result = $self->$method($objects->[$i], $objects->[$j]);
              if ($result) {
                if ($result == $objects->[$i]) {
                  $objects_to_delete{$objects->[$i]->dbID} = $objects->[$i];
                  last;
                }
                if ($result == $objects->[$j]) {
                  $objects_to_delete{$objects->[$j]->dbID} = $objects->[$j];
                  next;
                }
              }
              else {
                $self->warning($self->param('feature_type').' objects '.$objects->[$i]->display_id.' and '.$objects->[$j]->display_id.' have the same coordinates but does not seem to be the same objects');
              }
            }
          }
          else {
            last;
          }
        }
        else {
          $objects_to_delete{$objects->[$j]->dbID} = $objects->[$j];
        }
      }
    }
    else {
      $objects_to_delete{$objects->[$i]->dbID} = $objects->[$i];
    }
	}
	$self->output([values %objects_to_delete]);
}


=head2 write_output

 Arg [1]    : None
 Description: Delete the flagged objects from the database
 Returntype : None
 Exceptions : None

=cut

sub write_output {
	my ($self) = @_;

	foreach my $object (@{$self->output}) {
    $object->adaptor->remove($object) if ($object->dbID);
	}
}


=head2 process_gene

 Arg [1]    : Bio::EnsEMBL::Gene
 Arg [2]    : Bio::EnsEMBL::Gene
 Description: Try to find if one of the gene is broken or if the genes are duplicates
              It will delete the broken gene and Arg[2] if they are duplicates
 Returntype : Bio::EnsEMBL::Gene
 Exceptions : None

=cut

sub process_gene {
	my ($self, $gene_i, $gene_j) = @_;

	my $i_transcripts = $gene_i->get_all_Transcripts;
	my $j_transcripts = $gene_j->get_all_Transcripts;
	if (scalar(@$i_transcripts) and scalar(@$j_transcripts)) {
		for (my $m = 0; $m < @$i_transcripts; $m++) {
			my $hashkey_m = hashkey_Object($i_transcripts->[$m]->get_all_Exons);
			for (my $n = 0; $n < @$j_transcripts; $n++) {
				my $hashkey_n = hashkey_Object($j_transcripts->[$n]->get_all_Exons);
				if ($hashkey_n eq $hashkey_m) {
					if ($self->param('check_support')) {
						my $support_m = join(':', grep {$_->hit_name} sort {$a->hit_name cmp $b->hit_name} @{$i_transcripts->[$m]->get_all_supporting_evidences});
						my $support_n = join(':', grep {$_->hit_name} sort {$a->hit_name cmp $b->hit_name} @{$j_transcripts->[$n]->get_all_supporting_evidences});
						if ($support_n eq $support_m) {
							return $gene_j;
						}
					}
					else {
						return $gene_j;
					}
				}
			}
		}
	}
	else {
		return $gene_i unless (scalar(@$i_transcripts));
		return $gene_j unless (scalar(@$j_transcripts));
	}
	return; #Not sure it's needed but it should do the work
}


=head2 is_ok_gene

 Arg [1]    : Bio::EnsEMBL::Gene
 Description: Define if Arg[1] is complete by checing if it has transcripts,
              if the transcripts have exon where the start and end exons match
              with the transcript start and end. It also checks if the possible
              protein has a start and end exon.
 Returntype : Boolean
 Exceptions : None

=cut

sub is_ok_gene {
  my ($self, $gene) = @_;

  my $transcripts = $gene->get_all_Transcripts;
  if (scalar(@$transcripts)) {
    foreach my $t (@$transcripts) {
      return 0 unless (
          ($t->strand == 1 and $t->start_Exon  and $t->start_Exon->start == $t->start
          and $t->end_Exon and $t->end_Exon->end == $t->end)
        or ($t->strand == -1 and $t->start_Exon and $t->start_Exon->end == $t->end
          and $t->end_Exon and $t->end_Exon->start == $t->start));
      my $p = $t->translation;
      if ($p) {
        return 0 if (!$p->start_Exon or !$p->end_Exon);
      }
    }
  }
  else {
    return 0;
  }
  return 1;
}


=head2 process_prediction

 Arg [1]    : Bio::EnsEMBL::PredictionTranscript
 Arg [2]    : Bio::EnsEMBL::PredictionTranscript
 Description: Try to find if one of the transcript is broken or if the transcripts are duplicates
              It will delete the broken transcript and Arg[2] if they are duplicates
 Returntype : Bio::EnsEMBL::PredictionTranscript
 Exceptions : None

=cut

sub process_prediction {
	my ($self, $gene_i, $gene_j) = @_;

	my $exon_i = $gene_i->get_all_Exons;
	my $exon_j = $gene_j->get_all_Exons;
	if (scalar(@$exon_i) and scalar(@$exon_j)) {
		if (hashkey_Object($exon_i) eq hashkey_Object($exon_j)) {
			return $gene_j;
		}
	}
	else {
		return $gene_i unless (scalar(@$exon_i));
		return $gene_j unless (scalar(@$exon_j));
	}
}


=head2 is_ok_prediction

 Arg [1]    : Bio::EnsEMBL::PredictionTranscript
 Description: Define is Arg[1] is complete by checking if it has exons
              and if the start and end exon have matching start and end with
              the objectrg[1]
 Returntype : Boolean
 Exceptions : None

=cut

sub is_ok_prediction {
  my ($self, $transcript) = @_;

  if (scalar(@{$transcript->get_all_Exons})) {
    my $exon = $transcript->start_Exon;
    return 0 if (!$exon or $exon->start != $transcript->start);
    $exon = $transcript->end_Exon;
    return 0 if (!$exon or $exon->end != $transcript->end);
  }
  else {
    return 0;
  }
  return 1;
}
1;
