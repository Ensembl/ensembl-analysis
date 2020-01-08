=head1 LICENSE

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

Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds -

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds->new(
     -query_seqs     => \@q_seqs,
     -query_type     => 'dna',
     -target_seqs    => \@t_seqs,
     -options        => $options,
    );

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;

=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program,
to align clone sequences with genomic sequences.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEndsMapping;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Attribute;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature);

#
# Implementation of method in abstract superclass
#
sub parse_results {
  my ( $self, $fh ) = @_;

  my @features;
  my %misc_set;

  while (<$fh>) {

    print;
    next unless /^RESULT:/;

    chomp;

    my (
      $tag, $q_id, $q_start, $q_end, $q_strand,
      $t_id, $t_start, $t_end, $t_strand, $score,
      $perc_id, $q_length, $t_length, $gene_orientation,
      @vulgar_blocks
    ) = split;

    my $cigar_string='';
    while (@vulgar_blocks) {
      throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @vulgar_blocks . "] items left to process.")
      unless scalar @vulgar_blocks >= 3;

      my $match_type          = shift @vulgar_blocks;
      my $query_match_length  = shift @vulgar_blocks;
      my $target_match_length = shift @vulgar_blocks;


      if ($match_type eq "G"){
	if ($query_match_length == 0){
	    $match_type="D";
            $query_match_length = $target_match_length;
        }elsif ($target_match_length == 0){
            $match_type="I";
	}

      }

      $cigar_string .= $query_match_length.$match_type;

    }
      my ($clone_name, $query_id, $clone_direction, $clone_lib, $insert_size, $insert_dev, $trace_name) = split(':', $q_id);
    my $feature =
      $self->make_feature(
        $query_id, $q_length, $q_start, $q_end, $q_strand,
        $t_id, $t_length, $t_start, $t_end, $t_strand, $score,
        $perc_id, $cigar_string, 0, $t_strand
      );

    if($feature){
      push @features, $feature;
      $misc_set{$t_id}{$clone_name}->{lib} = $self->get_clone_set($clone_lib);
      push(@{$misc_set{$t_id}{$clone_name}{$clone_direction}}, { trace_name => $trace_name, insert_size => $insert_size, insert_dev => $insert_dev, feature => $feature });
    }else{
      warn "Clone end feature from probe :$q_id doesnt match well enough\n";
    }
  }
    my $misc_feat = $self->make_misc_feature(\%misc_set);
    foreach my $mf (@$misc_feat) {
        print STDERR $mf->seq_region_start, '-', $mf->seq_region_end, '^', $mf->seq_region_strand, "\n";
    }
    push(@features, @$misc_feat);

  return \@features;
}

sub make_misc_feature {
    my ($self, $misc_set_hash) = @_;

    my @misc_features;
    foreach my $key (keys %$misc_set_hash) {
        my $slice = $self->get_seq_region_slice($key);
        foreach my $clone_name (keys %{$misc_set_hash->{$key}}) {
            print STDERR $clone_name, "\n";
            my $clone_set = Bio::EnsEMBL::Attribute->new(
                -CODE => 'clone_name',
                -NAME => $clone_name,
                -DESCRIPTION => "$clone_name clone set",
                -VALUE => $clone_name,
            );
            my $clone_lib = $misc_set_hash->{$key}->{$clone_name}->{lib};
            for (my $f_index = 0; defined $misc_set_hash->{$key}{$clone_name}{F}[$f_index]; $f_index++) {
                my $low_limit = $misc_set_hash->{$key}{$clone_name}{F}[$f_index]->{insert_size}-$misc_set_hash->{$key}{$clone_name}{F}[$f_index]->{insert_dev};
                my $high_limit = $misc_set_hash->{$key}{$clone_name}{F}[$f_index]->{insert_size}+$misc_set_hash->{$key}{$clone_name}{F}[$f_index]->{insert_dev};
                for (my $r_index = 0; defined $misc_set_hash->{$key}{$clone_name}{R}[$r_index]; $r_index++) {
                    my $f_clone = $misc_set_hash->{$key}{$clone_name}{F}[$f_index];
                    my $r_clone = $misc_set_hash->{$key}{$clone_name}{R}[$r_index];
                    my $strand = $f_clone->{feature}->hstrand;
                    if ($f_clone->{feature}->hstrand != $r_clone->{feature}->hstrand) {
                        my $align_insert_size = 0;
                        my $start;
                        my $end;
                        if ($strand == 1) {
                            $start = $f_clone->{feature}->start-$f_clone->{feature}->hstart+1;
                            $end = $r_clone->{feature}->end+$r_clone->{feature}->hstart-1;
                            $align_insert_size = $r_clone->{feature}->start-$f_clone->{feature}->end+1;
                        }
                        else {
                            $start = $r_clone->{feature}->start-$r_clone->{feature}->hstart+1;
                            $end = $f_clone->{feature}->end+$f_clone->{feature}->hstart-1;
                            $align_insert_size = $f_clone->{feature}->start-$r_clone->{feature}->end+1;
                        }
                        if ($align_insert_size < 0) {
                            warning("NOT creating a MiscFeature for $clone_name: You have <==---------==> instead of ==>-------<==\n");
                        }
                        else {
                            if ($align_insert_size < $low_limit or $align_insert_size > $high_limit) {
                                warning("$align_insert_size is out of the allow limits: $low_limit < allowed < $high_limit for $clone_name\n");
                            }
                            my $misc_feat = Bio::EnsEMBL::MiscFeature->new(
                                -START => $start,
                                -END => $end,
                                -STRAND => $strand,
                                -SLICE => $slice,
                            );
                            $misc_feat->add_MiscSet($clone_lib);
                            $misc_feat->add_Attribute($clone_set);
                            foreach my $fr ('F', 'R') {
                                $misc_feat->add_Attribute(Bio::EnsEMBL::Attribute->new(
                                    -CODE => 'name',
                                    -VALUE => "$clone_name:$fr:".$misc_set_hash->{$key}{$clone_name}{F}[$f_index]->{insert_size},
                                ));
                            }
                            push(@misc_features, $misc_feat);
                        }
                    }
                    else {
                        warning("$clone_name mates are on the same strand ".$f_clone->{feature}->strand.'='.$r_clone->{feature}->strand);
                    }
                }
            }
        }
    }
    return \@misc_features;
}

sub get_seq_region_slice {
    my ($self, $slice_name) = @_;

    if (!exists $self->{_seq_region_slice}->{$slice_name}) {
        # We could ask for the slice but we don't want DB connection in a Runnable
        $slice_name =~ /^\w+:[^:]+:([^:]+):(\d+):(\d+)/;
        $self->{_seq_region_slice}->{$slice_name} = Bio::EnsEMBL::Slice->new(
            -START => $2,
            -END => $3,
            -SEQ_REGION_NAME => $1,
        );

    }
    return $self->{_seq_region_slice}->{$slice_name};
}

sub get_clone_set {
    my ($self, $clone_lib) = @_;

    if (!exists $self->{_clone_set}->{$clone_lib}) {
        $self->{_clone_set}->{$clone_lib} = Bio::EnsEMBL::MiscSet->new(
            -CODE => 'clone',
            -NAME => $clone_lib,
            -DESCRIPTION => "$clone_lib clone library",
            );
    }
    return $self->{_clone_set}->{$clone_lib};
}

1;
