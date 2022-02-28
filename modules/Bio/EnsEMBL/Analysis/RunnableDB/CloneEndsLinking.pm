=head1 LICENSE

=======
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

Bio::EnsEMBL::Analysis::RunnableDB::CloneEndsLinking -

=head1 SYNOPSIS

my $clonemap =
  Bio::EnsEMBL::Analysis::RunnableDB::CloneEndsLinking->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
    -database   => $EST_GENOMIC,
  );

$clonemap->fetch_input();
$clonemap->run();
$clonemap->write_output(); #writes to DB

=head1 DESCRIPTION

This object links clone ends sequences by using the
exonerate alignment run (using ExonerateAlignFeature)

=head1 METHODS

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::CloneEndsLinking;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::CloneEndsLinking qw(CLONE_END_LINKING_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(empty_Object);
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Attribute;


use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    $self->read_and_check_config($CLONE_END_LINKING_CONFIG_BY_LOGIC);
    return $self;
}

sub fetch_input {
    my( $self) = @_;

    # Don't use the disconect_when_inactive flag as the job is quick enough
    my $db = $self->get_dbadaptor($self->CLONE_ALIGNED_DB);
    my $dnaalign_adaptor = $db->get_DnaAlignFeatureAdaptor;
    my @features;
    $self->output_db($self->get_dbadaptor($self->OUTDB));
    $self->get_clone_information;
    foreach my $trace_id (keys %{$self->links}) {
        my $raw_features = $dnaalign_adaptor->fetch_all_by_hit_name($trace_id, $self->CLONE_LOGIC_NAME);
        $self->output($raw_features) if ($self->STORE_DNAALIGNFEATURES);
        foreach my $feature (@{$raw_features}) {
            push(@features, $feature) if ($feature->hcoverage > 50);
        }
    }
    info('There is '.scalar(@features).' for '.$self->input_id);
    $self->input_is_void(1) unless (scalar(@features) or scalar(@{$self->output}));
    $self->clones(\@features);
}

sub run {
    my $self = shift @_;

    if (@{$self->clones} > 1) {
        my $links = $self->links;
        my ($trace_id) = each %$links;
        my %misc_features;

        foreach my $feature (@{$self->clones}) {
            push(@{$misc_features{$links->{$feature->hseqname}{direction}}}, {
                    insert_size => $links->{$feature->hseqname}{insert_size},
                    insert_stdev => $links->{$feature->hseqname}{insert_stdev},
                    feature => $feature
                    }
                );
        }
        $self->make_misc_feature(\%misc_features, $self->misc_set) if (keys %misc_features == 2);
    }
}

sub write_output {
    my $self = shift @_;

    my $db = $self->output_db;
    my $daf_adaptor = $self->output_db->get_DnaAlignFeatureAdaptor;
    my $mf_adaptor = $self->output_db->get_MiscFeatureAdaptor;
    my $analysis_adaptor = $self->output_db->get_AnalysisAdaptor;
    my $analysis = $analysis_adaptor->fetch_by_logic_name($self->CLONE_LOGIC_NAME);
    foreach my $feature (@{$self->output}) {
        if ($self->STORE_DNAALIGNFEATURES and $feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature')) {
            empty_Object($feature, 0, $analysis);
            $daf_adaptor->store($feature);
        }
        else {
            $mf_adaptor->store($feature);
        }
    }
}

sub make_misc_feature {
    my ($self, $misc_set_hash, $clone_lib) = @_;

    my @misc_features;
    my $clone_name = $self->input_id;
    my $library = $clone_lib->code;
    my ($well_name) = $clone_name =~ /$library-?(\w+)/;
    my $multiplier = $self->LONG_INSERT_SIZE_MULTIPLIER;
    my $clone_set = Bio::EnsEMBL::Attribute->new(
            -CODE => 'clone_name',
            -VALUE => $library,
            );
    my $synonym = Bio::EnsEMBL::Attribute->new(
            -CODE => 'synonym',
            -VALUE => $clone_name,
            );
    my $highlimit_status = Bio::EnsEMBL::Attribute->new(
            -CODE => 'state',
            -VALUE => 'LongInsert',
            );
    my $consistent_status = Bio::EnsEMBL::Attribute->new(
            -CODE => 'state',
            -VALUE => 'Consistent',
            );
    my %other_seq_region;
    my $is_problematic = 0;
    my $length_f_clones = scalar(@{$misc_set_hash->{F}});
    my $length_r_clones = scalar(@{$misc_set_hash->{R}});
    my $multiple_status;
    if ($length_f_clones > 1 or $length_r_clones > 1) {
        $multiple_status = Bio::EnsEMBL::Attribute->new(
                -CODE => 'state',
                -VALUE => 'MultipleHits',
                );
        foreach my $clone (@{$misc_set_hash->{F}}, @{$misc_set_hash->{R}}) {
                print STDERR 'DEBUG: ', $clone->{feature}->seq_region_name, "\n";
                push(@{$other_seq_region{$clone->{feature}->hseqname}}, $clone->{feature}->seq_region_name);
        }
        $is_problematic = 1;
    }
    for (my $f_index = 0; $f_index < $length_f_clones; $f_index++) {
        my $f_clone = $misc_set_hash->{F}[$f_index];
        my $low_limit = $f_clone->{insert_size}-$f_clone->{insert_stdev};
        my $high_limit = $f_clone->{insert_size}+$f_clone->{insert_stdev};
        for (my $r_index = 0; $r_index < $length_r_clones; $r_index++) {
            my $r_clone = $misc_set_hash->{R}[$r_index];
            my $strand = $f_clone->{feature}->hstrand;
            if ($f_clone->{feature}->seq_region_name eq $r_clone->{feature}->seq_region_name) {
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
                        warning("WRONGORIENTATION: NOT creating a MiscFeature for $clone_name: You have <==---------==> instead of ==>-------<==");
                        $is_problematic = 1;
                    }
                    else {
                        my $multiple = 0;
                        if ($align_insert_size < $low_limit or $align_insert_size > $high_limit) {
                            warning("OUTOFLIMIT $align_insert_size is out of the allow limits: $low_limit < allowed < $high_limit for $clone_name");
                        }
                        my $misc_feat = Bio::EnsEMBL::MiscFeature->new(
                                -START => $start,
                                -END => $end,
                                -STRAND => $strand,
                                -SLICE => $f_clone->{feature}->slice,
                                );
                        $misc_feat->add_MiscSet($clone_lib);
                        $misc_feat->add_Attribute($clone_set);
                        $misc_feat->add_Attribute($synonym);
                        if (!$multiple_status and $align_insert_size > $multiplier*$high_limit) {
                            $misc_feat->add_Attribute($highlimit_status);
                            $is_problematic = 1;
                        }
                        $misc_feat->add_Attribute(Bio::EnsEMBL::Attribute->new(
                                    -CODE => 'well_name',
                                    -VALUE => $well_name,
                                    ));
                        my @other_alignments;
                        foreach my $clone ($f_clone, $r_clone) {
                            warning("MULTIPLE: $clone_name has multiple good alignments");
                            if (exists $other_seq_region{$clone->{feature}->hseqname}) {
                                next if (@{$other_seq_region{$clone->{feature}->hseqname}} == 1);
                                push(@other_alignments, @{$other_seq_region{$clone->{feature}->hseqname}});
                            }
                            $multiple = 1;
                        }
                        if ($multiple) {
                            my $value = 'Same region';
                            if (@other_alignments) {
                                warning("MULTIPLEALN: $clone_name has multiple good alignments");
                                my @seq_region_name;
                                foreach my $name (@other_alignments) {
                                    push(@seq_region_name, $name) if ($name ne $misc_feat->slice->seq_region_name);
                                }
                                $value = 'Regions: '.join(',', @seq_region_name) if (@seq_region_name);
                            }
                            my $remarks = Bio::EnsEMBL::Attribute->new(
                                    -CODE  => 'remark',
                                    -NAME  => 'Remark',
                                    -VALUE => $value,
                                    );
                            $misc_feat->add_Attribute($remarks);
                            $is_problematic = 1;

                        }
                        push(@misc_features, $misc_feat);
                    }
                }
                else {
                    warning("WRONGSTRAND: $clone_name mates are on the same strand ".$f_clone->{feature}->strand.'='.$r_clone->{feature}->strand);
                    $is_problematic = 1;
                }
            }
            else {
                warning("WRONGSEQREGION: $clone_name mates are on different regions ".$f_clone->{feature}->seq_region_name.'='.$r_clone->{feature}->seq_region_name);
                $is_problematic = 1;
            }
        }
    }
    foreach my $misc_feat (@misc_features) {
        $misc_feat->add_Attribute($multiple_status) if ($multiple_status);
        $misc_feat->add_Attribute($consistent_status) unless ($is_problematic);

    }
    $self->output(\@misc_features);
}

sub get_clone_information {
    my $self = shift @_;

    my %links;
    my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($self->CLONE_LOGIC_NAME);
    my $sqlquery = 'SELECT * FROM clones WHERE clone_id = "'.$self->input_id.'" AND analysis_id = '.$analysis->dbID;
    my $sth = $self->db->dbc->prepare($sqlquery);
    $sth->execute();
    my $library;
    while (my $hashref = $sth->fetchrow_hashref()) {
        $links{$hashref->{trace_id}} = $hashref;
        $library = $hashref->{library};
    }
    $self->links(\%links);
    my $miscset_adaptor = $self->output_db->get_MiscSetAdaptor;
    my $misc_set = $miscset_adaptor->fetch_by_code($library);
    throw("The MiscSet $library does not exist in your database ".$self->OUTDB) unless ($misc_set);
    $self->misc_set($misc_set);
}

sub get_clone_state {
    my $self = shift;

    my %clones_state;
    open(my $fh, $self->CLONE_STATE) || die("Could not open".$self->CLONE_STATE);
    while (<$fh>) {
        # This is the header from CloneDB:
        #Gi  CloneName   Stdn    Chrom   Phase   CloneState  GCenter Accession   SeqLen  LibAbbr
        if(/^\d+\t(\S+)\t(\S+)\t\S*\t(\d)\t(\S+)\t\S*\t(\S+)\t\d+\t(\S+)/) {
            next unless ($2 eq 'Y');
            next unless ($4 eq 'fin');
            $clones_state{$1} = {phase => $3, accesion => $4, library => $5};
        }
    }
    close($fh) || die("Could not close ".$self->CLONE_STATE);
    return \%clones_state;
}

sub clones {
  my ($self, $arg) = @_;
  if($arg){
    $self->{_clones} = $arg;
  }
  return $self->{_clones}
}

sub links {
  my ($self, $arg) = @_;
  if($arg){
    $self->{_links} = $arg;
  }
  return $self->{_links}
}

sub misc_set {
  my ($self, $arg) = @_;
  if($arg){
    $self->{_misc_set} = $arg;
  }
  return $self->{_misc_set}
}

sub output_db {
  my ($self, $arg) = @_;
  if($arg){
    $self->{_output_db} = $arg;
  }
  return $self->{_output_db}
}

sub OUTDB {
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTDB} = $arg;
  }
  return $self->{OUTDB}
}

sub CLONE_LOGIC_NAME {
  my ($self, $arg) = @_;
  if($arg){
    $self->{CLONE_LOGIC_NAME} = $arg;
  }
  return $self->{CLONE_LOGIC_NAME}
}

sub CLONE_ALIGNED_DB {
  my ($self, $arg) = @_;
  if($arg){
    $self->{CLONE_ALIGNED_DB} = $arg;
  }
  return $self->{CLONE_ALIGNED_DB}
}

sub STORE_DNAALIGNFEATURES {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{STORE_DNAALIGNFEATURES} = $arg;
  }
  return $self->{STORE_DNAALIGNFEATURES}
}

sub LONG_INSERT_SIZE_MULTIPLIER {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{LONG_INSERT_SIZE_MULTIPLIER} = $arg;
  }
  return $self->{LONG_INSERT_SIZE_MULTIPLIER}
}

1;
