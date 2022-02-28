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


package Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Pipeline::Tools::MM_Taxonomy;
use Bio::EnsEMBL::SimpleFeature;

use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished';

sub fetch_input {
	my ($self) = @_;
	my $slice =
	  $self->fetch_sequence( $self->input_id, $self->db,
		$ANALYSIS_REPEAT_MASKING, $SOFT_MASKING );
	$self->query($slice);

}

sub run {
	my ($self) = @_;

    my %params = %{ $self->parameters_hash() };
    my $max_coverage     = $params{max_coverage} || 10;
    my $hits_to_keep     = $params{hits_to_keep} || 3;
    my $max_hits_per_meta_cluster = $params{max_hits_per_meta_cluster} || 40;
    my $percentid_cutoff = $params{percentid_cutoff} || 0.0;
	my $orig_analysis_name = $params{ori_analysis};
	my $hit_db = $params{hit_db};
	my $mode = $params{mode};
	my $taxon_id = $params{taxon};
	my $no_filter = $params{no_filter};

	# Get the blast db version from the raw analysis and save it
	my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
	my $ana = $analysis_adaptor->fetch_by_logic_name($orig_analysis_name);
	my $stateinfocontainer = $self->db->get_StateInfoContainer;
	my $db_version = $stateinfocontainer->fetch_db_version($self->input_id,$ana);
	$self->db_version_searched($db_version);

    my $slice = $self->query();

	my $orig_features = $self->get_original_features($slice,$orig_analysis_name,$hit_db,$mode,$taxon_id, $no_filter);

	my ( $filtered_features, $saturated_zones) =
        $self->depth_filter($orig_features, $slice, $max_coverage,
        					$percentid_cutoff, $hits_to_keep, $no_filter, 
        					$max_hits_per_meta_cluster, $hit_db);

	$self->output($filtered_features, $saturated_zones);
}

sub get_original_features {
	my ($self,$slice,$analysis,$hit_db,$mode,$taxon_id, $no_filter) = @_;
	my $hit_db_features = [];
    my $prot_feat_a = $self->db->get_ProteinAlignFeatureAdaptor;
    my $dna_feat_a  = $self->db->get_DnaAlignFeatureAdaptor;

    my $orig_features = $dna_feat_a->fetch_all_by_Slice( $slice, $analysis );
	$orig_features = $prot_feat_a->fetch_all_by_Slice( $slice, $analysis ) if (!(@$orig_features));

	# The block of code below discard old sequence versions from the original set
	# of features (should add "mode => single" in analysis parameters)
	if($mode && $mode eq "single") {
		print STDERR "DepthFilter: $mode mode\n";
		my $single_hash;
		my $old_seq = [];
		foreach my $feature (@$orig_features) {
			my $hit_name	= $feature->hseqname;
			my ($acc,$ver)	= $hit_name =~ /(\w+-?\w+)\.(\d+)/;
			my $key			= $acc;
			if($single_hash->{$key}) {
				my $stored_feature = $single_hash->{$key}->[0];
				my ($sacc,$sver)	= $stored_feature->hseqname =~ /(\w+-?\w+)\.(\d+)/;
				if($ver > $sver) {
					$single_hash->{$key} = [$feature];
					#print STDERR "DepthFilter: drop ".$stored_feature->hseqname." (".$feature->hseqname.")\n";
					push @$old_seq, $stored_feature;
				} elsif($ver == $sver) {
					push @{$single_hash->{$key}}, $feature;
				} else {
					#print STDERR "DepthFilter: drop ".$feature->hseqname." (".$stored_feature->hseqname.")\n";
					push @$old_seq, $feature;
				}
			} else {
				$single_hash->{$key} = [$feature];
			}
		}
		print STDERR "DepthFilter: drop ".scalar(@$old_seq)." old sequences\n" if @$old_seq;
		$orig_features = [];
		map ( push(@$orig_features,@$_) , values %$single_hash);
	}

    if($hit_db || $taxon_id || $no_filter){
    	print STDERR "DepthFilter: hit db is $hit_db\n" if $hit_db;
    	my $plus_taxon_ids  = [];
    	my $minus_taxon_ids = [];
    	if ($taxon_id){
    		my @ti = split(/\|/,$taxon_id);
    		print STDERR "DepthFilter: taxonomy ids are ".join(" ",@ti)."\n";
    		for my $t (@ti) {
    			if($t < 0) {
    				$t = abs($t);
    				push @$minus_taxon_ids, @{$self->get_taxon_id_child($t)};
    				push @$minus_taxon_ids, $t;
    			} else {
    				push @$plus_taxon_ids, @{$self->get_taxon_id_child($t)};
    				push @$plus_taxon_ids, $t;
    			}
    		}
    	}


    	my $hit_hash = {map {$_->hseqname, undef} @$orig_features};
    	$self->get_hit_description($hit_hash);

	    my $error;

	    foreach my $feat (@$orig_features) {
	    	my $tag = 1;
	        if (my $desc = $self->get_hit_description($feat->hseqname)) {
	        	my $hit_taxon_id = $desc->taxon_id;

	            if($hit_db) {
	            	$tag = 0 if $desc->db_name ne $hit_db;
	            }
	            if(@$plus_taxon_ids) {
	            	$tag *= 0 unless grep(/^$hit_taxon_id$/,@$plus_taxon_ids);
	            }
	            if(@$minus_taxon_ids) {
	            	$tag *= 0 if grep(/^$hit_taxon_id$/,@$minus_taxon_ids);
	            }

				push(@$hit_db_features, $feat) if $tag;
	        } else {
	        	$error->{$feat->hseqname} = 1;
	        }
	    }

		throw("Missing hit_description entry for the following sequence(s):\n".join(',',keys %$error)) if keys %$error;

		print STDERR "DepthFilter: use ".scalar(@$hit_db_features)." features out of ".scalar(@$orig_features)."\n";
    } else {
    	return $orig_features;
    }

    return $hit_db_features
}


sub depth_filter {

	my ($self, $orig_features, $slice, $max_coverage, $percentid_cutoff, $hits_to_keep,
	$no_filter, $max_hits_per_meta_cluster, $hit_db) = @_;

	print STDERR "DepthFilter: MaxCoverage=$max_coverage\n";
	print STDERR "DepthFilter: PercentIdCutoff=$percentid_cutoff\n";
	print STDERR "DepthFilter: ".scalar(@$orig_features)." features before filtering\n";

    my %grouped_byname = ();
    my %is_kept_by_name = ();
	my %is_kept_by_acc = ();

    for my $af (@$orig_features) {
        my ($score, $percentid) = ($af->score(), $af->percent_id());
        if($percentid < $percentid_cutoff) {
            next;
        }

        my $node = $grouped_byname{$af->hseqname()} ||= {};
        if(%$node) { # nonempty
            $node->{max_score} = $score if $score>$node->{max_score};

            $node->{max_percentid} = $percentid if $percentid>$node->{max_percentid};
        } else {
            $node->{max_score} = $score;
            $node->{max_percentid} = $percentid;
        }
        if($self->get_hit_description($af->hseqname())) {
			$node->{taxon_id} = $self->get_hit_description($af->hseqname())->taxon_id;
        }

        push @{$node->{features}}, $af;
    }

	print STDERR "DepthFilter: ".scalar(keys %grouped_byname)." unique hitnames\n";

    my @bisorted =
        sort { ($b->{max_score} <=> $a->{max_score})
            || ($b->{max_percentid} <=> $a->{max_percentid}) }
        values %grouped_byname;

    my @coverage_map = ();
    my @filtered_features = ();

    for my $node (@bisorted) {
        my $keep_node = 0;
        my $hit_name;
        for my $af (sort {$a->start() <=> $b->start()} @{$node->{features}}) {
            for my $position ($af->start()..$af->end()) {
                my $depth = $coverage_map[$position] ||= 0;
                if($depth < $max_coverage) {
                    $keep_node = 1;

                }
                $coverage_map[$position]++;
            }
            $hit_name = $af->hseqname;
        }

        # add code to keep the hits with "no_filter" taxon_id
		$keep_node = 1 if $no_filter && grep(
								/^$node->{taxon_id}$/,
								@{$self->get_taxon_id_child($no_filter)},$no_filter
		);

        if($keep_node) {
            for my $af (@{$node->{features}}) {
                $af->analysis( $self->analysis );
                $af->dbID(0);
                $af->{adaptor} = undef;
                push @filtered_features, $af;
            }
            $is_kept_by_name{$hit_name} = 1;
            # trim off the splice variant number (if any) and version number
            $hit_name =~ s/(-\d+)?\.\d+$//;
            $is_kept_by_acc{$hit_name} = 1;
        }
    }

	# free up some memory ?
    %grouped_byname = ();

    # recover discarded Swissprot variants that should be shown
    if( $hit_db eq 'Swissprot' ) {
    	foreach my $f (@$orig_features) {
    		my $hit_name = $f->hseqname();
    		next unless !$is_kept_by_name{$hit_name};
    		$hit_name =~ s/(-\d+)?\.\d+$//;
			if($is_kept_by_acc{$hit_name}){
				print STDERR "DepthFilter: recover SW hit ".$f->hseqname."\n";
				$f->analysis( $self->analysis );
                $f->dbID(0);
                $f->{adaptor} = undef;
                push @filtered_features, $f;
			}
    	}
    }

	print STDERR "DepthFilter: ".scalar(@filtered_features)." features after filtering\n";

    my @saturated_zones = ();
    my $zone_start = undef;
    my $zone_score = 0;
    my $slice_length = $slice->length();
    for(my $i=1; $i<=$slice_length; $i++) {
        my $n = $coverage_map[$i] || 0;
        if ($zone_start) {
            $zone_score += $n;
            if ($n < $max_coverage) {
                my $new_zone = Bio::EnsEMBL::SimpleFeature->new(
                    -start  => $zone_start,
                    -end    => $i - 1,
                    -strand => 0,       # we mix both strands here
                    -score  => $zone_score,
                    -display_label => sprintf("avg depth = %.2f", $zone_score/($i-$zone_start)),
                    -analysis => $self->analysis(),
                );
                push(@saturated_zones, $new_zone);
                $zone_start = undef;
                $zone_score = 0;
            }
        }
        elsif ($n >= $max_coverage) {
            $zone_start = $i;
            $zone_score = $n;
        }
    }
    if ($zone_start) { # Are saturated up to end of contig:
        my $new_zone = Bio::EnsEMBL::SimpleFeature->new(
            -start  => $zone_start,
            -end    => $slice_length,
            -strand => 0,       # we mix both strands here
            -score  => $zone_score,
            -display_label => sprintf("avg depth = %.2f", $zone_score/($slice_length-$zone_start+1)),
            -analysis => $self->analysis(),
        );
        push(@saturated_zones, $new_zone);
    }

	print STDERR "DepthFilter: ".scalar(@saturated_zones)." saturated zones found\n";

    return (\@filtered_features, \@saturated_zones);
}

=head2 db_version_searched

    Title   :  db_version_searched
               [ distinguished from Runnable::*::get_db_version() ]
    Useage  :  $self->db_version_searched('version string')
               $obj->db_version_searched()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished::Blast
               This is just a holding place for the string in this
               module
    Returns :  String or undef
    Args    :  String
    Caller  :  $self::run()
               Job::run_module()

=cut

sub db_version_searched {

	my ( $self, $arg ) = @_;
	$self->{'_db_version_searched'} = $arg if $arg;
	return $self->{'_db_version_searched'};

}

sub get_hit_description {
	my ($self,$hit) = @_;
	if(ref($hit) eq 'HASH') {
		$self->{_hit_desc} = $hit;
		my $hit_desc_a = $self->db->get_HitDescriptionAdaptor;
		$hit_desc_a->fetch_HitDescriptions_into_hash($hit);

		return 1;
	}

	return $self->{_hit_desc}->{$hit};
}

sub get_taxon_id_child {
	my ( $self, $taxon_id ) = @_;
	$self->{taxonomy} ||= Bio::EnsEMBL::Pipeline::Tools::MM_Taxonomy->new();
	if(!$self->{_taxon}->{$taxon_id}) {
		$self->{_taxon}->{$taxon_id} = $self->{taxonomy}->get_all_children_id($taxon_id);
	}

	return $self->{_taxon}->{$taxon_id};
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter

=head2 AUTHOR
James Gilbert B<email> jgrg@sanger.ac.uk    - original implementation

=head2 AUTHOR
Mustapha Larbaoui B<email> ml6@sanger.ac.uk - new pipeline

=head2 AUTHOR
Leo Gordon B<email> lg4@sanger.ac.uk        - porting

