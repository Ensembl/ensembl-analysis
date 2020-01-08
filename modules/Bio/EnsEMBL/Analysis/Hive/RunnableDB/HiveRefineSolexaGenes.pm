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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::RefineSolexaGenes

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::RefineSolexaGenes->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $refine_genes->fetch_input();
  $refine_genes->run();
  $refine_genes->write_output(); #writes to DB



=head1 DESCRIPTION

This module takes intron spanning dna_align_features and combines them with
rough transcript models to build refined genes with CDS. The module produces
transcripts representing all possible combinations of introns and exons which
are then filtered according to the specifications in the config.
The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes;

use warnings ;
use strict;
use feature qw(say) ;

use Bio::DB::HTS;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(convert_to_ucsc_name);
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Runnable::RefineSolexaGenes;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: It will fetch all the proto transcript from the region specified in 'iid' and fetch all the
              reads representing the introns. Multiple intron files can be given and they can contain non spliced
              reads if MIXED is higher than 0.
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
    my( $self) = @_;

    my $input_id = $self->input_id;
    # I want to be able to use either slices or gene stable ids
    my $genes_db = $self->get_database_by_name('input_db');
    $genes_db->dbc->disconnect_when_inactive(0);
    my $reference_db = $self->get_database_by_name('dna_db');
    $reference_db->dbc->disconnect_when_inactive(0);
    $genes_db->dnadb($reference_db);
    $self->gene_slice_adaptor($reference_db->get_SliceAdaptor);
    $self->hrdb_set_con($self->get_database_by_name('output_db'), 'output_db');
    $self->hrdb_get_con('output_db')->dbc->disconnect_if_idle if ($self->param('disconnect_jobs'));
    my @rough_genes;
    my $real_slice_start;
    my $real_slice_end;
    my $chr_slice;
    if ($self->is_slice_name($input_id)) {
        my $genes;
        my $slice = $self->fetch_sequence($input_id, $genes_db);
        $real_slice_start = $slice->start;
        $real_slice_end = $slice->end;
        $chr_slice = $reference_db->get_SliceAdaptor->fetch_by_region( 'toplevel', $slice->seq_region_name);

        if ( $self->param('model_ln') ) {
            $genes = $slice->get_all_Genes_by_type( undef,$self->param('model_ln') );
            print STDERR "Got " .  scalar(@$genes) . " genes with logic name " . $self->param('model_ln') ."\n";
        }
        else {
            $genes = $slice->get_all_Genes(undef, undef, 1);
            print STDERR "Got " .  scalar(@$genes) . "  genes  \n";
        }
        foreach my $gene ( @$genes ) {
            # put them on the chromosome
            $gene = $gene->transfer($chr_slice);
            # reject genes that are from a different slice that overlap our slice at the start or end
            # say the models has to be > 10% on the slice
            my $os = $slice->start;
            $os = $gene->start if $gene->start > $slice->start;
            my $oe = $slice->end;
            $oe = $gene->end if $gene->end < $slice->end;
            my $overlap = $oe - $os +1;
            my $gc = int(($overlap / $gene->length) * 1000) / 10;
            my $sc =  int(($overlap / $slice->length) * 1000) /10;
            if ( $gc <= 10 && $sc <= 10) {
                print "Gene ", $gene->display_id, " has $gc% overlap with the slice\nSlice has $sc% overlap with the gene\n  Rejecting\n";
                next;
            }
            $real_slice_start = $gene->seq_region_start < $real_slice_start ? $gene->seq_region_start : $real_slice_start;
            $real_slice_end = $gene->seq_region_end > $real_slice_end ? $gene->seq_region_end : $real_slice_end;
            push(@rough_genes,$gene);
        }
        print STDERR "Got " . scalar(@rough_genes) . "  genes after filtering boundary overlaps  \n";
    }
    else {
        my $gene = $genes_db->get_GeneAdaptor->fetch_by_stable_id($input_id);
        $chr_slice = $reference_db->get_SliceAdaptor->fetch_by_region( 'toplevel', $gene->slice->seq_region_name);
        $real_slice_start = $gene->seq_region_start;
        $real_slice_end = $gene->seq_region_end;
        push(@rough_genes, $gene);
    }
    $self->chr_slice($chr_slice);
    if (scalar(@rough_genes)) {
        $self->create_analysis;
        if ( $self->param('intron_bam_files') ) {
            foreach my $intron_files (@{ $self->param('intron_bam_files')} ) {
                my $sam = Bio::DB::HTS->new(
                        -bam => $intron_files->{file},
                        -expand_flags => 1,
                        );
                $self->throw("Bam file " . $intron_files->{file} . "  not found \n") unless ($sam);
                my $count = 0;
                my $seq_region_name = $self->param('wide_use_ucsc_naming') ? convert_to_ucsc_name($self->chr_slice->seq_region_name, $self->chr_slice) : $self->chr_slice->seq_region_name;
                my $segment = $sam->segment($seq_region_name, $real_slice_start, $real_slice_end);
                $self->throw("Bam file segment not found for slice ".$self->chr_slice->seq_region_name."\n") unless ($segment);
                # need to seamlessly merge here with the dna2simplefeatures code
                $self->bam_2_intron_features($segment,$intron_files);
            }
        }
        else {
            # pre fetch all the intron features
            $self->dna_2_intron_features($real_slice_start, $real_slice_end);
        }
        $genes_db->dbc->disconnect_when_inactive(1) if ($self->param('disconnect_jobs'));
        my $runnable = Bio::EnsEMBL::Analysis::Runnable::RefineSolexaGenes->new (
                -analysis     => $self->analysis,
                -retained_intron_penalty => $self->param('retained_intron_penalty'),
                -filter_on_overlap => $self->param('filter_on_overlap'),
                -min_intron_size => $self->param('min_intron_size'),
                -max_intron_size => $self->param('max_intron_size'),
                -single_exon_model => $self->param('single_exon_model'),
                -min_single_exon => $self->param('min_single_exon'),
                -single_exon_cds => $self->param('single_exon_cds'),
                -strict_internal_splice_sites => $self->param('strict_internal_splice_sites'),
                -strict_internal_end_exon_splice_sites => $self->param('strict_internal_end_exon_splice_sites'),
                -best_score => $self->param('best_score'),
                -other_isoforms => $self->param('other_isoforms'),
                -other_num => $self->param('other_num'),
                -max_num => $self->param('max_num'),
                -bad_models => $self->param('bad_models'),
                -trim_utr => $self->param('trim_utr'),
                -max_3prime_exons => $self->param('max_3prime_exons'),
                -max_3prime_length => $self->param('max_3prime_length'),
                -max_5prime_exons => $self->param('max_5prime_exons'),
                -max_5prime_length => $self->param('max_5prime_length'),
                -reject_intron_cutoff => $self->param('reject_intron_cutoff'),
                -max_recursions => $self->param('max_recursions'),
                -chr_slice => $self->chr_slice,
                -query => $self->chr_slice,
                -rough_models => \@rough_genes,
                -intron_features => $self->intron_features,
                -extra_exons => $self->extra_exons,
                );
        $self->runnable($runnable);
    }
    else {
        $self->input_job->autoflow(0);
        $self->complete_early('No genes to process');
    }
}


=head2 run

  Arg [1]   : None
  Function  : Overrides run as we want to be able to disconnect from the database if 'disconnect_jobs'
              is set to 1.
  Returntype: None
  Exceptions: None

=cut

sub run {
    my ($self) = @_;

    $self->throw("Can't run - no runnable objects") unless ( $self->runnable );
    $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
    my ($runnable) = @{$self->runnable};
    $runnable->run;
    $self->output($runnable->output);
    return 1;
}


=head2 write_output

  Arg [1]   : None
  Function  : It writes the genes in the database specified by 'output_db' and it write the
              clustered introns in the dna_align_feature table if 'write_introns' is 1.
  Returntype: None
  Exceptions: Throws if genes or introns are not all stored

=cut

sub write_output {
    my ($self) = @_;

    my $outdb = $self->hrdb_get_con('output_db');
    $outdb->dbc->disconnect_when_inactive(0);
    my $gene_adaptor = $outdb->get_GeneAdaptor;

    my $fails = 0;
    my $total = 0;
    my $analysis = $self->analysis;
    my $source = $analysis->logic_name;
    $source =~ s/_rnaseq$//;
    foreach my $gene (@{$self->output}) {
        $gene->analysis($analysis);
        $gene->source($source);
        foreach my $tran ( @{$gene->get_all_Transcripts} ) {
            $tran->analysis($self->analysis);
        }
        # filter single exon genes that may have been made through UTR trimming
        my @exons = @{$gene->get_all_Exons};
        if ( scalar(@exons == 1 )) {
            if ( $self->param('single_exon_model') ) {
                next unless ($exons[0]->length >= $self->param('min_single_exon'));
                $gene->biotype($self->param('single_exon_model'));
            }
            else {
                # dont store it
                next;
            }
        }

        eval {
            $gene_adaptor->store($gene);
        };
        if ($@){
            $self->warning("Unable to store gene!!\n$@");
            $fails++;
        }
        $total++;
    }
    if ($fails > 0) {
        $self->throw("Not all genes could be written successfully " .
                "($fails fails out of $total)");
    }
    if ($self->param('write_introns')) {
        my $intron_adaptor = $outdb->get_DnaAlignFeatureAdaptor;
        my $intron_fails = 0;
        my $intron_total = 0;
        my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $self->param('introns_logic_name'));
        foreach my $intron (@{$self->intron_features}) {
            $intron->start($intron->start+1);
            $intron->end($intron->end-1);
            $intron->analysis($analysis);
            eval {
                $intron_adaptor->store($intron);
            };
            if ($@){
                $self->warning("Unable to store DnaAlignFeature!!\n$@");
                $intron_fails++;
            }
            $intron_total++;
        }
        if ($intron_fails > 0) {
            $self->throw("Not all introns could be written successfully ($intron_fails fails out of $intron_total)");
        }
    }
    if ($total == 0) {
        $self->input_job->autoflow(0);
    }
}


=head2 bam_2_intron_features

  Arg [1]    : Bio::DB::HTS::segment
  Arg [2]    : Arrayref of hashref containing the information about the introns
  Description: Fetches all alignments from the bam file segment, collapses them down into a
               non redundant set and builds a Bio::EnsEMBL::DnaDnaAlignFeature to
               represent it, then stores it in $self->intron_features
               analyses splice sites for consensus and non consensus splices as this data is
               not stored in the BAM.
               Also checks for small exons defined by a single read splicing at least twice
               stores any additional exons found this way in $self->extra_exons
  Returntype : Integer, 1
  Exceptions : None

=cut

sub bam_2_intron_features {
    my ($self,$segment,$intron_files) = @_;
    my $slice_adaptor = $self->gene_slice_adaptor;
    my @ifs;
    my $extra_exons = $self->extra_exons;
    my %id_list;
    my %read_groups;
    if (  $intron_files->{groupname} && scalar(@{$intron_files->{groupname}} > 0 ) ) {
        my @groups = @{$intron_files->{groupname}};
        print "Limiting to read groups ";
        foreach my $group ( @groups ) {
            print " $group";
            $read_groups{$group} = 1;
        }
        print "\n";
    }
    my $iterator = $segment->features(-iterator=>1);
    while (my $read = $iterator->next_seq) {
        my $spliced;
# ignore unspliced reads if the bam file is a mixture of spliced and
# unspliced reads
        if ( $intron_files->{mixed_bam} ) {
            $spliced = $read->get_tag_values('XS');
            next unless $spliced;
        }
# filter by read group if needed

# need to recreate the ungapped features code as the
# auto splitting code does not seem to work with > 2 features
        if ( $intron_files->{groupname}  && scalar(@{$intron_files->{groupname}} > 0 )) {
            next unless ($read_groups{$read->get_tag_values('RG')}) ;
        }
        my @mates = sort { $a->[2] <=> $b->[2] } @{$self->ungapped_features($read)};

# if mates > 2 then we have a possibility of adding in some extra exons into our rough models
# as the read has spliced into and out of an exon
# lets make them unique
        if ( scalar(@mates) > 2 ) {
            my $string;
            for ( my $i = 0 ; $i <= $#mates  ; $i++ ) {

                my $start  = $mates[$i]->[2];
                my $end    = $mates[$i]->[3];
                my $hstrand = $read->strand;
                $string .= $start .":" if $i > 0 ;
                $string .= $end .":" if $i < $#mates ;
            }
            $extra_exons->{$string} ++;
        }
        my $strand = $read->target->strand;
        if   ($intron_files->{mixed_bam} ) {
            $strand = 1 if $spliced eq '+';
            $strand = -1 if $spliced eq '-';
        }
# print "\nREAD " . $read->cigar_str;
        my $offset;
        for ( my $i = 0 ; $i <= $#mates  ; $i++ ) {
#   print "\n";
# intron reads should be split according to the CIGAR line
# the default split function seems to ad
# we want the ungapped features to make our introns
            my $name   = $read->seq_id;
# we dont allow . in the seq region name as we use them to delimit our paths
            $name =~ s/\./*/g;
            my $start  = $mates[$i]->[2];
            my $end    = $mates[$i]->[3];
            my $cigar  = $mates[$i]->[4];
            my $hstrand = $read->strand;
            next if $i == $#mates;
            my $unique_id = $name . ":" .
            $mates[$i]->[3] . ":" .
            $mates[$i+1]->[2] . ":" .
            $strand ;
            $id_list{$unique_id} ++;
        }
    }
# collapse them down and make them into simple features
    foreach my $key ( keys %id_list ) {
# filter on score if appropriate
        if ( $intron_files->{depth} ) {
            if ( $intron_files->{depth} > $id_list{$key} ) {
#print "Rejecting on score " . $id_list{$key} ."\n";
                next;
            }
        }
        my @data = split(/:/,$key) ;
        my $length =  $data[2] - $data[1] -1;
        next unless $length > 0 ;
        my $name = $data[0]. ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":";

        my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new(
            -start => $data[1],
            -end => $data[2],
            -strand => $data[3],
            -hstart => 1,
            -hend => $length,
            -hstrand => 1,
            -slice => $self->chr_slice,
            -analysis => $self->analysis,
            -score =>  $id_list{$key},
            -hseqname => "$name",
            -cigar_string => $length ."M",
            -align_type => 'ensembl',
        );
        my $canonical = 1;
# figure out if its cannonical or not
        my $left_splice = $slice_adaptor->fetch_by_region('toplevel',
            $if->seq_region_name,
            $if->start+1,
            $if->start+2,
            $if->strand
        );
        my $right_splice = $slice_adaptor->fetch_by_region('toplevel',
            $if->seq_region_name,
            $if->end-2,
            $if->end-1,
            $if->strand
        );
#  print "KEY $key " . $if->score ."\n";;
#  print "LEFT  " . $left_splice->start ." " . $left_splice->end  ." " . $left_splice->strand ." " . $left_splice->seq . "\n";
#  print "RIGHT " . $right_splice->start ." " . $right_splice->end  ." " . $right_splice->strand  ." " . $right_splice->seq ."\n\n";


        if ( $left_splice->seq eq 'NN' && $right_splice->seq eq 'NN' ) {
            warn("Cannot find dna sequence for $key this is used in detecting non cannonical splices\n");
        }
        else {
# is it cannonical
            if ( $if->strand  == 1 ) {
#	print "Splice type " . $left_splice->seq ."-".  $right_splice->seq ." ";
# is it GTAG?
                unless ( $left_splice->seq eq 'GT' && $right_splice->seq eq 'AG' ) {
                    $canonical = 0;
                }
            }
            else {
#	print "Splice type " . $right_splice->seq ."-".  $left_splice->seq ." ";
# is it GTAG?
                unless ( $right_splice->seq eq 'GT' && $left_splice->seq eq 'AG' ) {
                    $canonical = 0;
                }
            }
        }
        if ( $canonical ) {
            $if->hseqname($if->hseqname."canonical");
        }
        else {
            $if->hseqname($if->hseqname."non canonical");
        }
#print "INTRONS ".  $if->hseqname ."\n";
        push @ifs , $if;
    }
# sort them
    @ifs = sort {$a->start <=> $b->start} @ifs;
    if ($self->param('filter_on_overlap')) {
        my @tmp_array;
        my $threshold = $self->param('filter_on_overlap');
        my $array_length = scalar(@ifs);
        if ($array_length > 1) {
            for (my $j = 0; $j < $array_length-1; $j++) {
                my $k = 0;
                my $count = 1;
                my $overlapped_support = 0;
                while () {
                    ++$k;
                    if ($count > $threshold) {
                        if ($overlapped_support < $ifs[$j]->score) {
#                          print STDERR "\t",$ifs[$j+$k]->hseqname, ': ', $ifs[$j+$k]->start, ':', $ifs[$j+$k]->end, "\n";
                            push (@tmp_array, $ifs[$j]);
                        }
                        else {
#                          print STDERR 'THROWING: ', $ifs[$j]->hseqname, ': ', $ifs[$j]->start, ':', $ifs[$j]->end, "\n";
                        }
                        last;
                    }
                    $overlapped_support += $ifs[$j+$k]->score;
                    if (($ifs[$j]->end < $ifs[$j+$k]->start) or (($j+$k) == $array_length-1)) {
#                      print STDERR "\t",$ifs[$j+$k]->hseqname, ': ', $ifs[$j+$k]->start, ':', $ifs[$j+$k]->end, "\n";
                        push (@tmp_array, $ifs[$j]);
                        last;
                    }
#                      print STDERR $ifs[$j+$k]->hseqname, "\n";
                    next unless ($ifs[$j]->strand == $ifs[$j+$k]->strand);
                    ++$count;
                }
            }
            @ifs = @tmp_array;
        }
    }
#  print STDERR 'RES: ', scalar(@ifs), "\n";
    $self->intron_features(\@ifs);
    $self->extra_exons($extra_exons);
    print STDERR "Got " . scalar(@ifs)  . " unique introns  " ;
    print STDERR " and " . scalar(keys %$extra_exons) . " potential novel exons from " . $intron_files->{file} . "\n";
    return;
}


=head2 ungapped_features

 Arg [1]    : Bio::DB::HTS::Alignment
 Description: Create introns based on the cigar line of the read
 Returntype : Arrayref of array
              0 -> read name
              1 -> sequence id
              2 -> start
              3 -> end
              4 -> length of the match
 Exceptions : Throws if it cannot parse the rad cigar line

=cut

sub ungapped_features {
    my ($self,$read) = @_;
    my @ugfs;
    my @tmp_ugfs;
    my $string = $read->cigar_str;
    my $start = $read->start;
    my $end = $read->end;
    my @pieces = ( $string =~ /(\d*[MDN])/g );
    for my $piece ( @pieces ) {
        my ($length) = ( $piece =~ /^(\d*)/ );
        if( $length eq "" ) { $length = 1 }
        if( $piece =~ /M$/ ) {
            #
            # MATCH
            #
            my ( $qstart, $qend);
            $qstart = $start;
            $qend = $start + $length - 1;
            $start = $qend + 1;

            my $ugf;
            $ugf->[0] = $read->query->name;
            $ugf->[1] = $read->seq_id;
            $ugf->[2] = $qstart;
            $ugf->[3] = $qend;
            $ugf->[4] = $length."M";
            push @tmp_ugfs, $ugf;
            #print "UNGAPPED " .$ugf->[2] .
            #" " . $ugf->[3] . " " . $ugf->[4] ."\n";
        }
        elsif( $piece =~ /N$/ ) {
            #
            # INSERT
            #
            $start += $length;
            push @tmp_ugfs,"intron";
        }
        elsif( $piece =~ /D$/ ) {
            #
            # DELETION
            #
            $start += $length;
            push @tmp_ugfs,"deletion";
        } else {
            $self->throw( "Illegal cigar line $string!" );
        }
    }
    # only return the UGFS either side of splices
    my %used_pieces;
    foreach ( my $i = 0 ; $i < scalar(@pieces); $i++ )  {
        my $piece = $pieces[$i];
        if ( $piece =~ /\d*N/) {
            # it's a splice push the Matches either side of it
            for ( my $j = $i-1 ; $j >= 0 ; $j-- ) {
                if ( $tmp_ugfs[$j] && $pieces[$j] =~ /\d*M/ )  {
                    unless ( $used_pieces{$j} ) {
                        my $ugf =  $tmp_ugfs[$j];
                        $self->throw("Cannot find ugf $j\n") unless $ugf;
                        push @ugfs, $ugf;
                        $used_pieces{$j} =1;
                        last ;
                    }
                }
            }
            for ( my $j = $i+1 ; $j < scalar(@pieces)  ; $j++ ) {
                if ( $tmp_ugfs[$j] && $pieces[$j] =~ /\d*M/ )  {
                    unless ( $used_pieces{$j} ) {
                        my $ugf =  $tmp_ugfs[$j];
                        $self->throw("Cannot find ugf $j\n") unless $ugf;
                        push @ugfs, $ugf;
                        $used_pieces{$j} =1;
                        last ;
                    }
                }
            }
        }
    }
    return \@ugfs;
}


=head2 dna_2_intron_features

 Arg [1]    : Integer start
 Arg [2]    : Integer end
 Description: Fetches all dna_align_features from the intron db that lie within
              the range determined by start and end, collapses them down into a
              non redundant set and builds a Bio::EnsEMBL::DnaAlignFeature to
              represent it, then stores it in $self->intron_features
              also checks for small exons defined by a single read splicing at least twice
              stores any additional exons found this way in $self->extra_exons
 Returntype : None
 Exceptions : None

=cut

sub dna_2_intron_features {
    my ($self, $start, $end) = @_;
    my @ifs;
    my %id_list;
    my $intron_slice_adaptor = $self->get_database_by_name('intron_output_db')->get_SliceAdaptor;
    my $intron_slice = $intron_slice_adaptor->fetch_by_region(
            'toplevel',
            $self->chr_slice->seq_region_name,
            $start,
            $end,
            );
    # featch all the dna_align_features for this slice by logic name
    my @reads;
    print STDERR  "Fetching reads with logic names: ";
    foreach my $logic_name ( @{$self->param('logicname')} ) {
        print STDERR "$logic_name ";
        push @reads, @{$intron_slice->get_all_DnaAlignFeatures($logic_name)};
    }
    print STDERR "\n";
    # fetch them all if no logic name is Supplied
    if (scalar( @{$self->param('logicname')} ) == 0 ) {
        @reads =  @{$intron_slice->get_all_DnaAlignFeatures()};
    }
    print STDERR "Got " . scalar(@reads) . " reads\n";

    while ( scalar @reads > 0 ) {
        my $read = pop(@reads);
        my $type = 'canonical';
        $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
        $read = $read->transfer($self->chr_slice);
        my @ugfs = $read->ungapped_features;
        @ugfs = sort { $a->start <=> $b->start } @ugfs;
        next if ( scalar(@ugfs) == 0 ) ;
        for ( my $i = 0 ; $i < scalar(@ugfs) - 1 ; $i++ ) {
        # one read can span several exons so make all the features
        # cache them by internal boundaries
        # we use . to deliminate entries in our paths so dont allow them in the seq_region_name or it wont work
            my $name = $read->seq_region_name;
            $name =~ s/\./*/g;
            my $unique_id = $name . ":" .
            $ugfs[$i]->end . ":" .
            $ugfs[$i+1]->start . ":" .
            $read->strand .":$type";
            $id_list{$unique_id} ++;
        }
    }
    print STDERR "Got " . scalar( keys %id_list ) . " collapsed introns\n";
    # collapse them down and make them into simple features
    foreach my $key ( keys %id_list ) {
        my @data = split(/:/,$key) ;
        my $length =  $data[2] - $data[1] -1;
        next unless $length > 0;
        my $name = $self->chr_slice->seq_region_name . ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":".$data[4];
        my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new(
            -start => $data[1],
            -end => $data[2],
            -strand => $data[3],
            -hstart => 1,
            -hend => $length,
            -hstrand => 1,
            -slice => $self->chr_slice,
            -analysis => $self->analysis,
            -score =>  $id_list{$key},
            -hseqname => $name,
            -cigar_string => $length ."M",
            -align_type => 'ensembl',
        );
        push @ifs , $if;
    }
# sort them
    @ifs = sort {$a->start <=> $b->start} @ifs;
    $self->intron_features(\@ifs);
    print STDERR "Got " . scalar(@ifs) . " intron features\n";
    return;
}




##################################################################
# Containers

=head2 gene_slice_adaptor

 Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::SliceAdaptor
 Description: Getter/setter for a SliceAdaptor object on the input database
 Returntype : Bio::EnsEMBL::DBSQL::SliceAdaptor
 Exceptions : Throws if Arg[1] is not a Bio::EnsEMBL::DBSQL::SliceAdaptor

=cut

sub gene_slice_adaptor {
    my ($self, $val) = @_;

    if (defined $val) {
      $self->throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::SliceAdaptor!')
        unless (ref($val) eq 'Bio::EnsEMBL::DBSQL::SliceAdaptor');
      $self->param('gene_slice_adaptor', $val);
    }

    return $self->param('gene_slice_adaptor');
}


=head2 chr_slice

 Arg [1]    : (optional) Bio::EnsEMBL::Slice representing the whole sequence
 Description: Getter/setter for the whole sequence given as input_id
 Returntype : Bio::EnsEMBL::Slice
 Exceptions : None

=cut

sub chr_slice {
    my ($self, $val) = @_;

    if (defined $val) {
        $self->param('_chr_slice', $val);
    }

    return $self->param('_chr_slice');
}


=head2 intron_features

 Arg [1]    : (optional) Arrayref of Bio::EnsEMBL::DnaDnaAlignFeature representing the introns
 Description: Getter/setter for the introns. It make sure they are sorted based on the start
 Returntype : Arrayref of Bio::EnsEMBL::DnaDnaAlignFeature
 Exceptions : None

=cut

sub intron_features {
    my ($self, $val) = @_;

    if (!$self->param_is_defined('_introns')) {
        $self->param('_introns', []);
    }
    if (defined $val) {
#       make sure it is still sorted
        my @introns = sort { $a->start <=> $b->start } @{$self->param('_introns')}, @$val;
        $self->param('_introns', \@introns);
    }
    return $self->param('_introns');
}


=head2 extra_exons

 Arg [1]    : (optional) Hashref of String representing extra exons
 Description: Getter/setter for the extra exons.
 Returntype : Arrayref of String
 Exceptions : None

=cut

sub extra_exons {
    my ($self, $val) = @_;

    if (!$self->param_is_defined('extra_exons')) {
        $self->param('extra_exons', {});
    }
    if (defined $val) {
        $self->param('extra_exons', $val);
    }

    return $self->param('extra_exons');
}

1;
