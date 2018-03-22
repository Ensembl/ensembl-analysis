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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(convert_to_genes Transcript_info set_start_codon set_stop_codon list_evidence attach_Slice_to_Transcript);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info print_Gene attach_Analysis_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::KillList::KillList;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    rejected_set => [],
    exon_mask_list => [],
    gene_mask_list => [],
    id_pool_bin => 0,
    id_pool_index => 0,
  }
}

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Function  : This processes the input id, fetches both the sequence and the
  align features to run with and filters the align features where appropriate on
  the basis of a kill list and comparision to existing gene models
  Returntype: 1
  Exceptions:
  Example   :

=cut


sub fetch_input{
  my ($self) = @_;

  $self->create_analysis;
  my $dna_db = $self->get_database_by_name('dna_db');
  $self->hrdb_set_con($dna_db, 'dna_db');
  $self->query($self->fetch_sequence($self->parse_input_id, $self->get_database_by_name('source_db', $dna_db), $self->REPEATMASKING, $self->SOFTMASKING));
  my %hit_list;

  my $killed_count = 0;
  my $feature_count = 0;
  my %protein_count;
  foreach my $logic_name(@{$self->PAF_LOGICNAMES}){
    #print "LOGIC NAME : ",$logic_name, ' $ ', $self->paf_slice->get_seq_region_id, "\n";
    my $features = $self->query->get_all_ProteinAlignFeatures($logic_name, $self->PAF_MIN_SCORE_THRESHOLD);
    my %unique;
    foreach my $feature(@$features){
      $unique{$feature->hseqname} = 1;
    }
    #print "****HAVE ".@$features." features with ".$logic_name." and min score ".$self->PAF_MIN_SCORE_THRESHOLD."  with ".keys(%unique)." unique hit names from ".$self->paf_slice->adaptor->dbc->dbname."*****\n";
    logger_info("HAVE ".@$features." with ".$logic_name." and min score ".$self->PAF_MIN_SCORE_THRESHOLD);
    $feature_count += scalar(@$features);
    my $ids_to_ignore = {};
    $ids_to_ignore = %{$self->generate_ids_to_ignore($features)} if (@{$self->BIOTYPES_TO_MASK});
    #print "HAVE ".keys(%ids_to_ignore)." ids to ignore\n";
    logger_info("HAVE ".keys(%$ids_to_ignore)." ids to ignore");

     my %kill_list ;
     if ( scalar(@$features) > 0 ) {
       %kill_list = %{$self->kill_list} if($self->USE_KILL_LIST);
     }

  FEATURE:foreach my $feature(@$features){

      #print $feature->hseqname." is in this ID with paf_id ".$feature->dbID."\n";

      $protein_count{$feature->hseqname} = 1;
      my $temp_id = $feature->hseqname;
      $temp_id =~ s/\.\d+//;
      next FEATURE if($self->overlaps_fiveprime_end_of_slice($feature, $self->query));
      $killed_count++ if($kill_list{$temp_id});
      logger_info("SKIPPING ".$feature->hseqname." on kill list")
        if ($kill_list{$temp_id});
      next FEATURE if($kill_list{$temp_id});
      logger_info("SKIPPING ".$feature->hseqname." as greater than upper threshold")
        if($self->PAF_UPPER_SCORE_THRESHOLD and
           $feature->score < $self->PAF_UPPER_SCORE_THRESHOLD);
      next FEATURE if($self->PAF_UPPER_SCORE_THRESHOLD and
                      $feature->score < $self->PAF_UPPER_SCORE_THRESHOLD);
      logger_info("SKIPPING ".$feature->hseqname." on the mask list")
        if($self->PRE_GENEWISE_MASK &&
           exists $ids_to_ignore->{$feature->hseqname}) ;
      next FEATURE if($self->PRE_GENEWISE_MASK &&
                      exists $ids_to_ignore->{$feature->hseqname});
      if($self->use_id){
        if($feature->hseqname eq $self->use_id){
          push(@{$hit_list{$feature->hseqname}}, $feature);
        }
      }else{
        push(@{$hit_list{$feature->hseqname}}, $feature);
      }
    }
  }
  if($self->use_id && !(keys(%hit_list))){
    my $feature = $self->feature_factory->create_feature_pair($self->query->start,
                                                              $self->query->end,
                                                              1, 1, 1, 1, 1,
                                                              $self->use_id);
    $hit_list{$self->use_id} = [$feature];
  }
  logger_info("HAVE ".keys(%protein_count)." unique protein ids");
  logger_info("HAVE ".$self->id_pool_bin." and ".$self->id_pool_index) if($self->id_pool_bin && $self->id_pool_index);
  if($self->id_pool_bin and $self->id_pool_index){
    my @local_ids = keys %hit_list;
    @local_ids = sort {$a cmp $b} @local_ids;
    my (@restricted_list);
    for (my $i = $self->id_pool_index - 1;
         $i < @local_ids; $i += $self->id_pool_bin) {
      push @restricted_list, $local_ids[$i];
    }
    %hit_list = map { $_ => $hit_list{$_} } @restricted_list;
    #print "Hit name in list\n" if(exists $hit_list{'Q2PHF0.1'});

  }
  if($killed_count && $feature_count && $killed_count == $feature_count){
    $self->warning("RunnableDB:BlastMiniGenewise Killed all the potential ".
            "ids in this run there are NO IDS TO RUN");
    return;
  }
  if(scalar(keys(%hit_list)) == 0 && $feature_count >= 1){
    $self->warning("RunnableDB:BlastMiniGenewise not all ids were killed but you ".
            "still have NO IDS TO RUN with");
    return;
  }
  $self->create_bmg_runnables(\%hit_list);
  return 1;
}


=head2 generate_ids_to_ignore

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Features
  Function  : to produce a list of ids which shouldnt be considered as their features
  overlap with existing genes
  Returntype: hashref
  Exceptions:
  Example   :

=cut



sub generate_ids_to_ignore{
  my ($self, $features) = @_;
  my $exon_mask_list = $self->create_mask_list;
  my @exonmask_regions = @{$exon_mask_list};
  my %features;
  foreach my $f (@$features) {
    push @{$features{$f->hseqname}}, $f;
  }
  my %ids_to_ignore;
 SEQID:foreach my $sid (keys %features) {
    my $ex_idx = 0;
  FEAT: foreach my $f (sort {$a->start <=> $b->start}
                       @{$features{$sid}}) {
      for( ; $ex_idx < @exonmask_regions; ) {
        my $mask_exon = $exonmask_regions[$ex_idx];
        if ($mask_exon->{'start'} > $f->end) {
          # no exons will overlap this feature
          next FEAT;
        } elsif ( $mask_exon->{'end'} >= $f->start) {
          # overlap
          $ids_to_ignore{$f->hseqname} = 1;
          next SEQID;
        }  else {
          $ex_idx++;
        }
      }
    }
  }
  return \%ids_to_ignore;
}


=head2 post_mask

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Genes
  Function  : This method filters the prediction genes on the basis of existing
  predictions removing those which overlap either at the exon or the transcript
  level depending on the configuration
  Returntype: arrayref of Genes
  Exceptions:
  Example   :

=cut



sub post_mask{
  my ($self, $genes) = @_;
  my $exon_mask_list = $self->exon_mask_list;
  my @masked_genes;
  $exon_mask_list = $self->gene_mask_list if($self->GENE_BASED_MASKING);
  my @exonmask_regions = @{$exon_mask_list};

  foreach my $gene(@$genes){
    my $keep_gene = 1;
    my $mask_reg_idx = 0;
    my @test_regions;
    my @exons = sort {$a->start <=> $b->start} ( @{$gene->get_all_Exons});
    if($self->GENE_BASED_MASKING){
      @test_regions = ({start => $exons[0]->start, end => $exons[-1]->end});
    }else{
      @test_regions = map { { start => $_->start, end => $_->end } } @exons;
    }
  FEAT:foreach my $f(@test_regions){
      for( ; $mask_reg_idx < @exonmask_regions ;){
        my $mask_reg = $exonmask_regions[$mask_reg_idx];
        if ($mask_reg->{'start'} > $f->{'end'}) {
          # no mask regions will overlap this feature
          next FEAT;
        }
        elsif ( $mask_reg->{'end'} >= $f->{'start'}) {
          # overlap
          $keep_gene = 0;
          last FEAT;
        }
        else {
          $mask_reg_idx++;
        }
      }
    }
    push(@masked_genes, $gene) if($keep_gene);
    logger_info("SKIPPING ".Gene_info($gene)." masked out") if(!$keep_gene);
  }
  return \@masked_genes;
}



=head2 parse_input_id

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : String
  Function  : parse the input id to find all the appropriate information needed
  to fetch the input
  Returntype:
  Exceptions:
  Example   :

=cut



sub parse_input_id{
  my ($self, $input_id) = @_;

  $input_id = $self->input_id if(!$input_id);

  my ($slice_name, $first, $second) =
    ($input_id =~ /^([^:]*:[^:]*:[^:]*:[^:]*:[^:]*:[^:]*)(?::([^:]+):([^:]+))?$/);

  $self->throw("RunnableDB::BlastMiniGenewise Your input id ".$input_id." doesn't ".
        "match the regex and seems to be missing a slice name in the form ".
        "coordsystem:version:seq region name:start:end:strand") if(!$slice_name);

  if($first && $second){
    if($first =~ /^\d+$/ and
       $second =~ /^\d+$/){

      $self->id_pool_bin($first);
      $self->id_pool_index($second);

      if($self->id_pool_index > $self->id_pool_bin or
         $self->id_pool_index < 1){

        $self->warning("RunnableDB::BlastMiniGenewise could not get sensible values for ".
                "id pool bin ".$self->id_pool_bin." or id pool index ".
                $self->id_pool_index." will us complete set of proteins");
        $self->id_pool_bin(0);
        $self->id_pool_index(0);
      }

    }else{
      #assumimg this is a protein id and a logic name
      $self->PAF_LOGICNAMES([$first]);
      $self->use_id($second);

    }
  }

  return $slice_name;
}


=head2 create_bmg_runnables

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : hashref of hit names pointing at arrays of features
  Arg [3]   : Seqfetcher object
  Arg [4]   : Bio::EnsEMBL::Analysis object
  Function  : to create BlastMiniGenewise runnables on the basis of the features
  given
  Returntype:
  Exceptions:
  Example   :

=cut



sub create_bmg_runnables{
  my ($self, $hits_hash, $seqfetcher, $analysis) = @_;
  $self->throw("RunnableDB::BlastMiniGenewise Can't create BlastMiniGenewise ".
        "runnables without a hash of hits") unless($hits_hash);
  $seqfetcher = $self->seqfetcher if(!$seqfetcher);
  $analysis = $self->analysis if(!$analysis);
  my %bmg_params = %{$self->BLASTMINIGENEWISE_PARAMETERS};
  my %params_hash = %{$self->parameters_hash};
  if($self->LIMIT_TO_FEATURE_RANGE){
    foreach my $id(keys(%$hits_hash)){
      my $features = $hits_hash->{$id};
      foreach my $feature(@$features){
        #This should be made a non overlapping range
        my $start = $feature->seq_region_start - $self->FEATURE_RANGE_PADDING;
        my $end = $feature->seq_region_end + $self->FEATURE_RANGE_PADDING;
        $start = 1 if($start < 1);
        $end = $self->query->seq_region_length if($end > $self->query->seq_region_length);
        my $name = ($self->query->coord_system->name.":".
                    $self->query->coord_system->version.":".
                    $self->query->seq_region_name.":".
                    $start.":".$end.":".
                    $self->query->strand);
        my $query = $self->fetch_sequence($name, $self->get_database_by_name('source_db', $self->hrdb_get_con('dna_db')), $self->REPEATMASKING, $self->SOFTMASKING);
        logger_info("Creating BlastMiniGenewise Runnable over a limited range with ".$id." and ".$seqfetcher." to run on ".$query->name);
        my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise->
          new(
              -query => $query,
              -analysis => $analysis,
              -seqfetcher => $seqfetcher,
              -minigenewise_options => $self->MINIGENEWISE_PARAMETERS,
              -genewise_options => $self->GENEWISE_PARAMETERS,
              -multiminigenewise_options => $self->MULTIMINIGENEWISE_PARAMETERS,
              -ids => [$id],
              %bmg_params,
              %params_hash,
             );
        $self->runnable($runnable);
      }
    }
  }else{
    my @ids = sort keys(%$hits_hash);
    #print "***Creating BlastMiniGenewise Runnable with ".@ids." ids and ".$seqfetcher." to run on ".$self->query->name."***\n";
    logger_info("Creating BlastMiniGenewise Runnable with ".@ids." ids and ".$seqfetcher." to run on ".$self->query->name);
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastMiniGenewise->
      new(
          -query => $self->query,
          -analysis => $analysis,
          -seqfetcher => $seqfetcher,
          -minigenewise_options => $self->MINIGENEWISE_PARAMETERS,
          -genewise_options => $self->GENEWISE_PARAMETERS,
          -multiminigenewise_options => $self->MULTIMINIGENEWISE_PARAMETERS,
          -ids => \@ids,
          %bmg_params,
          %params_hash,
         );
    $self->runnable($runnable);
  }
}


=head2 create_mask_list

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : to generate the list of coordinates to compare predictions to
  to see if they should be ignored
  Returntype: arrayref of exon regions
  Exceptions: throws if both types of masking are turned on as they are mutually
  exclusive concepts
  Example   :

=cut



sub create_mask_list{
  my ($self, $slice) = @_;

  return [] if(!$self->GENE_BASED_MASKING && !$self->EXON_BASED_MASKING);
  $slice = $self->gene_slice if(!$slice);

  my (@mask_gene_regions, @mask_exon_regions);
  foreach my $biotype(@{$self->BIOTYPES_TO_MASK}){
    GENE:foreach my $gene(@{$slice->get_all_Genes_by_type($biotype)}){
			my @mask_exons = grep { $_->seqname eq $slice->id }
        (sort {$a->start <=> $b->start} @{$gene->get_all_Exons});
      next GENE unless (@mask_exons && scalar(@mask_exons) > 0);
			push @mask_gene_regions, { start => $mask_exons[0]->start,
                                 end   => $mask_exons[-1]->end };

      foreach my $mask_exon (@mask_exons) {
        push @mask_exon_regions, { start => $mask_exon->start,
                                   end   => $mask_exon->end };
      }
    }
  }

  my (@nr_mask_exon_regions, @nr_mask_gene_regions);
  foreach my $mask_exon_reg
    (sort {$a->{'start'} <=> $b->{'start'}} @mask_exon_regions) {
      if (@nr_mask_exon_regions and
          $nr_mask_exon_regions[-1]->{'end'} > $mask_exon_reg->{'start'}) {
        if ($mask_exon_reg->{'end'} > $nr_mask_exon_regions[-1]->{'end'}) {
          $nr_mask_exon_regions[-1]->{'end'} = $mask_exon_reg->{'end'};
        }
      } else {
        push @nr_mask_exon_regions, $mask_exon_reg;
      }
    }
  foreach my $mask_gene_reg
    (sort {$a->{'start'} <=> $b->{'start'}} @mask_gene_regions) {
      if (@nr_mask_gene_regions and
          $nr_mask_gene_regions[-1]->{'end'} > $mask_gene_reg->{'start'}) {
        if ($mask_gene_reg->{'end'} > $nr_mask_gene_regions[-1]->{'end'}) {
          $nr_mask_gene_regions[-1]->{'end'} = $mask_gene_reg->{'end'};
        }
      } else {
        push @nr_mask_gene_regions, $mask_gene_reg;
      }
    }
  $self->throw("RunnableDB:BlastMiniGenewise can't specify masking both Gene and Exon ".
        "level masking ") if($self->GENE_BASED_MASKING &&
                             $self->EXON_BASED_MASKING);
  $self->gene_mask_list(\@nr_mask_gene_regions);
  $self->exon_mask_list(\@nr_mask_exon_regions);
  return \@nr_mask_exon_regions;
}



=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Function  : This runs the created runnables, attached unmasked slices to the
  transcripts, converts the transcripts into genes, sets start and stop codons if
  possible then filters the results
  Returntype: arrayref of Bio::EnsEMBL::Gene objects
  Exceptions: throws if the runnable fails to run
  Example   :

=cut



sub run{
  my ($self) = @_;
  my @transcripts;
  #print "HAVE ".@{$self->runnable}." runnables to run\n";
  foreach my $runnable(@{$self->runnable}){
    my $output;
    $runnable->run;
    $output = $runnable->output;
    push(@transcripts, @$output) if($output);
  }

  my $complete_transcripts = $self->process_transcripts(\@transcripts);
  my @genes = @{convert_to_genes($complete_transcripts, $self->analysis,
                                 $self->OUTPUT_BIOTYPE)};
  #my $hash = $self->group_genes_by_id(\@genes);
  #foreach my $name(keys(%{$hash})){
  #  print "FILTER ".$name." has ".$hash->{$name}." genes before filter\n";
  #}
  my @masked_genes;
  logger_info("HAVE ".@genes." genes before masking") if($self->POST_GENEWISE_MASK);
  if($self->POST_GENEWISE_MASK){
    @masked_genes = @{$self->post_mask(\@genes)};
  }else{
    @masked_genes = @genes;
  }
  logger_info("HAVE ".@masked_genes." genesafter masking") if($self->POST_GENEWISE_MASK);

  $self->filter_genes(\@masked_genes);
  #$hash = $self->group_genes_by_id($self->output);
  #foreach my $name(keys(%{$hash})){
  #  print "FILTER ".$name." has ".$hash->{$name}." genes after filter\n";
  #}
  logger_info("HAVE ".@{$self->output}." genes to return");
  return $self->output;
}


=head2 filter_genes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : arrayref of Genes
  Function  : filters the gene sets using the defined filter
  Returntype: n/a, stores the genes in the output array
  Exceptions:
  Example   :

=cut


sub filter_genes{
  my ($self, $genes) = @_;
  if($self->filter_object){
    my ($filtered_set, $rejected_set) =
      $self->filter_object->filter_genes($genes);
    $self->output($filtered_set);
    $self->rejected_set($rejected_set);
  }else{
   $self->output($genes);
  }
}

=head2 process_transcripts

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Transcripts
  Function  : attaches unmasked slices, sets start and stop codons on transcripts
  Returntype: arrayref of Bio::EnsEMBL::Transcripts
  Exceptions:
  Example   :

=cut



sub process_transcripts{
  my ($self, $transcripts) = @_;
  my @complete_transcripts;
  foreach my $transcript(@$transcripts){
    my $unmasked_slice = $self->fetch_sequence($transcript->slice->name,
                                               $self->output_db);
    attach_Slice_to_Transcript($transcript, $unmasked_slice);
    my $start_t = set_start_codon($transcript);
    my $end_t = set_stop_codon($start_t);
    push(@complete_transcripts, $end_t);
  }
  return \@complete_transcripts
}


=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Function  : writes the accepted and if desired rejected genes to the specified
  database
  Returntype: 1
  Exceptions: throws if failed to store any genes, otherwise warnings are issued
  Example   :

=cut



sub write_output{

  my ($self) = @_;
  my $ga = $self->get_adaptor;
  my $sucessful_count = 0;
  logger_info("WRITE OUTPUT have ".@{$self->output}." genes to write");

  foreach my $gene(@{$self->output}){

    my $attach = 0;
    if(!$gene->analysis){
      my $attach = 1;
      attach_Analysis_to_Gene($gene, $self->analysis);
    }
    if($attach == 0){
    TRANSCRIPT:foreach my $transcript(@{$gene->get_all_Transcripts}){
        if(!$transcript->analysis){
          attach_Analysis_to_Gene($gene, $self->analysis);
          last TRANSCRIPT;
        }
      }
    }
    eval{
      $ga->store($gene);
    };
    if($@){
      $self->warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
    }else{
      $sucessful_count++;
      logger_info("STORED GENE ".$gene->dbID);
    }
  }
  if($sucessful_count != @{$self->output}){
    $self->throw("Failed to write some genes");
  }
  if($self->WRITE_REJECTED){
    $sucessful_count = 0;
    my $biotype = $self->REJECTED_BIOTYPE;
    $biotype = "reject_".$self->analysis->logic_name if(!$biotype);

    foreach my $gene(@{$self->rejected_set}){
      TRANSCRIPT:foreach my $transcript(@{$gene->get_all_Transcripts}){
          $transcript->biotype($biotype);
          if(!$gene->analysis||!$transcript->analysis){
            attach_Analysis_to_Gene($gene, $self->analysis);
            last TRANSCRIPT;
          }
        }

      $gene->biotype($biotype);

      eval{
        $ga->store($gene);
      };
      if($@){
        $self->warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
      }else{
        $sucessful_count++;
      }
    }
    if($sucessful_count != @{$self->rejected_set}){
      $self->throw("Failed to write some rejected genes");
    }
  }
  return 1;
}



=head2 exon_mask_list

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : various, String, Int, Arrayref
  Function  : various accessor methods for holding variables during the
  run. Some will create and return the approprite object if not already holding
  it
  Returntype: as given
  Exceptions: some thrown if given the wrong type
  Example   :

=cut


sub exon_mask_list{
  my ($self, $arg) = @_;

  if($arg){
    $self->param('exon_mask_list', $arg);
  }
  return $self->param('exon_mask_list');
}

=head2 gene_mask_list

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : various, String, Int, Arrayref
  Function  : various accessor methods for holding variables during the
  run. Some will create and return the approprite object if not already holding
  it
  Returntype: as given
  Exceptions: some thrown if given the wrong type
  Example   :

=cut
sub gene_mask_list{
  my ($self, $arg) = @_;

  if($arg){
    $self->param('gene_mask_list', $arg);
  }
  return $self->param('gene_mask_list');
}


sub id_pool_bin{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('id_pool_bin', $arg);
  }
  return $self->param('id_pool_bin');
}

sub id_pool_index{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('id_pool_index', $arg);
  }
  return $self->param('id_pool_index');
}

sub use_id{
  my ($self, $id) = @_;
  if($id){
    $self->param('use_id', $id);
  }
  return $self->param('use_id');
}

sub get_dbadaptor {
  my ($self, $dba) = @_;

  return $self->hrdb_get_dba($dba);
}

sub paf_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->param('paf_source_db', $db);
  }
  if(!$self->param_is_defined('paf_source_db')){
    my $db = $self->get_dbadaptor($self->PAF_SOURCE_DB);
    if ( $db->dnadb ) {
       $db->dnadb->disconnect_when_inactive(1);
    }
    $self->param('paf_source_db', $db);
  }
  return $self->param('paf_source_db');
}

sub gene_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->param('gene_source_db', $db);
  }
  if(!$self->param_is_defined('gene_source_db')){
    if (@{$self->BIOTYPES_TO_MASK}) {
      my $db = $self->get_dbadaptor($self->GENE_SOURCE_DB);
      if ( $db->dnadb ) {
         $db->dnadb->disconnect_when_inactive(1);
      }
      $self->param('gene_source_db', $db);
    }
  }
  return $self->param('gene_source_db');
}

sub output_db{
  my ($self, $db) = @_;
  if($db){
    $self->param('output_db', $db);
  }
  if (!$self->param_is_defined('output_db')){
    my $db = $self->get_dbadaptor($self->OUTPUT_DB);
    $self->param('output_db', $db);
  }
  return $self->param('output_db');
}

sub get_adaptor{
  my ($self) = @_;
  return $self->output_db->get_GeneAdaptor;
}

sub paf_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->param('paf_slice', $slice);
  }
  if (!$self->param_is_defined('paf_slice')){
    my $slice = $self->fetch_sequence($self->quey->name, $self->paf_source_db,
                                      $self->REPEATMASKING);
    $self->param('paf_slice', $slice);
  }
  return $self->param('paf_slice');
}

sub gene_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->param('gene_slice', $slice);
  }
  if (!$self->param_is_defined('gene_slice')){
    my $slice = $self->fetch_sequence($self->quey->name, $self->gene_source_db,
                                      $self->REPEATMASKING);
    $self->param('gene_slice', $slice);
  }
  return $self->param('gene_slice');
}

sub output_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->param('output_slice', $slice);
  }
  if (!$self->param_is_defined('output_slice')){
    my $slice = $self->fetch_sequence($self->quey->name, $self->output_db);
    $self->param('output_slice', $slice);
  }
  return $self->param('output_slice');
}

sub kill_list{
  my ($self, $arg) = @_;

  if($arg){
    $self->param('killlist_db', $arg);
  }
  if (!$self->param_is_defined('killlist_db')){
    my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => 'protein');
    $self->param('killlist_db', $kill_list_object->get_kill_list);
  }
  return $self->param('killlist_db');
}

sub rejected_set{
  my ($self, $arg) = @_;

  if($arg){
    $self->param('rejected_set', $arg);
  }
  return $self->param('rejected_set');
}

sub filter_object{
  my ($self, $arg) = @_;
  if($arg){
    $self->throw("RunnableDB::BlastMiniGenewise ".
          $arg." must have a method filter_genes")
      unless($arg->can("filter_genes"));
    $self->param('filter_object', $arg);
  }
  if(!$self->param_is_defined('filter_object') && $self->FILTER_OBJECT){
    $self->require_module($self->FILTER_OBJECT);
    my %params = %{$self->FILTER_PARAMS};
    $arg = $self->FILTER_OBJECT->new(
                                     -slice => $self->query,
                                     -seqfetcher => $self->seqfetcher,
                                     %params,
                                    );
    $self->param('filter_object', $arg);
  }
  return $self->param('filter_object');
}

sub seqfetcher{
  my ($self, $arg) = @_;

  if($arg){
    $self->throw("RunnableDB::BlastMiniGenewise ".
          $arg." must have a method get_Seq_by_acc")
      unless($arg->can("get_Seq_by_acc"));
    $self->param('seqfetcher', $arg);
  }
  if(!$self->param_is_defined('seqfetcher')){
    $self->require_module($self->SEQFETCHER_OBJECT);
    my %params = %{$self->SEQFETCHER_PARAMS};
    print $params{-db}->[0], "\n";
    $arg = $self->SEQFETCHER_OBJECT->new(
                                    %params,
                                   );
    $self->param('seqfetcher', $arg);
  }
  return $self->param('seqfetcher');
}


=head2 overlaps_fiveprime_end_of_slice

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Bio::EnsEMBL::Feature
  Arg [3]   : Bio::EnsEMBL::Slice
  Function  : returns 1 if the features starts of the 5 end of the slice
  Returntype: boolean
  Exceptions:
  Example   :

=cut



sub overlaps_fiveprime_end_of_slice{
  my ($self, $feature, $slice) = @_;
  if($feature->start < 1){
    return 1;
  }
  return 0;
}

#CONFIG METHODS


=head2 PAF_LOGICNAMES

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions:
  Example   :

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise


sub PAF_LOGICNAMES{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PAF_LOGICNAMES', $arg);
  }
  return $self->param('PAF_LOGICNAMES');
}

sub PAF_MIN_SCORE_THRESHOLD {
  my ($self, $arg) = @_;
  if( defined($arg) ) {
    $self->param('PAF_MIN_SCORE_THRESHOLD', $arg);
  }
  return $self->param('PAF_MIN_SCORE_THRESHOLD');
}

sub PAF_UPPER_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PAF_UPPER_SCORE_THRESHOLD', $arg);
  }
  return $self->param('PAF_UPPER_SCORE_THRESHOLD');
}

sub PAF_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PAF_SOURCE_DB', $arg);
  }
  return $self->param('PAF_SOURCE_DB');
}

sub GENE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('GENE_SOURCE_DB', $arg);
  }
  return $self->param('GENE_SOURCE_DB');
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OUTPUT_DB', $arg);
  }
  return $self->param('OUTPUT_DB');
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OUTPUT_BIOTYPE', $arg);
  }
  return $self->param('OUTPUT_BIOTYPE');
}

sub GENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('GENEWISE_PARAMETERS', $arg);
  }
  if ($self->param_is_defined('GENEWISE_PARAMETERS')) {
    return $self->param('GENEWISE_PARAMETERS');
  }
  else {
    return {};
  }
}

sub MINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MINIGENEWISE_PARAMETERS', $arg);
  }
  if ($self->param_is_defined('MINIGENEWISE_PARAMETERS')) {
    return $self->param('MINIGENEWISE_PARAMETERS');
  }
  else {
    return {};
  }
}

sub MULTIMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MULTIMINIGENEWISE_PARAMETERS', $arg);
  }
  if ($self->param_is_defined('MULTIMINIGENEWISE_PARAMETERS')) {
    return $self->param('MULTIMINIGENEWISE_PARAMETERS');
  }
  else {
    return {};
  }
}

sub BLASTMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('BLASTMINIGENEWISE_PARAMETERS', $arg);
  }
  if ($self->param_is_defined('BLASTMINIGENEWISE_PARAMETERS')) {
    return $self->param('BLASTMINIGENEWISE_PARAMETERS');
  }
  else {
    return {};
  }
}


sub EXONERATE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('EXONERATE_PARAMETERS', $arg);
  }
  if ($self->param_is_defined('EXONERATE_PARAMETERS')) {
    return $self->param('EXONERATE_PARAMETERS');
  }
  else {
    return {};
  }
}


sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FILTER_PARAMS', $arg);
  }
  return $self->param('FILTER_PARAMS');
}

sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FILTER_OBJECT', $arg);
  }
  return $self->param('FILTER_OBJECT');
}

sub BIOTYPES_TO_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('BIOTYPES_TO_MASK', $arg);
  }
  return $self->param('BIOTYPES_TO_MASK');
}

sub EXON_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('EXON_BASED_MASKING', $arg);
  }
  return $self->param('EXON_BASED_MASKING');
}

sub GENE_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('GENE_BASED_MASKING', $arg);
  }
  return $self->param('GENE_BASED_MASKING');
}


sub POST_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('POST_GENEWISE_MASK', $arg);
  }
  return $self->param('POST_GENEWISE_MASK');
}

sub PRE_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('PRE_GENEWISE_MASK', $arg);
  }
  return $self->param('PRE_GENEWISE_MASK');
}

sub REPEATMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('REPEATMASKING', $arg);
  }
  return $self->param('REPEATMASKING');
}

sub SEQFETCHER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('SEQFETCHER_OBJECT', $arg);
  }
  return $self->param('SEQFETCHER_OBJECT');
}

sub SEQFETCHER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('SEQFETCHER_PARAMS', $arg);
  }
  return $self->param('SEQFETCHER_PARAMS');
}

sub USE_KILL_LIST{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('USE_KILL_LIST', $arg);
  }
  return $self->param('USE_KILL_LIST');
}

sub LIMIT_TO_FEATURE_RANGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('LIMIT_TO_FEATURE_RANGE', $arg);
  }
  return $self->param('LIMIT_TO_FEATURE_RANGE');
}


sub FEATURE_RANGE_PADDING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FEATURE_RANGE_PADDING', $arg);
  }
  return $self->param('FEATURE_RANGE_PADDING');
}

sub WRITE_REJECTED{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->param('WRITE_REJECTED', $arg);
  }
  return $self->param('WRITE_REJECTED');
}

sub REJECTED_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('REJECTED_BIOTYPE', $arg);
  }
  return $self->param('REJECTED_BIOTYPE');
}

sub SOFTMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('SOFTMASKING', $arg);
  }
  return $self->param('SOFTMASKING');
}


sub MAKE_SIMGW_INPUT_ID_PARMAMS {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MAKE_SIMGW_INPUT_ID_PARMAMS', $arg);
  }
  return $self->param('MAKE_SIMGW_INPUT_ID_PARMAMS');
}

=head2 group_genes_by_id

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Genes
  Function  : groups the genes on the basis of transcript supporting features
  Returntype: hashref
  Exceptions: throws if a transcript has no supporting features
  Example   :

=cut


sub group_genes_by_id{
  my ($self, $genes) = @_;
  my %hash;
  foreach my $gene(@$genes){
    my $transcripts = $gene->get_all_Transcripts;
    my $transcript = $transcripts->[0];
    my $sfs = $transcript->get_all_supporting_features;
    $self->throw(id($transcript)." appears to have no supporting features") if(!@$sfs);
    my $id = $sfs->[0]->hseqname;
    if(!$hash{$id}){
      $hash{$id} = 1;
    }else{
      $hash{$id}++;
    }
  }
  return \%hash;
}

sub OPTIMAL_LENGTH {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OPTIMAL_LENGTH', $arg);
  }
  return $self->param('OPTIMAL_LENGTH');
}

1;

