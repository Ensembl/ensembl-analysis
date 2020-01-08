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

Bio::EnsEMBL::Analysis::RunnableDB::MakeSimilarityInputIDs - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::MakeSimilarityInputIDs;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise qw(GENEWISE_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (parse_config);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($output_logicname, $protein_count, $bmg_logicname, $optimal_length) = rearrange
    (['OUTPUT_LOGICNAME', 'PROTEIN_COUNT', 'BMG_LOGICNAME', 'OPTIMAL_LENGTH'], @args); 
  
  # The constructor arg like output_logicname, protein_count and bmg_logicname 
  # can be retrieved from analysis.parameters column in compara-style-format, ie.
  # " -protein_count => 10 , -output_logicname => bla, -bmg_logicname => bla
  #

   my $create_analysis =  ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{creation_analysis};
   my $s_regex = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{submission_logic_name};
   my $bmg_regex = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{bmg_logic_name};
   my $def_protein_count = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{protein_count} ; 
   my $max_padding = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{max_padding} ; 
   my $max_intron_size = ${$GENEWISE_CONFIG_BY_LOGIC}{DEFAULT}{MAKE_SIMGW_INPUT_ID_PARMAMS}{max_intron_size} ; 

   if ( $self->analysis->logic_name =~m/$create_analysis/i ) {     

     my $submission_logic_name = $self->analysis->logic_name;  # Create_simgw_vert_other 
     $submission_logic_name =~m/$s_regex/i;                    # 'Create_' 
     $submission_logic_name = $1 ;                             # simgw_vert_other 
     $submission_logic_name="Submit_".$submission_logic_name ; # Submission will be : Submit_simgw_vert_other  

     my $bmg  = $self->analysis->logic_name;   # Create_simgw_vert_other 
     $bmg =~m/$bmg_regex/i;                    # simgw_vert_other 
     $bmg = $1 ;                    
   
     # 
     # Create_simgw_vert_other is running MakeSimilarityInputIDS and will create input_ids for analysis Submit_simgw_vert_mammal 
     # thes will then run simgw_vert_mammal 
     # 
     
     print "BlastMiniGenewise - logic name : $bmg\n" ; 
     print "SUBMISSION ANAL   - logic name : $submission_logic_name\n" ;  

     $self->bmg_logicname($bmg);
     $self->output_logicname($submission_logic_name);  
 
     if ($def_protein_count) { 
       $self->protein_count($def_protein_count); 
     }
  } 

 
  ### Defaults are over-ridden by parameters given in analysis table...
  my $ph = $self->parameters_hash;
  $self->protein_count($ph->{-protein_count}) if(exists $ph->{-protein_count});
  $self->output_logicname($ph->{-output_logicname}) if(exists $ph->{-output_logicname});
  $self->bmg_logicname($ph->{-bmg_logicname}) if(exists $ph->{-bmg_logicname});
  $self->max_padding($ph->{-max_padding}) if(exists $ph->{-max_padding});
  $self->max_intron_size($ph->{-max_intron_size}) if(exists $ph->{-max_intron_size});

  ### ...which are over-ridden by constructor arguments. 
  $self->protein_count($protein_count);
  $self->output_logicname($output_logicname);
  $self->bmg_logicname($bmg_logicname);
  $self->max_padding($max_padding);
  $self->max_intron_size($max_intron_size);

  throw("Need an output logicname ".$self->output_logicname." and a bmg logicname ".
        $self->bmg_logicname." defined") if(!$self->output_logicname ||
                                            !$self->bmg_logicname);

  parse_config($self, $GENEWISE_CONFIG_BY_LOGIC, $self->bmg_logicname);
  $self->optimal_length($self->OPTIMAL_LENGTH);

  return $self;
}


sub fetch_input{
  my ($self) = @_;

  print "\n\n***Fetching sequence from ".$self->db->dbname."***\n\n";
  $self->query($self->fetch_sequence($self->input_id, $self->db, $self->REPEATMASKING)); 

  
  my $dba_pip = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( 
                                                      -host => $self->db->host , 
                                                      -dbname => $self->db->dbname , 
                                                      -user => $self->db->username , 
                                                      -pass => $self->db->password , 
                                                      -port => $self->db->port ,  
                                                       );      
  $self->pipeline_adaptor($dba_pip); 
  my $output_analysis = $dba_pip->get_AnalysisAdaptor->fetch_by_logic_name($self->output_logicname); 

  if(!$output_analysis) { 
    throw("Failed to find an analysis with the logic_name ". $self->output_logicname." Cannot continue") ; 
  }
  $self->output_analysis($output_analysis); 
}

sub run{
  my ($self) = @_;

  my %kill_list =  %{$self->kill_list} if($self->USE_KILL_LIST);
  my @mask_exons;
  my @a_iids;

  # remove masked and killed hits as will be done in the build itself
  foreach my $type ( @{$self->BIOTYPES_TO_MASK} ) {  
    print "\nmasking Gene-type : $type\n\n" ; 
    foreach my $mask_genes ( @{ $self->gene_slice->get_all_Genes_by_type($type) } ) {
      foreach my $mask_exon ( @{ $mask_genes->get_all_Exons } ) {
        if ( $mask_exon->seqname eq $self->gene_slice->id ) {
          push @mask_exons, $mask_exon;
        }
      }
    }
  }

  # make the mask list non-redundant. Much faster when checking against features
  my @mask_regions;
  foreach my $mask_exon ( sort { $a->start <=> $b->start } @mask_exons ) {
    if ( @mask_regions and $mask_regions[-1]->{'end'} > $mask_exon->start ) {
      if ( $mask_exon->end > $mask_regions[-1]->{'end'} ) {
        $mask_regions[-1]->{'end'} = $mask_exon->end;
      }
    } else {
      push @mask_regions, { start => $mask_exon->start, end => $mask_exon->end }

    }
  }  

  my $num_seeds = 0;  
  print STDERR 'masked : ' , scalar(@mask_regions), "\n";
  # If the length of the query is smaller than the length wanted, we use the "old" method
  if ($self->optimal_length >= $self->query->length) {

      foreach my $logicname(@{$self->PAF_LOGICNAMES}) {
          print 'FETCHING FEATURES FOR :'.$logicname."\n";
          my %features = %{$self->get_good_features($logicname, \@mask_regions, \%kill_list)};
          $num_seeds += scalar( keys %features );
      }

      # rule of thumb; split data so that each job constitutes one piece of
      # genomic DNA against ~20 proteins.
      #
      return () if($num_seeds == 0);
      @a_iids = @{$self->create_input_ids(undef, $num_seeds, $self->query->name)};
  }
  else {
      my @a_genes;

      foreach my $logicname(@{$self->PAF_LOGICNAMES}) {
          print 'FETCHING FEATURES FOR :'.$logicname."\n";
          my %features = %{$self->get_good_features($logicname, \@mask_regions, \%kill_list)};
          push @a_genes, @{$self->features2genes(\%features)};
      }

      return () unless (scalar(@a_genes));
      if (@a_genes > $self->protein_count) {
          my %h_types;
          $h_types{'good'} = ['protein_coding'];
          # We cluster without taking the strand into account and without checking the exons
          my ($ra_clustered, $ra_unclustered) = cluster_Genes(\@a_genes, \%h_types, 0, 1, 1);
          my @a_sorted_genes = sort {$a->end <=> $b->end} (@$ra_clustered, @$ra_unclustered);
          my $optimal_length = $self->optimal_length;
          my $max_padding = $self->max_padding;
          my $range_start = $self->query->start;
          my @a_input_genes;

          my $last_seq_region_start = $self->query->end-$optimal_length-1;
          for (my $index = $optimal_length+$self->query->start; $index < $last_seq_region_start; $index += $optimal_length) {
              while (my $gene = shift @a_sorted_genes) {
                  if ($gene->seq_region_end <= $index) {
                      push @a_input_genes, $gene;
                  }
                  else {
                      my $range_end;
                      if (scalar(@a_input_genes) == 0) {
                          if ($gene->seq_region_start < $index) {
                              $range_end = $gene->seq_region_end+$max_padding;
                          }
                          else {
                              @a_input_genes = ($gene);
                              $range_start = $gene->seq_region_start-$max_padding+1;
                              $range_start = 1 if ($range_start < 1);
                              $num_seeds = 0;
                              last;
                          }
                      }
                      else {
                          $range_end = $a_input_genes[$#a_input_genes]->seq_region_end+$max_padding;
                          my ($query_name, $strand) = $self->query->name =~ /(.*):\d+:\d+:(.+)$/;
                          $query_name .= ':'.$range_start.':'.$range_end.':'.$strand;
                          push(@a_iids, @{$self->create_input_ids(\@a_input_genes, $num_seeds, $query_name)});
                      }
                      @a_input_genes = ($gene);
                      $range_start = $gene->seq_region_start-$max_padding+1;
                      $range_start = 1 if ($range_start < 1);
                      $num_seeds = 0;
                      if ($gene->seq_region_start > $index) {
                          $index = int($gene->seq_region_start/$optimal_length)*$optimal_length;
                      }
                      last;
                  }
              }
          }
          my ($query_name, $range_end, $strand) = $self->query->name =~ /(.*):\d+:(\d+):(.+)$/;
          my $array_index = $#a_input_genes > 0 ? $#a_input_genes : 0;
          $range_end = $a_input_genes[$array_index]->seq_region_end+$max_padding if ($a_input_genes[$array_index] and $range_end > $a_input_genes[$array_index]->seq_region_end+$max_padding);
          $query_name .= ':'.$range_start.':'.$range_end.':'.$strand;
          push(@a_iids, @{$self->create_input_ids(\@a_sorted_genes, scalar(@a_sorted_genes), $query_name)});
      }
      else {
          @a_iids = @{$self->create_input_ids(undef, scalar(@a_genes), $self->query->name)};
      }
  }
  $self->output(\@a_iids);
}

sub create_input_ids {
      my ($self, $ra_input_genes, $num_seeds, $query_name) = @_;

      my @a_iids;
      print STDOUT "Inside\n";

      foreach my $gene (@{$ra_input_genes}) {
          $num_seeds += $gene->get_Gene_Count;
          foreach my $g (@{$gene->get_Genes}) {
              print 'cluster: ', $g->display_id,"\t",'start: ',$g->start."\t".'end: ',$g->end,"\n";
          }
      }

      my $num_chunks = int( $num_seeds / $self->protein_count ) + 1;

      for ( my $x = 1 ; $x <= $num_chunks ; $x++ ) {

          #generate input id : $chr_name.1-$chr_length:$num_chunks:$x
          my $new_iid = $query_name . ":$num_chunks:$x";
          push @a_iids, $new_iid;
      }

      print "HAVE ".@a_iids." to write to the ref database\n"; 

      for ( @a_iids ) {  
          print "$_\n" ; 
      } 

      return \@a_iids;
 }

sub get_good_features {
    my ($self, $logicname, $ra_mask_regions, $kill_list) = @_;

    my @a_features = @{$self->paf_slice->get_all_ProteinAlignFeatures($logicname, $self->PAF_MIN_SCORE_THRESHOLD-1)};
    print 'HAVE '.@a_features.' protein-align-features for ',$logicname,"\n";
    my %h_features;

    foreach my $f(@a_features){
#        next if($self->PAF_MIN_SCORE_THRESHOLD && $f->score < $self->PAF_MIN_SCORE_THRESHOLD);
        next if($self->PAF_UPPER_SCORE_THRESHOLD && $f->score > $self->PAF_UPPER_SCORE_THRESHOLD);
        push(@{$h_features{$f->hseqname}}, $f);
    }

    my @a_ids_to_ignore;
    if ($ra_mask_regions) {
#        print STDERR "entering masking\n";

SEQID:  foreach my $sid ( keys %h_features ) {
            my $ex_idx = 0;
            print STDERR "Looking at $sid\n";
FEAT:       foreach my $f ( sort { $a->start <=> $b->start } @{ $h_features{$sid} } ) {
#               printf STDERR "Feature: %d %d\n", $f->start, $f->end;
                for ( ; $ex_idx < @{$ra_mask_regions} ; ) {
                    my $mask_exon = $ra_mask_regions->[$ex_idx];
#                   printf STDERR " Mask exon %d %d\n", $mask_exon->{'start'}, $mask_exon->{'end'};
                    if ( $mask_exon->{'start'} > $f->end ) {
                        print "no exons will overlap this feature \n" ;  
                        # no exons will overlap this feature
                        next FEAT;
                    }
                    elsif ( $mask_exon->{'end'} >= $f->start ) {
                        # overlap
                        push @a_ids_to_ignore, $f->hseqname;
                        print 'Ignoring : ' . $f->hseqname . "\n" ;  
                        next SEQID;
                    }
                    else {
                        $ex_idx++;
                    }
                }
            }
        } 

    } 
    print 'Ignoring '.@a_ids_to_ignore." features\n";

    foreach my $dud_id ( @a_ids_to_ignore, keys %$kill_list ) {
        if ( exists $h_features{$dud_id} ) {
            delete $h_features{$dud_id};
        }
    }

    return \%h_features;
}

sub features2genes {
    my ($self, $rh_features) = @_;

    my @a_genes;
    my $transcripts_count = 1;

    foreach my $feature (keys %$rh_features) {
        my @plus = grep {$_->strand == 1} @{$rh_features->{$feature}};
        my @minus = grep {$_->strand == -1} @{$rh_features->{$feature}};
        if (@plus) {
            push @a_genes, @{$self->create_gene_from_feature(\@plus)};
        }
        if (@minus) {
            push @a_genes, @{$self->create_gene_from_feature(\@minus)};
        }
        ++$transcripts_count;
    }

    undef $rh_features;

    return \@a_genes;
}

sub create_gene_from_feature {
    my ($self, $r_features) = @_;

    my @a_sorted_features = sort { $a->end <=> $b->end }  @{$r_features};
    my $end = $a_sorted_features[$#a_sorted_features]->end;
    @a_sorted_features = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end }  @{$r_features};
    my $start = $a_sorted_features[0]->start;
    my @a_exons;
    my $last_end = $start-1;

    my $index = 0;
    foreach my $feature (@a_sorted_features) {
        my $exon = Bio::EnsEMBL::Exon->new(
            -START  => $feature->start,
            -END    => $feature->end,
            -STRAND => $feature->strand,
            -SLICE  => $self->paf_slice,
            -PHASE  => 0,
            -END_PHASE => 0,
            );
        if($last_end < $exon->start) {
            ++$index if ($exon->start-$last_end > $self->max_intron_size);
            push(@{$a_exons[$index]}, $exon);
        }
        elsif ($last_end < $exon->end) {
            $a_exons[$index][-1]->end($exon->end);
        }
        elsif ($last_end >= $exon->end) {
        }
        else {
            warning('This shouldn\'t happen');
        }
        $last_end = $a_exons[$index][-1]->end;
    }

    my @a_genes;
    foreach my $ra_exons (@a_exons) {
        my $transcript = Bio::EnsEMBL::Transcript->new(
                -START  => $start,
                -END    => $end,
                -STRAND => $a_sorted_features[0]->strand,
                -SLICE  => $self->paf_slice,
                -EXONS  => $ra_exons,
                );
        my $gene = Bio::EnsEMBL::Gene->new(
                -START => $start,
                -END => $end,
                -SLICE => $self->paf_slice,
                -DBID => $a_sorted_features[0]->hseqname,
                );
        $gene->add_Transcript($transcript);
        push (@a_genes, $gene);
    }

    return \@a_genes;
}



sub write_output{
  my ($self) = @_;  

  my $sic = $self->pipeline_adaptor->get_StateInfoContainer;  

#  my $dba_pip = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( 
#                                                      -host => $self->db->host , 
#                                                      -dbname => $self->db->dbname , 
#                                                      -user => $self->db->username , 
#                                                      -pass => $self->db->password , 
#                                                      -port => $self->db->port ,  
#                                                       );    
#  my $sic = $dba_pip->get_StateInfoContainer; 

  print "output analysis is : " . $self->output_analysis."\n" ;
  foreach my $iid(@{$self->output}){ 
    print "try to store input_id : $iid\n" ; 
    eval{
      $sic->store_input_id_analysis($iid, 
                                    $self->output_analysis, 
                                    '');
    };
    throw("Failed to store ".$iid." $@") if($@);
    logger_info("Stored ".$iid);
  }
}


sub kill_list{
  my ($self, $arg) = @_;

  if($arg){
    $self->{kill_list} = $arg;
  }
  if(!$self->{kill_list}){
    my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => 'protein');
    $self->{kill_list} = $kill_list_object->get_kill_list;
  }
  return $self->{kill_list};
}


sub output_analysis{
  my ($self, $value) = @_;
  if($value && 
     $value->isa('Bio::EnsEMBL::Pipeline::Analysis')){
    $self->{output_analysis} = $value;
  }
  return $self->{'output_analysis'};
} 

sub pipeline_adaptor{
  my ($self, $value) = @_;
  if($value && 
     $value->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
    $self->{pipeline_adaptor} = $value;
  }
  return $self->{'pipeline_adaptor'};
}

sub protein_count{
  my ($self, $value) = @_;
  if($value){
    $self->{'protein_count'} = $value;
  }
  return $self->{'protein_count'};
}

sub output_logicname{
  my ($self, $value) = @_;
  if($value){
    $self->{'output_logicname'} = $value;
  }
  return $self->{'output_logicname'};
}

sub optimal_length{
  my ($self, $value) = @_;
  if($value){
    $self->{'optimal_length'} = $value;
  }
  return $self->{'optimal_length'};
}

sub max_padding{
  my ($self, $value) = @_;
  if($value){
    $self->{'max_padding'} = $value;
  }
  return $self->{'max_padding'};
}

sub max_intron_size{
  my ($self, $value) = @_;
  if($value){
    $self->{'max_intron_size'} = $value;
  }
  return $self->{'max_intron_size'};
}

sub bmg_logicname{
  my ($self, $value) = @_;
  if($value){
    $self->{'bmg_logicname'} = $value;
  }
  return $self->{'bmg_logicname'};
}

sub paf_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->{paf_slice} = $slice;
  }
  if(!$self->{paf_slice}){
    my $slice = $self->fetch_sequence($self->input_id, $self->paf_source_db, $self->REPEATMASKING);
    $self->{paf_slice} = $slice;
  }
  return $self->{paf_slice};
}

sub gene_slice{
  my ($self, $slice) = @_;

  if($slice){
    $self->{gene_slice} = $slice;
  }
  if(!$self->{gene_slice}){
    my $slice = $self->fetch_sequence($self->input_id, $self->gene_source_db, $self->REPEATMASKING);
    $self->{gene_slice} = $slice;
  }
  return $self->{gene_slice};
}


sub paf_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->{paf_source_db} = $db;
  }
  if(!$self->{paf_source_db}){
    my $db = $self->get_dbadaptor($self->PAF_SOURCE_DB);
    $self->{paf_source_db} = $db;
  }
  return $self->{paf_source_db};
}

sub gene_source_db{
  my ($self, $db) = @_;
  if($db){
    $self->{gene_source_db} = $db;
  }
  if(!$self->{gene_source_db}){
    my $db = $self->get_dbadaptor($self->GENE_SOURCE_DB);
    $self->{gene_source_db} = $db;
  }
  return $self->{gene_source_db};
}


=head2 PAF_LOGICNAMES

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: 
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise


sub PAF_LOGICNAMES{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_LOGICNAMES} = $arg;
  }
  return $self->{PAF_LOGICNAMES}
}

sub PAF_MIN_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_MIN_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_MIN_SCORE_THRESHOLD}
}

sub PAF_UPPER_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_UPPER_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_UPPER_SCORE_THRESHOLD}
}



sub PAF_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_SOURCE_DB} = $arg;
  }
  return $self->{PAF_SOURCE_DB}
}

sub GENE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_SOURCE_DB} = $arg;
  }
  return $self->{GENE_SOURCE_DB}
}




sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_DB} = $arg;
  }
  return $self->{OUTPUT_DB}
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_BIOTYPE} = $arg;
  }
  return $self->{OUTPUT_BIOTYPE}
}

sub GENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENEWISE_PARAMETERS} = $arg;
  }
  return $self->{GENEWISE_PARAMETERS}
}

sub MINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MINIGENEWISE_PARAMETERS}
}

sub MULTIMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MULTIMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MULTIMINIGENEWISE_PARAMETERS}
}

sub BLASTMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BLASTMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{BLASTMINIGENEWISE_PARAMETERS}
}



sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_PARAMETERS} = $arg;
  }
  return $self->{FILTER_PARAMETERS}
}



sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_OBJECT} = $arg;
  }
  return $self->{FILTER_OBJECT}
}


sub BIOTYPES_TO_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BIOTYPES_TO_MASK} = $arg;
  }
  return $self->{BIOTYPES_TO_MASK}
}


sub EXON_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXON_BASED_MASKING} = $arg;
  }
  return $self->{EXON_BASED_MASKING}
}


sub GENE_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_BASED_MASKING} = $arg;
  }
  return $self->{GENE_BASED_MASKING}
}


sub POST_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{POST_GENEWISE_MASK} = $arg;
  }
  return $self->{POST_GENEWISE_MASK}
}

sub PRE_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PRE_GENEWISE_MASK} = $arg;
  }
  return $self->{PRE_GENEWISE_MASK}
}

sub REPEATMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REPEATMASKING} = $arg;
  }
  return $self->{REPEATMASKING}
}

sub SEQFETCHER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_OBJECT} = $arg;
  }
  return $self->{SEQFETCHER_OBJECT}
}

sub SEQFETCHER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_PARAMS} = $arg;
  }
  return $self->{SEQFETCHER_PARAMS}
}

sub USE_KILL_LIST{
  my ($self, $arg) = @_;
  if($arg){
    $self->{USE_KILL_LIST} = $arg;
  }
  return $self->{USE_KILL_LIST}
}


sub LIMIT_TO_FEATURE_RANGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{LIMIT_TO_FEATURE_RANGE} = $arg;
  }
  return $self->{LIMIT_TO_FEATURE_RANGE}
}


sub FEATURE_RANGE_PADDING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FEATURE_RANGE_PADDING} = $arg;
  }
  return $self->{FEATURE_RANGE_PADDING}
}

sub WRITE_REJECTED{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{WRITE_REJECTED} = $arg;
  }
  return $self->{WRITE_REJECTED};
}

sub REJECTED_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REJECTED_BIOTYPE} = $arg;
  }
  return $self->{REJECTED_BIOTYPE};
}

sub SOFTMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SOFTMASKING} = $arg;
  }
  return $self->{SOFTMASKING}
}

sub EXONERATE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXONERATE_PARAMETERS} = $arg;
  }
  return $self->{EXONERATE_PARAMETERS}
}

sub MAKE_SIMGW_INPUT_ID_PARMAMS { 
  my ($self, $arg) = @_;
  if($arg){
    $self->{MAKE_SIMGW_INPUT_ID_PARMAMS} = $arg;
  }
  return $self->{MAKE_SIMGW_INPUT_ID_PARMAMS}
}   

sub OPTIMAL_LENGTH { 
  my ($self, $arg) = @_;
  if($arg){
    $self->{OPTIMAL_LENGTH} = $arg;
  }
  return $self->{OPTIMAL_LENGTH}
}   
1;
