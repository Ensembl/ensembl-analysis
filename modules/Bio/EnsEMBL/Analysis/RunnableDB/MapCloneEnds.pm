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

Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds - 

=head1 SYNOPSIS

my $clonemap = 
  Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
    -database   => $EST_GENOMIC,
  );

$clonemap->fetch_input();
$clonemap->run();
$clonemap->write_output(); #writes to DB

=head1 DESCRIPTION

This object maps clone sequences to a genome by using the
exonerate alignment program,and write the results as Dna Align Features.
It relies on the following modules:
Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds;
Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds;

In the analysis table, you have to add 3 analysis, one that calls
MapCloneEnds, and two for the different exonerate steps:
MapCloneEnds (module -> MapCloneEnds),
EXONERATE_CLONE_ENDS (module -> MapCloneEnds),
REFINE_CLONE_ENDS(module -> MapCloneEnds)

MapCloneEnds reads the input_id_chunks from a file where each line 
contains a number of IDs separated by ":". In order to be able to correctly
parse the input ids, they should be given in a very specific format. The file 
can be automaticaly generated from an XML using chunker.pl or manually generated 
where each id should be in the format:
Clone_ID,length_of_clone,standard_deviation_of_length,Clone_end_ID,Direction_of_Clone
(i.e. CH243-307D10,184000.0,36800.0,1098421033278,F).

=head1 METHODS

=cut

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds;
use Bio::EnsEMBL::Analysis::Config::ExonerateCloneEnds;


use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

##########################################################################

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->{_refined_results}=[] ;
  $self->read_and_check_config($CLONE_CONFIG);
  return $self;
}

##########################################################################

sub fetch_input {
 
# It does nothing so no really code implementation here.
  my ($self) = @_;

}

##########################################################################

sub run {
  my ($self) = @_;

  my $trace_name = '';  # CloneEnd id as in database
  my $clone_id = '';# Clone id
  my $insert_size = '';# Length of the inserted sequence in clone
  my $insert_stdev = '';# Standard deviation of the clone insert. 
  # The two insert values will be used in the filter process to check quality of alignments

  # You need to run  perl ensembl-personal/jb16/scripts/general/chunker.pl in order to create the 
  # input_analysis_ids and create this file. To change the path to the list of chunks, you should also do it in
  # the chunker script
 
  my $chunks_list = $self->CHUNKSLIST;

  my %clones  = (); # hash where each clone object will be stored

  my %cloneEnd_ids = ();

  # Open the file with the list of seq_ids per chunk and load it into an array to be used by fetch_input
  open (INFILE,"<$chunks_list");
  my @listOfIDs = ();

  while (my $line = <INFILE>){

    chomp($line);
 
    # Split by ":" to get information of each clone as an independent element in an array. 
    my @clone = split (/:/,$line);
    my $listOfClones = '';

    foreach my $clone (@clone){

      # Split all the information of one cloneEnd and load each element into an array
      my @clone_data = split (/,/,$clone);

      if (!$clones{$clone_data[0]}){
        $clones{$clone_data[0]}=[];
      }
    
      # Group information of each cloneEnd using clone_id as reference (will be used in the filter step)  
      push (@{$clones{$clone_data[0]}}, $clone_data[3]); # Add cloneEnd_name
      push (@{$clones{$clone_data[0]}}, $clone_data[1]); # Add insert length
      push (@{$clones{$clone_data[0]}}, $clone_data[2]); # Add insert length standard deviation
      push (@{$clones{$clone_data[0]}}, $clone_data[4]); # Add cloneEnd direction (R or F)
      $cloneEnd_ids{$clone_data[3]} = $clone_data[0];

      if ($listOfClones eq ''){
	$listOfClones.= $clone_data[3];
      }else{
        $listOfClones.= ":".$clone_data[3];
      }
    }

    # Creates the array that contains the list of cloneEnd_ids that will go in each chunks
    push (@listOfIDs,$listOfClones);
  }

  close INFILE;

  # Run ExonerateClones to get a first idea on the aligned/location of the clone Ends 
    my $exonerate = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds->new(
                    -DB          => $self->db,
                    -INPUT_ID    => $self->input_id,
                    -ANALYSIS    => $self->fetch_analysis("EXONERATE_CLONE_ENDS"),
                    );

    $exonerate->fetch_input(\@listOfIDs);

    $exonerate->run();
  
  # Run exonerate and get the output.
  my $clone_alignments = $exonerate->output();

  my @selected_alignments = @{$self->filter_alignments($clone_alignments, \%clones, \%cloneEnd_ids)};

  foreach my $selected_alignment(@selected_alignments){
    
    if ($selected_alignment ne ''){
    
      my $clone_id = $selected_alignment->hseqname;
      my $chr_id   = $selected_alignment->seqname;
      my $start    = ($selected_alignment->start)-2000;
      my $end      = ($selected_alignment->end)+2000;

      my @chr_name = split (/:/, $chr_id);
  
      # Rebuild the input_id to send the coordinates for the target slice and the query sequence.
      my $input_id = $chr_name[0].":".$chr_name[1].":".$chr_name[2].":".$start.":".$end.":".$chr_name[5].":".$clone_id;

      my $refine = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds->new(
                   -DB          => $self->db,
                   -INPUT_ID    => $input_id,
                   -ANALYSIS    => $self->fetch_analysis("REFINE_CLONE_ENDS"),
                   );
    
      $refine ->fetch_input();
      $refine ->run();
      $self->refined_results($refine);
    }
  }
}

##########################################################################

sub write_output {

  my ( $self, @output ) = @_;
   
  foreach my $refine_object( @{ $self->refined_results }) {
      
      $refine_object->write_output();
  }
}


############################################################################################

sub fetch_analysis{
  my ($self, $logic_name) = @_;
 
  my  $db = $self->db;

  my $analysis_adaptor= $db->get_AnalysisAdaptor;

  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);
  
  return $analysis;
}

sub refined_results{
  my ($self, $refine_result) = @_;

 # if(!$self->{'refined_results'}){
 #   $self->{'refined_results'}=[];
 # }
 
  if ($refine_result) { 
   push @{$self->{_refined_results}}, $refine_result ; 
  }  
  
  return $self->{_refined_results};

}

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($CLONE_CONFIG);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (
    qw(
      CHUNKSLIST
    )
  ){ 
    if ( not defined $self->$config_var ){
      throw("You must define $config_var in config for logic '$logic'");
    }
  }
}

sub QUERYSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYSEQS'} ) ) {
    return $self->{'_CONFIG_QUERYSEQS'};
  } else {
    return undef;
  }
}

sub QUERYTYPE {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYTYPE'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYTYPE'} ) ) {
    return $self->{'_CONFIG_QUERYTYPE'};
  } else {
    return undef;
  }
}

sub GENOMICSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_GENOMICSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_GENOMICSEQS'} ) ) {
    return $self->{'_CONFIG_GENOMICSEQS'};
  } else {
    return undef;
  }
}

sub IIDREGEXP {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_IIDREGEXP'} = $value;
  }

  if ( exists( $self->{'_CONFIG_IIDREGEXP'} ) ) {
    return $self->{'_CONFIG_IIDREGEXP'};
  } else {
    return undef;
  }
}

sub OUTDB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OUTDB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OUTDB'} ) ) {
    return $self->{'_CONFIG_OUTDB'};
  } else {
    return undef;
  }
}

sub DNADB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_DNADB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_DNADB'} ) ) {
    return $self->{'_CONFIG_DNADB'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OPTIONS'} ) ) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}

sub CHUNKSLIST {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_CHUNKSLIST'} = $value;
  }

  if ( exists( $self->{'_CONFIG_CHUNKSLIST'} ) ) {
    return $self->{'_CONFIG_CHUNKSLIST'};
  } else {
    return undef;
  }
}

sub filter_alignments{

  my ($self, $clone_alignments, $clones_ref, $cloneEnd_ids_ref ) = @_;

  my %clone_cluster = ();
  my %aligned_clones = (); 
  my %clones = %{$clones_ref};
  my %cloneEnd_ids = %{$cloneEnd_ids_ref};
  my @selected_alignments;

  # Group alignments by clone name
  foreach my $clone_alignment(@{$clone_alignments}){
   
    # Get the cloneEnds name
    my $clone_name = $clone_alignment->hseqname();

    # Get the clone id to which the cloneEnds belongs
    my $complete_clone_name = $cloneEnd_ids{$clone_name};

    # Store the clone name in a hash, so only selected clones will be checked in the next step.
    # the object stored in the hash is not important, what we really want to get is a list of 
    # unique clones(complete_clone_name) with no duplicate entries
    $aligned_clones{$complete_clone_name} = $clone_name;


    if(!$clone_cluster{$clone_name}){
      $clone_cluster{$clone_name} = [];
    }
    push (@{$clone_cluster{$clone_name}}, $clone_alignment);
  }


  foreach my $pair (keys %aligned_clones){
    
    # Check if the clone has two cloneEnds
    if ((scalar @{$clones{$pair}})/4 > 1){

      # Check if at least one of the cluster pair alignments is selected.
      # in case none is selected it will send cloneEnds to the single_filter
      my $cluster_selected = 0;

      my %clean_cloneEnds = ();
     
      my %status = ();

      my %clone_dir = ();

      # Get the length of the clone and the name of the two cloneEnds
      my $clone_length = $clones{$pair}[1]+$clones{$pair}[2]+200000;
      for (my $numberOfEnd=0;$numberOfEnd < scalar @{$clones{$pair}};$numberOfEnd+=4){
       
        my $cloneEnd =$clones{$pair}[$numberOfEnd];

        # Clone_dir used to avoid pairing of two cloneEnds in the same end as for some clones
        # there are more than two cloneEnds where more than one correspond to the same end that
        # was sequenced more than once.
        $clone_dir{$numberOfEnd} = $clones{$pair}[$numberOfEnd+3];
	
        my @cloneEnd_clean  = @{$self->clean_clusters($clone_cluster{$cloneEnd})};
	if (!$clean_cloneEnds{$numberOfEnd}){
          $clean_cloneEnds{$numberOfEnd} = \@cloneEnd_clean;
        }
       
        # Set initial status of cloneEnd to 0 which means that none of its alignments pair with the other cloneEnd
        if (!$status{$numberOfEnd}){
	  $status{$numberOfEnd} = 0;
        }
      }

      foreach my $clean_cloneEnd1 (keys %clean_cloneEnds){

        foreach my $clean_cloneEnd2 (keys %clean_cloneEnds){
          # first check that one cloneEnd is not compared with itself. Then check that if we have pair
          # A->B don't use again B->A and finally get only pairs of F + R cloneEnds
      	  if ($clean_cloneEnd1 != $clean_cloneEnd2 && $clean_cloneEnd1 < $clean_cloneEnd2 
              && $clone_dir{$clean_cloneEnd1} ne $clone_dir{$clean_cloneEnd2}){

            # Use this variable to check where an alignment occur and avoid duplicate alignments because of
            # short consecutive sequences that align near by and are both selected for the next exonerate step
	    my $first_prev_chr_start = 0;
            my $first_chr_hit = 'Null';

            foreach my $fa (@{$clean_cloneEnds{$clean_cloneEnd1}}){

              # Check if two alignments don't belong to a cluster or they align in different chromosomes.
              # This is made to avoid getting duplicate alignments in the next exonerate step when very close
              # coordinates are selected for the target sequence
	      if (($fa->seqname() eq $first_chr_hit && ($fa->start()-$first_prev_chr_start) > 4000)
                   || $fa->seqname() ne $first_chr_hit){

	        my $second_prev_chr_start = 0;

                foreach my $sa (@{$clean_cloneEnds{$clean_cloneEnd2}}){
	
                  if ($fa->seqname() eq $sa->seqname()){

                    my $diff = ($fa->start())-($sa->end());
                    my $abs_diff = abs($diff);

                    # Check that the two cloneEnds are close enough as to be considered as paired and that
                    # the second cloneEnd don't belong to a cluster that was previously selected
                    # This is made to avoid the same cloneEnds to be paired more than once when there is more 
                    # than one alignment in a short region.
                    if ($abs_diff <= $clone_length  && (($sa->start()-$second_prev_chr_start)>4000)){
             
                      # In case two alignments are selected store them for the next exonerate step
                      push (@selected_alignments, $fa);
                      push (@selected_alignments, $sa); 
                    
                      $status{$clean_cloneEnd1} = 1;
                      $status{$clean_cloneEnd2} = 1;
                      $cluster_selected = 1;
                      $first_prev_chr_start = $fa->start();
                      $second_prev_chr_start = $sa->start();
                      $first_chr_hit= $fa->seqname();
     	            }
	          }
	        }
              }
            }
          }
        }
      }

      foreach my $clone_status (keys %status){
	if ($status{$clone_status} == 0){

          my $cloneEnd =$clones{$pair}[$clone_status];
          my $selected_align = $self->single_filter_alignments($clone_cluster{$cloneEnd});
          push (@selected_alignments, ${$selected_align});
        
        }
      }
    }else{
      # if the clone has only one cloneEnd do:
      # get the name of the cloneEnds,
      my $first_cloneEnd =$clones{$pair}[0];
      # send all the alignments for that cloneEnd stored in clone_cluster to the single_filter_alignment 
      my $selected_align = $self->single_filter_alignments($clone_cluster{$first_cloneEnd});
      # add the result of filtering to the selected alignments array.
      push (@selected_alignments, ${$selected_align});
      
    }
  }
  return \@selected_alignments;

}

sub single_filter_alignments{

  my ($self, $single_clone_alignments) = @_;
 
  # Cluster alignments obtained for each clone
  my %chr_cluster = ();
  
  foreach my $alignment(@{$single_clone_alignments}){
   
    my $hit_name= $alignment->hseqname;
    my $chr_name= $alignment->seqname;
    
    # Generate a unique key that will be found only when only clone aligns more than
    # once against the same chromosomal sequence
    my $cluster_key = $chr_name."-".$hit_name;

    # check if there is any sequence aligned from the same clone against the chromosome
    # In case there isn't, start the count of alignments for that clone.
    if(!$chr_cluster{$cluster_key}){
      $chr_cluster{$cluster_key} = [];
    }
    push (@{$chr_cluster{$cluster_key}}, $alignment);
  }

  my $selected_alignment = "";
  my $cluster_size = 0;
  my $biggest_score = 0;

  # For each clone check how many clusters were obtained.
  # Check the clusters and selected the best alignment cluster
  # to be used as template for the more exhaustive exonerate alignment.

  foreach my $chr_cluster_key (keys %chr_cluster){
    # If the cluster contains more than one sequence, order the sequences and check if they are
    # separated by less than 1000 bases. If the separation is bigger don't count them as part of the cluster
    if (scalar(@{$chr_cluster{$chr_cluster_key}})>1){

      # counter to get the real number of alignment within the cluster.
      my $size_counter = 0;

      # get the first aligment of the cluster to be used in the exonerate step
      my $first_selected_alignment;

      # order the alignments in the clster to be able to check the real distance among them
      my  @ordered_aligns = sort{$a->hstart() <=> $b->hstart()} @{$chr_cluster{$chr_cluster_key}};    

      my $previous_end = 0;
      my $previous_alignment;

      my $end_of_cluster = 0;
     
      my %cluster_alignments = ();

      # Use this value to avoid the comparison of the first alignment against null, that in case the start
      # point of the alignment is lower than the cutoff which lead to the clustering of the first alignment
      # with a NULL alignment
      my $first_alignment = 0;
      
      foreach my $ordered_align(@ordered_aligns){
        
        # if the two condiditions are acomplished check if the alignments form a cluster
        if ($first_alignment!=0 && $end_of_cluster == 0){
	  my $cloneEnd_distance = abs($ordered_align->start)-($previous_end);
          # if the distance is sorter that the threshold the alignments cluster
          if($cloneEnd_distance<4000){
          
            if(!$cluster_alignments{$ordered_align->hseqname}){
              $cluster_alignments{$ordered_align->hseqname} = [];
              push (@{$cluster_alignments{$ordered_align->hseqname}}, $previous_alignment);
              $size_counter++;  
              $first_selected_alignment = $previous_alignment;
            }
            push (@{$cluster_alignments{$ordered_align->hseqname}}, $ordered_align);
            $size_counter++;
         
          }elsif($cluster_alignments{$ordered_align->hseqname}){
            $end_of_cluster = 1;
            
	  }
          $previous_alignment = $ordered_align;
          $previous_end = $ordered_align->end;
        }

        # Change this value after the first alignment is checked.
        if ($first_alignment == 0){
          $first_alignment = 1;
          $previous_alignment = $ordered_align;
          $previous_end = $ordered_align->end;
        }
      }
      if ($size_counter > $cluster_size){
        $cluster_size= $size_counter;
        $selected_alignment = $first_selected_alignment;

      } 
    }
    elsif($cluster_size==0 && ${$chr_cluster{$chr_cluster_key}}[0]->score()> $biggest_score){
      $biggest_score = ${$chr_cluster{$chr_cluster_key}}[0]->score();
      $selected_alignment = $chr_cluster{$chr_cluster_key}[0];
    }
  }
  return \$selected_alignment; 
}

sub clean_clusters{

  my ($self, $clone_cluster_alignments) = @_;
 
  # Cluster alignments obtained for each clone
  my %chr_cluster = ();
  
  foreach my $alignment(@{$clone_cluster_alignments}){
     
    my $hit_name= $alignment->hseqname;
    my $chr_name= $alignment->seqname;
  
    # Generate a unique key that will be found only when a cloneEnd aligns more than
    # once against the same chromosomal sequence
    my $cluster_key = $chr_name."-".$hit_name;

    # check if there is any sequence aligned from the same clone against the chromosome
    # In case there isn't, start the count of alignments for that clone.
    if(!$chr_cluster{$cluster_key}){
      $chr_cluster{$cluster_key} = [];
    }
    push (@{$chr_cluster{$cluster_key}}, $alignment);
  }

  my @selection = ();

  # For each clone check how many clusters were obtained.
  # Check the clusters and selected the best alignment cluster
  # to be used as template for the more exhaustive exonerate alignment.

  foreach my $chr_cluster_key (keys %chr_cluster){
    # If the cluster contains more than one sequence, order the sequences and check if they are
    # separated by less than 4000 bases. If the separation is bigger don't count them as part of the cluster

    if (scalar(@{$chr_cluster{$chr_cluster_key}})>1){
      # get the first aligment of the cluster to be used in the exonerate step
      my $first_selected_alignment;

      # order the alignments in the cluster to be able to check the real distance among them
      my  @ordered_aligns = sort{$a->hstart() <=> $b->hstart()} @{$chr_cluster{$chr_cluster_key}};    

      my $initial_end = 0;
      my $initial_alignment;
    
      my %cluster_alignments = ();

      # Use this value to avoid the comparison of the first alignment against null, that in case the start
      # point of the alignment is lower than the cutoff which lead to the clustering of the first alignment
      # with a NULL alignment
      my $first_alignment = 0;
      my $alignment_added = 0;
      
      foreach my $ordered_align(@ordered_aligns){
        
        # if the two condiditions are acomplished check if the alignments form a cluster
        if ($first_alignment!=0){
     
	  my $cloneEnd_distance = abs(($ordered_align->start)-($initial_end));
          # if the distance is sorter that the threshold the alignments cluster
          if($cloneEnd_distance < 4000){
 
            $initial_alignment = $ordered_align;
            $initial_end = $ordered_align->end;
            if ($alignment_added==0){
              push (@selection, $initial_alignment);
	      $alignment_added = 1;
	    }
          }else{
            push (@selection, $ordered_align);
            $initial_alignment = $ordered_align;
            $initial_end = $ordered_align->end;         
          }
        }
        # Change this value after the first alignment is checked.
        if ($first_alignment == 0){
          $first_alignment = 1;
          $initial_alignment = $ordered_align;
          $initial_end = $ordered_align->end;
        }
      }

    }else{
      push (@selection, $chr_cluster{$chr_cluster_key}[0]);
    }
  }
  return \@selection; 
}
