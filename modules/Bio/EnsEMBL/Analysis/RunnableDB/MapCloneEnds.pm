=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds;

=head1 SYNOPSIS


=head1 DESCRIPTION

This object maps clone sequences to a genome,
and writing the results as Dna Align Features. 

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
#use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateClones;
#use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateRefinedCloneEnds;
use Bio::EnsEMBL::Analysis::Config::MapCloneEnds;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);



# Run ExonerateClones to get a first idea on the aligned/location of the clone Ends 
 my $exonerate = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateClones->new(
  );

$exonerate->fetch_input();

## Run exonerat and get the output.
@clone_alignments = @{$exonerate->run()};


my  @ordered_alignments = sort{$a->hseqname() <=> $b->hseqname()} @clone_alignments;


### Filter the results of the previous step and prepare data to send to the refine step.


my %clone_cluster = ();

# Group alignments by clone name
foreach my $clone_alignment(@ordered_alignments){

  my $clone_name = $clone_alignment->hseqname();

 if(!$clone_cluster{$clone_name}){
    $clone_cluster{$clone_name} = [];
  }
 push (@{$clone_cluster{$clone_name}}, $clone_alignment);
}

my @selected_alignments = [];

# Cluster alignments obtained for each clone
foreach my $clone_key(keys %clone_cluster){  
  my @alignments = @{$clone_cluster{$clone_key}};

  foreach my $alignment(@alignments){
 
  my %chr_cluster = ();
  my $hit_name= $alignment->hseq_name;
  my $chr_name= $alignment->seq_name;
  
# print $seqname,"\n";

   
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
     

      # Use this value to avoid the comparison of the first alignment against null, that in case the start
      # point of the alignment is lower than the cutoff which lead to the clustering of the first alignment
      # with a NULL alignment
      my $first_alignment = 0;
    
      
      foreach my $ordered_align(@ordered_aligns){
        
        # if the two condiditions are acomplished check if the alignments form a cluster
        if ($first_alignment!=0 && $end_of_cluster == 0){
        
          # if the distance is sorter that the threshold the alignments cluster
          if((($ordered_align->start)-($previous_end))<1000 && (($ordered_align->start)-($previous_end))>0){
            
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
        $selected_alignment = $

      } 
    }
    elsif($cluster_size==0 && {$chr_cluster{$chr_cluster_key}}[0]->score()> $biggest_score){
     $biggest_score = {$chr_cluster{$chr_cluster_key}}[0]->score();
     $selected_alignment = $chr_cluster{$chr_cluster_key}[0];
    }
  }
  push (@selected_alignments, $selected_alignment); 
}
 
#### Prepare to run ExonerateRefineCloneEnds.pm

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'ensro',
        -dbname => 'homo_sapiens_core_36_35i',
        -host   =>  'ecs2',
        -port   =>  '3364',
        -driver => 'mysql');

my  $slice_adaptor = $db->get_SliceAdaptor();


####
# foreach $selected_alignment(@selectted_alignemnts){
#  
# my $clone_id=$selected_alignment->hseqname;
# my $chr_id = $selected_alignment->seqname;
# my $start = ($selected_alignment->start)-1000;
# my $end = ($selected_alignment->end)+1000;


#It has to be run for one sequence at a time so it might go inside a foreach loop??

my $slices = $slice_adaptor->fetch_by_region('chromosome',$chr, $start, $end);
#my $slices = $slice_adaptor->fetch_by_region('supercontig',$chr, $start, $end);

my $logic = $self->analysis->logic_name;

 my $refine = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateRefinedCloneEnds->new(

  );

$refine ->fetch_input();
$refine ->run();
$refine ->write_output();

#}
