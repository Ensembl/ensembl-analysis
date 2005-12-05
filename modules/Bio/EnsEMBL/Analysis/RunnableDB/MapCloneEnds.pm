=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds;

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

This object maps clone sequences to a genome by using
exonerate alignment program,and write the results as 
Dna Align Features.
It needs to have installed the following modules:
Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds;
Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds;

=head1 CONTACT

Post general queries to <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::MapCloneEnds;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

##########################################################################

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->{_refined_results}=[] ; 
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

# Run ExonerateClones to get a first idea on the aligned/location of the clone Ends 
  my $exonerate = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds->new(
    -DB          => $self->db,
    -INPUT_ID    => $self->input_id,
    -ANALYSIS    => $self->fetch_analysis("EXONERATE_CLONE_ENDS"),
  );

  $exonerate->fetch_input();

  $exonerate->run();

  ## Run exonerat and get the output.
  my $clone_alignments = $exonerate->output();

  my ($pointer, @selected_alignments) = @{$self->filter_alignments($clone_alignments)};

 # foreach my $test( @selected_alignments){
 #   print "This is : ",$test,"\n";
 # }

  my @refine_output;
  
  foreach my $selected_alignment(@selected_alignments){
    if ($selected_alignment ne ''){
    print "Test here: ", $selected_alignment,"\n"; 
    my $clone_id=$selected_alignment->hseqname;
    my $chr_id = $selected_alignment->seqname;
    my $start = ($selected_alignment->start)-1000;
    my $end = ($selected_alignment->end)+1000;

 my $clone_id=$selected_alignment->hseqname;
    my @chr_name = split (/:/, $chr_id);
  
    my $input_id = $chr_name[0].":".$chr_name[1].":".$chr_name[2].":".$start.":".$end.":".$chr_name[5].":".$clone_id;
    print $input_id,"\n";
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

sub filter_alignments{
  my ($self, $clone_alignments) = @_;
 
  my  @ordered_alignments = sort{$a->hseqname() <=> $b->hseqname()} @{$clone_alignments};


  ### Filter the results of the previous step and prepare data to send to the refine step.

  my %clone_cluster = ();

  # Group alignments by clone name
  foreach my $clone_alignment(@ordered_alignments){
   
    my $clone_name = $clone_alignment->hseqname();

  #  print $clone_name,"\n";

    if(!$clone_cluster{$clone_name}){
      $clone_cluster{$clone_name} = [];
    }
    push (@{$clone_cluster{$clone_name}}, $clone_alignment);
  }

  my @selected_alignments = [];

  # Cluster alignments obtained for each clone
  foreach my $clone_key(keys %clone_cluster){  
    my @alignments = @{$clone_cluster{$clone_key}};
    my %chr_cluster = ();
  
    foreach my $alignment(@alignments){
   
      my $hit_name= $alignment->hseqname;
      my $chr_name= $alignment->seqname;
  
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
     
        my %cluster_alignments = ();

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
          $selected_alignment = $first_selected_alignment;

        } 
      }
      elsif($cluster_size==0 && ${$chr_cluster{$chr_cluster_key}}[0]->score()> $biggest_score){
        $biggest_score = ${$chr_cluster{$chr_cluster_key}}[0]->score();
        $selected_alignment = $chr_cluster{$chr_cluster_key}[0];
      }
    }
   # print $selected_alignment->seqname,"\t",$selected_alignment->hseqname,"\t",$selected_alignment->start,"\n";
    push (@selected_alignments, $selected_alignment); 
  }
 # print @selected_alignments,"\n";
  return \@selected_alignments;
}

