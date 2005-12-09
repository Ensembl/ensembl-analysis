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

#  my @cloneEnd = ();
  my $trace_name = '';  # ClonEnd id as in database
  my $clone_id = '';# Clone id
  my $insert_size = '';# Length of the inserted sequence in clone
  my $insert_stdev = '';# Standard deviation of the clone insert. 
  # The two insert values will be used in the filter process to check quality of alignments


  ### MOVE THIS TO A CONFIG FILE ???###
  my $xmlfile = '/ecs2/scratch2/jb16/sheep/traceinfo/sheep_trace_info.xml'; # Allowed to be a file or a dir of xml files

  # You need to run  perl ~/cvs_checkout/ensembl-personal/jb16/scripts/general/chunker.pl in order to create the 
  # input_analysis_ids and create this file. To change the path to the list of chunks, you should also do it in
  # the chunker script
  my $chunks_list = '/ecs2/scratch2/jb16/sheep/listOfChunks.txt';
  ##################################
  my $single_entry = 0;

  my %clones  = (); # hash where each clone object will be stored

  my %cloneEnd_ids = ();

  print "I'm in step 1\n";

  open (IN,"<$xmlfile");

  # parses the XML file to get the information of which cloneEnds belong to the same Clone and which
  # is the average size of the clone insert
  while(my $line =<IN>){
  
    if ($line =~ /^\s+<trace>/){
      $single_entry = 1;
    }

    if ($single_entry == 1 && ($line =~ /^\s+<trace_name>/)){
      $line =~ /^\s+<trace_name>(\d+)<\/trace_name>/;
      $trace_name = $1;
    }

    if ($single_entry == 1 && ($line =~ /^\s+<clone_id>/)){
      $line =~ /^\s+<clone_id>(\w+[-_]\w+)<\/clone_id>/;
      $clone_id = $1;
    }

    if ($single_entry == 1 && ($line =~ /^\s+<insert_size>/)){
      $line =~ /^\s+<insert_size>(\d+.\d*)<\/insert_size>/;
      $insert_size = $1;
    }

    if ($single_entry == 1 && ($line =~ /^\s+<insert_stdev>/)){
      $line =~ /^\s+<insert_stdev>(\d+.\d*)<\/insert_stdev>/;
      $insert_stdev = $1;
    }
  
    if (($line =~ /^\s+<\/trace>/) && $single_entry == 1) {

      # Create the clone object that will store information of corresponding cloneEnds
      if (!$clones{$clone_id}){
        $clones{$clone_id}=[];
      }
      
      # Add the information of a cloneEnd to the corresponding clone
      push (@{$clones{$clone_id}}, $trace_name);
      push (@{$clones{$clone_id}}, $insert_size);  
      push (@{$clones{$clone_id}}, $insert_stdev);

      # Hash where with relation between clone and cloneEnd, will be used in the filter step      
      $cloneEnd_ids{$trace_name} = $clone_id;

      $single_entry = 0;
    }
     
  }  

  close IN;
 
  print "I'm in step 2\n";

  # Open the file with the list of seq_ids per chunks and load it into an array to be used by fetch_input
  open (INFILE,"<$chunks_list");

  my @chunks_list = <INFILE>;

  close INFILE;

  # Run ExonerateClones to get a first idea on the aligned/location of the clone Ends 
    my $exonerate = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds->new(
                    -DB          => $self->db,
                    -INPUT_ID    => $self->input_id,
                    -ANALYSIS    => $self->fetch_analysis("EXONERATE_CLONE_ENDS"),
                    );

    $exonerate->fetch_input(\@chunks_list);

    $exonerate->run();
  
  # Run exonerate and get the output.
  my $clone_alignments = $exonerate->output();

  print "I'm going to filter the alignments\n";
  my @selected_alignments = @{$self->filter_alignments($clone_alignments, \%clones, \%cloneEnd_ids)};

  foreach my $selected_alignment(@selected_alignments){
    
    if ($selected_alignment ne ''){
    
      #print "Test here: ", $selected_alignment,"\n"; 
      my $clone_id=$selected_alignment->hseqname;
      my $chr_id = $selected_alignment->seqname;
      my $start = ($selected_alignment->start)-1000;
      my $end = ($selected_alignment->end)+1000;
      print "Clone id: ",$clone_id,"\n";
      my @chr_name = split (/:/, $chr_id);
  
      my $input_id = $chr_name[0].":".$chr_name[1].":".$chr_name[2].":".$start.":".$end.":".$chr_name[5].":".$clone_id;
      #print $input_id,"\n";

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

    # Get the clone to which the cloneEnds belongs
    my $complete_clone_name = $cloneEnd_ids{$clone_name};

    # Store the clone name in a hash, so only selected clones will be checked in the next step
    # the object stored in the hash is not important, what we really want to get is a list of 
    # unique clones(complete_clone_name) with no duplicate entries
    $aligned_clones{$complete_clone_name} = $clone_name;

  #  print $clone_name,"\n";


    if(!$clone_cluster{$clone_name}){
      $clone_cluster{$clone_name} = [];
    }
    push (@{$clone_cluster{$clone_name}}, $clone_alignment);
  }


  foreach my $pair (keys %aligned_clones){
    
    # Check if the clone has two cloneEnds
    if ($clones{$pair}[3]){

      # Chenk if at least one of the cluster pair alignments is selected.
      # in case none is selected it will send cloneEnds to the single_filter
      my $cluster_selected = 0;
      # Get the length of the clone and the name of the two cloneEnds
      my $clone_length = $clones{$pair}[1]+$clones{$pair}[2]+200000;
      my $first_cloneEnd =$clones{$pair}[0];
      my $second_cloneEnd =$clones{$pair}[3];

     # print "This is the first: ",$clone_cluster{$first_cloneEnd},"\n";
      my @first_cloneEnd_clean  = @{$self->clean_clusters($clone_cluster{$first_cloneEnd})};
      my @second_cloneEnd_clean = @{$self->clean_clusters($clone_cluster{$second_cloneEnd})};

      # Compare all the alignments in both cloneEnds to find those which pair according to the
      # position where they align  and the length of the clone (distance between both cloneEnds
      foreach my $fa (@first_cloneEnd_clean){
	
        foreach my $sa (@second_cloneEnd_clean){
	
          my $chr_first = $fa->seqname();
          my $chr_second = $sa->seqname();
   #       print  "First name: ",$chr_first,"  Second name: ",$chr_second,"\n";          
          if ($chr_first eq $chr_second){
                           
            my $diff = ($fa->start())-($sa->end());
            my $abs_diff = abs($diff);
           # print "Absolute difference: ",$abs_diff,"\n";
            if ($abs_diff <= $clone_length){
              # In case two alignments are selected store them for the next exonerate step
              push (@selected_alignments, $fa);
              push (@selected_alignments, $sa); 
              $cluster_selected = 1;
	    }
          }
        }
      }
      if ($cluster_selected == 0){
        # send all the alignments for that cloneEnd stored in clone_cluster to the single_filter_alignment 
        my $first_selected_align = $self->single_filter_alignments($clone_cluster{$first_cloneEnd});
        my $second_selected_align = $self->single_filter_alignments($clone_cluster{$second_cloneEnd});
       # add the result of filtering to the selected alignments array.
        push (@selected_alignments, ${$first_selected_align});
        push (@selected_alignments, ${$second_selected_align});
      }
    }else{
	"I enter to the single_cloneEnd\n";
      # if the clone has only one cloneEnd do:
      # get the name of the cloneEnds,
      my $first_cloneEnd =$clones{$pair}[0];
      # send all the alignments for that cloneEnd stored in clone_cluster to the single_filter_alignment 
      my $selected_align = $self->single_filter_alignments($clone_cluster{$first_cloneEnd});
      # add the result of filtering to the selected alignments array.
      push (@selected_alignments, $selected_align);
      
    }
  }
  print "selected in filter alignments: ",@selected_alignments,"\n";
  return \@selected_alignments;

}

sub single_filter_alignments{

  my ($self, $single_clone_alignments) = @_;
 
  # Cluster alignments obtained for each clone
  my %chr_cluster = ();
  
#  print "Selected single filtering\n";
  foreach my $alignment(@{$single_clone_alignments}){
   
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
  # print $selected_alignment->seqname,"\t",$selected_alignment->hseqname,"\t",$selected_alignment->start,"\n";
  print "Selected in single filtering: ",$selected_alignment,"\n";
  return \$selected_alignment; 
}

sub clean_clusters{

  my ($self, $clone_cluster_alignments) = @_;
 
  # Cluster alignments obtained for each clone
  my %chr_cluster = ();
  
  print "I'm in here:", $clone_cluster_alignments,"\n";
  foreach my $alignment(@{$clone_cluster_alignments}){
 #    print "Alignment: ",$alignment,"\n";
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

  my @selection = ();

  # For each clone check how many clusters were obtained.
  # Check the clusters and selected the best alignment cluster
  # to be used as template for the more exhaustive exonerate alignment.

  foreach my $chr_cluster_key (keys %chr_cluster){
    # If the cluster contains more than one sequence, order the sequences and check if they are
    # separated by less than 1000 bases. If the separation is bigger don't count them as part of the cluster
    if (scalar(@{$chr_cluster{$chr_cluster_key}})>1){
	print "In the same chromosome\n";
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
      
      foreach my $ordered_align(@ordered_aligns){
        
        # if the two condiditions are acomplished check if the alignments form a cluster
        if ($first_alignment!=0){
     
	my $cloneEnd_distance = abs(($ordered_align->start)-($initial_end));
	print "Clean distance: ",$cloneEnd_distance,"\n";
          # if the distance is sorter that the threshold the alignments cluster
          if($cloneEnd_distance < 4000){
 
            $initial_alignment = $ordered_align;
            $initial_end = $ordered_align->end;
  
          }else{
            push (@selection, $initial_alignment);
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
  # print $selected_alignment->seqname,"\t",$selected_alignment->hseqname,"\t",$selected_alignment->start,"\n";
  print"Selected in clean Clusters: ", scalar(@selection),"\n";
  return \@selection; 
}
