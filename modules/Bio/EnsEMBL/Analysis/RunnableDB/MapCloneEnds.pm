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
      # Empty the array where data is going to be stored
#      @cloneEnd = ();

      # Add cloneEnd, clone length and clone standard deviation to the array
#      push (@cloneEnd, $trace_name);
#      push (@cloneEnd, $insert_size);  
#      push (@cloneEnd, $insert_stdev);
   
      # Create the clone object that will store information of corresponding cloneEnds
      if (!$clones{$clone_id}){
        $clones{$clone_id}=[];
      }
      
      # Add the information of a cloneEnd to the corresponding clone
 #     push (@{$clones{$clone_id}}, @cloneEnd);
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
  open (INFILE,"<$chunks_list");
  my @chunks_list = <INFILE>;
  close INFILE;
#  $self->chunks_list(\@chunks_list);

#  my %chunks = %{$self->create_chunks(\%clones)};

  #foreach chunk run the more generic exonerate analysis. 
  #foreach my $chunk_id (keys %chunks){

  # Run ExonerateClones to get a first idea on the aligned/location of the clone Ends 
    my $exonerate = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateCloneEnds->new(
                    -DB          => $self->db,
                    -INPUT_ID    => $self->input_id,
                    -ANALYSIS    => $self->fetch_analysis("EXONERATE_CLONE_ENDS"),
                    );

    $exonerate->fetch_input(\@chunks_list);

    $exonerate->run();
  
  ## Run exonerate and get the output.
  my $clone_alignments = $exonerate->output();
#  }
  print "I'm in step 3\n";
  my ($pointer, @selected_alignments) = @{$self->filter_alignments($clone_alignments, \%clones, \%cloneEnd_ids)};

 # foreach my $test( @selected_alignments){
 #   print "This is : ",$test,"\n";
 # }

 # my @refine_output;
  
    foreach my $selected_alignment(@selected_alignments){
    
      if ($selected_alignment ne ''){
    
        #print "Test here: ", $selected_alignment,"\n"; 
        my $clone_id=$selected_alignment->hseqname;
        my $chr_id = $selected_alignment->seqname;
        my $start = ($selected_alignment->start)-1000;
        my $end = ($selected_alignment->end)+1000;

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
#  }
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

#sub chunks_list{
#  my ($self, $chunks_list) = @_;
#
#  if ($chunks_list) { 
#   push @{$self->{_chunks_list}}, $chunks_list; 
#  }  
#
#  return $self->{_chunks_list};
#
#}

sub create_chunks{
  my ($self, $clones2chunk)= @_;
 
  my %clones = %{$clones2chunk};
  my %chunks; # hash where the arrays of seuqence chunks will be stored
  my @clone_pairs = ();
  my @single_clones = ();

  #Counter of files added to a chunk, once reacher a certain number of files, let's say 30,
  # the chunk is send to the pipeline to run exonerate.
  my $files_in_chunk = 0;
  my $single_files_chunk = 0;

  my $chunk_number = 0;   

  foreach my $clone (keys %clones){
    
    # Check if the clone has two clone ends or it's lacking one of the ends
    if ($clones{$clone}[3]){

      # check the number of files in the chunk
      if ($files_in_chunk < 30){

        # Add the two cloneEnd ids to the list of files in the same chunk
	push (@clone_pairs, $clones{$clone}[0]);
	push (@clone_pairs, $clones{$clone}[3]);

        # Plus 2 becuase two files are added each time in the case of cloneEnd pairs
        $files_in_chunk += 2;
 
      # If files in chunk is 30 the load the chunk to the hash and start a new chunk     
      }elsif ($files_in_chunk == 30){
        
        # Load the chunk to the hash
	$chunks{$chunk_number}= @clone_pairs;
     
        # Empty the existing hash
        @clone_pairs = ();

        # Add the first two elements to the hash
      	push (@clone_pairs, $clones{$clone}[0]);
	push (@clone_pairs, $clones{$clone}[3]);

        # Files in chunk is set to 2 becuase two files were already added when starting new chunk
        $files_in_chunk = 2;
 
        # Add one to change the chunk number in the hash.
        $chunk_number++;    
      }
    }else{
      # Hash where you add the sequences where there are no CloneEnd pairs
      if ($single_files_chunk < 30){
        push (@single_clones, $clones{$clone}[0]);
        $single_files_chunk++;
      }else{
        # Load the chunk to the hash
        $chunks{$chunk_number}= @single_clones;

        # Empty the existing hash
        @single_clones = ();

        # Add the first element to the hash
        push (@single_clones, $clones{$clone}[0]);

        $single_files_chunk = 0;
        # Add one to change the chunk number in the hash.
        $chunk_number++;         

      }
    }
  }

  # Add the last chunk of pairs to the hash (this is needed in case the number of sequences is not
  # a multiple of 30)
  if ($files_in_chunk > 2 && $files_in_chunk < 30){
    $chunks{$chunk_number}= @clone_pairs;

    # add one to chunk number becuase single_cloneEnds array was not added yet.
    $chunk_number++;
  }

  # Add the last chunk of single_cloneEnds to the hash
  if ($single_files_chunk < 30){

  # Add single cloneEnd array.
  $chunks{$chunk_number} = @single_clones;
  }

  return \%chunks;
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
   
    my $clone_name = $clone_alignment->hseqname();

    # Get the clone to which the cloneEnds belongs
    my $complete_clone_name = $cloneEnd_ids{$clone_name};

    # Store the clone name in a hash, so only selected clones will be checked in the next step
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
      
      # Get the length of the clone and the name of the two cloneEnds
      my $clone_length = $clones{$pair}[1]+$clones{$pair}[2];
      my $first_cloneEnd =$clones{$pair}[0];
      my $second_cloneEnd =$clones{$pair}[3];

      # Compare all the alignments in both cloneEnds to find those which pair according to the
      # position where they align  and the length of the clone (distance between both cloneEnds
      foreach my $fa (@{$clone_cluster{$first_cloneEnd}}){
        foreach my $sa (@{$clone_cluster{$second_cloneEnd}}){
          my $chr_first = $fa->seqname();
          my $chr_second = $sa->seqname();
                         
          if ($chr_first eq $chr_second){
                           
            my $diff = $fa->start()-$sa->end();
            my $abs_diff = abs($diff);
            if ($abs_diff <= $clone_length){
              # In case two alignments are selected store them for the next exonerate step
              push (@selected_alignments, $fa);
              push (@selected_alignments, $sa); 
	    }
          }
        }
      }
    }else{
      # if the clone has only one cloneEnds do:
      # get the name of the cloneEnds,
      my $first_cloneEnd =$clones{$pair}[0];
      # send all the alignments for that cloneEnd stored in clone_cluster to the single_filter_alignment 
      my $selected_align = single_filter_alignment(\$clone_cluster{$first_cloneEnd});
      # add the result of filtering to the selected alignments array.
      push (@selected_alignments, $selected_align);
      
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
  return $selected_alignment; 
}

