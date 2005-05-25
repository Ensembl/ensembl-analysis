# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Blast
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::BlastRfam

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastRfam->
  new(
      -query => $slice,
      -program => 'wublastn',
      -database => 'embl_vertrna',
      -options => 'hitdist=40 -cpus=1',
      -parser => $bplitewrapper,
      -filter => $featurefilter,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

Modified blast runnable for specific use with RFAMSEQ.
Use for running BLASTN of genomic vs RFAMSEQ prior to 
ncRNA analysis using Infernal.
Keeps the coverage in the dna_align_feature score field.
Also clusters overlapping hits and picks the one with the lowest 
evalue to represent that cluster.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk


=cut


package Bio::EnsEMBL::Analysis::Runnable::BlastRfam;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Blast);

=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string, directory
  Function  : a generic run method. This checks the directory specifed
  to run it, write the query sequence to file, marks the query sequence
  file and results file for deletion, runs the analysis parses the 
  results and deletes any files
  Returntype: 1
  Exceptions: throws if no query sequence is specified
  Example   : 

=cut


sub run{
  my ($self, $dir) = @_;
  my $coverage = 70; # Coverage filter for very common blast hits
  $self->workdir($dir) if($dir);
  throw("Can't run ".$self." without a query sequence") 
    unless($self->query);
  $self->checkdir();
  my $filename = $self->write_seq_file();
  $self->files_to_delete($filename);
  $self->files_to_delete($self->resultsfile);
  $self->run_ncbi_analysis("-W7 -F F -b 1000000 -v 1000000","/ecs2/work2/sw4/RFAM/Partition/low_copy.fasta");  
  $self->parse_results(20);
  # evil way of emptying the result array!
  $self->{'results_files'} = ();
  $self->database("/ecs2/work2/sw4/RFAM/Partition/high_copy.fasta");
  # NOTE TO SELF need to set the analysis to use wublastn or it will have a benny unless you want to hard code it?
  $self->run_analysis;
  $self->parse_results($coverage);
  $self->delete_files;
  return 1;
}


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Blast
  Function  : override results parsing to allow for 
coverage calculations
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub parse_results{
  my ($self,$coverage_cutoff) = @_;
  my $results = $self->results_files;
  my @daf_coverage_results;
  my $filtered_output;
  my $bplite = $self->parser->get_parsers($results);
  foreach my $blast (@{$bplite}){
      while( my $subject = $blast->nextSbjct){
	 while (my $hsp = $subject->nextHSP) {
	   my @daf_results;
	   my $hsp_length = $hsp->length."\t";
	   my $subject_length = $subject->{'LENGTH'};
	   $subject_length = 1 unless ($subject_length);
	   my $coverage = $hsp_length/$subject_length*100;
	   $coverage =~ s/\.\d+//;
	   if ($coverage_cutoff){
	     next unless($coverage > $coverage_cutoff);
	   }
	   $subject->name =~ /^(\S+)\/\S+\s+(\w+);\w+/;
	   my $name = $2."-".$1;
	   push @daf_results, $self->parser->split_hsp($hsp,$name);
	   # add coverage into daf score?
	   foreach my $daf(@daf_results){
	     # swaps strands over if blast aligns +strand genomic to -ve strand in 
	     # RFAM file, Thing is RFAM file sequences are the correct orientation
	     # Whereas the genomic can be either, with unspliced DNA it becomes
	     # impossibe to tell automatically
	     if ($daf->hstrand == -1 && $daf->strand == 1){
	       $daf->strand(-1);
	       $daf->hstrand(1);
	     }
	     $daf->score($coverage);
	     push  @daf_coverage_results, $daf;
	   }
	 }
      }
    }
#  print "Before clustering ".scalar( @daf_coverage_results)."\n";
  return undef unless  ( @daf_coverage_results);
  my $output = $self->cluster(\@daf_coverage_results);
  $self->output($output);
}

# foreach hit, look at all the other hits and push into an array anything that overlaps with it
# then sort the array and take the top scoring hits for heach familly that lives in it.
# you have 3 instances of a variable called RFAM in the same method all of which hold different values

sub cluster{
  my ($self,$dafs_ref)=@_;
  my @dafs = @$dafs_ref;
  my $start =0;
  my @representative_sequences;
 DAFS: foreach my $daf (@dafs){
    my %family_cluster;
    $start ++;
    next DAFS unless($daf);	
    my $RFAM = substr($daf->hseqname,0,7);
    push @{$family_cluster{$RFAM}},$daf;
  MATCHES:  for (my $index = $start; $index <= $#dafs ; $index ++){
      next MATCHES unless ($dafs[$index]);
      if ($daf->overlaps($dafs[$index])){
	my $familly = substr($dafs[$index]->hseqname,0,7);
	push @{$family_cluster{$familly}},$dafs[$index];
	$dafs[$index] = undef;
      }
    }
    foreach my $domain (keys %family_cluster){
      my @temp_array = @{$family_cluster{$domain}};
      # get highest scoring alignment for each family to be representative
      @temp_array = sort{$a->p_value <=> $b->p_value} @temp_array;
      push @representative_sequences,shift @temp_array;
    }
  }
  return \@representative_sequences;
}

#####################################################
# USE THIS TO RUN NCBI BLAST IF YOU ARE SO INCLINED #
#####################################################


=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Blast
  Function  : gets a list of databases to run against and constructs
  commandlines against each one and runs them
  Returntype: none
  Exceptions: throws if there is a problem opening the commandline or
  if blast produces an error
  Example   : 

=cut



sub run_ncbi_analysis {
  my ($self,$options,$database) = @_;

  # This routine expands the database name into $db-1 etc for
  # split databases
  
  #allow system call to adapt to using ncbi blastall. 
  #defaults to WU blast
  my $command  = "/usr/local/ensembl/bin/blastall";
  my $filename = $self->queryfile;
  my $results_file = $self->create_filename("Infernal", 'blast.out');
  $self->files_to_delete($results_file);
  $self->results_files($results_file);
  $command .= " -p blastn -d $database -i $filename ";
  $command .= "$options 2>&1 > ".$results_file;
  print "Running blast ".$command."\n";
  open(my $fh, "$command |") || 
    throw("Error opening Blast cmd <$command>." .
	  " Returned error $? BLAST EXIT: '" . 
	  ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) . 
	  "', There was " . ($? & 128 ? 'a' : 'no') . 
	  " core dump");
  # this loop reads the STDERR from the blast command
  # checking for FATAL: messages (wublast) [what does ncbi blast say?]
    # N.B. using simple die() to make it easier for RunnableDB to parse.
  while(<$fh>){
    if(/FATAL:(.+)/){
      my $match = $1;
      print $match;
      if($match =~ /no valid contexts/){
	die qq{"VOID"\n}; # hack instead
      }elsif($match =~ /Bus Error signal received/){
	die qq{"BUS_ERROR"\n}; # can we work out which host?
      }elsif($match =~ /Segmentation Violation signal received./){
	die qq{"SEGMENTATION_FAULT"\n}; # can we work out which host?
      }elsif($match =~ /Out of memory;(.+)/){
	# (.+) will be something like "1050704 bytes were last 
	#requested."
	die qq{"OUT_OF_MEMORY"\n}; 
	# resenD to big mem machine by rulemanager
      }elsif($match =~ /the query sequence is shorter 
	     than the word length/){
	#no valid context 
          die qq{"VOID"\n}; # hack instead
        }else{
          warning("Something FATAL happened to BLAST we've not ".
                  "seen before, please add it to Package: " 
                  . __PACKAGE__ . ", File: " . __FILE__);
          die ($self->unknown_error_string."\n"); 
          # send appropriate string 
          #as standard this will be failed so job can be retried
          #when in pipeline
        }
    }elsif(/WARNING:(.+)/){
      # ONLY a warning usually something like hspmax=xxxx was exceeded
      # skip ...
    }elsif(/^\s{10}(.+)/){ # ten spaces
      # Continuation of a WARNING: message
      # Hope this doesn't catch more than these.
      # skip ...
    }
  }
  unless(close $fh){
    # checking for failures when closing.
    # we should't get here but if we do then $? is translated 
    #below see man perlvar
    warning("Error running Blast cmd <$command>. Returned ".
	    "error $? BLAST EXIT: '" . ($? >> 8) . 
	    "', SIGNAL '" . ($? & 127) . "', There was " . 
              ($? & 128 ? 'a' : 'no') . " core dump");
    die ($self->unknown_error_string."\n"); 
    }
}




1;
