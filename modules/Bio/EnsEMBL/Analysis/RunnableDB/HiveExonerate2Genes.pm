=head1 LICENSE

# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes - 

=head1 SYNOPSIS

my $exonerate2genes = Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -input_id => $chunk_file_name
			     );

$exonerate2genes->fetch_input();
$exonerate2genes->run();
$exonerate2genes->output();
$exonerate2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
It is meant to provide the interface for mapping ESTs to the genome
sequence and writing the results as genes. By the way Exonerate is run
we do not cluster transcripts into genes and only write one transcript per gene.
we then create a dbadaptor for the target database.


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::HiveExonerate2Genes;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::KillList::KillList;
use Bio::SeqIO;

use Data::Dumper;
use feature 'say';
use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable',
            'Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild');

sub fetch_input {
  my($self) = @_;

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                            );
  $self->analysis($analysis);

  $self->db('dna_db',$self->get_dba($self->param('input_db')));
  $self->db('output_db',$self->get_dba($self->param('output_db'),undef,'dna_db'));
  ##########################################
  # set up the target (genome)
  ##########################################

  my @db_files;
  my @target_list = $self->param('GENOMICSEQS');


  foreach my $target (@target_list){ 

    if(ref $target eq 'ARRAY'){
      #check to see if we have multiple files or directories:

      my $dir = 0;
      foreach my $alt_target (@$target){
        if (-d $alt_target){
            $dir = 1;
            last;
        }
      }

      # genome is in multiple directories;  the order of directories determines
      # which file is used in case of duplicates. New versions should therefore
      # be in the directory listed first.

      if ($dir) {

        foreach my $chr_name ($self->get_chr_names) {
          my $found = 0;
          DIRCHECK:
          foreach my $alt_target (@$target){
            if (-s "$alt_target/$chr_name.fa") {
              push @db_files, "$alt_target/$chr_name.fa";
              $found = 1;
              last DIRCHECK;
            }
          }
          if(!$found){
	        warning( "Could not find fasta file for '$chr_name' in directories:\n".
            join("\n\t", @$target)."\n");
          }
        }
      }else{
        foreach my $alt_target (@$target){
          if (-s $alt_target){
             push @db_files, $alt_target; 
          }
        }
      }
    } # // end target is a directory 
    else{
      $target =~s/^\s+//;  
      if (-e $target and -d $target) {
        # genome is in a directory; the directory must contain the complete
        # genome else we cannot do best-in-genome filtering. 
        # 
        # We would like to use exonerate's ability to accept a directory as
        # target (because bestn then works), but we must check that the directory 
        # contains only toplevel sequence files

        my %dir_contents;
        opendir DIR, $target;
        while(my $entry = readdir DIR) {
          if ($entry ne '.' and $entry ne '..') {
            $dir_contents{$entry} = 0;
          }
        }
        closedir(DIR);

        foreach my $chr ($self->get_chr_names) {
          my $seq_fname = "$chr.fa"; 
          if (-s "$target/$seq_fname") {
            $dir_contents{$seq_fname}++;
            push @db_files, "$target/$seq_fname";
          } else {
            warning( "Could not find fasta file for '$chr' in '$target'\n");
          }
        }

        # if all files in dir were expected, we can revert to having
        # the whole directory as target
        if (not grep { $dir_contents{$_} == 0 } keys %dir_contents) {
          @db_files = ($target);
        }
      }
      elsif (-e $target and -s $target) {
        # genome sequence is in a single file
        @db_files = ($target);
      } else {
        throw("'$target' refers to something that could not be made sense of");
      }
    }
  }

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################

  my ($query_file, $chunk_number, $chunk_total);

  my $query = $self->param('QUERYSEQS');
  say "QUERY: ".$query;

  if (-e $query and -d $query) {
    # query seqs is a directory; input id will be a file in that directory
    # As input_id returns a string, I've made it parse out the file name. I don't
    # like this solution but it is the quickest for the moment
    my $query_string = $self->input_id;
    unless($query_string =~ /.+\=\>.+\"(.+)\"/) {
      die "Could not find the chunk file in the input id. Input id:\n".$query_string;
    }

    $query_string = $1;
    $query_file = "$query/" . $query_string;
    if (not -e $query_file) {
      throw( "Query file '$query_file' does not exist'\n");
    }
    if ($self->param('USE_KILL_LIST')) {
      $query_file = filter_killed_entries($query_file, $self->param('KILL_TYPE'), $self->input_id);
      $self->filtered_query_file($query_file);
    }
  }
  elsif (-e $query and -s $query) {
    # query seqs is a single file; input id will correspond to a chunk number
    $query_file = $query;
    my $iid_regexp = $self->param('IIDREGEXP');

    throw("When your input ids are not filenames, you must define ".
          "IIDREGEXP in config to enable inference of chunk number and total")
        if not defined $iid_regexp;

    ($chunk_number, $chunk_total) = $self->input_id =~ /$iid_regexp/;

    ###
    ### DO THE KILL LIST FILTER FOR QUERY FILE. AGAIN THE FILE CAN CONTAIN MULTIPLE ENTIRES
    ###
    if ($self->param('USE_KILL_LIST')) {
      $query_file = filter_killed_entries($query_file, $self->param('KILL_TYPE'));
    }
  } else {
    throw("'$query' refers to something that could not be made sense of\n");
  }

  ##########################################
  # Annotation file with CDS positions
  ##########################################


  ##########################################
  # setup the runnables
  ##########################################

  my %parameters = %{$self->parameters_hash};
  if (not exists($parameters{-options}) and
      defined $self->param('OPTIONS')) {
    $parameters{-options} = $self->param('OPTIONS');
  }
  if (not exists($parameters{-coverage_by_aligned}) and
      defined $self->param('COVERAGE_BY_ALIGNED')) {
    $parameters{-coverage_by_aligned} = $self->param('COVERAGE_BY_ALIGNED');
  }

  if (defined $self->param('PROGRAM') && defined $self->analysis->program_file) {
    if ($self->param('PROGRAM') ne $self->analysis->program_file) {
# I'm just warning because for debugging it's easier to change just the PROGRAM parameters...
      warning("CONFLICT: You have defined -program in your config file and ".
            "-program_file in your analysis table.");
    }
  }

  foreach my $database ( @db_files ){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
        ->new(
              -program  => $self->param('PROGRAM') ? $self->param('PROGRAM') : $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_file    => $database,
              -query_type     => $self->param('QUERYTYPE'),
              -query_file     => $query_file,
              -annotation_file => $self->param('QUERYANNOTATION') ? $self->param('QUERYANNOTATION') : undef,
              -query_chunk_number => $chunk_number ? $chunk_number : undef,
              -query_chunk_total => $chunk_total ? $chunk_total : undef,
              %parameters,
              );
    $self->runnable($runnable);
  }

}

############################################################

sub run {
  my ($self) = @_;
  my @results;

  throw("Can't run - no runnable objects") unless ($self->runnable);

  foreach my $runnable (@{$self->runnable}){
    $runnable->run;
    push ( @results, @{$runnable->output} );
  }
  if ($self->param('USE_KILL_LIST')) {
    unlink $self->filtered_query_file;
    # print "Removed temporary query file ".$self->filtered_query_file."\n";
  }
  if ($self->param('filter')) {
    my $filtered_transcripts = $self->filter->filter_results(\@results);
    @results = @$filtered_transcripts;
  }

  my @genes = $self->make_genes(@results);
  $self->{'output_genes'} = \@genes;
}


############################################################

#sub write_output_jobs {
#  my ($self) = @_;
#  return 1;
#}

sub write_output {
  my ($self,@output) = @_;

  my $dna_db = $self->db('dna_db');
  my $outdb = $self->db('output_db');

  my $gene_adaptor = $outdb->get_GeneAdaptor;

  say Dumper(@{$self->{'output_genes'}});

  unless (@output){
    @output = @{$self->{'output_genes'}};
  }

  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output){
    say "FM2 ATTEMPTING TO STORE GENE!!!!!!!!";
    empty_Gene($gene);
    eval {
      $gene_adaptor->store($gene);
    };
    if ($@){
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ($fails > 0) {
    throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }
}

############################################################

sub db {
  my ($self,$adaptor_name,$value) = @_;

  if($value){
    $self->{$adaptor_name} = $value;
  }
  return $self->{$adaptor_name};
}

sub get_dba {
   my ($self,$connection_info, $non_standard_db_adaptor, $dna_db_name) = @_;
   my $dba;

   if(defined $non_standard_db_adaptor) {
     if($non_standard_db_adaptor eq 'compara') {
       $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
                                                            %$connection_info
                                                          );
       print "Not attaching a dna db to: ".$dba->dbname."\n";
     }
   }

   else {
     $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                     %$connection_info
                                              );

     if($dna_db_name) {
            my $dnadb = $self->db($dna_db_name);

            # try to get default asm+ species name for OTHER db - does not work
            # for comapra database
            my $core_db_asm = $dba->get_MetaContainer->get_default_assembly();
            my $core_db_species =
            $dba->get_MetaContainer->get_common_name();

            # get the same for dna-db
            my $dna_db_asm = $dnadb->get_MetaContainer->get_default_assembly();
            my $dna_db_species =
            $dnadb->get_MetaContainer()->get_common_name();

            my $dna_db_and_core_db_are_compatible = 1;

            unless ( $core_db_species eq $dna_db_species ) {
            warning( "You try to add a DNA_DB with species ".$dna_db_species." to "
                . "a core database with species: '" .$core_db_species . "' - this does not work. \n"
                . "Check that you are using the correct DNA_DB and that the species.common_name values in the meta tables match\n"
            );
            $dna_db_and_core_db_are_compatible = 0;

          }

            if ($dna_db_and_core_db_are_compatible) {
              $dba->dnadb($dnadb);
              print "\nAttaching DNA_DB "
              . $dnadb->dbname . " to "
              . $dba->dbname . "\n";
            }
          }

     else {
            print "Not attaching a dna db to: ".$dba->dbname."\n";
          }

   }

  $dba->dbc->disconnect_when_inactive(1) ;
  return $dba;

 }



sub make_genes{
  my ($self,@transcripts) = @_;

  my (@genes);

  my $slice_adaptor = $self->db('dna_db')->get_SliceAdaptor;
  my %genome_slices;

  foreach my $tran ( @transcripts ){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);

    ############################################################
    # put a slice on the transcript
    my $slice_id = $tran->start_Exon->seqname;
    if (not exists $genome_slices{$slice_id}) {
      # assumes genome seqs were named in the Ensembl API Slice naming
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};

    foreach my $exon (@{$tran->get_all_Exons}){
      $exon->slice($slice);
      foreach my $evi (@{$exon->get_all_supporting_features}){
        $evi->slice($slice);
        $evi->analysis($self->analysis);
      }
    }
    foreach my $evi (@{$tran->get_all_supporting_features}) {
      $evi->slice($slice);
      $evi->analysis($self->analysis);
    }

    if (!$slice){
        my ($sf);

        if (@{$tran->get_all_supporting_features}) {
          ($sf) = @{$tran->get_all_supporting_features};
        } else {
          my @exons = @{$tran->get_all_Exons};
          ($sf) = @{$exons[0]->get_all_supporting_features};    
        }
        print $sf->hseqname."\t$slice_id\n";
    }
 
    throw("Have no slice") if(!$slice);
    $tran->slice($slice);
    $gene->add_Transcript($tran);
    push( @genes, $gene);
  }
  return @genes;
}

############################################################

sub get_chr_names{
  my ($self) = @_;
  my @chr_names;
  my @chromosomes;

  my $chr_adaptor = $self->db->get_SliceAdaptor;
  #also fetching non-reference regions like DR52 for human by default.
  #specify in Exonerate2Genes config-file.
  if(defined($self->param('NONREF_REGIONS'))){
    @chromosomes = @{$chr_adaptor->fetch_all('toplevel', undef, 1)};
  }
  else{
    @chromosomes = @{$chr_adaptor->fetch_all('toplevel')};
  }

  foreach my $chromosome ( @chromosomes ){
    push( @chr_names, $chromosome->seq_region_name );
  }

  return @chr_names;
}


############################################################

sub get_output_db {
  my ($self) = @_;

  my $outdb;

  if ($self->OUTDB) {
    if ( ref($self->OUTDB)=~m/HASH/) {

      $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$self->OUTDB},
                                                -dnadb => $self->db);
    }else{
      $outdb = $self->get_dbadaptor($self->OUTDB);
    }
  } else {
    $outdb = $self->db;
  }
  $self->db->dbc->disconnect_when_inactive(1) ;
  $outdb->dbc->disconnect_when_inactive(1) ;
  return $outdb;
}


############################################################
#
# get/set methods
#
############################################################


sub query_seqs {
  my ($self, @seqs) = @_;
  if( @seqs ) {
    unless ($seqs[0]->isa("Bio::PrimarySeqI") || $seqs[0]->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    push( @{$self->{_query_seqs}}, @seqs);
  }
  return @{$self->{_query_seqs}};
}

############################################################

sub genomic {
  my ($self, $seq) = @_;
  if ($seq){
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_genomic} = $seq ;
  }
  return $self->{_genomic};
}


############################################################

sub database {
  my ($self, $database) = @_;
  if ($database) {
    $self->{_database} = $database;
  }
  return $self->{_database};
}

############################################################

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->{_transcript_filter} = $val;
  }
  return $self->{_transcript_filter};
}

############################################################

sub filtered_query_file {
  my ($self, $val) = @_;
  if ($val) { 
    $self->{_filtered_query_file} = $val;
  }
  return $self->{_filtered_query_file};
}

#############################################################
# Declare and set up config variables
#############################################################

#sub read_and_check_config {
#  my $self = shift;

#  $self->SUPER::read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  ##########
  # CHECKS
  ##########
#  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
#  foreach my $config_var (qw(QUERYSEQS 
#                             QUERYTYPE
#                             GENOMICSEQS)) {

#   throw("You must define $config_var in config for logic '$logic'")
#        if not defined $self->$config_var;
#  }
  
#  throw("QUERYANNOTATION '" . $self->QUERYANNOTATION . "' in config must be readable")
#      if $self->QUERYANNOTATION and not -e $self->QUERYANNOTATION;

  # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
#  if ($self->FILTER) {
#    if (not ref($self->FILTER) eq "HASH" or
#        not exists($self->FILTER->{OBJECT}) or
#        not exists($self->FILTER->{PARAMETERS})) {
          
#      throw("FILTER in config fo '$logic' must be a hash ref with elements:\n" . 
#            "  OBJECT : qualified name of the filter module;\n" .
#            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
#    } else {
#      my $module = $self->FILTER->{OBJECT};
#      my $pars   = $self->FILTER->{PARAMETERS};
      
#      (my $class = $module) =~ s/::/\//g;
#      eval{
#        require "$class.pm";
#      };
#      throw("Couldn't require ".$class." Exonerate2Genes:require_module $@") if($@);
#    
#      $self->filter($module->new(%{$pars}));
#    }
#  }
#}



###############################################
###     end of config
###############################################

sub filter_killed_entries {
  my ($orig_query_filename, $mol_type, $inputID) = @_;
  my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => $mol_type);
  my %kill_list = %{ $kill_list_object->get_kill_list() };

  my $seqin  = new Bio::SeqIO(-file   => "<$orig_query_filename",
                            -format => "Fasta",
                          );

  my $filtered_seqout_filename = "/tmp/$inputID"."_filtered";
  print "Filename for my filtered sequence: $filtered_seqout_filename.\n";

  my $seqout = new Bio::SeqIO(-file   => ">$filtered_seqout_filename",
                              -format => "Fasta"
                           );

  while( my $query_entry = $seqin->next_seq ){
    my $display_id  = $query_entry->display_id;
    my $no_ver_id;
    # Depending on the display ID's format, strip off the
    # version number because the kill_list hash keys are
    # without version numbers
  
    if ($display_id =~/\w+\.\d/) {
      ($no_ver_id) = $display_id =~/(\w+)\.\d/;
    } elsif ($display_id =~/\w+\-\d/) {
      ($no_ver_id) = $display_id =~/(\w+)\-\d/;
    } elsif ($display_id =~/\w+/ ) {
      ($no_ver_id) = $display_id;
    }
    if ( !$kill_list{$no_ver_id} ) {
      $seqout->write_seq($query_entry);
    } elsif ( $kill_list{$no_ver_id} ) {
      print "$mol_type $display_id is in the kill_list. Discarded from analysis.\n";
    }
  }
  return $filtered_seqout_filename; 
}    

1;
