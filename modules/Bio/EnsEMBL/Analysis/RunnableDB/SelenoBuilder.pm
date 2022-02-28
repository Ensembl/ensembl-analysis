=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::SelenoBuilder - 

=head1 SYNOPSIS

my $exonerate4selenos = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateForSelenos->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -input_id => $selenos_file_name
			     );

$exonerate4selenos->fetch_input();
$exonerate4selenos->run();
$exonerate4selenos->output();
$exonerate4selenos->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
It is meant to provide the interface for aligning cdnas of selenocysteine 
containing transcripts to the genome sequence and writing the results as genes. 
By the way Exonerate is run we do not cluster transcripts into genes and only 
write one transcript per gene.
we then create a dbadaptor for the target database.


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::SelenoBuilder;

use warnings ;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name write_seqfile);


use Bio::EnsEMBL::Analysis::Config::SelenoBuild;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


############################################################
sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  return $self;
}


sub fetch_input {
  my( $self) = @_;

  my $logic = $self->analysis->logic_name;

  $self->throw("No input id") unless defined($self->input_id);
    
  $self->fetch_sequence();

#  my $discarded_db = $self->get_dbadaptor("DISCARDED_DB");

#  print "DISCARDED GENE DB: ", $discarded_db->dbname,"\n";

  # database where the genebuild produced genes are
  my $seleno_db = $self->get_dbadaptor("SELENO_DB") ;

  print "ENSEMBL DB : ",  $seleno_db->dbname,"\n";
  
  my $ref_db = $self->get_dbadaptor("REFERENCE_DB");
  
  print $self->input_id,"\n";

  #@input_id = split("-",$self->input_id);

  my $slice = $ref_db->get_SliceAdaptor->fetch_by_name($self->input_id);

  my $gene_slice = $seleno_db->get_SliceAdaptor->fetch_by_name($self->input_id);

  print $slice,"\n";
   
  my @genes = @{$gene_slice->get_all_Genes};
 

  ##########################################
  # set up the target (genome)
  ##########################################

  my $target_file = $self->QUERYSEQS."/".$slice->name;
  
  my $targetseqobj = Bio::Seq->new( -display_id =>$slice->name,
                                    -seq => $slice->seq);

  #write_seqfile($seqobj,
  #             create_file_name($dna->display_name(), "fa",
  #             "/tmp/"));


  my $seqout = Bio::SeqIO->new(
                               -file => ">".$target_file,
                               -format => 'fasta',
                               );

  $seqout->write_seq($targetseqobj);

  print "Fetching files:\n ";
  print $target_file,"\n",

  ##########################################
  # setup the runnables
  ##########################################

  my %parameters = %{$self->parameters_hash};
  if (not exists($parameters{-options}) and
      defined $self->OPTIONS) {
    $parameters{-options} = $self->OPTIONS
  }
  if (not exists($parameters{-coverage_by_aligned}) and
      defined $self->COVERAGE_BY_ALIGNED) {
    $parameters{-coverage_by_aligned} = $self->COVERAGE_BY_ALIGNED;
  }

  if (defined $self->PROGRAM && defined $self->analysis->program_file) {
    if ($self->PROGRAM ne $self->analysis->program_file) {
      throw("CONFLICT: You have defined -program in your config file and ".
            "-program_file in your analysis table.");
    }
  }

 # my @runnables;

  foreach my $gene(@genes){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      foreach my $evidence (@{ $transcript->get_all_supporting_features }){
        print "THIS IS YOUR EVIDENCE: ",$evidence->hseqname,"\n";

        my $query_file = $self->QUERYSEQS."/".$evidence->hseqname;

        #print "DOING A NEW EXONERATE ALIGNMENT\n";

        
        my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
            ->new(
                  -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
                  -analysis => $self->analysis,
                  -target_file    => $target_file,
                  -query_type     => $self->QUERYTYPE,
                  -query_file     => $query_file,
                  -annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
                  %parameters,
                  );
        $self->runnable($runnable);
      }
    }
  }

}

############################################################

sub run{
  my ($self) = @_;
  my @results;

  #print "THIS IS YOUR RUNNABLE: ",join(" : ",@{$self->runnable}),"\n";
  
  throw("Can't run - no runnable objects") unless ($self->runnable);
  
  foreach my $runnable (@{$self->runnable}){
    
    $runnable->run;     
    
   # print "Storing output: ",@{$runnable->output},"\n";
    push ( @results, @{$runnable->output} );
    
    
    if ($self->filter) {
      my $filtered_transcripts = $self->filter->filter_results(\@results);
      @results = @$filtered_transcripts;
    }
  }
  my @genes = $self->make_genes(@results);
  print "YOU HAVE ",scalar(@genes)," genes\n";

  $self->output(\@genes);
}


############################################################

sub write_output{
  my ($self,@output) = @_;

  my $outdb = $self->get_output_db;
  my $gene_adaptor = $outdb->get_GeneAdaptor;  

 # unless (@output){
 #   @output = @{$self->output};
 # }
  
  my @genes = $self->filter_redundant_transcript;

  print "this is genes: ",@genes, " with scalar: ",scalar(@genes),"\n";

  my $fails = 0;
  my $total = 0;
  
  foreach my $gene (@genes){
    if($gene == 0){
      print "No genes output so nothing has been written\n"; 
    }else{
      
      print "this is your gene: ", $gene,"\n";
      print "number of transcripts is : ",scalar(@{$gene->get_all_Transcripts}),"\n";
      
      
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
}

############################################################
sub make_genes{
  my ($self,@transcripts) = @_;
  
  my (@genes);

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  
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

    ########################################
    # Add seleno location attrib

    my $sel_location = 0;
    
    my @seq_aas = split(//,$tran->translate->seq);
      
    foreach my $locat (@seq_aas){
      $sel_location++;
      
      if($locat eq "*"){
        
        print "YOU HAVE A SELENO IN LOCATION: ",$sel_location,"\n";
        my $seq_edit = 
            Bio::EnsEMBL::SeqEdit->new(
                                       -CODE    => '_selenocysteine',
                                       -NAME    => 'Selenocysteine',
                                       -DESC    => 'Selenocysteine',
                                       -START   => $sel_location,
                                       -END     => $sel_location,
                                       -ALT_SEQ => 'U'
                                       );

        my $attribute = $seq_edit->get_Attribute();
        
        ##my $translation = $tran->translation();
        
        $tran->add_Attributes($attribute);
        
      }
    }
        
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
  if(defined($self->NONREF_REGIONS)){
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
sub filter_redundant_transcript{
  my ($self) = @_;

  my @all_genes = @{$self->output};

  if (scalar(@all_genes) == 0){
    return 0;
  }

  my @all_transcripts;
  
  my @transcripts;
  
  foreach my $all_gene (@all_genes){
    foreach my $all_transcript(@{$all_gene->get_all_Transcripts}){
      push(@all_transcripts, $all_transcript);
    }
  }

  @all_transcripts = sort { $b->length <=> $a->length} @all_transcripts;

  foreach my $transcript (@all_transcripts){
    my $transcript_exists = 0;

    if(@transcripts && scalar(@transcripts) > 0){
      NEW: foreach my $new_trans(@transcripts){
        my @exons = @{$transcript->get_all_Exons};
        my @new_exons = @{$new_trans->get_all_Exons};
      
        my @c_exons = @{$transcript->get_all_translateable_Exons};
        my @new_c_exons = @{$new_trans->get_all_translateable_Exons};

        print "Number of new coding exons: ",scalar(@new_c_exons),"\n";
        print "Number of test coding exons: ",scalar(@c_exons),"\n";

        # First we check that the transcripts have the exact same coding structure
        next NEW unless (scalar(@c_exons) == scalar(@new_c_exons));
        
        for (my $i = 0; $i < scalar(@c_exons); $i++){
          next NEW unless ($c_exons[$i]->start == $new_c_exons[$i]->start &&
                           $c_exons[$i]->end == $new_c_exons[$i]->end &&
                           $c_exons[$i]->strand == $new_c_exons[$i]->strand);
          
        }

        # Non we want to check that the rest of the non-coding exons are the same apart from the terminal exons

        next NEW unless (scalar(@exons) == scalar(@new_exons));
        print "YOUR EXON STRAND: ",$exons[0]->strand,"\n";

        if($exons[0]->strand == 1){

          next NEW unless(#$exons[0]->start == $new_exons[0]->start &&
                          $exons[0]->end == $new_exons[0]->end &&
                          $exons[0]->strand == $new_exons[0]->strand &&
                          $exons[-1]->start == $new_exons[-1]->start &&
                          #$exons[-1]->end == $new_exons[-1]->end &&
                          $exons[-1]->strand == $new_exons[-1]->strand);
          
          if (scalar(@exons) > 2){
            for (my $i = 1; $i < scalar(@exons)-1; $i++){
              next NEW unless ($exons[$i]->start == $new_exons[$i]->start &&
                               $exons[$i]->end == $new_exons[$i]->end &&
                               $exons[$i]->strand == $new_exons[$i]->strand);
              
            }
          } 
        }else{
          next NEW unless($exons[0]->start == $new_exons[0]->start &&
                          #$exons[0]->end == $new_exons[0]->end &&
                          $exons[0]->strand == $new_exons[0]->strand &&
                          #$exons[-1]->start == $new_exons[-1]->start &&
                          $exons[-1]->end == $new_exons[-1]->end &&
                          $exons[-1]->strand == $new_exons[-1]->strand);
          
          if (scalar(@exons) > 2){
            for (my $i = 1; $i < scalar(@exons)-1; $i++){
              next NEW unless ($exons[$i]->start == $new_exons[$i]->start &&
                               $exons[$i]->end == $new_exons[$i]->end &&
                               $exons[$i]->strand == $new_exons[$i]->strand);
              
            }
          } 
        }
        # If you reach here it means that both your transcripts share the same coding structure
        $transcript_exists = 1;

  
        # keep track of features already transferred, so that we do not duplicate
        #Transfer exons supporting features
        for (my $i = 0; $i < scalar(@exons); $i++){

          my %unique_evidence;
          my %hold_evidence;
          SOURCE_FEAT:
          foreach my $feat ( @{$exons[$i]->get_all_supporting_features}){
            next SOURCE_FEAT unless $feat->isa("Bio::EnsEMBL::FeaturePair");
            
            # skip duplicated evidence objects
            next SOURCE_FEAT if ( $unique_evidence{ $feat } );
            
            # skip duplicated evidence 
            if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } ){
              #print STDERR "Skipping duplicated evidence\n";
              next SOURCE_FEAT;
            }
            
          TARGET_FEAT:
            foreach my $tsf (@{$new_exons[$i]->get_all_supporting_features}){
              next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");
              
              if($feat->start    == $tsf->start &&
                 $feat->end      == $tsf->end &&
                 $feat->strand   == $tsf->strand &&
                 $feat->hseqname eq $tsf->hseqname &&
                 $feat->hstart   == $tsf->hstart &&
                 $feat->hend     == $tsf->hend){
                
                #print STDERR "feature already in target exon\n";
                next SOURCE_FEAT;
              }
            }
            #print STDERR "from ".$source_exon->dbID." to ".$target_exon->dbID."\n";
            #$self->print_FeaturePair($feat);
            # I may need to add a paranoid check to see that no exons longer than the current one are transferred 
            $new_exons[$i]->add_supporting_features($feat);
            $unique_evidence{ $feat } = 1;
            $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } = 1;
          }
        }  

        #Transfer transcript supporting features
        my %t_unique_evidence;
        my %t_hold_evidence;
        T_SOURCE_FEAT:
        foreach my $t_feat ( @{$transcript->get_all_supporting_features}){
          next T_SOURCE_FEAT unless $t_feat->isa("Bio::EnsEMBL::FeaturePair");
          
          # skip duplicated evidence objects
          next T_SOURCE_FEAT if ( $t_unique_evidence{ $t_feat } );
          
          # skip duplicated evidence 
          if ( $t_hold_evidence{ $t_feat->hseqname }{ $t_feat->start }{ $t_feat->end }{ $t_feat->hstart }{ $t_feat->hend } ){
              #print STDERR "Skipping duplicated evidence\n";
            next T_SOURCE_FEAT;
          }
          
          T_TARGET_FEAT:
          foreach my $t_tsf (@{$new_trans->get_all_supporting_features}){
            next T_TARGET_FEAT unless $t_tsf->isa("Bio::EnsEMBL::FeaturePair");
            
            if($t_feat->start    == $t_tsf->start &&
               $t_feat->end      == $t_tsf->end &&
               $t_feat->strand   == $t_tsf->strand &&
               $t_feat->hseqname eq $t_tsf->hseqname &&
               $t_feat->hstart   == $t_tsf->hstart &&
               $t_feat->hend     == $t_tsf->hend){
              
              #print STDERR "feature already in target exon\n";
              next T_SOURCE_FEAT;
            }
          }
          #print STDERR "from ".$source_exon->dbID." to ".$target_exon->dbID."\n";
            #$self->print_FeaturePair($feat);
          # I may need to add a paranoid check to see that no exons longer than the current one are transferred 
          $new_trans->add_supporting_features($t_feat);
          $t_unique_evidence{ $t_feat } = 1;
          $t_hold_evidence{ $t_feat->hseqname }{ $t_feat->start }{ $t_feat->end }{ $t_feat->hstart }{ $t_feat->hend } = 1;
        }
      }
      if($transcript_exists == 0){
        push(@transcripts,$transcript);
      }

    }else{
      push(@transcripts,$transcript);
    }
    
  }
 
  my $gene = new Bio::EnsEMBL::Gene;
  $gene->analysis($self->analysis);
  $gene->biotype($self->analysis->logic_name);

  foreach my $transcript (@transcripts){
    foreach my $se(@{$transcript->get_all_supporting_features}){
      $se->slice($transcript->slice);
     # print "YOU SP LOOKS LIKE: ",join(" - ",%{$se}),"\n";
    }

    #print "Transcript Stable ID: ",$transcript->dbID,"\n";
    $gene->add_Transcript($transcript);
  }
  
  return $gene;
  
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


#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (qw(QUERYSEQS 
                             QUERYTYPE
                             GENOMICSEQS)) {

   throw("You must define $config_var in config for logic '$logic'")
        if not defined $self->$config_var;
  }
  
  throw("QUERYANNOTATION '" . $self->QUERYANNOTATION . "' in config must be readable")
      if $self->QUERYANNOTATION and not -e $self->QUERYANNOTATION;

  # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
  if ($self->FILTER) {
    if (not ref($self->FILTER) eq "HASH" or
        not exists($self->FILTER->{OBJECT}) or
        not exists($self->FILTER->{PARAMETERS})) {
          
      throw("FILTER in config fo '$logic' must be a hash ref with elements:\n" . 
            "  OBJECT : qualified name of the filter module;\n" .
            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
    } else {
      my $module = $self->FILTER->{OBJECT};
      my $pars   = $self->FILTER->{PARAMETERS};
      
      (my $class = $module) =~ s/::/\//g;
      eval{
        require "$class.pm";
      };
      throw("Couldn't require ".$class." Exonerate2Genes:require_module $@") if($@);
    
      $self->filter($module->new(%{$pars}));
    }
  }

}

sub QUERYSEQS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_QUERYSEQS'} = $value;
  }

  if (exists($self->{'_CONFIG_QUERYSEQS'})) {
    return $self->{'_CONFIG_QUERYSEQS'};
  } else {
    return undef;
  }
}

sub QUERYTYPE {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_QUERYTYPE'} = $value;
  }

  if (exists($self->{'_CONFIG_QUERYTYPE'})) {
    return $self->{'_CONFIG_QUERYTYPE'};
  } else {
    return undef;
  }
}


sub QUERYANNOTATION {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_QUERYANNOTATION'} = $value;
  }

  if (exists($self->{'_CONFIG_QUERYANNOTATION'})) {
    return $self->{'_CONFIG_QUERYANNOTATION'};
  } else {
    return undef;
  }
}



sub GENOMICSEQS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_GENOMICSEQS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_GENOMICSEQS'})) {
    return $self->{'_CONFIG_GENOMICSEQS'};
  } else {
    return undef;
  }
}


sub IIDREGEXP {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_IIDREGEXP'} = $value;
  }

  if (exists($self->{'_CONFIG_IIDREGEXP'})) {
    return $self->{'_CONFIG_IIDREGEXP'};
  } else {
    return undef;
  }
}

sub OUTDB {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_OUTDB'} = $value;
  }

  if (exists($self->{'_CONFIG_OUTDB'})) {
    return $self->{'_CONFIG_OUTDB'};

  } else {
    return undef;
  }
}

sub COVERAGE_BY_ALIGNED {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_COVERAGE'} = $value;
  }

  if (exists($self->{'_CONFIG_COVERAGE'})) {
    return $self->{'_CONFIG_COVERAGE'};
  } else {
    return undef;
  }
}


sub FILTER {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_FILTER'} = $value;
  }
  
  if (exists($self->{'_CONFIG_FILTER'})) {
    return $self->{'_CONFIG_FILTER'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if (exists($self->{'_CONFIG_OPTIONS'})) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}

sub NONREF_REGIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_NONREF_REGIONS'} = $value;
  }

  if (exists($self->{'_CONFIG_NONREF_REGIONS'})) {
    return $self->{'_CONFIG_NONREF_REGIONS'};
  } else {
    return undef;
  }
}

sub PROGRAM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PROGRAM'} = $value;
  }

  if (exists($self->{'_CONFIG_PROGRAM'})) {
    return $self->{'_CONFIG_PROGRAM'};
  } else {
    return undef;
  }
}


###############################################
###     end of config
###############################################




1;
