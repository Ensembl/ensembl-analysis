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

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateForSeleno - 

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

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateForSeleno;

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

  ##########################################
  # set up the target (genome)
  ##########################################

  my $target_file = $self->GENOMICSEQS;
  

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################


  my $query_file = $self->QUERYSEQS."/".$self->input_id;


  my @query_seqs = $self->process_query_file($query_file);
  ##########################################
  # Annotation file with CDS positions
  ##########################################

  open (IN, ">",$self->QUERYANNOTATION);

  my $annotation_file = $self->QUERYANNOTATION;

  print "Fetching files:\n ";
  print $target_file,"\n",$query_file,"\n",$annotation_file,"\n";

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

  foreach my $query ( @query_seqs ){
    #print "DOING A NEW EXONERATE ALIGNMENT\n";

    print "YOUR QUERY IS: ",$query,"\n";

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
        ->new(
              -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_file    => $target_file,
              -query_type     => $self->QUERYTYPE,
              -query_file     => $query,
              #-annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
              %parameters,
              );
$self->runnable($runnable);
   # push (@runnables, $runnable);
  }
  #print "Runnables: ",join(" : ",@runnables),"\n";

 # $self->runnable(@runnables);
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

  unless (@output){
    @output = @{$self->output};
  }
  
  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output){

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

sub process_query_file{

  my ($self, $query_file) = @_;

  my $dnafile = Bio::SeqIO->new(
                                '-file'   => "< $query_file",
                                '-format' => 'embl',
                                );

  my @sequences;

  while((my $dna = $dnafile->next_seq())) {
    #print "DNA: ",$dna,"\n";
    
    my $contains_seleno = 0;
    my $seq;
    
    my $cds_start = 1;
    my $cds_end = 1;
    my %intron_lengths;
    
    my $gene_start = 1; 
    my $gene_end   = $dna->length;
    
    #print "LOOKING AT ENTRY: ",$dna->display_name,"\n";
    
    my $prot_seq = "";
    
    foreach my $feat ($dna->get_SeqFeatures){
      #We have to decompose the sequence and put it back together so we get rid of retained introns.
      
      if ($feat->primary_tag eq "CDS"){
        # We get all the locations which will take into account the retained introns
        my $locations = $feat->location();
        my @locats = $locations->each_Location();
        $cds_start = $locats[0]->start;
       # print "This is where your CDS begins: ",$cds_start,"\n";
       # print "Your primary tag is: ",$feat->primary_tag,"\n";
 
        # First we add the first bit of the non coding sequence
        unless ($cds_start == $gene_start) {$seq .= $dna->subseq($gene_start,$cds_start-1);}
        
        my $prev_locat_end = 0;
        foreach my $locat(@locats){
          # Now we add all the coding sequence
          my $intron_length = $locat->start - ($prev_locat_end+1);
          $prev_locat_end = $locat->end;
          $intron_lengths{$locat->start} = $intron_length;
          #print $locat->start," - ",$locat->end,"\n";
          # print "your Spliced seq is: ",$dna->subseq($locat->start,$locat->end),"\n";
          $seq .= $dna->subseq($locat->start,$locat->end);
          
          # Set the cds end
          if ($locat->end > $cds_end){
            $cds_end = $locat->end;
          }
        }
        #And we add the last bit of non-coding sequence
        unless ($cds_end == $gene_end) {$seq .= $dna->subseq($cds_end + 1,$dna->length);}
        
        #And we recalculate the length of the CDS as we may have removed retained introns
        my $intron;
        foreach my $lengths (keys %intron_lengths){
          $intron += $intron_lengths{$lengths};
        }
        $cds_end -= $intron;
        #print "your seq is: ",$seq,"\n";
      }
      
      foreach my $tag ( $feat->all_tags() ) {
        #print "TAG: ",$tag,"\n";
        if ($tag eq "codon_start"){
          # Some etries have problems as they don't start in the ATG. What should we do?
          foreach my $value ($feat->each_tag_value($tag)){
            #print "VALUE: ",$value,"\n";
            $cds_start += $value-1;
          }
        }
        if ($tag eq "note"){
          foreach my $value ($feat->each_tag_value($tag)){
            if ($value =~ /seleno/){
              $contains_seleno = 1;
             # print "Feature has tag ", $tag, " with values ",$value, "\n";
            }
          }
        }elsif($tag eq "transl_except"){
          foreach my $value ($feat->each_tag_value($tag)){
            #print "Feature has tag ", $tag, " with values ",$value, "\n";
            # This is a double check as sometimes the entry don't have the note with value "seleno"
            if ($value =~ /sec/i){
              $contains_seleno = 1;
            }

          }
        }
      }
    }
    if($contains_seleno == 1){
      my $seqobj = Bio::Seq->new( -display_id =>$dna->display_name() ,
                                  -seq => $seq);

     

      my $sel_output_file = $self->QUERYSEQS."/".$dna->display_name;

      #write_seqfile($seqobj,
      #                               create_file_name($dna->display_name(), "fa",
      #                                                "/tmp/"));


      my $seqout = Bio::SeqIO->new(
                               -file => ">".$sel_output_file,
                               -format => 'fasta',
                                   );

      $seqout->write_seq($seqobj);


      my $annot_output_file = $self->QUERYANNOTATION;
  
      open(ANNOT, ">>$annot_output_file") or die 
          "can open output file: $annot_output_file for writting\n";
      
      #print "ANNOTATION: ",$dna->display_name()," + ",$cds_start," ",$cds_end,"\n";
      
      print ANNOT $dna->display_name()," + ",$cds_start," ",$cds_end,"\n";
      
      close ANNOT;
      
      push (@sequences, $sel_output_file);
       
    }
  }
  return @sequences;
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
