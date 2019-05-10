=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::MakeCdna2GenomeRegionInputIDs - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::MakeCdna2GenomeRegionInputIDs;
use strict;
use warnings;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::MakeCdna2GenomeRegionInputIDs qw(CDNA2GENOME_REGION_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (parse_config);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
#Needed for subroutine &Bio::DB::Flat::OBDAIndex::FileHandle
#called at bioperl-live/Bio/DB/Flat/OBDAIndex.pm
use FileHandle;
use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $submit_logic_name, $gene_db, $pipe_db, $expansion, $annotation, $seqfetcher_object, 
       $seqfetcher_params, $gene_biotypes, $gene_logic_names ) =
    rearrange( [ 'SUBMIT_LOGIC_NAME', 'GENE_DB',
                 'PIPE_DB',           'EXPANSION',
                 'ANNOTATION',        'SEQFETCHER_OBJECT',
                 'SEQFETCHER_PARAMS', 'GENE_BIOTYPES',
                 'GENE_LOGIC_NAMES',
               ], @args );

  # read default config entries
  # actually will read default twice by giving default as logic_name
  # but also does some nice checking and I want the config precedence
  # to be default, param hash, constructor and finally logic_name
  # config over-riding others
  parse_config($self, $CDNA2GENOME_REGION_CONFIG_BY_LOGIC, 'DEFAULT');

  # Defaults are over-ridden by parameters given in analysis table...
  my $ph = $self->parameters_hash;
  $self->SUBMIT_LOGIC_NAME($ph->{-submit_logic_name})   if $ph->{-submit_logic_name};
  $self->GENE_DB($ph->{-gene_db})                       if $ph->{-gene_db};
  $self->PIPE_DB($ph->{-pipe_db})                       if $ph->{-pipe_db};
  $self->EXPANSION($ph->{-expansion})                   if $ph->{-expansion};
  $self->ANNOTATION($ph->{-annotation})                 if $ph->{-annotation};
  $self->SEQFETCHER_OBJECT($ph->{-seqfetcher_object})   if $ph->{-seqfetcher_object};
  $self->SEQFETCHER_PARAMS($ph->{-seqfetcher_params})   if $ph->{-seqfetcher_params};
  $self->GENE_BIOTYPES($ph->{-gene_biotypes})           if $ph->{-gene_biotypes};
  $self->GENE_LOGIC_NAMES($ph->{-gene_logic_names})     if $ph->{-gene_logic_names};

  #...which are over-ridden by constructor arguments. 
  $self->SUBMIT_LOGIC_NAME($submit_logic_name);
  $self->GENE_DB($gene_db);
  $self->PIPE_DB($pipe_db);
  $self->EXPANSION($expansion);
  $self->ANNOTATION($annotation);
  $self->SEQFETCHER_OBJECT($seqfetcher_object);
  $self->SEQFETCHER_PARAMS($seqfetcher_params);
  $self->GENE_BIOTYPES($gene_biotypes);
  $self->GENE_LOGIC_NAMES($gene_logic_names);

  # Finally, analysis specific config
  # use uc as parse_config call above switches logic name to upper case
  $self->SUBMIT_LOGIC_NAME(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{SUBMIT_LOGIC_NAME});
  $self->GENE_DB(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{GENE_DB});
  $self->PIPE_DB(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{PIPE_DB});
  $self->EXPANSION(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{EXPANSION});
  $self->ANNOTATION(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{ANNOTATION});
  $self->SEQFETCHER_OBJECT(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{SEQFETCHER_OBJECT});
  $self->SEQFETCHER_PARAMS(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{SEQFETCHER_PARAMS});
  $self->GENE_BIOTYPES(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{GENE_BIOTYPES});
  $self->GENE_LOGIC_NAMES(${$CDNA2GENOME_REGION_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{GENE_LOGIC_NAMES});

  #throw if something vital is missing
  if ( (    !$self->SUBMIT_LOGIC_NAME
         || !$self->GENE_DB
         || !$self->PIPE_DB
         || !$self->EXPANSION ) )
  {
    throw(   "Need to specify SUBMIT_LOGIC_NAME, GENE_DB, "
           . "PIPE_DB and EXPANSION." );
  }

  print "Annotation file is:".$self->ANNOTATION."\n";

  if ( defined $self->ANNOTATION and not -e $self->ANNOTATION ) {
    throw( "Annotation file " . $self->ANNOTATION . " is not readable." );
  }
  #check that if these options have been defined they are arrays
  if( defined $self->GENE_BIOTYPES and not (ref($self->GENE_BIOTYPES) eq "ARRAY")){
    throw ("Gene biotypes are not in an array.\n");
  }
  if( defined $self->GENE_LOGIC_NAMES and not (ref($self->GENE_LOGIC_NAMES) eq "ARRAY")){
    throw ("Gene logic names are not in an array.\n");
  }


  return $self;
}

sub fetch_input {
  my ($self) = @_;

  # Check submit analysis exists in pipe db, where input ids of that
  # type will be written
  my $dba_pip = $self->get_pipe_dba;
  my $submit_analysis =
    $dba_pip->get_AnalysisAdaptor->fetch_by_logic_name( $self->SUBMIT_LOGIC_NAME );

  if ( !$submit_analysis ) {
    throw(   "Failed to find an analysis with the logic_name "
           . $self->SUBMIT_LOGIC_NAME
           . " in the pipeline db. Cannot continue." );
  }

  $self->submit_analysis($submit_analysis);

  # get genes
  my $gene_dba = $self->get_gene_dba;
  my @genes;
  my @logic_names = @{$self->GENE_LOGIC_NAMES};
  my @biotypes = @{$self->GENE_BIOTYPES};

  #if gene logic names and biotypes have been specified
  if($biotypes[0] and $logic_names[0]){
    foreach my $logic_name(@logic_names){
      foreach my $biotype(@biotypes){
        push (@genes, @{ $gene_dba->get_SliceAdaptor->fetch_by_name( $self->input_id )->
        get_all_Genes($logic_name, undef, undef, undef, $biotype) });
      }
    }
  }#only biotype
  elsif($biotypes[0]){
    foreach my $biotype(@biotypes){
      push (@genes, @{ $gene_dba->get_SliceAdaptor->fetch_by_name( $self->input_id )->
      get_all_Genes(undef, undef, undef, undef, $biotype) });
    }
  }#only logic name
  elsif($logic_names[0]){
    foreach my $logic_name(@logic_names){
      push (@genes, @{ $gene_dba->get_SliceAdaptor->fetch_by_name( $self->input_id )->
      get_all_Genes($logic_name) });
    }
  }#neither specified, so get all
  else{
    @genes = @{ $gene_dba->get_SliceAdaptor->fetch_by_name( $self->input_id )->get_all_Genes() };
  }
  print "There are " . scalar(@genes) . " genes.\n";
  $self->genes( \@genes );
}

sub run {
  my ($self) = @_;

  my $expansion = $self->EXPANSION;
  my @iids;

  GENE: foreach my $gene ( @{ $self->genes } ) {
    # Check that there is only one piece of supporting evidence per gene
    my @transcripts = @{ $gene->get_all_Transcripts };
    throw "Multi-transcript gene." if ( scalar(@transcripts) > 1 );

    my $transcript = $transcripts[0];
    my @evidence   = @{ $transcript->get_all_supporting_features };
    throw "Multiple supporting features for transcript." if ( scalar(@evidence) > 1 );

    # query
    my $hit_name = $evidence[0]->hseqname();

    #if dir, index has been specified and will check that we have the query seq
    #otherwise index hasn't been specified so assume user doesn't want check at this stage
    my %index_params = %{$self->SEQFETCHER_PARAMS};
    if(-d $index_params{-db}->[0]){
      my $seqfetcher = $self->seqfetcher;
      my $query_seq  = $seqfetcher->get_Seq_by_acc($hit_name);
      if(!$query_seq){
        print "WARNING: ".$hit_name."  had no entry in query sequence index.\n";
        next GENE;
      }
    }

    if ( $self->ANNOTATION ) {
      open(F, $self->ANNOTATION) or throw("Could not open supplied annotation file for reading");
      my $unmatched = 1;
      LINE: while (<F>) {
        my @fields = split;
        $unmatched = 0 if $hit_name eq $fields[0];
        last LINE if !$unmatched;
      }
      close(F);
      print "WARNING: ".$hit_name." had no entry in annotation file.\n" if $unmatched;
      next GENE if $unmatched;
    }
    # genomic
    my $genomic_slice = $self->get_gene_dba->get_SliceAdaptor->
    fetch_by_transcript_id($transcript->dbID, $expansion);

    #check expansion hasn't extended slice beyond seq_region
    if(($genomic_slice->start < 1) or ($genomic_slice->end > $genomic_slice->seq_region_length)){

      if(($genomic_slice->start < 1) and ($genomic_slice->end > $genomic_slice->seq_region_length)){
        $genomic_slice = $self->get_gene_dba->get_SliceAdaptor->fetch_by_region($genomic_slice->coord_system_name,$genomic_slice->seq_region_name, 1, $genomic_slice->seq_region_length, $genomic_slice->strand);
      }
      elsif($genomic_slice->start < 1){
        $genomic_slice = $self->get_gene_dba->get_SliceAdaptor->fetch_by_region($genomic_slice->coord_system_name,$genomic_slice->seq_region_name, 1, $genomic_slice->end,$genomic_slice->strand);
      }
      else{
        $genomic_slice = $self->get_gene_dba->get_SliceAdaptor->fetch_by_region($genomic_slice->coord_system_name,$genomic_slice->seq_region_name, $genomic_slice->start, $genomic_slice->seq_region_length,$genomic_slice->strand);
      }

    }

    #use :_: to separate the two parts of the output input ID
    my $id = $genomic_slice->name . ':_:' . $hit_name;

    push @iids, $id;
  }#end GENE loop

  #Get rid of any duplicate ids that may have been generated
  my $num_ids = scalar(@iids);
  my %unique_ids;
  @unique_ids{@iids} = ();
  @iids = keys %unique_ids;
  if(scalar(@iids) < $num_ids){
    print "WARNING: ".scalar(@iids)." unique ids from ".$num_ids." generated ids.\n"
  }
  $self->output( \@iids );
}

sub write_output {
  my ($self) = @_;
  my $sic = $self->get_pipe_dba->get_StateInfoContainer;

  print "output submit analysis is : " . $self->SUBMIT_LOGIC_NAME . "\n";

  foreach my $iid ( @{ $self->output } ) {
    print "try to store input_id : $iid\n";
    eval {
      $sic->store_input_id_analysis( $iid, $self->submit_analysis, '' );
    };
    throw( "Failed to store " . $iid . " $@" ) if ($@);
    logger_info( "Stored " . $iid );
  }
}

sub get_gene_dba {
  my ($self) = @_;
  my $genedba;

  if ( $self->GENE_DB ) {
    if ( ref( $self->GENE_DB ) =~ m/HASH/ ) {
      $genedba = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $self->GENE_DB },
                                                    -dnadb => $self->DNA_DB );
    } else {
      $genedba = $self->get_dbadaptor( $self->GENE_DB );
    }
  }
  return $genedba;
}

sub get_pipe_dba {
  my ($self) = @_;
  my $pipedba;

  if ( $self->PIPE_DB ) {
    if ( ref( $self->PIPE_DB ) =~ m/HASH/ ) {
      $pipedba =
        new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( %{ $self->PIPE_DB },
                                                    -dnadb => $self->DNA_DB );
    } else {
      $pipedba = $self->get_dbadaptor( $self->PIPE_DB, 'pipeline' );
    }
  }
  return $pipedba;
}

sub GENE_BIOTYPES {
  my ($self, $value) = @_;
  if($value){
    $self->{'GENE_BIOTYPES'} = $value;
  }
  return $self->{'GENE_BIOTYPES'};
}

sub GENE_LOGIC_NAMES {
  my ($self, $value) = @_;
  if($value){
    $self->{'GENE_LOGIC_NAMES'} = $value;
  }
  return $self->{'GENE_LOGIC_NAMES'};
}

sub EXPANSION {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'EXPANSION'} = $value;
  }
  return $self->{'EXPANSION'};
}

sub SUBMIT_LOGIC_NAME {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'SUBMIT_LOGIC_NAME'} = $value;
  }
  return $self->{'SUBMIT_LOGIC_NAME'};
}

sub ANNOTATION {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'ANNOTATION'} = $value;
  }
  return $self->{'ANNOTATION'};
}


sub GENE_DB {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'GENE_DB'} = $value;
  }
  return $self->{'GENE_DB'};
}

sub PIPE_DB {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'PIPE_DB'} = $value;
  }
  return $self->{'PIPE_DB'};
}

sub seqfetcher {
  my ( $self, $arg ) = @_;

  if ($arg) {
    throw(   "RunnableDB::Exonerate2GenesRegion " . $arg
           . " must have a method get_Seq_by_acc" )
      unless ( $arg->can("get_Seq_by_acc") );
    $self->{seqfetcher} = $arg;
  }
  if ( !$self->{seqfetcher} ) {
    $self->require_module( $self->SEQFETCHER_OBJECT );
    my %params = %{ $self->SEQFETCHER_PARAMS };
    print $params{-db}->[0], "\n";
    $self->{seqfetcher} = $self->SEQFETCHER_OBJECT->new( %params, );
  }
  return $self->{seqfetcher};
}

sub submit_analysis {
  my ( $self, $value ) = @_;
  if (    $value
       && $value->isa('Bio::EnsEMBL::Pipeline::Analysis') )
  {
    $self->{submit_analysis} = $value;
  }
  return $self->{'submit_analysis'};
}

sub genes {
  my ( $self, $value ) = @_;
  if ( $value && ( ref($value) eq "ARRAY" ) ) {
    $self->{genes} = $value;
  }
  return $self->{'genes'};
}

sub pipeline_adaptor {
  my ( $self, $value ) = @_;
  if (    $value
       && $value->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor') )
  {
    $self->{pipeline_adaptor} = $value;
  }
  return $self->{'pipeline_adaptor'};
}

sub protein_count {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'protein_count'} = $value;
  }
  return $self->{'protein_count'};
}

sub output_logicname {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'output_logicname'} = $value;
  }
  return $self->{'output_logicname'};
}

sub bmg_logicname {
  my ( $self, $value ) = @_;
  if ($value) {
    $self->{'bmg_logicname'} = $value;
  }
  return $self->{'bmg_logicname'};
}

sub paf_slice {
  my ( $self, $slice ) = @_;

  if ($slice) {
    $self->{paf_slice} = $slice;
  }
  if ( !$self->{paf_slice} ) {
    my $slice = $self->fetch_sequence( $self->input_id,
                                       $self->paf_source_db,
                                       $self->REPEATMASKING );
    $self->{paf_slice} = $slice;
  }
  return $self->{paf_slice};
}

sub gene_slice {
  my ( $self, $slice ) = @_;

  if ($slice) {
    $self->{gene_slice} = $slice;
  }
  if ( !$self->{gene_slice} ) {
    my $slice = $self->fetch_sequence( $self->input_id,
                                       $self->gene_source_db,
                                       $self->REPEATMASKING );
    $self->{gene_slice} = $slice;
  }
  return $self->{gene_slice};
}

sub paf_source_db {
  my ( $self, $db ) = @_;
  if ($db) {
    $self->{paf_source_db} = $db;
  }
  if ( !$self->{paf_source_db} ) {
    my $db = $self->get_dbadaptor( $self->PAF_SOURCE_DB );
    $self->{paf_source_db} = $db;
  }
  return $self->{paf_source_db};
}

sub gene_source_db {
  my ( $self, $db ) = @_;
  if ($db) {
    $self->{gene_source_db} = $db;
  }
  if ( !$self->{gene_source_db} ) {
    my $db = $self->get_dbadaptor( $self->GENE_SOURCE_DB );
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

# Note the function of these variables is better described in the
# config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise

sub PAF_LOGICNAMES {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{PAF_LOGICNAMES} = $arg;
  }
  return $self->{PAF_LOGICNAMES};
}

sub PAF_MIN_SCORE_THRESHOLD {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{PAF_MIN_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_MIN_SCORE_THRESHOLD};
}

sub PAF_UPPER_SCORE_THRESHOLD {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{PAF_UPPER_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_UPPER_SCORE_THRESHOLD};
}

sub SEQFETCHER_OBJECT {
  my ($self, $value) = @_;
  if($value){
    $self->{'SEQFETCHER_OBJECT'} = $value;
  }
  return $self->{'SEQFETCHER_OBJECT'};
}

sub SEQFETCHER_PARAMS {
  my ($self, $value) = @_;
  if($value){
    $self->{'SEQFETCHER_PARAMS'} = $value;
  }
  return $self->{'SEQFETCHER_PARAMS'};
}

sub PAF_SOURCE_DB {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{PAF_SOURCE_DB} = $arg;
  }
  return $self->{PAF_SOURCE_DB};
}

sub GENE_SOURCE_DB {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{GENE_SOURCE_DB} = $arg;
  }
  return $self->{GENE_SOURCE_DB};
}

sub OUTPUT_DB {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{OUTPUT_DB} = $arg;
  }
  return $self->{OUTPUT_DB};
}

sub OUTPUT_BIOTYPE {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{OUTPUT_BIOTYPE} = $arg;
  }
  return $self->{OUTPUT_BIOTYPE};
}

sub GENEWISE_PARAMETERS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{GENEWISE_PARAMETERS} = $arg;
  }
  return $self->{GENEWISE_PARAMETERS};
}

sub MINIGENEWISE_PARAMETERS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{MINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MINIGENEWISE_PARAMETERS};
}

sub MULTIMINIGENEWISE_PARAMETERS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{MULTIMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MULTIMINIGENEWISE_PARAMETERS};
}

sub BLASTMINIGENEWISE_PARAMETERS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{BLASTMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{BLASTMINIGENEWISE_PARAMETERS};
}

sub FILTER_PARAMS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{FILTER_PARAMETERS} = $arg;
  }
  return $self->{FILTER_PARAMETERS};
}

sub FILTER_OBJECT {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{FILTER_OBJECT} = $arg;
  }
  return $self->{FILTER_OBJECT};
}

sub BIOTYPES_TO_MASK {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{BIOTYPES_TO_MASK} = $arg;
  }
  return $self->{BIOTYPES_TO_MASK};
}

sub EXON_BASED_MASKING {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{EXON_BASED_MASKING} = $arg;
  }
  return $self->{EXON_BASED_MASKING};
}

sub GENE_BASED_MASKING {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{GENE_BASED_MASKING} = $arg;
  }
  return $self->{GENE_BASED_MASKING};
}

sub POST_GENEWISE_MASK {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{POST_GENEWISE_MASK} = $arg;
  }
  return $self->{POST_GENEWISE_MASK};
}

sub PRE_GENEWISE_MASK {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{PRE_GENEWISE_MASK} = $arg;
  }
  return $self->{PRE_GENEWISE_MASK};
}

sub REPEATMASKING {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{REPEATMASKING} = $arg;
  }
  return $self->{REPEATMASKING};
}

sub USE_KILL_LIST {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{USE_KILL_LIST} = $arg;
  }
  return $self->{USE_KILL_LIST};
}

sub LIMIT_TO_FEATURE_RANGE {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{LIMIT_TO_FEATURE_RANGE} = $arg;
  }
  return $self->{LIMIT_TO_FEATURE_RANGE};
}

sub FEATURE_RANGE_PADDING {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{FEATURE_RANGE_PADDING} = $arg;
  }
  return $self->{FEATURE_RANGE_PADDING};
}

sub WRITE_REJECTED {
  my ( $self, $arg ) = @_;
  if ( defined($arg) ) {
    $self->{WRITE_REJECTED} = $arg;
  }
  return $self->{WRITE_REJECTED};
}

sub REJECTED_BIOTYPE {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{REJECTED_BIOTYPE} = $arg;
  }
  return $self->{REJECTED_BIOTYPE};
}

sub SOFTMASKING {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{SOFTMASKING} = $arg;
  }
  return $self->{SOFTMASKING};
}

sub EXONERATE_PARAMETERS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{EXONERATE_PARAMETERS} = $arg;
  }
  return $self->{EXONERATE_PARAMETERS};
}

sub MAKE_SIMGW_INPUT_ID_PARMAMS {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{MAKE_SIMGW_INPUT_ID_PARMAMS} = $arg;
  }
  return $self->{MAKE_SIMGW_INPUT_ID_PARMAMS};
}
1;
