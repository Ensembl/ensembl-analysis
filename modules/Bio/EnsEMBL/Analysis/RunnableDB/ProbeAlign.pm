
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ProbeAlign;




=head1 SYNOPSIS

my $affy = 
  Bio::EnsEMBL::Analysis::RunnableDB::ProbeAlign->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
    -database   => $EST_GENOMIC,
    -query_seqs => \@sequences,
  );

$affy->fetch_input();
$affy->run();
$affy->write_output(); #writes to DB

=head1 DESCRIPTION

This object maps probes to a genome,
and writing the results as Funcgen::ProbeFeatures. You must FIRST
have created and persisted all the necessary Arrays and Probe
objects using the ImportArrays module.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProbeAlign;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBEntry;

use Bio::EnsEMBL::Analysis::Config::ProbeAlign;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

############################################################


#TO DO
#1. We need a sanity check somewhere that we have loaded the transcript coordinate system.

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  
  

  #Because we dont know whether the sort of adapter (compara, hive, core)
  #we'll be passed, it's better to just remake the core dbadaptor and
  #keep it on ourselves as the dna db - this has to be the intent of the
  #dataadaptor input through the rulemanager / beekeeper / whatever

  $self->read_and_check_config;


  #Set up and check DBs before we do the mapping
  $self->outdb->dbc->db_handle;
  $self->outdb->dnadb->dbc->db_handle;

  $self->arrays({});
  $self->probes({});
  #$self->features({});#??
  #$self->{'_unmapped_objects'} = [];

 
  my $logic = $self->analysis->logic_name;
  my $schema_build = $self->outdb->_get_schema_build($self->outdb->dnadb);
  my $species = $self->outdb->species;
  throw('Must provide a -species to the OUTDB in Bio::EnsEMBL::Analysis::Config::ProbeAlign') if ! $species;
	 
  my ($db_name, $display_name);

  if($logic =~ /Transcript/){
	 $self->{'mapping_type'} = 'transcript';
	 
	 #We need to check for ensembl_Core_Transcript external db here
	 #This should also be linked to the release
	 #Should we store here? May over-write each other?
	 #No problem if it fails once as the pipeline will resubmit the job
	 
	 #Will this cause problems with maintaining a master external_dbs file
	 #As this may will grow *3 for every release.
	 
	 #But once we allow GENE, TRANSCRIPT, PROTEIN in xref.info_type
	 #Then this would only be one extra db per release
	 
	 #Would be DBEntry::ADaptor->get_all_external_dbs_by_name

	 #This is tricky as we may be aligning against an old DB
	 #i.e. if there was a new array but not a new build
	 #Current DBEntry queries are naive to external_db version

	 $db_name      = $species.'_core_Transcript';
	 $display_name = 'EnsemblTranscript';
  }
  else{
	 $self->{'mapping_type'} = 'genomic';
	 $db_name      = $species.'_core_Genome';
	 $display_name = 'EnsemblGenome';	 

	 #This is really just a place holder for UnmappedObjects against the entire genome
	 #Should we redo this as EnsemblSpecies and actually have xref entries for each species name?
  }


  my $sql = "select external_db_id from external_db where db_name='$db_name' and db_release='$schema_build'";
  my ($extdb_id) = $self->outdb->db_handle->selectrow_array($sql);
	
 
  if(! $extdb_id){
	warn 'No external_db found for '.$self->{'mapping_type'}." mapping, inserting $db_name $schema_build";
	
	#status is dubious here as this should really be on object_xref
	my $insert_sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type)'.
	  " values('$db_name', '$schema_build', 'KNOWNXREF', 1, 5, '$display_name', 'MISC')";

	$self->outdb->db_handle->do($insert_sql);	
	($extdb_id) = $self->outdb->db_handle->selectrow_array($sql);
  }

  $self->{'_external_db_id'} = $extdb_id;

  return $self;
}


sub mapping_type{
  return $_[0]->{'mapping_type'};
}

sub fetch_input {
  my ($self) = @_;


  ##########################################
  # set up the target (genome)
  ##########################################

  my $target = $self->TARGETSEQS;

  if ( -e $target ){
    if(-d $target ) {
      warn ("Target $target is a directory of files\n");
    }elsif (-s $target){
      warn ("Target $target is a whole-genome file\n");
    }else{
      throw("'$target' isn't a file or a directory?");
    }
  } else {
    throw("'$target' could not be found");
  }

  ##########################################
  # set up the query probe seq
  ##########################################

  my ($query_file, $chunk_number, $chunk_total);

  my $query = $self->QUERYSEQS;

  if ( -e $query and -s $query ) {

    # query seqs is a single file; input id will correspond to a chunk number
    $query_file = $query;
    my $iid_regexp = $self->IIDREGEXP;

    if (not defined $iid_regexp){
      throw("You must define IIDREGEXP in config to enable inference of chunk number and total from your single fasta file" )
    }

    ( $chunk_number, $chunk_total ) = $self->input_id =~ /$iid_regexp/;
    if(!$chunk_number || !$chunk_total){
      throw "I can't make sense of your input id  using the IIDREGEXP in the config!\n";
    }

    #store this for reference later
    $self->query_file($query_file);

  } else {

    throw("'$query'  must refer to a single fasta file with all probe sequences referenced by affy_probe_id\n");

  }

  ##########################################
  # setup the runnable
  ##########################################

  my %parameters = %{ $self->parameters_hash };
 
  #Parameter hash is picked up from the DB analysis table
  #else the config -options value is used
  

  if (not exists $parameters{-options}) {
    if ($self->OPTIONS) {
      $parameters{-options} = $self->OPTIONS;
    } else {
	  throw('You have not options stored in the DB or accessible from the analysis config hash for '.$self->analysis->logic_name);
      #$parameters{-options} = "";
    }
  }
  
  print STDERR "PROGRAM FILE: ".$self->analysis->program_file."\n";


  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe->new
	(
	 -program            => $self->analysis->program_file,
	 -analysis           => $self->analysis,
	 -target_file        => $target,
	 -query_type         => $self->QUERYTYPE,#dna
	 -query_file         => $query_file,
	 -query_chunk_number => $chunk_number,
	 -query_chunk_total  => $chunk_total,
	 -max_mismatches     => $self->MAX_MISMATCHES,
	 -mapping_type       => $self->mapping_type,
	 #-filter_method      => $self->FILTER_METHOD,
	 %parameters,
  );
  
  $self->runnable($runnable);

}

############################################################
#Overrides RunnableDB::run

sub run {
  my ($self) = @_;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my $runnable = @{ $self->runnable }[0];
  
  $runnable->run;

  my $features = $self->filter_features($runnable->output);

  #Why are we setting this in two attrs?
  #These must be used in the rule_manager.pl
  #Which is where write_output must also be called ffrom
  #$self->output($features);#Moved this after filter/set_probe so we get the write output reported
  $self->features($features);
}

############################################################
#overrides RunnableDB::write_output

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb           = $self->outdb;
  my $probe_adaptor   = $outdb->get_ProbeAdaptor;
  my $feature_adaptor = $outdb->get_ProbeFeatureAdaptor;
  my $dbe_adaptor     = $outdb->get_DBEntryAdaptor;
  
  #Add analysis, slices to features, and make
  #sure they're pointing at the persistent array / probe instances
  #instead of the fake arrays & probes they were created with
  my $features = $self->set_probe_and_slice($self->features);

  #Now set correct output
  #This is only used to report the number of features
  $self->output($features);#$self->features);
  
  foreach my $feature_xref(@{$features}){
	my ($feature, $xref) = @$feature_xref;


	if(! defined $feature->slice()){
	  warn "Skipping feature storage for Probe, no slice available for ".$feature->seqname();
	  next;
	}

    eval{ $feature_adaptor->store($feature)};

    if ($@) {

	  warn $feature->slice->name;
      $self->throw('Unable to store ProbeFeature for probe '.$feature->probe_id." !\n $@");
    }

	if($xref){

	  #Remove ignore release flag so we store on the correct schema build
	  eval{ $dbe_adaptor->store($xref, $feature->dbID, 'ProbeFeature') };

	  if ($@) {
		$self->throw('Unable to store ProbeFeature DBEntry for probe '.$feature->dbID." !\n $@");
	  }
	}
  }
}

############################################################

sub filter_features {
  my ($self, $features) = @_;

  my (%hits_by_probe, @kept_hits);
  my $analysis     = $self->analysis;
  my $logic_name   = $analysis->logic_name;
  my $mapping_type = $self->mapping_type;
  my $max_hits     = $self->HIT_SATURATION_LEVEL;#default is 100
  my $uo_adaptor   = $self->outdb->get_UnmappedObjectAdaptor;

  #We don't really want one for every logic name, 
  #so let's simplify this to just ProbeAlign and ProbeTranscriptAlign
  my $uo_type;
  ($uo_type = $logic_name) =~ s/^.*_//;




  #This is tricky
  #We may map something >100 times to the genome but only a few times to transcripts.
  #Do we need a clean up step which removes all transcript mapping if we have an 
  #UnmappedObject due to exceeding the genomic hit saturation???
  #Or keep the transcript mapping and warn in the annotation
  #Do we need a HIT_SATURATION_LEVEL explcitly for transcripts?
  #I think yes and yes, yes we keep and warn and yes we need a more stringent saturation level,
  #probably needs to be number of genes mapped rather than simply number of mappings as we will always
  #get multiple mappings across alt trans
  #In fact this is why we ignored it for transcript mapping.
  #But we are throwing away non-gapped alignments for transcripts!
  #So it is very unlikely that we will ever reach this level!
  #We are also not throwing ungapped away before we do this test
  #So must not filter here as we may throw away a probe which matches many times ungapped
  #but only once gapped...but we're stil matching transcripts rather than the whole genome
  #here so this is still valid.
  

  foreach my $hit (@$features) {


	#we need to capture first an dlast transcript ID's to
	#enable cleaning of potential duplicates

	push @{$hits_by_probe{$hit->probe_id}}, $hit;
  }

  foreach my $probe_id (keys %hits_by_probe) {
    my @hits = @{$hits_by_probe{$probe_id}};
	my $num_hits = scalar(@hits);
	
	#if($logic_name eq 'ProbeAlign'){
	#Only ProbeAlign as we do no use HIT_SATURATION_LEVEL for ProbeTranscriptAlign

	if (scalar(@hits) > $max_hits) {
	  @hits = grep { $_->mismatchcount == 0 } @hits;

	  $num_hits = scalar(@hits);
	  
	  if (scalar(@hits) <= $max_hits) {
		warn "Keeping only perfect matches to $probe_id\n";
		push @kept_hits, @hits;
	  } else {
		warn "Too many hits to $probe_id so rejecting all hits\n";
			
		#So we have issues with what we consider for transcript mapping.
		#If it has failed genomic mapping then we currently do not consider if for xrefing.
		#But it may only map to one transcript, but all over the genome
		#We don't store the genomic mappings, but maybe we should store the transcript mapping which we can then warn against using the unmapped object info from the genomic mapping.
		#Okay, so now we don't have interdepedancy of the mapping steps
		
		

		#We need to grab the external_db_id here if there is only one
		#Leave blank if this exists in multiple DBs e.g. affy 3' UTR designs
		#This assumes that the same ID in two different DBs are the same entity
		#This is true for affy but maybe not for other arrays
		#The simple way to do this is to pull back the array names using the probe_id and direct sql
		#This means yet another call for each unmapped probe
		

		#push @{$self->unmapped_objects},

	
		$uo_adaptor->store(Bio::EnsEMBL::UnmappedObject->new
						   (
							-type       => $uo_type,#Currently get's set to NULL as can only have xref or probe2transcript
							-analysis   => $analysis,
							-ensembl_id => $probe_id,
							-ensembl_object_type => 'Probe',
							-external_db_id => $self->{'_external_db_id'},
							-identifier     => $mapping_type,
							-summary    => 'Promiscuous probe',
							-full_desc  => "Probe exceeded maximum allowed number of $mapping_type mappings(${num_hits}/${max_hits})"
						   ));
	  }
	} 
	#}
	
	#Both ProbeAlign and ProbeTranscriptAlign
	if($num_hits == 0){
 
	  #Very unlikely, and probably only ProbeTranscriptAlign
	  #push @{$self->unmapped_objects},	

	  #We want to an external_db_id here for if this is Trancsript
	  #What about genomic mapping?
	  #Yes we need to reimplement species DB?

	  $uo_adaptor->store(Bio::EnsEMBL::UnmappedObject->new
						 (
						  -type       => $uo_type,#Currently get's set to NULL as can only have xref or probe2transcript
						  -analysis   => $analysis,
						  -ensembl_id => $probe_id,
						  -ensembl_object_type => 'Probe',
						  -external_db_id => $self->{'_external_db_id'},
						  -identifier     => $mapping_type,
						  -summary    => 'Unmapped probe',
						  -full_desc  => "Probe has no $mapping_type mappings",
						 ));
	}
	else {
      push @kept_hits, @hits;
    }
  }

  return \@kept_hits;
}

############################################################

=head2 get_display_name_by_stable_id

  Args [0]   : stable ID from core DB.
  Args [1]   : stable feature type e.g. gene, transcript, translation
  Example    : $self->validate_and_store_feature_types;
  Description: Builds a cache of stable ID to display names.
  Returntype : string - display name
  Exceptions : Throws is type is not valid.
  Caller     : General
  Status     : At risk

=cut

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
#Need to update cache if we're doing more than one 'type' at a time
# as it will never get loaded for the new type!

sub get_display_name_by_stable_id{
  my ($self, $stable_id, $type) = @_;

  $type = lc($type);

  if($type !~ /(gene|transcript|translation)/){
	throw("Cannot get display_name for stable_id $stable_id with type $type");
  }
  
  if(! exists $self->{'display_name_cache'}->{$stable_id}){
	#warn "Generating $type display_name cache\n";
	($self->{'display_name_cache'}->{$stable_id}) = $self->outdb->dnadb->dbc->db_handle->selectrow_array("SELECT x.display_label FROM ${type}_stable_id s, $type t, xref x where t.display_xref_id=x.xref_id and s.${type}_id=t.gene_id and s.stable_id='${stable_id}'");
  }

  return $self->{'display_name_cache'}->{$stable_id};
}


sub set_probe_and_slice {
  my ( $self, $features ) = @_;

  my $db = $self->outdb;
  my $slice_adaptor = $db->get_SliceAdaptor;
  my $probe_adaptor = $db->get_ProbeAdaptor;
  my $trans_adaptor = $db->dnadb->get_TranscriptAdaptor;
  my $mapping_type  = $self->mapping_type;
  my $gene_adaptor  = $db->dnadb->get_GeneAdaptor;
  #my $dbe_adaptor   = $db->get_DBEntryAdaptor;
  my $analysis      = $self->analysis;
  my $schema_build  = $self->outdb->_get_schema_build($self->outdb->dnadb);
  my $edb_name      = $self->outdb->species.'_core_Transcript';

	

  my (%slices, $slice_id, $trans_mapper, $align_type, $align_length, @tmp, $gap_length);
  my (@trans_cigar_line, @genomic_blocks, $genomic_start, $genomic_end, $cigar_line, $gap_start, $block_end);
  my ($block_start, $slice, @gaps, @features, %gene_hits, $gene, @stranded_cigar_line, $gene_sid);
  my ($query_perc, $display_name, $gene_hit_key, %transcript_cache);

  foreach my $feature (@$features) {
    # get the slice based on the seqname stamped on in the runnable
    my $seq_id = $feature->seqname;
	my $probe_id = $feature->probe_id;
	my $load_feature = 1;
	my $xref;

	#Get the slice
	if($mapping_type eq 'transcript'){

	  if(! exists $transcript_cache{$seq_id}){
		$transcript_cache{$seq_id} = $trans_adaptor->fetch_by_stable_id($seq_id);
	  }

	  
	  $slice_id  =  $transcript_cache{$seq_id}->seq_region_name;

	  if ( not exists $slices{$slice_id} ) {
		$slices{$slice_id} = 	$transcript_cache{$seq_id}->slice;
	  }
	}
	else{
	  $slice_id= $seq_id;
	  
	  if ( not exists $slices{$slice_id} ) {
		# assumes genome seqs were named in the Ensembl API Slice naming
		# convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
		$slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
	  }
	}

    $slice = $slices{$slice_id};
	

	#Now project onto genomic coords if we are on cdna
	if($mapping_type eq 'transcript'){
	  #This will basically need to change the start and stops accordingly
	  #and insert D's the size of introns
	  #reject id no D's found as we will already have this as genomic mapping.

	  #The cigar line should always be displayed with reference to the +ve strand of the
	  #target feature i.e. the -ve strand if the transcript is on the -ve strand
	  #These are local stranded start and ends with respect to transcript strand
	  my $transcript_start  = $feature->start;
	  my $transcript_end    = $feature->end;
	  my $transcript_strand = $transcript_cache{$seq_id}->feature_Slice->strand;
	
	  #warn "trans($seq_id) start end strand $transcript_start $transcript_end $transcript_strand" if $warn;


	  $trans_mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript_cache{$seq_id});
	  @genomic_blocks = $trans_mapper->cdna2genomic($transcript_start, $transcript_end);
	

	  #Next feature is this is an ungapped alignment (1 block)
	  #or just representing a flank seq overhang
	  #which will have been caught by the genomic mapping (2 blocks)

	  #print "Found genomic blocks @genomic_blocks\n";

	  next if(! (scalar(@genomic_blocks) >2));
 


	  #So these appear to be 1bp shorter than they should be?
	  #But M add up to 50, is this just and partial intron match?
	  #And also in reverse for -ve strand genes
	  #Are these being reversed incorrectly by the ExonerateProbe
	  #module?
	  #This all depends on what sequence is dumped for the transcript?
	  #Surely this should be feature slice?
	  


	  #else alter the start stop values and rebuild the cigarline
	  my $cigar_line = '';
	  @trans_cigar_line = split/:/, $feature->cigar_line;
		 
	  if($transcript_strand == -1){
		#Always report start end and cigar line wrt the +ve strand of the target seq.
		#Even if the match is against the -ve strand. This is what we do in compara.
		#This is always the genome even for transcript mapping, as we a re projecting back.
		#Need to flip genomic blocks around so that we're still dealing with +ve order
		#This is counterintuitive to how the rest of the API deals with stranded features
		#Normally returns lowest start first wrt to +ve strand even if you use a -ve strand slice, no?
		@genomic_blocks      = reverse(@genomic_blocks);
		@stranded_cigar_line = reverse(@trans_cigar_line);

		#This only makes a difference if there are any m's (seq mismatches)
		#print "Transcript is -ve strand so reverse cigar is @stranded_cigar_line\n";

	  }
	  else{
		@stranded_cigar_line = @trans_cigar_line;
	  }

	  $genomic_start = $genomic_blocks[0]->start;
	  $genomic_end   = $genomic_blocks[$#genomic_blocks]->end;
	  #genomic end is sometimes returning the loci of a gap in transcript coords
	  #This should not be possible for the last block as this would denote a mismatch, but we are
	  #mapping cdna, so it should not have a gap at the end, only in the middle?
	  #Is this because we are actually have a flank seq mismatch
	  #So we actually want the number of genomic blocks to be 3!
	  #warn "genomic $genomic_start $genomic_end $transcript_strand" if $warn;

	  #Only reports match coordinate blocks
	  #And gaps as gaps with reference to transcript
	  #i.e. overhangs of cDNA
	  #Very unlikely to have a Gap at the start or end but we need
	  #to model this in the case of very short exons, where a probe
	  #may span an extire 5' or 3' exon, part of a neighbouring exon 
	  #and also overhang the cDNA
	  
	  #There is no way of knowing whether this overhang is a seq match 
	  #against the corresponding dna seq, so we can't build an accurate 
	  #genomic cigar line for these

	  

	  #Generate start and length values for gaps
	  foreach my $block(@genomic_blocks){
		#warn $block.' '.$block->start.' '.$block->end;#.' '.$block->strand;

		#So calculate gap coordinates
		if(@gaps){
		  push @gaps, ($block->start - $gaps[$#gaps]);
		}

		push @gaps, ($block->end + 1);
	  }

	  #remove last value as this is the end of the match and not a gap start
	  pop @gaps;
	  
	  #Insert intron gaps into cDNA cigarline
	  $gap_start  = shift @gaps;
	  $gap_length = shift @gaps;

	  #Set this to a hypothetical previous block end
	  #So the initial block_start calculation will be valid
	  $block_end  = ($genomic_start - 1);
	  my $last_gap       = 0;	
	  my $match_count    = 0;
	  my $mismatch_count = 0;
	  my $score          = 0;
	  my $q_length       = 0;




	  foreach my $block(@stranded_cigar_line){
		@tmp = split//, $block;
		$align_type = pop @tmp;
		($align_length = $block) =~ s/$align_type//;

		$q_length += $align_length;

		if($align_type eq 'M'){
		  $match_count += $align_length;
		  $score       += ($align_length * 5);
		}
		else{# 'm'  mismatch
		  $mismatch_count += $align_length;
		  $score          -= ($align_length * 4);
		}


		#We need to calculate the genomic end of the current block
		#given the previous block end or the genomic_start
		$block_start = ($block_end + 1);
		$block_end += $align_length;

		if($block_end >= $gap_start){

		  #Could have multiple deletions
		  while(($block_end >= $gap_start) && ! $last_gap){
			#Insert the match first
			$align_length = ($gap_start - $block_start);
			$cigar_line  .= $align_length.$align_type if $align_length;

			#Deletion wrt probe i.e. intron
			$cigar_line .= $gap_length.'D';

			#Now redefine start and end values
			$block_start += $align_length + $gap_length;
			$block_end   += $gap_length;
		
			#Now grab the next gap
			if(@gaps){
			  $gap_start  = shift @gaps;
			  $gap_length =  shift @gaps;
			}
			else{
			  $last_gap = 1;
			}
		  }
		  
		  #We have reached the end of the gaps in this block
		  #so just redefine block here
		  $align_length = ($block_end - $block_start + 1);
		  $block = ($align_length) ? $align_length.$align_type : '';
		}

		$cigar_line .= $block;


		#print "Epxanded cigarline is $cigar_line\n";

	  }

	  #We could assign the start end directly
	  $feature->start($genomic_start);
	  $feature->end($genomic_end);
	  $feature->strand($transcript_strand);
	  $feature->cigar_line($cigar_line);
	  #warn "Final Feature $genomic_start $genomic_end $cigar_line" if $warn;


	  #Test if we have already seen this alignment
	  $gene       = $gene_adaptor->fetch_by_transcript_stable_id($seq_id);
	  #The only way of doing this is to test the genomic_start/end and the genomic cigarline with the gene_stable_id and the probe_id
	  $gene_sid = $gene->stable_id;
	  $gene_hit_key = "${gene_sid}:${probe_id}:${genomic_start}:${genomic_end}:${cigar_line}";


	  if(exists $gene_hits{$gene_hit_key}){
		$load_feature = 0;
	  }else{
		#No need to count hits here
		$gene_hits{$gene_hit_key} = undef;
	  }

	  #Now store the IDXref for this probe transcript hit
	  #This will mean we don't have recalculate this during the probe2transcript mapping step


	  #% ID over aligned region
	  #Which can be different if we have I|Ds
	  #As these will give different seq lengths
	  #But this is the ungapped alignment
	  #So ignore Target ID?
	  #cigar_line is reported wrt to strand of target
	  #we filter out all -ve hits by this point so don't need to account for here.
	  $query_perc = ($match_count/($match_count + $mismatch_count)) * 100;
	  $display_name = $self->get_display_name_by_stable_id($seq_id, 'transcript');

	  #$id_xref = Bio::EnsEMBL::IdentityXref->new
	  $xref = Bio::EnsEMBL::DBEntry->new
		(
		 #-XREF_IDENTITY => $query_perc,
		 #-TARGET_IDENTITY => 90.1,
		 #-SCORE => $score,
		 #-EVALUE => 12,
		 #-CIGAR_LINE => $cigar_line,
		 #-XREF_START => 1,#We are currently padding with mismatches to full length of query
		 #-XREF_END => $q_length,
		 #-ENSEMBL_START => $transcript_start,#target/hit_start
		 #-ENSEMBL_END => $transcript_end,#target/hit_end
		 #-ANALYSIS => $analysis,
		 -PRIMARY_ID => $seq_id,
		 -DISPLAY_ID => $display_name,
		 -DBNAME  => $edb_name,
		 -release => $schema_build,
		 -info_type => 'MISC',
		 -info_text => 'TRANSCRIPT',
		 -linkage_annotation => "ProbeTranscriptAlign",#Add query_perc here when we have analysis
		 #-info_text => , #? What is this for? Is used in unique key so we get duplicated if null!!!
		 #-version => , #version of transcript sid?

		);
	  #No strand here! Always +ve?!
	  #$dbe_adaptor->store($id_xref, $probe_id, 'Probe', 1);#Do we need ignore release flag here?


	  ###This cannot be done until we have ox.analysis_id in place v54?
	  #Yes we can, just use the linkage_annotation for now!


	  #No store this as a ProbeFeature DBEntry in line with how the individual probe xrefs
	  #are stored in probe2transcript
	  #we don't really need this extra translation and score info here
	  #and we are storing the cigar_line as part of the probe_feature
	  #so we can reconstitute the alignment if required

	  #This will mean we can separate the individual feature xrefs from the Probe/ProbeSet level full xref annotation
	  #Will mean some duplication for single probe arrays
	  #But will allow different sets of rules between ProbeAlign and probe2transcript
	  #i.e. we could relax the mapping rules but keep the xreffin rules more stringent

	  #We don't know what the ProbeFeature dbID is yet so we will have to pass this to store with the feature

	  #$dbe_adaptor->store($xref, $probe_id, 'Probe', 1);#Do we need ignore release flag here?
	}
	else{#No introns - reformat cigar line
	  $feature->cigar_line(join('', (split/:/, $feature->cigar_line)));
	}


	if($load_feature){

	  #Reset start ends for non-ref slices

	  if($slice->start != 1){

		#Need to reslice, but we don't want to affect
		#the slice in the cahse as this will screw up
		#further feature start/end tranforms

		my ($level, undef, $name) = split /:/, $slice->name;

		#warn "Resetting slice for ".$slice->name;
		my $start = $feature->start + $slice->start - 1;
		my $end   = $feature->end   + $slice->start - 1;
		$feature->start($start);
		$feature->end($end);
		$slice = $slice_adaptor->fetch_by_region($level, $name, 1, $end);
	  }

	  

	  $feature->slice($slice);

	  # Recover the probe from the cache or DB
	  my $real_probe = $self->probes->{$probe_id};
	  
	  if(!$real_probe){
		$real_probe = $probe_adaptor->fetch_by_dbID($probe_id);
      
		if (!$real_probe){
		  throw "Trying to clean up features for persistence: cant find probe with dbID $probe_id in database!\n";
		}
		
		$self->probes->{$probe_id} = $real_probe;
	  }

	  $feature->probe($real_probe);
	  $feature->analysis($analysis);
	  push @features, [ $feature, $xref ];
	}
  }

  return $self->features(\@features);
}


sub outdb {
  my ($self) = @_;

  my $outdb;
  my $dnadb;

  #Do we need to alter this???????????????????????????????????????????????????????????????????????????????????????????????????
  #Look at ImportArrays
  #There is method duplication here which we could move to a ProbeDB.pm?

  if(! defined $self->{'efg_db'}){

	#This DNADB testing is a work around to avoid having to edit Config::ProbeAlign
	my $dnadb;

	#not defined as an empty env var will give the defined null string

	if($self->DNADB->{-dbname}){
	  $dnadb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(%{ $self->DNADB });
	}

	$self->{'efg_db'} = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
	  (
	   %{ $self->OUTDB }, 
	   -dnadb => $dnadb,
	  );
   
	if(! $self->DNADB->{-dbname}){
	  print "WARNING: Using default DNADB ". $self->{'efg_db'}->dnadb->dbname."\n";
	}

	#We are now forcing use of OUTDB/DNADB config
	#else {
	#  #Historical, but sensible?
	#  $self->{'efg_db'} = $self->SUPER::db;
	#}
  }
  
  return $self->{'efg_db'};
}

sub query_file {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_query_file'} = $value;
  }

  if ( exists( $self->{'_query_file'} ) ) {
    return $self->{'_query_file'};
  } else {
    return undef;
  }
}


#############################################################

sub unmapped_objects{
  my ( $self, $value ) = @_;

  $self->{'_unmnapped_objects'} = $value if defined $value;
  
  return $self->{'_unmapped_objects'};
  
}

#############################################################


#############################################################

sub arrays {
  my ( $self, $value ) = @_;

  $self->{'_arrays'} = $value if defined $value;
  
  return $self->{'_arrays'};
  
}

#############################################################

sub probes {
  my ( $self, $value ) = @_;

  $self->{'_probes'} = $value if  defined $value;
  return $self->{'_probes'};
}

#############################################################

sub features {
  my ( $self, $value ) = @_;

  $self->{'_features'} = $value  if defined $value;
  return $self->{'_features'};#do we need to test and return undef?
}

#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($PROBE_CONFIG);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  #We need to set the QUERYSEQS here dependant on the config
  #We are running all array formats from one cat'd file so we don't need to know any special format specific config

  
  # check that compulsory options have values
  # Remove DNADB??
  # Removed ARRAY_DESIGN

  foreach my $config_var (
						  qw(
							 QUERYSEQS
							 QUERYTYPE
							 TARGETSEQS
							 OUTDB
							 DNADB
							 OPTIONS
							 HIT_SATURATION_LEVEL
							 MAX_MISMATCHES
							)
						 ){
	#We could make MAX_MISMATCHES optional and have % ID here instead?
	# FILTER_METHOD

    if ( not defined $self->$config_var ){
      throw("You must define $config_var in config for logic '$logic'");
    }
  }

  #Test DB config hashes
  for my $db(qw(OUTDB DNADB)){
	
	if (ref( $self->$db ) ne "HASH" || ! defined $self->$db->{-dbname}) {
	  throw("$db in config for '$logic' must be a hash ref of db connection parameters.");
	}
  }

}

#sub ARRAY_DESIGN{
#  my ( $self, $value ) = @_;

#  $self->{'ARRAY_DESIGN'} = $value if defined $value;
  
#  return $self->{'ARRAY_DESIGN'};
#}

#sub FILTER_METHOD{
#  my ( $self, $value ) = @_;

#  $self->{'_CONFIG_FILTER_METHOD'} = $value if defined $value;
  
#  return $self->{'_CONFIG_FILTER_METHOD'} || undef;
#}

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

sub TARGETSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_TARGETSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_TARGETSEQS'} ) ) {
    return $self->{'_CONFIG_TARGETSEQS'};
  } else {
    return undef;
  }
}

sub MAX_MISMATCHES {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_MAX_MISMATCHES'} = $value;
  }

  if ( exists( $self->{'_CONFIG_MAX_MISMATCHES'} ) ) {
    return $self->{'_CONFIG_MAX_MISMATCHES'};
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

sub HIT_SATURATION_LEVEL {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{'_CONFIG_HIT_SATURATION_LEVEL'} = $val;
  }

  if (exists $self->{'_CONFIG_HIT_SATURATION_LEVEL'}) {
    return $self->{'_CONFIG_HIT_SATURATION_LEVEL'};
  } else {
    # default to 100
	# This is not used for ProbeTranscriptAlign
    return 100;
  }
}

###############################################
###     end of config
###############################################

1;
