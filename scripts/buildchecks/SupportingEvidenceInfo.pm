package SupportingEvidenceInfo;


=pod 

=head1 NAME

SupportingEvidenceInfo

=head1 SYNOPSIS

A module to identify supporting evidence associated with transcripts and
genes on the basis of both gene ids and evidence ids

=head1 DESCRIPTION

This module can be given an evidence id (either protein or cdna) or a gene
or transcript id then find either the gene ids which it supports of the
evidence which supports it

=head1 CONTACT

ensembl dev mailing list <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

use vars qw(@ISA $AUTOLOAD);
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw( Bio::EnsEMBL::Root );




=head2 Container methods

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : value of variable
  Function  : container for specified variable. This pod refers to the
  six methods below, db, verbose, info, evidence_type, id_type and 
  have_pfetch. These are simple containers which dont do more than hold and 
  return an given value
  Returntype: value of variable
  Exceptions: some throw if type is incorrect
  Example   : my $db = SupportingEvidenceInfo->db;

=cut


sub db{
  my ($self, $db) = @_;
  if($db){
    if(!($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
      throw("Must pass SupportingEvidenceInfo:db a DBAdaptor not a $db");
    }
    $self->{'db'} = $db;
  }
  return $self->{'db'};
}


sub verbose{
  my ($self, $value) = @_;
  if(defined($value)){
    $self->{'verbose'} = $value;
  }
  return $self->{'verbose'};
}

sub info{
  my ($self, $value) = @_;
  if(defined($value)){
    $self->{'info'} = $value;
  }
  return $self->{'info'};
}

sub evidence_type{
   my ($self, $value) = @_;
  if($value){
    $self->{'evidence_type'} = $value;
  }
  return $self->{'evidence_type'} || 'protein_align_feature';
}

sub id_type{
   my ($self, $value) = @_;
  if($value){
    $self->{'id_type'} = $value;
  }
  return $self->{'id_type'} || 'gene';
}

sub have_pfetch{
  my ($self, $value) = @_;
  if(defined($value)){
    $self->{'have_pfetch'} = $value;
  }
  return $self->{'have_pfetch'};
}




=head2 evidence_ids_from_feature_id

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : int, transcript or gene dbID
  Arg [3]   : string, table name, gene or transcript
  Arg [4]   : string, evidence type dna_align_feature or protein_align_feature
  Arg [5]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : find which evidence ids of the type specified are associated
  with this dbID
  Returntype: none
  Exceptions: none
  Example   : $supportingevidenceinfo->(1);

=cut



sub evidence_ids_from_feature_id{
  my ($self, $transcript_id, $id_type, $evidence_type, $db) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$evidence_type){
    $evidence_type = $self->evidence_type;
  }
  if(!$id_type){
    $id_type = $self->id_type;
  }
  my @protein_ids;
  my $sql = ("select distinct(hit_name) ".
             "from $evidence_type, supporting_feature, ".
             "exon_transcript, exon, transcript ".
             "where ".$evidence_type."_id = feature_id ".
             "and feature_type = '$evidence_type'". 
             "and supporting_feature.exon_id = exon_transcript.exon_id ". 
             "and exon_transcript.transcript_id = ".
             "transcript.transcript_id ".
             "and transcript.".$id_type."_id = $transcript_id");
  print "SQL: ".$sql."\n" if($self->verbose);
  my $sth = $db->dbc->prepare($sql);
  $sth->execute();
  while(my ($protein_id) = $sth->fetchrow){
    print"Have protein id ".$protein_id."\n" if($self->verbose);
    push(@protein_ids, $protein_id);
  }
  $self->fetch_descriptions(\@protein_ids);
}



=head2 feature_ids_from_evidence_id

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : string, evidence identifier
  Arg [3]   : string, gene type
  Arg [4]   : string, table name, gene or transcript
  Arg [5]   : string, evidence type dna_align_feature or 
              protein_align_feature
  Arg [6]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : get the dbIDs of any genes or transcripts which have this
  as part of their supporting evidence
  Returntype: none
  Exceptions: none
  Example   : $supportingevidenceinfo->
  feature_ids_from_evidence_id('Q8K1F0');

=cut



sub feature_ids_from_evidence_id{
  my ($self, $evidence_id, $gene_type, $evidence_type, $id_type, $db) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$evidence_type){
    $evidence_type = $self->evidence_type;
  }
  if(!$id_type){
    $id_type = $self->id_type;
  }
  my $sql = ("SELECT distinct(transcript.".$id_type."_id) ".
             "FROM transcript, exon_transcript, supporting_feature, ".
             "protein_align_feature, gene ".
             "WHERE hit_name = '$evidence_id' ".
             "AND ".$evidence_type."_id = feature_id ".
             "AND feature_type = '$evidence_type' ".
             "AND supporting_feature.exon_id = exon_transcript.exon_id ".
             "AND exon_transcript.transcript_id = ".
             "transcript.transcript_id ".
             "AND transcript.gene_id = gene.gene_id");
  $sql .= " AND type = '".$gene_type."'" if($gene_type);
  print "SQL:".$sql."\n" if($self->verbose);
  my $sth = $db->dbc->prepare($sql);
  $sth->execute();
  while(my ($feature_id) = $sth->fetchrow){
    if($self->info){
      print $evidence_id." ";
      my $adaptor = $self->get_adaptor($id_type, $db);
      my $feature = $adaptor->fetch_by_dbID($feature_id);
      $self->print_feature_info($feature);
    }else{
      print $evidence_id ." ".$feature_id."\n";
    }
  }
}



=head2 fetch_descriptions

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : Arrayref (array of ids for pfetch)
  Function  : fetch descriptions for given identifiers from pfetch
  Returntype: none
  Exceptions: throws if faileds to open or close the pmatch command
  Example   : $self->fetch_descriptions(\@protein_ids);

=cut


sub fetch_descriptions{
  my ($self, $protein_ids) = @_;
  
 ID:foreach my $id(@$protein_ids){
    if($id =~ /(S+)-\1/){
      $id = 1;
    }
    if($self->have_pfetch){
      my $command = "pfetch -D ".$id."\n";
      print $command."\n" if($self->verbose);
      open(FH, $command." | ") or throw("failed to open ".$command);
    LINE:while(<FH>){
        if(/no\s+match/){
          print $id." has no match\n";
          next LINE;
        }
        print;
      }
      close(FH) or throw("Failed to close $command");
    }else{
      print $id."\n";
    }
  }
}



=head2 print_feature_info

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : Bio::EnsEMBL::Feature
  Function  : prints location information about feature give
  Returntype: none
  Exceptions: none
  Example   : $self->print_feature_info($transcript);

=cut



sub print_feature_info{
  my ($self, $feature) = @_;
  throw("Must pass print_feature_info a feature") if(!$feature);
  my $info_string = $feature->dbID." ";
  $info_string .= $feature->stable_id if($feature->stable_id);
  $info_string .= $feature->slice->seq_region_name." ".
    $feature->slice->coord_system_name." ".$feature->start." ".
      $feature->end." ".$feature->strand."\n";
  print $info_string."\n";
}


=head2 get_adaptor

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : string, id type gene or transcript
  Arg [3]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : get an appropriate buisness adaptor for the features
  wanted
  Returntype: Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor 
  Exceptions: throws if not passed a recognised type
  Example   : 

=cut



sub get_adaptor{
  my ($self, $id_type, $db) = @_;
  $db = $self->db if(!$db);
  $id_type = $id_type if(!$id_type);
  if($id_type eq 'gene'){
    return $db->get_GeneAdaptor;
  }elsif($id_type eq 'transcript'){
    return $db->get_TranscriptAdaptor;
  }
  throw("Failed to get adaptor for feature type ".$id_type);
}


=head2 read_id_file

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : string, filename
  Function  : read file to produce a list of ids, takes the first word
  from each line
  Returntype: arrayref
  Exceptions: throws if fails to open or close file
  Example   : my $ids = SupportingEvidenceInfo->read_id_file($file)

=cut



sub read_id_file{
  my ($self, $file) = @_;
  open(FH, $file) or throw("Failed to open ".$file);
  my @ids;
  while(<FH>){
    chomp;
    my ($id) = split;
    push(@ids, $id);
  }
  close(FH) or throw("Failed to close ".$file);
  return \@ids;
}
