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

ensembl dev mailing list <http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

use vars qw(@ISA $AUTOLOAD);
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw( );


sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;
  return $self;
}



=head2 Container methods

  Arg [1]   : SupportingEvidenceInfo
  Arg [2]   : value of variable
  Function  : container for specified variable. This pod refers to the
  seven methods below, db, verbose, info, evidence_type, id_type
  have_pfetch and primary_evidence. These are simple containers which do no
  more than hold and return an given value
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

sub primary_evidence{
  my ($self, $value) = @_;
  if(defined($value)){
    $self->{'primary_evidence'} = $value;
  }
  return $self->{'primary_evidence'};
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
  my ($self, $transcript_id, $id_type, $evidence_type, $primary, $db) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$evidence_type){
    $evidence_type = $self->evidence_type;
  }
  if(!$id_type){
    $id_type = $self->id_type;
  }
  if(!$primary){
    $primary = $self->primary_evidence;
  }
  my @protein_ids;
  my $sql;
  if($primary){
    $sql = $self->primary_evidence_from_feature_id_sql($transcript_id,
                                                       $evidence_type,
                                                       $id_type);
  }else{
    $sql = $self->evidence_from_feature_id_sql($transcript_id,
                                               $evidence_type,
                                               $id_type);
  }
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
  my ($self, $evidence_id, $gene_type, $evidence_type, $id_type,
      $primary, $db) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$evidence_type){
    $evidence_type = $self->evidence_type;
  }
  if(!$id_type){
    $id_type = $self->id_type;
  }
  if(!$primary){
    $primary = $self->primary_evidence;
  }
  my $sql;
  if($primary){
    $sql = $self->feature_id_from_primary_evidence_sql($evidence_id,
                                                       $evidence_type,
                                                       $id_type,
                                                       $gene_type);
  }else{
    $sql = $self->feature_id_from_evidence_sql($evidence_id,
                                               $evidence_type,
                                               $id_type,
                                               $gene_type);
  }
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



sub dbID_from_stable_id{
  my ($self, $stable_id, $id_type, $db) = @_;
  $db = $self->db if(!$db);
  $id_type = $self->id_type if(!$id_type);
  my $sql = "select ".$id_type."_id ".
    "from ".$id_type."_stable_id ".
      "where stable_id = ?";
  print "SQL: $sql \n" if($self->verbose);
  my $sth = $db->dbc->prepare($sql);
  $sth->execute($stable_id);
  my ($dbID) = $sth->fetchrow;
  return $dbID;
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


#Methods to return the SQL


=head2 sql methods

  Arg [1]   : SupportingEvidenceInfo
  Function  : These methods all take arguments to be used to produce
  an sql string to be used in order to fetch either evidence or feature
  ids. What is different is if they fetch on the basis of primary evidence
  associated with whole transcripts or all the evidence associated with each
  exon
  Returntype: string
  Exceptions: none
  Example   :

=cut



sub primary_evidence_from_feature_id_sql{
  my ($self, $feature_id, $evidence_type, $id_type) = @_;
  my $sql = ("SELECT distinct(hit_name) ".
             "FROM transcript, ".
             "transcript_supporting_feature, ".
             "$evidence_type ".
             "WHERE transcript.".$id_type."_id = ".$feature_id." ".
             "and transcript.transcript_id = ".
             "transcript_supporting_feature.transcript_id ".
             "and feature_type = '".$evidence_type."' ".
             "and feature_id = ".$evidence_type."_id ");
  return $sql;
}

sub evidence_from_feature_id_sql{
  my ($self, $feature_id, $evidence_type, $id_type) = @_;
  my $sql = ("select distinct(hit_name) ".
             "from $evidence_type, supporting_feature, ".
             "exon_transcript, exon, transcript ".
             "where ".$evidence_type."_id = feature_id ".
             "and feature_type = '$evidence_type'".
             "and supporting_feature.exon_id = exon_transcript.exon_id ".
             "and exon_transcript.transcript_id = ".
             "transcript.transcript_id ".
             "and transcript.".$id_type."_id = $feature_id");
  return $sql;
}

sub feature_id_from_primary_evidence_sql{
  my ($self, $evidence_id, $evidence_type, $id_type, $gene_type) = @_;
  my $sql = ("SELECT distinct(transcript.".$id_type."_id) ".
             "FROM transcript, transcript_supporting_feature, ".
             "gene, $evidence_type ".
             "WHERE transcript.gene_id = gene.gene_id ".
             "and transcript.transcript_id = ".
             "transcript_supporting_feature.transcript_id ".
             "and feature_type  = '$evidence_type' ".
             "and feature_id = ".$evidence_type."_id ".
             "and hit_name = '$evidence_id' ");
  $sql .= "and biotype = '".$gene_type."' " if($gene_type);
  return $sql;
}

sub feature_id_from_evidence_sql{
  my ($self, $evidence_id, $evidence_type, $id_type, $gene_type) = @_;
  my $sql = ("SELECT distinct(transcript.".$id_type."_id) ".
             "FROM transcript, exon_transcript, supporting_feature, ".
             "$evidence_type, gene ".
             "WHERE hit_name = '$evidence_id' ".
             "AND ".$evidence_type."_id = feature_id ".
             "AND feature_type = '$evidence_type' ".
             "AND supporting_feature.exon_id = exon_transcript.exon_id ".
             "AND exon_transcript.transcript_id = ".
             "transcript.transcript_id ".
             "AND transcript.gene_id = gene.gene_id");
  $sql .= " AND biotype = '".$gene_type."'" if($gene_type);
  return $sql;
}


1;
