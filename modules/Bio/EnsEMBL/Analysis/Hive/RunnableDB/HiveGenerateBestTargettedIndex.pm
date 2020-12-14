=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateBestTargettedIndex

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateBestTargettedIndex;

use strict;
use warnings;

use Bio::EnsEMBL::IO::Parser::Genbank;
use Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::Seq;
use Bio::SeqIO;

use LWP::UserAgent;

use parent qw(Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB);


=head2 fetch_input

 Arg [1]    : None
 Description: Parse the genbank file if present as it contains the cDNAs and
              check that the protein files containing selenoproteins and PE1
              and PE2 proteins from Uniprot. It produces a warning if any of
              the protein files does not exist.
              If there is no cDNAs or proteins, the job finishes early.
 Returntype : None
 Exceptions : Throws if there are cDNAs but no protein files exist

=cut

sub fetch_input {
  my ($self) = @_;

  if (-e $self->param_required('genbank_file')) {
    my $genbank_parser = Bio::EnsEMBL::IO::Parser::Genbank->open($self->param_required('genbank_file'));

    $self->param('seqfetcher', Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher->new( -db => $self->param_required('seqfetcher_index')));
    $self->param_required('fasta_filename');
    my $dba = $self->get_database_by_name('source_db');
    my $adaptor = $dba->get_DnaAlignFeatureAdaptor;
    my %accessions;
    foreach my $f (@{$adaptor->fetch_all}) {
      next if (exists $accessions{$f->hseqname});
      $accessions{$f->hseqname} = 1;
    }
    my @seqs;
    while($genbank_parser->next) {
      my $cdna_accession = $genbank_parser->get_sequence_name;
      next unless (exists $accessions{$cdna_accession});
      foreach my $feature (@{$genbank_parser->get_features}) {
        if ($feature->{header} eq 'CDS' and exists $feature->{translation}) {
          push(@seqs, Bio::Seq->new(-id => $cdna_accession, -desc => $feature->{protein_id}->[0], -seq => join('', @{$feature->{translation}})));
          last;
        }
      }
    }
    $genbank_parser->close;
    if (@seqs) {
      $self->output(\@seqs);
    }
    else {
      $self->warning('Strange, you do not seem to have any full length cdna');
    }
  }
  my $count = 0;
  if ($self->param_is_defined('protein_files')) {
    foreach my $file (@{$self->param('protein_files')}) {
      if (-e $file) {
        ++$count;
      }
      else {
        # It would be better to throw because we expect to have PE1 and PE2 proteins if we have cDNAs.
        # it is possible that there is no selenoprotein for these species as they usually need manual curation.
        # However, it is also possible that there is only PE3 selenoproteins which we still want as there should
        # be better than the Genblast models
        $self->warning('Could not find '.$file);
      }
    }
  }
  if (@{$self->output} and !$count) {
    $self->throw('You have a cDNA file but none of your protein files exist') unless ($count);
  }
  elsif (!@{$self->output} and !$count) {
    $self->input_job->autoflow(0);
    $self->complete_early("Could not find any cdnas");
  }
}


=head2 run

 Arg [1]    : None
 Description: Retrieve the UniProt accession corresponding to the RefSeq/INSDC accession
              in order to avoid creating models from the same sequence
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $seqfetcher = $self->param('seqfetcher');
  my @embl_ids;
  my @refseq_ids;
  foreach my $sequence (@{$self->output}) {
    my $protein_accession = $sequence->desc;
    my $seq = $seqfetcher->get_entry_by_acc($protein_accession);
    if (!$seq) {
      if ($sequence->desc =~ /[AN]P_/) {
        push(@refseq_ids, $sequence);
      }
      else {
        push(@embl_ids, $sequence);
      }
    }
  }
  my $params = {
    to => 'ACC',
    format => 'tab',
    columns => 'id,version(sequence)',
  };
  if (@embl_ids) {
    $params->{from} = 'EMBL';
    $self->_get_uniprot_accession($params, \@embl_ids);
  }
  if (@refseq_ids) {
    $params->{from} = 'P_REFSEQ_AC';
    $self->_get_uniprot_accession($params, \@refseq_ids);
  }
}


=head2 _get_uniprot_accession

 Arg [1]    : Hashref containing the parameter for the REST query to UniProt
 Arg [2]    : Array ref of Bio::Seq
 Description: Query UniProt to find the Uniprot accession of Refseq or INSDC
              protein accessions
              It updates the description field when a UniProt accession is found
 Returntype : None
 Exceptions : Throws if the REST query failed

=cut

sub _get_uniprot_accession {
  my ($self, $params, $seqs) = @_;

  $params->{query} = join(' ', map {$_->desc} @$seqs);
  my $query_url = 'https://www.uniprot.org/uploadlists/';
  my $ua = LWP::UserAgent->new(agent => 'libwwww-perl '.$self->param('email'));
  $ua->env_proxy();
  push(@{$ua->requests_redirectable}, 'POST');
  my %missing;
  my $response = $ua->post($query_url, $params);
  while (my $wait = $response->header('Retry-After')) {
    sleep $wait;
    $response = $ua->get($response->base);
  }
  if ($response->is_success ) {
    if ($response->content =~ /^Entry/) {
      my $result = $response->content;
      while($result =~ /(\w+)\s+(\d+)\s+(\S+)/mgc) {
        foreach my $acc (split(',', $3)) {
          $missing{$acc} = "$1.$2";
        }
      }
    }
    else {
      $response = $ua->post($query_url, $params);
      my $url_edit = $response->request->uri; 
      $url_edit =~ s/https:\/\/www.uniprot.org\/uniprot/https:\/\/www.uniprot.org\/uniparc/; 
      my $new_response = $ua->post($url_edit); 
      while (my $wait = $new_response->header('Retry-After')) {
        sleep $wait;
        $new_response = $ua->get($new_response->base);
      }
      if ($new_response->content =~ /^Entry/) {
        my $result = $new_response->content;
        while($result =~ /(\w+)\s+(\d+)\s+(\S+)/mgc) {
          foreach my $acc (split(',', $3)) {
            $missing{$acc} = "$1.$2";
          }
        }
      }  
      else {
        # I tried hard to get this id but I failed. Check what is the issue. 
      	$self->throw('Cannot find ids. The response is: '. $response->status_line.' for '.$response->request->uri.
      	  ' and the content \n'. $response->content . '\n Also check ' . $url_edit . 'but nothing');
      }
    }
  }
  foreach my $seq (@$seqs) {
    $seq->desc($missing{$seq->desc}) if (exists $missing{$seq->desc});
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Write a new FASTA file containing all possible sequences downloaded previously
              When the seqeunce has multiple accession they are all in the header
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my %seen;
  my $fasta_file = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->param('fasta_filename'));
  foreach my $seq (@{$self->output}) {
    $fasta_file->write_seq($seq);
    $seen{$seq->id} = 1;
  }
  if ($self->param_is_defined('protein_files')) {
    foreach my $file (@{$self->param('protein_files')}) {
      my $protein_file = Bio::SeqIO->new(-format => 'fasta', -file => $file);
      while (my $seq = $protein_file->next_seq) {
        my $id = $seq->id;
        if ($id =~ s/\w{2}\|(\w+)\|\w+.*/$1/) {
          $seq->desc =~ /SV=(\d+)/;
          $id .= '.'.$1;
          $seq->id($id);
          $seq->desc('');
        }
        next if (exists $seen{$seq->id});
        $fasta_file->write_seq($seq);
        $seen{$seq->id} = 1;
      }
    }
  }
  $fasta_file->close;
}

1;
