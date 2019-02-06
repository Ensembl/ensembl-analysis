#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion;

use strict;
use warnings;

use File::Fetch;
use File::Path qw(make_path);

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::SeqEdit;

use Bio::EnsEMBL::IO::Parser::Genbank;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               logic_name => 'mt_genbank_import',
               force_loading => 0,
               source => 'RefSeq',
               biotype => 'protein_coding',
               biotype_prefix => 'Mt_',
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    logic_name => 'mt_genbank_import',
    force_loading => 0,
    source => 'RefSeq',
    biotype => 'protein_coding',
    biotype_prefix => 'Mt_',
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Prepare the options for the load mitochondria script
              If 'mt_accession' is set, it will check the name of
              the toplevel coord system and download the genbank file
              Otherwise it completes early
              Creates 'output_path' if it does not exist
              If 'skip_analysis' is set to 1, it will complete early
 Returntype : None
 Exceptions : Throws if 'enscode_root_dir' is not set
              Throws if 'target_db' is not set
              Throws if 'output_path' is not set

=cut

sub fetch_input {
  my $self = shift;

  if($self->param_is_defined('skip_analysis') && $self->param('skip_analysis')) {
    $self->complete_early('The skip_analysis flag has been set to 1, skipping loading');
  }

  my $enscode_dir = $self->param_required('enscode_root_dir');
  my $target_dba = $self->get_database_by_name('target_db');
  my $output_path = $self->param_required('output_path');
  if (!-d $output_path) {
    make_path($output_path);
  }

  my $mt_filename;
  my $mt_accession;
  if ($self->param_is_defined('mt_accession')) {
    chdir($output_path);
    my $fetcher = File::Fetch->new(uri => 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmod=text&rettype=gb&id='.$self->param('mt_accession'));
    $mt_filename = $fetcher->fetch();
    $self->throw('Could not fetch '.$fetcher->output_file.' from '.$fetcher->uri) unless ($mt_filename);
    $mt_accession = $self->param('mt_accession');
  }
  else {
    $self->complete_early('No mitochondria for this species');
  }

  my $genbank_parser = Bio::EnsEMBL::IO::Parser::Genbank->open($mt_filename);
  $genbank_parser->next;
  my $meta_container = $target_dba->get_MetaContainer();
  if ($genbank_parser->get_taxon_id != $meta_container->get_taxonomy_id) {
    $self->throw('Your taxon ids differ: '.$meta_container->get_taxonomy_id.' and '.$genbank_parser->get_taxon_id);
  }
  my $chromosomes = $target_dba->get_SliceAdaptor->fetch_all_karyotype;
  if (@$chromosomes) {
    $self->param('mt_karyotype', $chromosomes->[-1]->get_all_Attributes('karyotype_rank')->[0]+1);
  }
  my $cs_adaptor = $target_dba->get_CoordSystemAdaptor;
  my %required_cs;
  my %slices;
  my $insdc_db_id = $target_dba->get_DBEntryAdaptor->get_external_db_id('INSDC');
  my $refseq_db_id = $target_dba->get_DBEntryAdaptor->get_external_db_id('RefSeq_genomic');
  my @coord_systems = sort {$a->rank <=> $b->rank} @{$cs_adaptor->fetch_all_by_attrib('default_version')};
  foreach my $cs (@coord_systems)  {
    my $name = $mt_accession;
    my @synonyms;
    if ($cs->is_sequence_level) {
      $required_cs{sequence_level} = $cs;
      my $comments = join(' ', @{$genbank_parser->get_raw_comment});
      if ($comments =~ /reference\s+sequence\s+[a-z ]+([A-Z]{1,2}\d+)./) {
        $name = $1;
      }
      else {
        $self->throw("could not get contig accession from $comments");
      }
    }
    if ($cs->rank == 1) {
      $required_cs{toplevel} = $cs;
      if ($cs->name eq 'primary_assembly') {
        push(@synonyms, [$mt_accession, $refseq_db_id]);
        if ($self->param_is_defined('mt_karyotype')) {
          push(@synonyms, [$name, $insdc_db_id]);
          $name = 'MT';
        }
        else {
          push(@synonyms, ['MT']);
        }
      }
      elsif ($cs->name ne 'chromosome') {
        push(@synonyms, ['MT']);
      }
      else {
        $name = 'MT';
      }
    }
    $slices{$cs->name} = Bio::EnsEMBL::Slice->new(
      -coord_system      => $cs,
      -start             => 1,
      -end               => $genbank_parser->get_length,
      -strand            => 1,
      -seq_region_name   => $name,
      -seq_region_length => $genbank_parser->get_length,
    );
    foreach my $synonym (@synonyms) {
      $slices{$cs->name}->add_synonym(@$synonym);
    }
  }
  $self->throw('Could not find a sequence level coordinate system') unless (exists $required_cs{sequence_level});
  $self->throw('Could not find a toplevel coordinate system') unless (exists $required_cs{toplevel});
  my $slice = $target_dba->get_SliceAdaptor->fetch_by_region('toplevel', 'MT');
  if ($slice) {
    my $genes = $slice->get_all_Genes;
    if (@$genes) {
      $self->warning(scalar(@$genes).' are already on the mitochondria');
      if ($self->param('force_loading')) {
        $self->warning('The genes will be deleted');
      }
      else {
#        $self->complete_early('Genes are already loaded');
      }
    }
    if ($slice->length != $genbank_parser->get_length) {
      $self->throw('Length are different between the loaded sequence '.$slice->length.' and the sequence from the file '.$genbank_parser->get_length);
    }
    $self->param('slices', [$slice]);
  }
  else {
    $self->param('slices', [sort {$a->coord_system->rank <=> $b->coord_system->rank} keys %slices]);
  $self->param('parser', $genbank_parser);

  my $analysis = $target_dba->get_AnalysisAdaptor->fetch_by_logic_name($self->param('logic_name'));
  if ($analysis) {
    $self->analysis($analysis);
  }
  else {
    $self->create_analysis(1, {-db_file => $mt_filename});
  }
  $self->hrdb_set_con($target_dba, 'target_db');
}


=head2 run

 Arg [1]    : None
 Description: Execute the load mitochondria script with the parameters set
              in fetch_input like the name of the toplevel coord system
 Returntype : None
 Exceptions : Throws if the command fails

=cut

sub run {
  my $self = shift;

  my @genes;
  my %codon_table;
  my $analysis = $self->analysis;
  my $genbank_parser = $self->param('parser');
  my $slice = $self->param('slices')->[0];
  foreach my $feature (@{$genbank_parser->get_features}) {
    next unless ($feature->{header} eq 'CDS' or $feature->{header} =~ /[tr]RNA/);
    if ($feature->{position} =~ tr/././ > 2) {
      if ($feature->{note}->[0] =~ /frameshift/i) {
        $self->warning('There is a frameshift in '.$feature->{product}->[0].", the gene has more than one exon!\n");
      }
      else {
        $self->throw($feature->{product}->[0]." has more than one exon!\n");
      }
    }
    if (exists $feature->{note}->[0] and $feature->{note}->[0] =~ /tRNAscan-SE/) {
      $self->warning('Skipping '.$feature->{product}->[0].' : '.$feature->{note}->[0]."\n");
      next;
    }
    my $transcript = new Bio::EnsEMBL::Transcript;
    my $exon_number = 0;
    my $start_exon;
    my $end_exon;
    my $strand = 1;
    $strand = -1 if ($feature->{position} =~ /complement/);
    while ($feature->{position} =~ /(\d+)\.\.(\d+)/gc) {
      my $exon = new Bio::EnsEMBL::Exon;
      $exon->start($1);
      $exon->end($2);
      $exon->strand($strand);
      $exon->slice($slice);

      if ($feature->{header} eq 'CDS') {
        $exon->phase(0);
        $exon->end_phase(($exon->end - $exon->start + 1)%3);
      } else {
        $exon->phase(-1);
        $exon->end_phase(-1);
      }

      $transcript->add_Exon($exon);
      $start_exon = $exon if ($exon_number == 0);
      $end_exon = $exon;
      ++$exon_number;
    }
    $transcript->start_Exon($start_exon);
    $transcript->end_Exon  ( $end_exon );
    my $type;
    my $gene = new Bio::EnsEMBL::Gene;
    if ($feature->{header} =~ /(\w)RNA/) {
      if ($1 eq 't' or $1 eq 'r') {
        $type = $self->param('biotype_prefix').$feature->{header};
      }
      else {
        $type = 'UNKNOWN';
        $self->warning('Unknow type for '.$feature->{header}."\n");
      }
    }
    elsif ($feature->{header} eq 'CDS') {
      $codon_table{$feature->{transl_table}->[0]} = $feature->{transl_table}->[0];
      $self->throw('Translation table is 1') if ($feature->{transl_table}->[0] == 1);
      $self->throw('Too many translation table '.$feature->{transl_table}->[0]) if (scalar(keys %codon_table) > 1);

      my $translation = new  Bio::EnsEMBL::Translation(
          -START_EXON => $start_exon,
          -END_EXON   => $end_exon,
          -SEQ_START  => 1,
          -SEQ_END    => $end_exon->length,
          );

      if (exists $feature->{note}->[0] and $feature->{note}->[0] =~ /frameshift/) {
        $self->warning('There is a frameshift in: '.$feature->{gene}->[0]."\nNote: ".$feature->{note}->[0]);
      }
      if (exists $feature->{transl_except}) {
        if ($feature->{transl_except}->[0] =~ /pos:(\d+),aa:(\w+)/) {
          my $alt_seq;
          my $pos = $1;
          if ($2 eq 'TERM') {
            $alt_seq = 'AA' if ($feature->{note}->[0] =~ /TAA/);
            $self->warning("Adding SeqEdit for TAA stop codon completion by 3' AA residues addition.");
            my $seq_edit = Bio::EnsEMBL::SeqEdit->new(
                -CODE    => '_rna_edit',
                -START   => $pos,
                -END     => $pos+1,
                -ALT_SEQ => $alt_seq
                );
            $transcript->add_Attributes($seq_edit->get_Attribute());
            $feature->{position} =~ /(\d+)\.\.(\d+)/;
            my $length = $2-$1+2;
            $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(
              -code => '_transl_end',
              -value => $length,
            ));
          }
        }
      }
      if (exists $feature->{fragment}) {
        $self->warning('The gene '.$feature->{gene}->[0]." is fragmented, a methionine will be added!\n");
      }
      $transcript->translation($translation);
      my $protein = $transcript->translate()->seq();
      if ($protein =~ /\*/) {
        $self->throw("Stop codon found in translation ".$protein);
      }
      if ($protein =~ /^[^M]/) {
        $self->warning("Adding SeqEdit for non-methionine start codon in translation ".$protein);
        my $seqedit = Bio::EnsEMBL::SeqEdit->new(
            -CODE    => 'amino_acid_sub',
            -START   => 1,
            -END     => 1,
            -ALT_SEQ => 'M'
            );
        $transcript->translation()->add_Attributes($seqedit->get_Attribute());
      }
      $type = $self->param('biotype');
    }
    $gene->biotype($type);
    $gene->analysis($analysis);
    $gene->description($feature->{product}->[0]);
    $transcript->biotype($type);
    $transcript->source($self->param('source'));
    $transcript->analysis($analysis);
    $gene->add_Transcript($transcript);
    $gene->canonical_transcript($transcript);
    $gene->source($self->param('source'));
    push(@genes, $gene);
  }
  if (@genes == 37) {
    my ($codon_table_id) = keys %codon_table;
    $self->param('codon_table', $codon_table_id);
    $self->output(\@genes);
  }
  else {
    $self->throw('You have '.scalar(@genes).' which could be a problem');
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Nothing
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  my $target_db = $self->hrdb_get_con('target_db');
  my $slices = $self->param('slices');
  my $toplevel = $slices->[0];
  if (!$toplevel->is_stored) {
    my $slice_adaptor = $target_db->get_SliceAdaptor;
    my $attribute_adaptor = $target_db->get_AttributeAdaptor;
    foreach my $slice (@$slices) {
      $slice_adaptor->store($slice);
      if ($slice->coord_system->rank == 1) {
        $attribute_adaptor->store_on_Slice($slice, Bio::EnsEMBL::Attribute->new(
          -code => 'toplevel',
          -value => 1,
        ));
        if ($self->param_is_defined('karyotype_rank')) {
          $attribute_adaptor->store_on_Slice($slice, Bio::EnsEMBL::Attribute->new(
            -code => 'karyotype_rank',
            -value => $self->param('karyotype_rank'),
          ));
        }
      }
      $attribute_adaptor->store_on_Slice($slice, Bio::EnsEMBL::Attribute->new(
        -code => 'codon_table',
        -value => $self->param('codon_table'),
      ));
    }
  }
  my $gene_adaptor = $target_db->get_GeneAdaptor;
  if ($self->param('force_loading')) {
    foreach my $gene (@{$gene_adaptor->fetch_all_by_Slice($toplevel)}) {
      $gene_adaptor->remove($gene);
    }
  }
  foreach my $gene (@{$self->output}) {
    $gene_adaptor->store($gene);
  }
}

1;
