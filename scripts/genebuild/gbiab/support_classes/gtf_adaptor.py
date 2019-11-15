# Copyright [2019] EMBL-European Bioinformatics Institute
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

import os
import re

from gene import Gene
from transcript import Transcript
from exon import Exon
from sequence import Sequence

class GeneAdaptorGTF:

  def __init__(self, gtf_file, sequence_fasta_file=None):

    if not os.path.exists(gtf_file):
      raise Exception("GTF file does not exist on specified path. Path specified:\n" + gtf_file)

    self.gtf_file = gtf_file
    self.sequence_fasta_file = sequence_fasta_file


  def load_all_genes(self, location_name_constraints=None):
    # location_name_constraints will be a list and the code will check if a location name is in the list
    gtf_in = open(self.gtf_file)
    gtf_lines = gtf_in.readlines()

    processed_genes = []
    current_gene_record = []
    current_gene_id = None
    for line in gtf_lines:
      if re.search(r'^#', line):
        continue

      line = line.rstrip()
      elements = re.split(r'\t+', line)
      location_name = elements[0]
      if location_name_constraints is not None and location_name not in location_name_constraints:
        continue

      unparsed_attributes = elements[8]
      parsed_attributes = self.set_attributes(unparsed_attributes)
      del elements[8]

      gene_id = parsed_attributes['gene_id']
      if current_gene_id is None:
        current_gene_id = gene_id

      if gene_id != current_gene_id:
        gene = self.create_gene(current_gene_record)
        if gene:
          processed_genes.append(gene)

        current_gene_record = []
        current_gene_record.append([elements, parsed_attributes])
        current_gene_id = gene_id     
      else:
        current_gene_record.append([elements, parsed_attributes])

    # Process the final record
    gene = self.create_gene(current_gene_record)
    if gene:
      processed_genes.append(gene)

    return processed_genes


  def create_gene(self, gene_record):
    gene_data = None
    transcript_data = {}
    exon_data = []
    cds_data = []

    for record_entry in gene_record:
      unit_type = record_entry[0][2]
      if unit_type == 'gene':
        gene_data = record_entry
      elif unit_type == 'transcript':
        transcript_data[record_entry[1]['transcript_id']] = record_entry
      elif unit_type == 'exon':
        exon_data.append(record_entry)
      elif unit_type == 'CDS':
        cds_data.append(record_entry)

    grouped_exons = self.build_and_group_exons(exon_data)
    grouped_cds_exons = self.build_and_group_exons(cds_data)
    transcripts = []
    gene_id = None
    for transcript_id, exons in grouped_exons.items():
      print("Building transcript: " + transcript_id)
      transcript = Transcript(exons)
      transcript.public_identifier = transcript_id
      transcripts.append(transcript)
      transcript.attributes = transcript_data[transcript_id][1]

      if transcript_id in grouped_cds_exons:
        GeneAdaptorGTF.attach_cds(transcript, grouped_cds_exons[transcript_id])

      # This is a little clunky, but will do for now. Many gtf files don't have gene entries
      # so getting the gene id from the exon is the safest method
      gene_id = exons[0].attributes['gene_id']

    gene = Gene(transcripts)
    gene.public_identifier = gene_id
    if gene_data is not None:
      gene.attributes = gene_data[1]

    return gene


  def build_and_group_exons(self, exon_data):
    grouped_exons = {}
    exon_unit_types = ['exon', 'CDS']
    for record_entry in exon_data:
      elements = record_entry[0]
      attributes = record_entry[1]
      unit_type = elements[2]
      if unit_type not in exon_unit_types:
        continue
      
      location_name = elements[0]
      source = elements[1]
      start = elements[3]
      end = elements[4]
      strand = elements[6]
      if strand not in ['+', '-']:
        strand = '+'
      phase = elements[7]

      transcript_id = attributes['transcript_id']
      if transcript_id not in grouped_exons:
        grouped_exons[transcript_id] = []

      exon = Exon(int(start), int(end), strand, location_name, self.sequence_fasta_file)
      exon.attributes = attributes
      if 'exon_id' in exon.attributes and  exon.attributes['exon_id'] is not None:
        exon.public_identifier = exon.attributes['exon_id']

      grouped_exons[transcript_id].append(exon)

    return grouped_exons


  def attach_cds(transcript, cds_exons):
    # Note this sort doesn't reverse on the negative strand, getting the sequence will be
    # forward strand based regardless, it just makes more sense to sort on ends for the
    # reverse strand in case there are some oddities in the GTF with overlapping exons
    if transcript.strand == '+':
      cds_exons.sort(key=lambda x: x.start)
      transcript.cds_genomic_start = cds_exons[0].start
      transcript.cds_genomic_end = cds_exons[-1].end
    else:
      cds_exons.sort(key=lambda x: x.end)
      transcript.cds_genomic_start = cds_exons[0].start
      transcript.cds_genomic_end = cds_exons[-1].end

# Build a temp transcript with cds exons
# fetch sequence => cds seq
# Transcript should have a cds start/end exon index, cds/translation just goes over that range
# The genomic start/end could be forward strand based since get_sequence would automatically
# reverse. Regardless the genomic start and end reported back should be in logical order I think..
# Not sure how to force translate to compute a direct translation. I assume it can

#    temp_transcript = Transcript(cds_exons)
#    transcript.cds_sequence = temp_transcript.get_sequence()
#    print("GTF CDS: " + transcript.cds_sequence)
#    transcript.translation_sequence = Transcript.local_translate(transcript.cds_sequence.upper())

    # Some code to get stop?
    # One method of getting the stop would be to take the genome end, convert to seq pos and
    # then get the next three bases (if they exist)

  def set_attributes(self, unparsed_attributes):
    final_attribute_dict = {}
    split_attributes = re.split(r';', unparsed_attributes)
    for split_attribute in split_attributes:
      split_attribute = re.sub('^ +', '', split_attribute)
      split_attribute = re.sub(' +$', '', split_attribute)
      attribute_pair = re.split(r' ', split_attribute)
   
      if len(attribute_pair) < 2:
        continue
      attribute_pair[1] = re.sub(r'^"', '', attribute_pair[1])
      attribute_pair[1] = re.sub(r'"$', '', attribute_pair[1])
      final_attribute_dict[attribute_pair[0]] = attribute_pair[1]

    return final_attribute_dict


  @staticmethod
  def write_genes_to_file(genes, output_file):
    gtf_out = open(output_file, "w")
    for gene_idx, gene in  enumerate(genes):
      gene_attributes = {}
      if hasattr(gene, 'attributes'):
        gene_attributes = gene.attributes

      if gene.public_identifier is not None:
        gene_attributes['gene_id'] = gene.public_identifier
      else:
        gene_attributes['gene_id'] = gene.internal_identifier

      gene_attribute_string = GeneAdaptorGTF.create_attribute_string(gene_attributes)
      gtf_out.write(GeneAdaptorGTF.create_feature_entry(gene, 'gene', gene_attribute_string))

      transcripts = gene.transcripts
      for transcript_idx, transcript in enumerate(transcripts):
        transcript_attributes = {}
        if hasattr(transcript, 'attributes'):
          transcript_attributes = transcript.attributes
 
        transcript_attributes['gene_id'] = gene_attributes['gene_id']

        if transcript.public_identifier is not None:
          transcript_attributes['transcript_id'] = transcript.public_identifier
        else:
          transcript_attributes['transcript_id'] = transcript.internal_identifier
    
        transcript_attribute_string = GeneAdaptorGTF.create_attribute_string(transcript_attributes)
        gtf_out.write(GeneAdaptorGTF.create_feature_entry(transcript, 'transcript', transcript_attribute_string))

        exons = transcript.exons
        for exon_idx, exon in enumerate(exons):
          exon_attributes = {}
          if hasattr(exon, 'attributes'):
            exon_attributes = exon.attributes

          exon_attributes['gene_id'] = gene_attributes['gene_id']
          exon_attributes['transcript_id'] = transcript_attributes['transcript_id']
          exon_attributes['exon_number'] = exon_idx + 1

          if exon.public_identifier is not None:
            exon_attributes['exon_id'] = exon.public_identifier
          else:
            exon_attributes['exon_id'] = exon.internal_identifier

          exon_attribute_string = GeneAdaptorGTF.create_attribute_string(exon_attributes)
          gtf_out.write(GeneAdaptorGTF.create_feature_entry(exon, 'exon', exon_attribute_string))

    gtf_out.close()


  @staticmethod
  def create_attribute_string(attributes):
    attribute_string = ''
    for key, val in attributes.items():
      attribute_string = attribute_string + str(key) + ' "' + str(val) + '"; '
      re.sub(r' $', '', attribute_string)
    return attribute_string


  @staticmethod
  def create_feature_entry(feature, unit_type, attribute_string):
    location_name = feature.location_name
    source = 'ensembl'
    if hasattr(feature, source) and feature.source is not None:
      source = feature.source

    start = feature.start
    end = feature.end
    col_5 = '.'
    strand = feature.strand
    phase = '.'
    if hasattr(feature, phase) and feature.phase is not None:
      phase = feature.phase

    feature_set = [location_name, source, unit_type, str(start), str(end), col_5, strand, phase, attribute_string]

    return "\t".join(feature_set) + "\n"
