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


  def create_gene(self, gene_record):
    grouped_exons = self.build_and_group_exons(gene_record)
    transcripts = []
    for transcript_id, exons in grouped_exons.items():
      print("Building transcript: " + transcript_id)
      transcript = Transcript(exons)
      transcripts.append(transcript)

    gene = Gene(transcripts)
    print(gene.gene_string())

    return None


  def build_and_group_exons(self, gene_record):
    grouped_exons = {}
    for record_entry in gene_record:
      elements = record_entry[0]
      attributes = record_entry[1]
      unit_type = elements[2]
      if unit_type != 'exon':
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
      print("Exon: " + exon.exon_string())
      exon.attributes = attributes

      grouped_exons[transcript_id].append(exon)

    return grouped_exons


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
