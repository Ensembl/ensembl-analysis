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

from exon import Exon

class Transcript:

  def __init__(self, exons, fasta_file=None, internal_identifier=None, public_identifier=None):

    self.exons = exons
    self.build_transcript(exons)
    self.sequence = None
    self.fasta_file = fasta_file      
    self.internal_identifier = internal_identifier
    self.public_identifier = public_identifier 


  def build_transcript(self, exons):
    # Check the integrity of the exons
    strand = exons[0].strand
    location_name = exons[0].location_name
    for exon in exons:
      if exon.strand != strand:
        raise Exception("Inconsistent strands on the exons. Transplicing not supported")
      if exon.location_name != location_name:
        raise Exception("Inconsistent location names for the exons. Exons should belong to the same parent sequence")

    if strand == '+':
      exons.sort(key=lambda x: x.start)
      self.start = exons[0].start
      self.end = exons[-1].end
    else:
      exons.sort(key=lambda x: x.start, reverse=True)
      self.end = exons[0].start
      self.start = exons[-1].end

    self.strand = strand
    self.location_name = location_name


  def add_exons(self, exons):
    # Add a list of exons onto the existing set of exons. Rebuild transcript and
    # set any sequence to None, since that will need to be re-calculated
    self.exons = self.exons + exons
    self.build_transcript(self.exons)
    self.sequence = None


  def get_sequence(self):
    if self.sequence is None:
      sequence = ''
      for exon in self.exons:
        sequence = sequence + exon.get_sequence()

    self.sequence = sequence
    return self.sequence
