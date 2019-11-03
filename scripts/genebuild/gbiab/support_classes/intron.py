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

from sequence import Sequence

class Intron:

  canonical_splice_sites = ['GTAG', 'ATAC', 'GCAG']

  def __init__(self, exons, fasta_file=None, internal_identifier=None, public_identifier=None, exon_start_phase=None, exon_end_phase=None):
    self.build_intron(exons)
    self.fasta_file = fasta_file
    self.sequence = None
    self.internal_identifier = internal_identifier
    self.public_identifier = public_identifier 


  def build_intron(self, exons):
    if exons[0].start > exons[1].start:
      exons.sort(key=lambda x: x.start)
      print("Left exon start coord > right exon start coord, will swap")

    exon_left = exons[0]
    exon_right = exons[1]

    self.start = exon_left.end + 1
    self.end = exon_right.start - 1
    self.length = self.end - self.start + 1
    self.strand = exon_left.strand
    self.location_name = exon_left.location_name
    if hasattr(exon_left, 'fasta_file'):
      self.fasta_file = exon_left.fasta_file


  def get_sequence(self):
    if self.sequence is None:
      sequence = Sequence(self.start, self.end, self.strand, self.location_name, self.fasta_file)
      self.sequence = sequence

    if self.sequence.sequence is not None:
      return self.sequence.sequence
    else:
      sequence_string = sequence.get_sequence()
      return sequence_string


  def is_splice_canonical(self):
    sequence = self.get_sequence()
    donor = sequence[:2]
    acceptor = sequence[-2:]
    splice_site = donor + acceptor
    if splice_site in Intron.canonical_splice_sites:
      return True
    else:
      return False


  def intron_string(self, verbose=None):
    start = self.start
    end = self.end
    if self.strand == '-':
      start = self.end
      end = self.start

    intron_string = "<" + str(start) + ".." + str(end) + ">"

    if verbose:
      intron_string = intron_string+ ":" + self.strand + ":" + self.location_name

    return intron_string
