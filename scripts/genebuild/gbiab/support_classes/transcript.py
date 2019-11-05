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

import tempfile
import subprocess
import io
import os
import re

from exon import Exon
from intron import Intron

class Transcript:

  internal_identifier = 1
  translate_path = 'translate'

  def __init__(self, exons, fasta_file=None, internal_identifier=None, public_identifier=None, translate_path=None):

    self.exons = exons
    self.build_transcript(exons)
    self.sequence = None
    self.fasta_file = fasta_file      
    if internal_identifier is not None:
      self.internal_identifier = internal_identifier
    else:
      self.internal_identifier = Transcript.internal_identifier
      Transcript.internal_identifier += 1
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
      exons.sort(key=lambda x: x.end, reverse=True)
      self.start = exons[-1].start
      self.end = exons[0].end

    # If we have multiple exons then calculate the introns
    if len(exons) > 1:
      introns = []
      for idx,exon in enumerate(exons[:-1]):
        intron = Intron([exon,exons[idx + 1]])
        introns.append(intron)
      if strand == '+':
        introns.sort(key=lambda x: x.start)
      else:
        introns.sort(key=lambda x: x.end, reverse=True)
      self.introns = introns


    if self.start >= self.end:
      raise Exception("Transcript start was >= end, this should not be")

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


  def compute_translation(self):
    translations_methonine_required = []
    translations_methonine_not_required = []
    translations_methonine_required = self.run_translate(1)
    translations_methonine_not_required = self.run_translate(0)

    if translations_methonine_required[0][2] < 100 and translations_methonine_not_required[0][2] > (2 * translations_methonine_required[0][2]):
      self.translations = translations_methonine_not_required
    else:
      self.translations = translations_methonine_required

    if len(self.translations) > 0:
      self.has_translation = 1
      self.primary_translation = self.translations[0]

    print(self.translations[0][3])


  def run_translate(self, require_methonine):
    translations = []
    translate_path = Transcript.translate_path
    sequence = self.get_sequence()
    with tempfile.NamedTemporaryFile(mode='w+t', delete=False) as sequence_temp_file:
      sequence_temp_file.writelines(">tempseq\n" + sequence + "\n")
      sequence_temp_file.close()

      translate_command = []
      if require_methonine:
        translate_command = [translate_path, '-m', '-l', '50', sequence_temp_file.name]
      else:
        translate_command = [translate_path, '-l', '50', sequence_temp_file.name]

      translate_output = subprocess.Popen(translate_command, stdout=subprocess.PIPE)
      decoded_translate_output = []
      longest_frame = None
      for idx,line in enumerate(io.TextIOWrapper(translate_output.stdout, encoding="utf-8")):
        decoded_translate_output.append(line.rstrip())

      for idx,line in enumerate(decoded_translate_output):
        if not re.search(r'^>', line):
          continue

        match = re.search(r'^>.+ nt (\d+)\.\.(\d+)', line)
        start = int(match.group(1))
        end = int(match.group(2))
        frame = start % 3

        # in this case the translation is on the opposite strand, so ignore
        if start >= end:
          continue

        if longest_frame is None:
          longest_frame = frame

        if frame != longest_frame:
          continue

        translations.append([start, end, len(decoded_translate_output[idx + 1]), decoded_translate_output[idx + 1]])

      os.remove(sequence_temp_file.name)

    return translations


  def transcript_string(self, verbose=None):
    transcript_string = "transcript; location='" + self.location_name + "'; strand='" + self.strand + "'; structure="
    intron_count  = len(self.exons) - 1
    for idx,exon in enumerate(self.exons):
      transcript_string = transcript_string + exon.exon_string()
      if idx < intron_count:
        transcript_string = transcript_string + self.introns[idx].intron_string()

    return transcript_string
