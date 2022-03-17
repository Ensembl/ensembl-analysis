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
import copy

from exon import Exon
from intron import Intron


class Transcript:

    internal_identifier = 1
    translate_path = "translate"

    def __init__(
        self,
        exons,
        fasta_file=None,
        internal_identifier=None,
        public_identifier=None,
        translate_path=None,
    ):

        self.exons = exons
        self.build_transcript(exons)
        self.sequence = None
        self.cds_genomic_start = None
        self.cds_genomic_end = None
        self.cds_sequence = None
        self.translation_sequence = None
        self.translations = []
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
                raise Exception(
                    "Inconsistent strands on the exons. Transplicing not supported"
                )
            if exon.location_name != location_name:
                raise Exception(
                    "Inconsistent location names for the exons. Exons should belong to the same parent sequence"
                )

        if strand == "+":
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
            for idx, exon in enumerate(exons[:-1]):
                intron = Intron([exon, exons[idx + 1]])
                introns.append(intron)
            if strand == "+":
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
            sequence = ""
            for exon in self.exons:
                sequence = sequence + exon.get_sequence()
            self.sequence = sequence

        return self.sequence

    def get_cds_sequence(self):
        if (
            self.cds_sequence is None
            and self.cds_genomic_start is not None
            and self.cds_genomic_end is not None
        ):
            self.construct_cds(
                self.cds_genomic_start, self.cds_genomic_end, self.strand, self.exons
            )

        return self.cds_sequence

    def get_translation_sequence(self):
        #    if self.translation_sequence is None and self.cds_genomic_start is not None and self.cds_genomic_end is not None:
        #      self.construct_translation(self.cds_genomic_start, self.cds_genomic_end, self.strand, self.exons)

        cds_sequence = self.get_cds_sequence()
        if cds_sequence is not None:
            self.construct_translation(cds_sequence)

        return self.translation_sequence

    def construct_cds(self, genomic_start, genomic_end, strand, exons):
        # Make a copy of the exons, based on start, then loop over the range of exons that
        # the cds start and end cover. Then edit the boundaries of the start and end exons
        # At the moment I've just made a temp transcript with the cds exons and directly
        # translates that
        forward_sorted_exons = copy.deepcopy(exons)
        forward_sorted_exons.sort(key=lambda x: x.start)
        start_exon_idx = Transcript.get_feature_index(genomic_start, forward_sorted_exons)
        end_exon_idx = Transcript.get_feature_index(genomic_end, forward_sorted_exons)

        if start_exon_idx is None or end_exon_idx is None:
            raise Exception(
                "Start or end exon index was not found based on genomic coordinate, this is wrong"
            )

        cds_exons = []
        for idx in range(start_exon_idx, end_exon_idx + 1):
            exon = forward_sorted_exons[idx]
            if idx == start_exon_idx and start_exon_idx == end_exon_idx:
                exon.sequence = None
                exon.start = genomic_start
                exon.end = genomic_end
                cds_exons.append(exon)
            elif idx == start_exon_idx:
                exon.sequence = None
                exon.start = genomic_start
                cds_exons.append(exon)
            elif idx == end_exon_idx:
                exon.sequence = None
                exon.end = genomic_end
                cds_exons.append(exon)
            else:
                cds_exons.append(exon)

        tmp_transcript = Transcript(cds_exons)
        cds_sequence = tmp_transcript.get_sequence().upper()

        self.cds_sequence = cds_sequence

    def construct_translation(self, cds_sequence):
        # Just does a direct translation of a cds sequence that has already been calculated
        self.translation_sequence = Transcript.local_translate(cds_sequence)

    def compute_translation(self):
        # First remove any existing cds/translation info
        self.cds_sequence = None
        self.cds_genomic_start = None
        self.cds_genomic_end = None
        self.translation_sequence = None
        translations_methonine_required = []
        translations_methonine_not_required = []
        translations_methonine_required = Transcript.run_translate(self.get_sequence(), 1)
        translations_methonine_not_required = Transcript.run_translate(
            self.get_sequence(), 0
        )

        best_translation_met = None
        best_translation_no_met = None
        primary_translation = None

        if len(translations_methonine_required) > 0:
            best_translation_met = translations_methonine_required[0]

        if len(translations_methonine_not_required) > 0:
            best_translation_no_met = translations_methonine_not_required[0]

        if best_translation_met and best_translation_no_met:
            if translations_methonine_required[0][
                2
            ] < 100 and translations_methonine_not_required[0][2] > (
                2 * translations_methonine_required[0][2]
            ):
                primary_translation = best_translation_no_met
            else:
                primary_translation = best_translation_met
        elif best_translation_met:
            primary_translation = best_translation_met
        elif best_translation_no_met:
            primary_translation = best_translation_no_met

        if primary_translation is not None:
            sequence_start = primary_translation[0]
            sequence_end = primary_translation[1]
            if self.strand == "+":
                self.cds_genomic_start = Transcript.sequence_to_genomic_coord(
                    sequence_start, self.exons
                )
                self.cds_genomic_end = Transcript.sequence_to_genomic_coord(
                    sequence_end, self.exons
                )
            else:
                self.cds_genomic_start = Transcript.sequence_to_genomic_coord(
                    sequence_end, self.exons
                )
                self.cds_genomic_end = Transcript.sequence_to_genomic_coord(
                    sequence_start, self.exons
                )

            # Now store the cds and translation seqeunce
            #      self.construct_cds(self.cds_genomic_start, self.cds_genomic_end, self.strand, self.exons)
            self.get_translation_sequence()

    #      self.translation_sequence = primary_translation[3]

    @staticmethod
    def run_translate(sequence, require_methonine=None, min_length=None):

        if require_methonine is None:
            require_methonine = 1
        if min_length is None:
            min_length = 50

        translations = []
        translate_path = Transcript.translate_path
        with tempfile.NamedTemporaryFile(mode="w+t", delete=False) as sequence_temp_file:
            sequence_temp_file.write(">tempseq\n" + sequence + "\n")
            sequence_temp_file.close()

            translate_command = [translate_path]
            if require_methonine:
                translate_command.append("-m")

            if min_length:
                translate_command.append("-l")
                translate_command.append(str(min_length))

            translate_command.append(sequence_temp_file.name)

            translate_output = subprocess.Popen(translate_command, stdout=subprocess.PIPE)

            longest_frame = None
            translate_output_string = ""
            for idx, line in enumerate(
                io.TextIOWrapper(translate_output.stdout, encoding="utf-8")
            ):
                translate_output_string += line

            fasta_regex = "^>.+ nt (\d+)\.\.(\d+)\n(([a-zA-Z\*]+\n)+)"
            match = re.search(fasta_regex, translate_output_string)
            while match is not None:
                start = int(match.group(1))
                end = int(match.group(2))
                translation_sequence = match.group(3)
                frame = start % 3

                translation_sequence = re.sub("\n", "", translation_sequence)

                # in this case the translation is on the opposite strand, so ignore
                if start >= end:
                    translate_output_string = re.sub(
                        fasta_regex, "", translate_output_string
                    )
                    match = re.search(fasta_regex, translate_output_string)
                    continue

                translations.append(
                    [start, end, len(translation_sequence), translation_sequence]
                )
                translate_output_string = re.sub(fasta_regex, "", translate_output_string)
                match = re.search(fasta_regex, translate_output_string)

            os.remove(sequence_temp_file.name)

        return translations

    @staticmethod
    def get_feature_index(genomic_position, features):
        for idx, feature in enumerate(features):
            if genomic_position >= feature.start and genomic_position <= feature.end:
                return idx
        return None

    @staticmethod
    def sequence_to_genomic_coord(
        sequence_position, features, feature_start_offset=None, feature_end_offset=None
    ):
        # This loops through a set features with an associated sequence to place a pair of sequence coords onto the genome
        # A couple of straightforward use cases are converting protein and transcript coords to genomic coords
        # Since the sequence might not cover all features, a feature start and end offset can be provided
        # For example a CDS sequence might start/end in the middle of an exon and thus the offset is needed
        combined_length = 0
        for feature in features:
            next_combined_length = combined_length + (feature.end - feature.start + 1)
            if next_combined_length >= sequence_position:
                remaining_offset = (sequence_position - combined_length) - 1
                if feature.strand == "+":
                    return feature.start + remaining_offset
                else:
                    return feature.end - remaining_offset

            combined_length = next_combined_length

        return None

    @staticmethod
    def local_translate(sequence):
        translation_table = {
            "ATA": "I",
            "ATC": "I",
            "ATT": "I",
            "ATG": "M",
            "ACA": "T",
            "ACC": "T",
            "ACG": "T",
            "ACT": "T",
            "AAC": "N",
            "AAT": "N",
            "AAA": "K",
            "AAG": "K",
            "AGC": "S",
            "AGT": "S",
            "AGA": "R",
            "AGG": "R",
            "CTA": "L",
            "CTC": "L",
            "CTG": "L",
            "CTT": "L",
            "CCA": "P",
            "CCC": "P",
            "CCG": "P",
            "CCT": "P",
            "CAC": "H",
            "CAT": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGA": "R",
            "CGC": "R",
            "CGG": "R",
            "CGT": "R",
            "GTA": "V",
            "GTC": "V",
            "GTG": "V",
            "GTT": "V",
            "GCA": "A",
            "GCC": "A",
            "GCG": "A",
            "GCT": "A",
            "GAC": "D",
            "GAT": "D",
            "GAA": "E",
            "GAG": "E",
            "GGA": "G",
            "GGC": "G",
            "GGG": "G",
            "GGT": "G",
            "TCA": "S",
            "TCC": "S",
            "TCG": "S",
            "TCT": "S",
            "TTC": "F",
            "TTT": "F",
            "TTA": "L",
            "TTG": "L",
            "TAC": "Y",
            "TAT": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGC": "C",
            "TGT": "C",
            "TGA": "*",
            "TGG": "W",
        }

        translation = ""
        if len(sequence) % 3 == 0:
            for i in range(0, len(sequence), 3):
                codon = sequence[i : i + 3]
                translation += translation_table[codon]

        else:
            raise Exception(
                "Sequence passed in for local translation was not zero mod three"
            )

        return translation

    def transcript_string(self, verbose=None):
        transcript_string = (
            "transcript; location='"
            + self.location_name
            + "'; strand='"
            + self.strand
            + "'; structure="
        )
        intron_count = len(self.exons) - 1
        for idx, exon in enumerate(self.exons):
            transcript_string = transcript_string + exon.exon_string()
            if idx < intron_count:
                transcript_string = transcript_string + self.introns[idx].intron_string()

        return transcript_string
