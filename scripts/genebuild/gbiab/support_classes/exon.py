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


class Exon:

    internal_identifier = 1

    def __init__(
        self,
        start,
        end,
        strand,
        location_name,
        fasta_file=None,
        sequence=None,
        internal_identifier=None,
        public_identifier=None,
        exon_start_phase=None,
        exon_end_phase=None,
    ):
        self.start = start
        self.end = end
        if self.start > self.end:
            raise Exception(
                "Exon start was >= end, this should not be.\nExon start: "
                + str(start)
                + "\nExon end: "
                + str(end)
            )

        self.strand = strand
        self.location_name = location_name
        self.fasta_file = fasta_file
        self.sequence = sequence
        if internal_identifier is not None:
            self.internal_identifier = internal_identifier
        else:
            self.internal_identifier = Exon.internal_identifier
            Exon.internal_identifier += 1
        self.public_identifier = public_identifier
        self.exon_start_phase = exon_start_phase
        self.exon_end_phase = exon_end_phase

    def get_sequence(self):
        if self.sequence is None:
            sequence = Sequence(
                self.start, self.end, self.strand, self.location_name, self.fasta_file
            )
            self.sequence = sequence

        if self.sequence.sequence is not None:
            return self.sequence.sequence
        else:
            sequence_string = sequence.get_sequence()
            return sequence_string

    def exon_string(self, verbose=None):

        start = self.start
        end = self.end
        if self.strand == "-":
            start = self.end
            end = self.start

        exon_string = "(" + str(start) + ".." + str(end) + ")"

        if verbose:
            exon_string = exon_string + ":" + self.strand + ":" + self.location_name

        return exon_string
