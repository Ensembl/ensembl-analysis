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

from transcript import Transcript


class Gene:

    internal_identifier = 1

    def __init__(
        self,
        transcripts,
        fasta_file=None,
        internal_identifier=None,
        public_identifier=None,
    ):

        self.transcripts = transcripts
        self.build_gene(transcripts)
        self.fasta_file = fasta_file
        if internal_identifier is not None:
            self.internal_identifier = internal_identifier
        else:
            self.internal_identifier = Gene.internal_identifier
            Gene.internal_identifier += 1
        self.public_identifier = public_identifier

    def build_gene(self, transcripts):
        # Check the integrity of the exons
        strand = transcripts[0].strand
        location_name = transcripts[0].location_name
        for transcript in transcripts:
            if transcript.strand != strand:
                raise Exception(
                    "Inconsistent strands on the transcripts. Transcripts should all reside on same strand"
                )
            if transcript.location_name != location_name:
                raise Exception(
                    "Inconsistent location names for the transcripts. Transcripts should belong to the same parent sequence"
                )

        # This is not needed really, but might be useful when clustering and doing thing like that
        if strand == "+":
            transcripts.sort(key=lambda x: x.start)
            self.start = transcripts[0].start
            self.end = transcripts[-1].end
        else:
            transcripts.sort(key=lambda x: x.end, reverse=True)
            self.start = transcripts[-1].start
            self.end = transcripts[0].end

        if self.start >= self.end:
            raise Exception("Gene start was >= end, this should not be")

        self.strand = strand
        self.location_name = location_name

    def add_transcripts(self, transcripts):
        # Add a list of transcript onto the existing set of transcripts.
        # Rebuild gene as a result to ensure everything is consistent
        self.transcripts = self.transcripts + transcripts
        self.build_gene(self.transcripts)

    def gene_string(self):
        start = self.start
        end = self.end
        if self.strand == "-":
            start = self.end
            end = self.start

        gene_string = (
            "gene; start="
            + str(start)
            + "; end="
            + str(end)
            + "; strand='"
            + self.strand
            + "'; location='"
            + self.location_name
            + "';"
        )
        for transcript in self.transcripts:
            gene_string = gene_string + "\n    " + transcript.transcript_string()

        return gene_string
