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


class Sequence:

    bedtools_path = "bedtools"
    samtools_path = "samtools"
    rev_matrix = str.maketrans("atgcATGC", "tacgTACG")

    def __init__(
        self,
        start=None,
        end=None,
        strand=None,
        location_name=None,
        fasta_file=None,
        sequence=None,
        bedtools_path=None,
        samtools_path=None,
    ):
        self.start = start
        self.end = end
        self.strand = strand
        self.location_name = location_name
        self.sequence = sequence
        self.fasta_file = fasta_file
        self.bedtools_path = bedtools_path

    def get_sequence(self, fasta_file=None):
        start = self.start - 1
        end = self.end
        length = end - start + 1
        strand = self.strand
        location_name = self.location_name
        if fasta_file is None:
            fasta_file = self.fasta_file

        bedtools_path = Sequence.bedtools_path

        # This creates a tempfile and writes the bed info to it based on whatever information
        # has been passed in about the sequence. Then runs bedtools getfasta. The fasta file
        # should have a faidx. This can be created with the create_faidx static method prior
        # to fetching sequence
        with tempfile.NamedTemporaryFile(mode="w+t", delete=False) as bed_temp_file:
            bed_temp_file.writelines(location_name + "\t" + str(start) + "\t" + str(end))
            bed_temp_file.close()

            bedtools_command = [
                bedtools_path,
                "getfasta",
                "-fi",
                fasta_file,
                "-bed",
                bed_temp_file.name,
            ]
            bedtools_output = subprocess.Popen(bedtools_command, stdout=subprocess.PIPE)
            for idx, line in enumerate(
                io.TextIOWrapper(bedtools_output.stdout, encoding="utf-8")
            ):
                if idx == 1:
                    if self.strand == "+":
                        self.sequence = line.rstrip()
                    else:
                        self.sequence = Sequence.reverse_complement(line.rstrip())

            os.remove(bed_temp_file.name)
            return self.sequence

    @staticmethod
    def reverse_complement(sequence):
        return sequence.translate(Sequence.rev_matrix)[::-1]

    @staticmethod
    def create_faidx(fasta_file=None, samtools_path=None):
        if samtools_path is None:
            samtools_path = Sequence.samtools_path

        faidx_command = [samtools_path, "faidx", fasta_file]
        subprocess.run(faidx_command)
