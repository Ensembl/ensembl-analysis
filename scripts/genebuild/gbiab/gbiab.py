# Copyright [2019-2021] EMBL-European Bioinformatics Institute
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

import argparse
import os
import errno
import shutil
import subprocess
import glob
import re
import multiprocessing
import random
import tempfile
import io
import gc
import logging
import math
from pathlib import Path

# filename=os.environ[work_dir] + "myLogFile.log"
# logging.basicConfig(filename, filemode='w', format='%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s',datefmt='%H:%M:%S', level=logging.DEBUG)
# logging.warning('This will get logged to a file')


def create_dir(main_output_dir, dir_name):

    if dir_name:
        target_dir = os.path.join(main_output_dir, dir_name)
    else:
        target_dir = main_output_dir

    if os.path.exists(target_dir):
        logging.warning("Directory already exists, will not create again")
        return target_dir

    logging.info("Attempting to create target dir: %s" % target_dir)

    try:
        os.mkdir(target_dir)

    except OSError:
        logging.error("Creation of the dir failed, path used: %s" % target_dir)
    else:
        logging.info(
            "Successfully created the dir on the following path: %s" % target_dir
        )

    return target_dir


def load_results_to_ensembl_db(
    main_script_dir,
    load_to_ensembl_db,
    genome_file,
    main_output_dir,
    db_details,
    num_threads,
):

    db_loading_script = os.path.join(
        main_script_dir, "support_scripts_perl", "load_gtf_ensembl.pl"
    )
    db_loading_dir = create_dir(main_output_dir, "db_loading")

    # Should collapse this into a function
    annotation_results_gtf_file = os.path.join(
        main_output_dir, "annotation_output", "annotation.gtf"
    )
    if os.path.exists(annotation_results_gtf_file):
        logging.info("Loading main geneset to db")
        batch_size = 200
        load_type = "gene"
        analysis_name = "ensembl"
        gtf_records = batch_gtf_records(
            annotation_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find the main gene annotation file, so not loading. Path checked:\n"
            + annotation_results_gtf_file
        )

    rfam_results_gtf_file = os.path.join(main_output_dir, "rfam_output", "annotation.gtf")
    if os.path.exists(rfam_results_gtf_file):
        logging.info("Loading Rfam-based sncRNA genes to db")
        batch_size = 500
        load_type = "gene"
        analysis_name = "ncrna"
        gtf_records = batch_gtf_records(
            rfam_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find an Rfam annotation file, so not loading. Path checked:\n"
            + rfam_results_gtf_file
        )

    trnascan_results_gtf_file = os.path.join(
        main_output_dir, "trnascan_output", "annotation.gtf"
    )
    if os.path.exists(trnascan_results_gtf_file):
        logging.info("Loading tRNAScan-SE tRNA genes to db")
        batch_size = 500
        load_type = "gene"
        analysis_name = "ncrna"
        gtf_records = batch_gtf_records(
            trnascan_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find an tRNAScan-SE annotation file, so not loading. Path checked:\n"
            + trnascan_results_gtf_file
        )

    dust_results_gtf_file = os.path.join(main_output_dir, "dust_output", "annotation.gtf")
    if os.path.exists(dust_results_gtf_file):
        logging.info("Loading Dust repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "dust"
        gtf_records = batch_gtf_records(
            dust_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find a Dust annotation file, so not loading. Path checked:\n"
            + dust_results_gtf_file
        )

    red_results_gtf_file = os.path.join(main_output_dir, "red_output", "annotation.gtf")
    if os.path.exists(red_results_gtf_file):
        logging.info("Loading Red repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "repeatdetector"
        gtf_records = batch_gtf_records(
            red_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find a Red annotation file, so not loading. Path checked:\n"
            + red_results_gtf_file
        )

    trf_results_gtf_file = os.path.join(main_output_dir, "trf_output", "annotation.gtf")
    if os.path.exists(trf_results_gtf_file):
        logging.info("Loading TRF repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "trf"
        gtf_records = batch_gtf_records(
            trf_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find a TRF annotation file, so not loading. Path checked:\n"
            + trf_results_gtf_file
        )

    cpg_results_gtf_file = os.path.join(main_output_dir, "cpg_output", "annotation.gtf")
    if os.path.exists(cpg_results_gtf_file):
        logging.info("Loading CpG islands to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "cpg"
        gtf_records = batch_gtf_records(
            cpg_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find a CpG annotation file, so not loading. Path checked:\n"
            + cpg_results_gtf_file
        )

    eponine_results_gtf_file = os.path.join(
        main_output_dir, "eponine_output", "annotation.gtf"
    )
    if os.path.exists(eponine_results_gtf_file):
        logging.info("Loading Eponine repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "eponine"
        gtf_records = batch_gtf_records(
            eponine_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logging.error(
            "Did not find an Eponine annotation file, so not loading. Path checked:\n"
            + eponine_results_gtf_file
        )

    logging.info("Finished loading records to db")


def generic_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file,
    db_details,
    db_loading_dir,
    load_type,
    analysis_name,
    gtf_records,
    num_threads,
):

    pool = multiprocessing.Pool(int(num_threads))
    for record_batch in gtf_records:
        pool.apply_async(
            multiprocess_load_records_to_ensembl_db,
            args=(
                load_to_ensembl_db,
                db_loading_script,
                genome_file,
                db_details,
                db_loading_dir,
                load_type,
                analysis_name,
                record_batch,
            ),
        )

    pool.close()
    pool.join()


def multiprocess_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file,
    db_details,
    output_dir,
    load_type,
    analysis_name,
    record_batch,
):

    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as gtf_temp_out:
        for line in record_batch:
            gtf_temp_out.writelines(line)
            gtf_temp_file_path = gtf_temp_out.name

    (db_name, db_host, db_port, db_user, db_pass) = db_details.split(",")

    loading_cmd = [
        "perl",
        db_loading_script,
        "-genome_file",
        genome_file,
        "-dbname",
        db_name,
        "-host",
        db_host,
        "-port",
        str(db_port),
        "-user",
        db_user,
        "-pass",
        db_pass,
        "-gtf_file",
        gtf_temp_file_path,
        "-analysis_name",
        analysis_name,
        "-load_type",
        load_type,
    ]

    if load_type == "gene" and analysis_name == "ensembl":
        loading_cmd.extend(
            [
                "-protein_coding_biotype",
                "gbiab_protein_coding",
                "-non_coding_biotype",
                "gbiab_lncRNA",
            ]
        )

        if load_to_ensembl_db == "single_transcript_genes":
            loading_cmd.append("-make_single_transcript_genes")

    logging.info(" ".join(loading_cmd))
    subprocess.run(loading_cmd)
    gtf_temp_out.close()
    os.remove(gtf_temp_file_path)  # doesn't seem to be working
    logging.info("Finished: " + gtf_temp_file_path)
    gc.collect()


def batch_gtf_records(input_gtf_file, batch_size, output_dir, record_type):

    records = []
    gtf_in = open(input_gtf_file)
    if record_type == "gene":

        # NOTE that the neverending variations on GTF reading/writing/merging is becoming very messy
        # need to create a set of utility methods outside of this script
        # This one assumes the file has unique ids for the parent features. It then batches them into
        # sets of records based on the batch size passed in
        record_counter = 0
        current_record_batch = []
        current_gene_id = ""
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)"', line)
            gene_id = match.group(1)

            if not current_gene_id:
                record_counter += 1
                current_gene_id = gene_id

            if not gene_id == current_gene_id:
                record_counter += 1
                if record_counter % batch_size == 0:
                    records.append(current_record_batch)
                    current_record_batch = []
                current_gene_id = gene_id

            current_record_batch.append(line)
            line = gtf_in.readline()

        records.append(current_record_batch)

    elif record_type == "single_line_feature":
        record_counter = 0
        current_record_batch = []
        current_gene_id = ""
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            record_counter += 1

            if record_counter % batch_size == 0:
                records.append(current_record_batch)
                current_record_batch = []

            current_record_batch.append(line)
            line = gtf_in.readline()

        records.append(current_record_batch)

    gtf_in.close()

    return records


def run_find_orfs(genome_file, main_output_dir):

    min_orf_length = 600

    orf_output_dir = create_dir(main_output_dir, "orf_output")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    for region_name in seq_region_lengths:
        region_length = seq_region_lengths[region_name]
        seq = get_sequence(region_name, 1, region_length, 1, genome_file, orf_output_dir)
        for phase in range(0, 6):
            find_orf_phased_region(
                region_name, seq, phase, min_orf_length, orf_output_dir
            )


def find_orf_phased_region(region_name, seq, phase, min_orf_length, orf_output_dir):

    current_index = phase
    orf_counter = 1
    if phase > 2:
        seq = reverse_complement(seq)
        current_index = current_index % 3

    orf_file_path = os.path.join(
        orf_output_dir, (region_name + ".phase" + str(phase) + ".orf.fa")
    )
    orf_out = open(orf_file_path, "w+")

    while current_index < len(seq):
        codon = seq[current_index : current_index + 3]
        if codon == "ATG":
            orf_seq = codon
            for j in range(current_index + 3, len(seq), 3):
                next_codon = seq[j : j + 3]
                if next_codon == "TAA" or next_codon == "TAG" or next_codon == "TGA":
                    orf_seq += next_codon
                    if len(orf_seq) >= min_orf_length:
                        orf_out.write(
                            ">"
                            + region_name
                            + "_phase"
                            + str(phase)
                            + "_orf"
                            + str(orf_counter)
                            + "\n"
                        )
                        orf_out.write(orf_seq + "\n")
                        orf_counter += 1
                        orf_seq = ""
                        break

                # If there's another met in phase, then put i to the start of the codon after j so that only the longest ORF is found
                if next_codon == "ATG":
                    current_index = j + 3
                orf_seq += next_codon
        current_index += 3
    orf_out.close()


def run_repeatmasker_regions(
    genome_file, repeatmasker_path, library, species, main_output_dir, num_threads
):

    if not repeatmasker_path:
        repeatmasker_path = "RepeatMasker"

    if not library:
        library = "homo"

    check_exe(repeatmasker_path)
    repeatmasker_output_dir = create_dir(main_output_dir, "repeatmasker_output")
    os.chdir(repeatmasker_output_dir)

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(repeatmasker_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logging.info("Repeatmasker gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    if not library:
        if not species:
            species = "homo"
            generic_repeatmasker_cmd = [
                repeatmasker_path,
                "-nolow",
                "-species",
                species,
                "-engine",
                "crossmatch",
                "-dir",
                repeatmasker_output_dir,
            ]
        else:
            generic_repeatmasker_cmd = [
                repeatmasker_path,
                "-nolow",
                "-species",
                species,
                "-engine",
                "crossmatch",
                "-dir",
                repeatmasker_output_dir,
            ]

    else:
        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-lib",
            library,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]

    logging.info("Running RepeatMasker")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_repeatmasker,
            args=(
                generic_repeatmasker_cmd,
                slice_id,
                genome_file,
                repeatmasker_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(repeatmasker_output_dir, ".rm.gtf", 1, "repeat_id", "repeatmask")


def multiprocess_repeatmasker(
    generic_repeatmasker_cmd, slice_id, genome_file, repeatmasker_output_dir
):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing slice to find repeats with RepeatMasker: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, repeatmasker_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(repeatmasker_output_dir, region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_path = region_fasta_file_path + ".rm.gtf"
    repeatmasker_output_file_path = region_fasta_file_path + ".out"
    repeatmasker_masked_file_path = region_fasta_file_path + ".masked"
    repeatmasker_tbl_file_path = region_fasta_file_path + ".tbl"
    repeatmasker_log_file_path = region_fasta_file_path + ".log"
    repeatmasker_cat_file_path = region_fasta_file_path + ".cat"

    repeatmasker_cmd = generic_repeatmasker_cmd.copy()
    repeatmasker_cmd.append(region_fasta_file_path)
    logging.info(" ".join(repeatmasker_cmd))
    subprocess.run(repeatmasker_cmd)

    create_repeatmasker_gtf(
        repeatmasker_output_file_path, region_results_file_path, region_name
    )

    os.remove(region_fasta_file_path)
    if os.path.exists(region_results_file_path):
        os.remove(region_results_file_path)
    if os.path.exists(repeatmasker_masked_file_path):
        os.remove(repeatmasker_masked_file_path)
    if os.path.exists(repeatmasker_tbl_file_path):
        os.remove(repeatmasker_tbl_file_path)
    if os.path.exists(repeatmasker_log_file_path):
        os.remove(repeatmasker_log_file_path)
    if os.path.exists(repeatmasker_cat_file_path):
        os.remove(repeatmasker_cat_file_path)


def create_repeatmasker_gtf(
    repeatmasker_output_file_path, region_results_file_path, region_name
):

    repeatmasker_in = open(repeatmasker_output_file_path, "r")
    repeatmasker_out = open(region_results_file_path, "w+")
    line = repeatmasker_in.readline()
    repeat_count = 1
    while line:
        result_match = re.search(r"^\s*\d+\s+", line)
        if result_match:
            results = line.split()
            if results[-1] == "*":
                results.pop()
            if not len(results) == 15:
                continue

            score = results[0]
            start = results[5]
            end = results[6]
            strand = results[8]
            repeat_name = results[9]
            repeat_class = results[10]
            if strand == "+":
                repeat_start = results[11]
                repeat_end = results[12]
            else:
                repeat_start = results[13]
                repeat_end = results[12]
                strand = "-"

            gtf_line = (
                region_name
                + "\tRepeatMasker\trepeat\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + strand
                + "\t.\t"
                + 'repeat_id "'
                + str(repeat_count)
                + '"; repeat_name "'
                + repeat_name
                + '"; repeat_class "'
                + repeat_class
                + '"; repeat_start "'
                + str(repeat_start)
                + '"; repeat_end "'
                + str(repeat_end)
                + '"; score "'
                + str(score)
                + '";\n'
            )
            repeatmasker_out.write(gtf_line)
            repeat_count += 1
        line = repeatmasker_in.readline()
    repeatmasker_in.close()
    repeatmasker_out.close()


def run_eponine_regions(
    genome_file, java_path, eponine_path, main_output_dir, num_threads
):

    if not java_path:
        java_path = "java"

    if not eponine_path:
        eponine_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/eponine/libexec/eponine-scan.jar"

    check_file(eponine_path)
    check_exe(java_path)

    eponine_output_dir = create_dir(main_output_dir, "eponine_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(eponine_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logging.info("Eponine gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    threshold = "0.999"
    generic_eponine_cmd = [
        java_path,
        "-jar",
        eponine_path,
        "-threshold",
        threshold,
        "-seq",
    ]
    logging.info("Running Eponine")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_eponine,
            args=(
                generic_eponine_cmd,
                slice_id,
                genome_file,
                eponine_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(eponine_output_dir, ".epo.gtf", 1, "feature_id", "eponine")


def multiprocess_eponine(generic_eponine_cmd, slice_id, genome_file, eponine_output_dir):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing slice to find transcription start sites with Eponine: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, eponine_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(eponine_output_dir, region_fasta_file_name)

    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".epo.gtf"
    region_results_file_path = os.path.join(eponine_output_dir, region_results_file_name)

    eponine_output_file_path = region_fasta_file_path + ".epo"
    eponine_out = open(eponine_output_file_path, "w+")

    eponine_cmd = generic_eponine_cmd.copy()
    eponine_cmd.append(region_fasta_file_path)

    logging.info(" ".join(eponine_cmd))
    subprocess.run(eponine_cmd, stdout=eponine_out)
    eponine_out.close()

    create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name)
    os.remove(eponine_output_file_path)
    os.remove(region_fasta_file_path)


def create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name):

    eponine_in = open(eponine_output_file_path, "r")
    eponine_out = open(region_results_file_path, "w+")

    line = eponine_in.readline()
    feature_count = 1
    while line:
        result_match = re.search(r"^" + region_name, line)

        if result_match:
            results = line.split()
            start = int(results[3])
            end = int(results[4])
            score = float(results[5])
            strand = results[6]

            # There's a one base offset on the reverse strand
            if strand == "-":
                start -= 1
                end -= 1

            gtf_line = (
                region_name
                + "\tEponine\tsimple_feature\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + strand
                + "\t.\t"
                + 'feature_id "'
                + str(feature_count)
                + '"; score "'
                + str(score)
                + '";\n'
            )
            eponine_out.write(gtf_line)
            feature_count += 1
        line = eponine_in.readline()
    eponine_in.close()
    eponine_out.close()


def run_cpg_regions(genome_file, cpg_path, main_output_dir, num_threads):

    if not cpg_path:
        cpg_path = "cpg_lh"

    check_exe(cpg_path)
    cpg_output_dir = create_dir(main_output_dir, "cpg_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(cpg_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logging.info("Cpg gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    logging.info("Running CpG")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_cpg,
            args=(
                cpg_path,
                slice_id,
                genome_file,
                cpg_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(cpg_output_dir, ".cpg.gtf", 1, "feature_id", "cpg")


def multiprocess_cpg(cpg_path, slice_id, genome_file, cpg_output_dir):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing slice to find CpG islands with cpg_lh: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, cpg_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(cpg_output_dir, region_fasta_file_name)

    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".cpg.gtf"
    region_results_file_path = os.path.join(cpg_output_dir, region_results_file_name)

    cpg_output_file_path = region_fasta_file_path + ".cpg"
    cpg_out = open(cpg_output_file_path, "w+")

    cpg_cmd = [cpg_path, region_fasta_file_path]
    logging.info(" ".join(cpg_cmd))
    subprocess.run(cpg_cmd, stdout=cpg_out)
    cpg_out.close()

    create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name)
    os.remove(cpg_output_file_path)
    os.remove(region_fasta_file_path)


def create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name):

    cpg_min_length = 400
    cpg_min_gc_content = 50
    cpg_min_oe = 0.6

    cpg_in = open(cpg_output_file_path, "r")
    cpg_out = open(region_results_file_path, "w+")
    line = cpg_in.readline()
    feature_count = 1
    while line:
        result_match = re.search(r"^" + region_name, line)
        if result_match:
            results = line.split()
            start = int(results[1])
            end = int(results[2])
            length = end - start + 1
            score = float(results[3])
            gc_content = float(results[6])
            oe = results[7]

            if oe == "-" or oe == "inf":
                oe = 0
            else:
                oe = float(oe)

            if (
                length >= cpg_min_length
                and gc_content >= cpg_min_gc_content
                and oe >= cpg_min_oe
            ):
                gtf_line = (
                    region_name
                    + "\tCpG\tsimple_feature\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t.\t+\t.\t"
                    + 'feature_id "'
                    + str(feature_count)
                    + '"; score "'
                    + str(score)
                    + '";\n'
                )
                cpg_out.write(gtf_line)
                feature_count += 1
        line = cpg_in.readline()
    cpg_in.close()
    cpg_out.close()


def run_trnascan_regions(
    genome_file, trnascan_path, trnascan_filter_path, main_output_dir, num_threads
):

    if not trnascan_path:
        trnascan_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/tRNAscan-SE"
    logging.info(trnascan_path)
    if not trnascan_filter_path:
        trnascan_filter_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/EukHighConfidenceFilter"
    logging.info(trnascan_filter_path)
    check_exe(trnascan_path)
    logging.info(trnascan_path)
    # check_exe(trnascan_filter_path)
    check_file(trnascan_filter_path)
    logging.info(trnascan_filter_path)

    trnascan_output_dir = create_dir(main_output_dir, "trnascan_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(trnascan_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Trnascan gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    generic_trnascan_cmd = [
        trnascan_path,
        None,
        "-o",
        None,
        "-f",
        None,
        "-H",
        "-q",
        "--detail",
        "-Q",
    ]
    logging.info("Running tRNAscan-SE")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trnascan,
            args=(
                generic_trnascan_cmd,
                slice_id,
                genome_file,
                trnascan_filter_path,
                trnascan_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(trnascan_output_dir, ".trna.gtf", 1, None, None)


def multiprocess_trnascan(
    generic_trnascan_cmd, slice_id, genome_file, trnascan_filter_path, trnascan_output_dir
):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing slice to find tRNAs using tRNAscan-SE: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, trnascan_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(trnascan_output_dir, region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".trna.gtf"
    region_results_file_path = os.path.join(trnascan_output_dir, region_results_file_name)

    trnascan_output_file_path = region_fasta_file_path + ".trna"
    trnascan_ss_output_file_path = trnascan_output_file_path + ".ss"

    # The filter takes an output dir and a prefix and then uses those to make a path to a .out file
    trnascan_filter_file_prefix = region_fasta_file_name + ".filt"
    trnascan_filter_file_name = trnascan_filter_file_prefix + ".out"
    trnascan_filter_log_file_name = trnascan_filter_file_prefix + ".log"
    trnascan_filter_ss_file_name = trnascan_filter_file_prefix + ".ss"
    trnascan_filter_file_path = os.path.join(
        trnascan_output_dir, trnascan_filter_file_name
    )
    trnascan_filter_log_file_path = os.path.join(
        trnascan_output_dir, trnascan_filter_log_file_name
    )
    trnascan_filter_ss_file_path = os.path.join(
        trnascan_output_dir, trnascan_filter_ss_file_name
    )

    trnascan_cmd = generic_trnascan_cmd.copy()
    trnascan_cmd[1] = region_fasta_file_path
    trnascan_cmd[3] = trnascan_output_file_path
    trnascan_cmd[5] = trnascan_ss_output_file_path

    logging.info("tRNAscan-SE command:")
    logging.info(" ".join(trnascan_cmd))
    subprocess.run(trnascan_cmd)

    # If we have a blank output file at this point we want to stop and remove whatever files
    # have been created instead of moving onto the filter
    if os.stat(trnascan_output_file_path).st_size == 0:
        os.remove(trnascan_output_file_path)
        os.remove(region_fasta_file_path)
        if os.path.exists(trnascan_ss_output_file_path):
            os.remove(trnascan_ss_output_file_path)
        return

    filter_cmd = [
        trnascan_filter_path,
        "--result",
        trnascan_output_file_path,
        "--ss",
        trnascan_ss_output_file_path,
        "--output",
        trnascan_output_dir,
        "--prefix",
        trnascan_filter_file_prefix,
    ]
    logging.info("tRNAscan-SE filter command:")
    logging.info(" ".join(filter_cmd))
    subprocess.run(filter_cmd)

    create_trnascan_gtf(region_results_file_path, trnascan_filter_file_path, region_name)
    if os.path.exists(trnascan_output_file_path):
        os.remove(trnascan_output_file_path)
    if os.path.exists(trnascan_ss_output_file_path):
        os.remove(trnascan_ss_output_file_path)
    if os.path.exists(trnascan_filter_file_path):
        os.remove(trnascan_filter_file_path)
    if os.path.exists(trnascan_filter_log_file_path):
        os.remove(trnascan_filter_log_file_path)
    if os.path.exists(trnascan_filter_ss_file_path):
        os.remove(trnascan_filter_ss_file_path)
    if os.path.exists(region_fasta_file_path):
        os.remove(region_fasta_file_path)


def create_trnascan_gtf(region_results_file_path, trnascan_filter_file_path, region_name):

    trna_in = open(trnascan_filter_file_path, "r")
    trna_out = open(region_results_file_path, "w+")
    line = trna_in.readline()
    gene_counter = 1
    while line:
        result_match = re.search(r"^" + region_name, line)
        if result_match:
            results = line.split()
            start = int(results[2])
            end = int(results[3])
            trna_type = results[4]
            score = results[8]

            strand = "+"
            if start > end:
                strand = "-"
                temp_end = start
                start = end
                end = temp_end

            biotype = "tRNA_pseudogene"
            high_confidence_match = re.search(r"high confidence set", line)
            if high_confidence_match:
                biotype = "tRNA"

            transcript_string = (
                region_name
                + "\ttRNAscan\ttranscript\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + strand
                + "\t.\t"
                + 'gene_id "'
                + str(gene_counter)
                + '"; transcript_id "'
                + str(gene_counter)
                + '"; biotype "'
                + biotype
                + '";\n'
            )
            exon_string = (
                region_name
                + "\ttRNAscan\texon\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + strand
                + "\t.\t"
                + 'gene_id "'
                + str(gene_counter)
                + '"; transcript_id "'
                + str(gene_counter)
                + '"; exon_number "1"; biotype "'
                + biotype
                + '";\n'
            )

            trna_out.write(transcript_string)
            trna_out.write(exon_string)
            gene_counter += 1
        line = trna_in.readline()
    trna_in.close()
    trna_out.close()


def run_dust_regions(genome_file, dust_path, main_output_dir, num_threads):

    if not dust_path:
        dust_path = "dustmasker"

    check_exe(dust_path)
    dust_output_dir = create_dir(main_output_dir, "dust_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(dust_output_dir, "annotation.gtf")
    logging.info(output_file)
    if Path(output_file).is_file():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logging.info("Dust gtf file exists")
            return 0
    else:
        logging.info("No gtf file, go on with the analysis")

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    generic_dust_cmd = [dust_path, "-in"]
    logging.info("Running Dust")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_dust,
            args=(
                generic_dust_cmd,
                slice_id,
                genome_file,
                dust_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(dust_output_dir, ".dust.gtf", 1, "repeat_id", "dust")


def multiprocess_dust(generic_dust_cmd, slice_id, genome_file, dust_output_dir):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing slice to find low complexity regions with Dust: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, dust_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(dust_output_dir, region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".dust.gtf"
    region_results_file_path = os.path.join(dust_output_dir, region_results_file_name)

    dust_output_file_path = region_fasta_file_path + ".dust"
    dust_out = open(dust_output_file_path, "w+")
    dust_cmd = generic_dust_cmd.copy()
    dust_cmd.append(region_fasta_file_path)
    logging.info(" ".join(dust_cmd))
    subprocess.run(dust_cmd, stdout=dust_out)
    dust_out.close()

    create_dust_gtf(dust_output_file_path, region_results_file_path, region_name)
    os.remove(dust_output_file_path)
    os.remove(region_fasta_file_path)


def create_dust_gtf(dust_output_file_path, region_results_file_path, region_name):

    dust_in = open(dust_output_file_path, "r")
    dust_out = open(region_results_file_path, "w+")
    line = dust_in.readline()
    repeat_count = 1
    while line:
        result_match = re.search(r"(\d+)\ - (\d+)", line)
        if result_match:
            start = int(result_match.group(1)) + 1
            end = int(result_match.group(2)) + 1
            gtf_line = (
                region_name
                + "\tDust\trepeat\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t+\t.\t"
                + 'repeat_id "'
                + str(repeat_count)
                + '";\n'
            )
            dust_out.write(gtf_line)
            repeat_count += 1
        line = dust_in.readline()
    dust_in.close()
    dust_out.close()


def run_trf_repeats(genome_file, trf_path, main_output_dir, num_threads):

    if not trf_path:
        trf_path = "trf"

    check_exe(trf_path)
    trf_output_dir = create_dir(main_output_dir, "trf_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(trf_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logging.info("Trf gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    match_score = 2
    mismatch_score = 5
    delta = 7
    pm = 80
    pi = 10
    minscore = 40
    maxperiod = 500

    trf_output_extension = (
        "."
        + str(match_score)
        + "."
        + str(mismatch_score)
        + "."
        + str(delta)
        + "."
        + str(pm)
        + "."
        + str(pi)
        + "."
        + str(minscore)
        + "."
        + str(maxperiod)
        + ".dat"
    )

    generic_trf_cmd = [
        trf_path,
        None,
        str(match_score),
        str(mismatch_score),
        str(delta),
        str(pm),
        str(pi),
        str(minscore),
        str(maxperiod),
        "-d",
        "-h",
    ]
    logging.info("Running TRF")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trf,
            args=(
                generic_trf_cmd,
                slice_id,
                genome_file,
                trf_output_dir,
                trf_output_extension,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(trf_output_dir, ".trf.gtf", 1, "repeat_id", "trf")


def multiprocess_trf(
    generic_trf_cmd, slice_id, genome_file, trf_output_dir, trf_output_extension
):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing slice to find tandem repeats with TRF: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, trf_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(trf_output_dir, region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".trf.gtf"
    region_results_file_path = os.path.join(trf_output_dir, region_results_file_name)

    # TRF writes to the current dir, so swtich to the output dir for it
    os.chdir(trf_output_dir)
    trf_output_file_path = region_fasta_file_path + trf_output_extension
    trf_cmd = generic_trf_cmd.copy()
    trf_cmd[1] = region_fasta_file_path
    logging.info(" ".join(trf_cmd))
    subprocess.run(trf_cmd)
    create_trf_gtf(trf_output_file_path, region_results_file_path, region_name)
    os.remove(trf_output_file_path)
    os.remove(region_fasta_file_path)


def create_trf_gtf(trf_output_file_path, region_results_file_path, region_name):

    trf_in = open(trf_output_file_path, "r")
    trf_out = open(region_results_file_path, "w+")
    line = trf_in.readline()
    repeat_count = 1
    while line:
        result_match = re.search(r"^\d+", line)
        if result_match:
            results = line.split()
            if not len(results) == 15:
                continue
            start = results[0]
            end = results[1]
            period = float(results[2])
            copy_number = float(results[3])
            percent_matches = float(results[5])
            score = float(results[7])
            repeat_consensus = results[13]
            if (
                score < 50 and percent_matches >= 80 and copy_number > 2 and period < 10
            ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                gtf_line = (
                    region_name
                    + "\tTRF\trepeat\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t.\t+\t.\t"
                    + 'repeat_id "'
                    + str(repeat_count)
                    + '"; score "'
                    + str(score)
                    + '"; repeat_consensus "'
                    + repeat_consensus
                    + '";\n'
                )
                trf_out.write(gtf_line)
                repeat_count += 1
        line = trf_in.readline()
    trf_in.close()
    trf_out.close()


def run_cmsearch_regions(
    genome_file,
    cmsearch_path,
    rfam_cm_db_path,
    rfam_seeds_file_path,
    rfam_accession_file,
    main_output_dir,
    num_threads,
):

    if not cmsearch_path:
        cmsearch_path = "cmsearch"

    check_exe(cmsearch_path)
    rfam_output_dir = create_dir(main_output_dir, "rfam_output")

    os.chdir(rfam_output_dir)

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(rfam_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Cmsearch gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    rfam_dbname = "Rfam"
    rfam_user = "rfamro"
    rfam_host = "mysql-rfam-public.ebi.ac.uk"
    rfam_port = "4497"
    #  rfam_accession_query_cmd = ["mysql -h", rfam_host,"-u",rfam_user,"-P",rfam_port,"-NB -e",rfam_dbname,"'select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE \'%" + clade + "%\' AND (f.type LIKE \'%snRNA%\' OR f.type LIKE \'%rRNA%\' OR LOWER(f.rfam_id) LIKE \'%rnase%\' OR LOWER(f.rfam_id) LIKE \'%vault%\' OR LOWER(f.rfam_id) LIKE \'%y_rna%\' OR f.rfam_id LIKE \'%Metazoa_SRP%\') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE \'%bacteria%\' AND rfam_id NOT LIKE \'%archaea%\' AND rfam_id NOT LIKE \'%microsporidia%\';'"]

    # mysql -hmysql-rfam-public.ebi.ac.uk -urfamro -P4497 Rfam -NB -e "select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE '%insect%' AND (f.type LIKE '%snRNA%' OR f.type LIKE '%rRNA%' OR LOWER(f.rfam_id) LIKE '%rnase%' OR LOWER(f.rfam_id) LIKE '%vault%' OR LOWER(f.rfam_id) LIKE '%y_rna%' OR f.rfam_id LIKE '%Metazoa_SRP%') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE '%bacteria%' AND rfam_id NOT LIKE '%archaea%' AND rfam_id NOT LIKE '%microsporidia%';"

    #  rfam_accession_file = '/hps/nobackup2/production/ensembl/fergal/production/test_runs/non_verts/butterfly/rfam_insect_ids.txt'
    # rfam_accession_file = os.path.join(main_output_dir,'rfam_accessions.txt')
    rfam_cm_db_path = (
        "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.cm"
    )
    rfam_seeds_file_path = (
        "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed"
    )
    rfam_selected_models_file = os.path.join(rfam_output_dir, "rfam_models.cm")
    with open(rfam_accession_file) as rfam_accessions_in:
        rfam_accessions = rfam_accessions_in.read().splitlines()

    with open(rfam_cm_db_path, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    rfam_cm_out = open(rfam_selected_models_file, "w+")

    for model in rfam_models:
        # The Rfam.cm file has INFERNAL and HMMR models, both are needed at this point
        # Later we just want the INFERNAL ones for looking at thresholds
        match = re.search(r"(RF\d+)", model)
        if match:
            model_accession = match.group(1)
            if model_accession in rfam_accessions:
                rfam_cm_out.write(model + "//\n")
    rfam_cm_out.close()

    seed_descriptions = get_rfam_seed_descriptions(rfam_seeds_file_path)
    cv_models = extract_rfam_metrics(rfam_selected_models_file)

    logging.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 100000, 0, 5000)

    generic_cmsearch_cmd = [
        cmsearch_path,
        "--rfam",
        "--cpu",
        "1",
        "--nohmmonly",
        "--cut_ga",
        "--tblout",
    ]
    logging.info("Running Rfam")
    pool = multiprocessing.Pool(int(num_threads))
    results = []
    failed_slice_ids = []
    memory_limit = 3 * 1024**3
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_cmsearch,
            args=(
                generic_cmsearch_cmd,
                slice_id,
                genome_file,
                rfam_output_dir,
                rfam_selected_models_file,
                cv_models,
                seed_descriptions,
                memory_limit,
            ),
        )
    pool.close()
    pool.join()

    # Need to figure something more automated out here. At the moment it's just limiting to 5 cores and 5GB vram
    # Ideally we could look at the amount of mem requested and put something like 10GB per core and then figure
    # out how many cores to use (obviously not using more than the amount specified)
    memory_limit = 5 * 1024**3
    if num_threads > 5:
        num_threads = 5
    pool = multiprocessing.Pool(num_threads)
    exception_files = glob.glob(rfam_output_dir + "/*.rfam.except")
    for exception_file_path in exception_files:
        exception_file_name = os.path.basename(exception_file_path)
        logging.info("Running himem job for failed region:\n" + exception_file_path)
        match = re.search(r"(.+)\.rs(\d+)\.re(\d+)\.", exception_file_name)
        if match:
            except_region = match.group(1)
            except_start = match.group(2)
            except_end = match.group(3)
            except_slice_id = [except_region, except_start, except_end]
            pool.apply_async(
                multiprocess_cmsearch,
                args=(
                    generic_cmsearch_cmd,
                    except_slice_id,
                    genome_file,
                    rfam_output_dir,
                    rfam_selected_models_file,
                    cv_models,
                    seed_descriptions,
                    memory_limit,
                ),
            )
    pool.close()
    pool.join()

    slice_output_to_gtf(rfam_output_dir, ".rfam.gtf", 1, None, None)


def prlimit_command(command_list, virtual_memory_limit):

    """
    Prepend memory limiting arguments to a command list to be run with subprocess.

    This method uses the `prlimit` program to set the memory limit.

    The `virtual_memory_limit` size is in bytes.

    prlimit arguments:
    -v, --as[=limits]
           Address space limit.
    """
    return ["prlimit", f"-v{virtual_memory_limit}"] + command_list


def multiprocess_cmsearch(
    generic_cmsearch_cmd,
    slice_id,
    genome_file,
    rfam_output_dir,
    rfam_selected_models_file,
    cv_models,
    seed_descriptions,
    memory_limit,
):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logging.info(
        "Processing Rfam data using cmsearch against slice: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, rfam_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(rfam_output_dir, region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_tblout_file_name = slice_file_name + ".tblout"
    region_tblout_file_path = os.path.join(rfam_output_dir, region_tblout_file_name)
    region_results_file_name = slice_file_name + ".rfam.gtf"
    region_results_file_path = os.path.join(rfam_output_dir, region_results_file_name)

    exception_results_file_name = slice_file_name + ".rfam.except"
    exception_results_file_path = os.path.join(
        rfam_output_dir, exception_results_file_name
    )

    cmsearch_cmd = generic_cmsearch_cmd.copy()
    cmsearch_cmd.append(region_tblout_file_path)
    cmsearch_cmd.append(rfam_selected_models_file)
    cmsearch_cmd.append(region_fasta_file_path)
    logging.info(" ".join(cmsearch_cmd))

    if memory_limit is not None:
        cmsearch_cmd = prlimit_command(cmsearch_cmd, memory_limit)

    logging.info(" ".join(cmsearch_cmd))

    return_value = None
    try:
        return_value = subprocess.check_output(cmsearch_cmd)
    except subprocess.CalledProcessError as ex:
        # Note that writing to file was the only option here. If return_value was passed back, eventually it would clog the
        # tiny pipe that is used by the workers to send info back. That would mean that it would eventually just be one
        # worked running at a time
        logging.error(
            "Issue processing the following region with cmsearch: "
            + region_name
            + " "
            + str(start)
            + "-"
            + str(end)
        )
        logging.error("Return value: " + str(return_value))
        exception_out = open(exception_results_file_path, "w+")
        exception_out.write(region_name + " " + str(start) + " " + str(end) + "\n")
        exception_out.close()
        os.remove(region_fasta_file_path)
        os.remove(region_tblout_file_path)
        return

    initial_table_results = parse_rfam_tblout(region_tblout_file_path, region_name)
    unique_table_results = remove_rfam_overlap(initial_table_results)
    filtered_table_results = filter_rfam_results(unique_table_results, cv_models)
    create_rfam_gtf(
        filtered_table_results,
        cv_models,
        seed_descriptions,
        region_name,
        region_results_file_path,
        genome_file,
        rfam_output_dir,
    )
    os.remove(region_fasta_file_path)
    os.remove(region_tblout_file_path)
    gc.collect()


def get_rfam_seed_descriptions(rfam_seeds_path):

    descriptions = {}
    rfam_seeds = []

    # NOTE: for some reason the decoder breaks on the seeds file, so I have made this ignore errors
    with open(rfam_seeds_path, encoding="utf-8", errors="ignore") as rfam_seeds_in:
        rfam_seeds = rfam_seeds_in.read().splitlines()

    for seed in rfam_seeds:
        domain_match = re.search("^\#=GF AC   (.+)", seed)
        if domain_match:
            domain = domain_match.group(1)
            descriptions[domain] = {}
            continue

        description_match = re.search("^\#=GF DE   (.+)", seed)
        if description_match:
            description = description_match.group(1)
            descriptions[domain]["description"] = description
            continue

        name_match = re.search("^\#=GF ID   (.+)", seed)
        if name_match:
            name = name_match.group(1)
            descriptions[domain]["name"] = name
            continue

        type_match = re.search("^\#=GF TP   Gene; (.+)", seed)
        if type_match:
            rfam_type = type_match.group(1)
            descriptions[domain]["type"] = rfam_type
            continue

    return descriptions


def extract_rfam_metrics(rfam_selected_models_file):

    with open(rfam_selected_models_file, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    parsed_cm_data = {}
    for model in rfam_models:
        temp = model.split("\n")
        model_name_match = re.search(r"NAME\s+(\S+)", model)
        match_infernal = re.search(r"INFERNAL", model)
        if model_name_match and match_infernal:
            model_name = model_name_match.group(1)
            parsed_cm_data[model_name] = {}
            for line in temp:
                name_match = re.search(r"^NAME\s+(\S+)", line)
                if name_match:
                    parsed_cm_data[model_name]["-name"] = name_match.group(1)
                    continue

                description_match = re.search(r"^DESC\s+(\S+)", line)
                if description_match:
                    parsed_cm_data[model_name]["-description"] = description_match.group(
                        1
                    )
                    continue

                length_match = re.search(r"^CLEN\s+(\d+)", line)
                if length_match:
                    parsed_cm_data[model_name]["-length"] = length_match.group(1)
                    continue

                max_length_match = re.search(r"^W\s+(\d+)", line)
                if max_length_match:
                    parsed_cm_data[model_name]["-maxlength"] = max_length_match.group(1)
                    continue

                accession_match = re.search(r"^ACC\s+(\S+)", line)
                if accession_match:
                    parsed_cm_data[model_name]["-accession"] = accession_match.group(1)
                    continue

                threshold_match = re.search(r"^GA\s+(\d+)", line)
                if threshold_match:
                    parsed_cm_data[model_name]["-threshold"] = threshold_match.group(1)
                    continue

    return parsed_cm_data


def parse_rfam_tblout(region_tblout_file_path, region_name):

    with open(region_tblout_file_path, "r") as rfam_tbl_in:
        rfam_tbl_data = rfam_tbl_in.read()

    tbl_results = rfam_tbl_data.split("\n")

    all_parsed_results = []
    for result in tbl_results:
        parsed_tbl_data = {}
        if not re.match(r"^" + region_name, result):
            continue

        hit = result.split()
        accession = hit[3]
        target_name = hit[0]
        query_name = hit[2]
        hstart = hit[5]
        hend = hit[6]
        start = hit[7]
        end = hit[8]
        if hit[9] == "+":
            strand = 1
        else:
            strand = -1
        evalue = hit[15]
        score = hit[14]

        parsed_tbl_data["accession"] = accession
        parsed_tbl_data["start"] = start
        parsed_tbl_data["end"] = end
        parsed_tbl_data["strand"] = strand
        parsed_tbl_data["query_name"] = query_name
        parsed_tbl_data["score"] = score
        all_parsed_results.append(parsed_tbl_data)
    return all_parsed_results


# NOTE some of the code above and the code commented out here is to do with creating
# a DAF. As we don't have a python concept of this I'm leaving it out for the moment
# but the code below is a reference
#    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
#      -slice          => $self->queries,
#      -start          => $strand == 1 ? $start : $end,
#      -end            => $strand == 1 ? $end : $start,
#      -strand         => $strand,
#      -hstart         => $hstart,
#      -hend           => $hend,
#      -hstrand        => $strand,
#      -score          => $score,
#      -hseqname       => length($target_name) > 39 ? substr($target_name, 0, 39) : $target_name,,
#      -p_value  => $evalue,
#      -align_type => 'ensembl',
#      -cigar_string  => abs($hend - $hstart) . "M",
#   );


def remove_rfam_overlap(parsed_tbl_data):

    excluded_structures = {}
    chosen_structures = []
    for structure_x in parsed_tbl_data:
        chosen_structure = structure_x
        structure_x_start = int(structure_x["start"])
        structure_x_end = int(structure_x["end"])
        structure_x_score = float(structure_x["score"])
        structure_x_accession = structure_x["accession"]
        structure_x_string = ":".join(
            [
                str(structure_x_start),
                str(structure_x_end),
                str(structure_x_score),
                structure_x_accession,
            ]
        )
        for structure_y in parsed_tbl_data:
            structure_y_start = int(structure_y["start"])
            structure_y_end = int(structure_y["end"])
            structure_y_score = float(structure_y["score"])
            structure_y_accession = structure_y["accession"]
            structure_y_string = ":".join(
                [
                    str(structure_y_start),
                    str(structure_y_end),
                    str(structure_y_score),
                    structure_y_accession,
                ]
            )
            if structure_y_string in excluded_structures:
                continue

            if (
                structure_x_start <= structure_y_end
                and structure_x_end >= structure_y_start
            ):
                if structure_x_score < structure_y_score:
                    chosen_structure = structure_y
                    excluded_structures[structure_x_string] = 1
                else:
                    excluded_structures[structure_y_string] = 1

        chosen_structures.append(chosen_structure)

    return chosen_structures


def filter_rfam_results(unfiltered_tbl_data, cv_models):

    filtered_results = []
    for structure in unfiltered_tbl_data:
        query = structure["query_name"]
        if query in cv_models:
            threshold = cv_models[query]["-threshold"]
            if query == "LSU_rRNA_eukarya":
                threshold = 1700

            elif query == "LSU_rRNA_archaea":
                continue

            elif query == "LSU_rRNA_bacteria":
                continue

            elif query == "SSU_rRNA_eukarya":
                threshold = 1600

            elif query == "5_8S_rRNA":
                threshold = 85

            elif query == "5S_rRNA":
                threshold = 75

            if threshold and float(structure["score"]) >= float(threshold):
                filtered_results.append(structure)

    return filtered_results


# NOTE: The below are notes from the perl code (which has extra code) about possible improvements that are not implemented there
# Although not included in RefSeq filters, additional filters that consider sizes and score_to_size ratios can be applied
# in future work to further exclude FPs
#
# my $is_valid_size = $mapping_length > $min_length && $mapping_length < $max_length ? 1 : 0;
# my $score_size_ratio = $result->{'score'} / $mapping_length;


def create_rfam_gtf(
    filtered_results,
    cm_models,
    descriptions,
    region_name,
    region_results_file_path,
    genome_file,
    rfam_output_dir,
):

    if not filtered_results:
        return

    rfam_gtf_out = open(region_results_file_path, "w+")
    gene_counter = 1
    for structure in filtered_results:
        query = structure["query_name"]
        accession = structure["accession"]
        if query in cm_models:
            model = cm_models[query]
            if accession in descriptions:
                description = descriptions[accession]
                if "type" in description:
                    rfam_type = description["type"]
                else:
                    description = None
                    rfam_type = "misc_RNA"
            domain = structure["query_name"]
            padding = model["-length"]
            #      rfam_type = description['type']
            gtf_strand = structure["strand"]
            rnafold_strand = structure["strand"]
            if gtf_strand == 1:
                start = structure["start"]
                end = structure["end"]
                gtf_strand = "+"
            else:
                start = structure["end"]
                end = structure["start"]
                score = structure["score"]
                gtf_strand = "-"
                rnafold_strand = -1

            biotype = "misc_RNA"
            if re.match(r"^snRNA;", rfam_type):
                biotype = "snRNA"
            if re.match(r"^snRNA; snoRNA", rfam_type):
                biotype = "snoRNA"
            if re.match(r"^snRNA; snoRNA; scaRNA;", rfam_type):
                biotype = "scaRNA"
            if re.match(r"rRNA;", rfam_type):
                biotype = "rRNA"
            if re.match(r"antisense;", rfam_type):
                biotype = "antisense"
            if re.match(r"antitoxin;", rfam_type):
                biotype = "antitoxin"
            if re.match(r"ribozyme;", rfam_type):
                biotype = "ribozyme"
            if re.match(r"" + domain, rfam_type):
                biotype = domain
            if re.match(r"" + domain, rfam_type):
                biotype = "Vault_RNA"
            if re.match(r"" + domain, rfam_type):
                biotype = "Y_RNA"

            rna_seq = get_sequence(
                region_name, start, end, rnafold_strand, genome_file, rfam_output_dir
            )
            valid_structure = check_rnafold_structure(rna_seq, rfam_output_dir)

            if not valid_structure:
                continue

            transcript_string = (
                region_name
                + "\tRfam\ttranscript\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + gtf_strand
                + "\t.\t"
                + 'gene_id "'
                + str(gene_counter)
                + '"; transcript_id "'
                + str(gene_counter)
                + '"; biotype "'
                + biotype
                + '";\n'
            )
            exon_string = (
                region_name
                + "\tRfam\texon\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + gtf_strand
                + "\t.\t"
                + 'gene_id "'
                + str(gene_counter)
                + '"; transcript_id "'
                + str(gene_counter)
                + '"; exon_number "1"; biotype "'
                + biotype
                + '";\n'
            )

            rfam_gtf_out.write(transcript_string)
            rfam_gtf_out.write(exon_string)
            gene_counter += 1
    rfam_gtf_out.close()


def check_rnafold_structure(seq, rfam_output_dir):

    # Note there's some extra code in the RNAfold Perl module for encoding the structure into an attrib
    # Could consider implementing this when running for loading into an Ensembl db
    structure = 0
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=rfam_output_dir
    ) as rna_temp_in:
        rna_temp_in.writelines(">seq1\n" + seq + "\n")
        rna_in_file_path = rna_temp_in.name

    rnafold_cmd = ["RNAfold", "--infile", rna_in_file_path]
    rnafold_output = subprocess.Popen(rnafold_cmd, stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(rnafold_output.stdout, encoding="utf-8"):
        match = re.search(r"([().]+)\s\(\s*(-*\d+.\d+)\)\n$", line)
        if match:
            structure = match.group(1)
            score = match.group(2)
            break
    rna_temp_in.close()
    os.remove(rna_in_file_path)

    return structure


def slice_output_to_gtf(
    output_dir, extension, unique_ids, feature_id_label, new_id_prefix
):

    # Note that this does not make unique ids at the moment
    # In many cases this is fine because the ids are unique by seq region, but in cases like batching it can cause problems
    # So will add in a helper method to make ids unique

    # This holds keys of the current slice details with the gene id to form unique keys. Each time a new key is added
    # the overall gene counter is incremented and the value of the key is set to the new gene id. Any subsequent
    # lines with the same region/gene id key will then just get the new id without incrementing the counter
    gene_id_index = {}
    gene_transcript_id_index = {}
    gene_counter = 1

    # Similar to the gene id index, this will have a key that is based on the slice details, gene id and transcript id. If there
    # is no existing entry, the transcript key will be added and the transcript counter is incremented. If there is a key then
    # the transcript id will be replaced with the new transcript id (which is based on the new gene id and transcript counter)
    # Example key KS8000.rs1.re1000000.gene_1.transcript_1 =
    transcript_id_count_index = {}

    feature_counter = 1

    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    if not extension:
        extension = ".gtf"
    gtf_files = glob.glob(output_dir + "/*" + extension)
    gtf_file_path = os.path.join(output_dir, "annotation.gtf")
    gtf_out = open(gtf_file_path, "w+")
    for gtf_file_path in gtf_files:
        if os.stat(gtf_file_path).st_size == 0:
            logging.info("File is empty, will skip:\n" + gtf_file_path)
            continue

        gtf_file_name = os.path.basename(gtf_file_path)
        match = re.search(r"\.rs(\d+)\.re(\d+)\.", gtf_file_name)
        start_offset = int(match.group(1))
        gtf_in = open(gtf_file_path, "r")
        line = gtf_in.readline()
        while line:
            values = line.split("\t")
            if len(values) == 9 and (values[2] in feature_types):
                values[3] = str(int(values[3]) + (start_offset - 1))
                values[4] = str(int(values[4]) + (start_offset - 1))
                if unique_ids:
                    # Maybe make a unique id based on the feature type
                    # Basically region/feature id should be unique at this point, so could use region_id and current_id is key, value is the unique id that is incremented
                    attribs = values[8]

                    # This bit assigns unique gene/transcript ids if the line contains gene_id/transcript_id
                    match_gene_type = re.search(
                        r'(gene_id +"([^"]+)").+(transcript_id +"([^"]+)")', line
                    )
                    if match_gene_type:
                        full_gene_id_string = match_gene_type.group(1)
                        current_gene_id = match_gene_type.group(2)
                        full_transcript_id_string = match_gene_type.group(3)
                        current_transcript_id = match_gene_type.group(4)
                        gene_id_key = gtf_file_name + "." + str(current_gene_id)
                        transcript_id_key = gene_id_key + "." + str(current_transcript_id)
                        if gene_id_key not in gene_id_index:
                            new_gene_id = "gene" + str(gene_counter)
                            gene_id_index[gene_id_key] = new_gene_id
                            attribs = re.sub(
                                full_gene_id_string,
                                'gene_id "' + new_gene_id + '"',
                                attribs,
                            )
                            transcript_id_count_index[gene_id_key] = 1
                            gene_counter += 1
                        else:
                            new_gene_id = gene_id_index[gene_id_key]
                            attribs = re.sub(
                                full_gene_id_string,
                                'gene_id "' + new_gene_id + '"',
                                attribs,
                            )
                        if transcript_id_key not in gene_transcript_id_index:
                            new_transcript_id = (
                                gene_id_index[gene_id_key]
                                + ".t"
                                + str(transcript_id_count_index[gene_id_key])
                            )
                            gene_transcript_id_index[
                                transcript_id_key
                            ] = new_transcript_id
                            attribs = re.sub(
                                full_transcript_id_string,
                                'transcript_id "' + new_transcript_id + '"',
                                attribs,
                            )
                            transcript_id_count_index[gene_id_key] += 1
                        else:
                            new_transcript_id = gene_transcript_id_index[
                                transcript_id_key
                            ]
                            attribs = re.sub(
                                full_transcript_id_string,
                                'transcript_id "' + new_transcript_id + '"',
                                attribs,
                            )
                        values[8] = attribs

                    # If you don't match a gene line, try a feature line
                    else:
                        match_feature_type = re.search(
                            r"(" + feature_id_label + ' +"([^"]+)")', line
                        )
                        if match_feature_type:
                            full_feature_id_string = match_feature_type.group(1)
                            current_feature_id = match_feature_type.group(2)
                            new_feature_id = new_id_prefix + str(feature_counter)
                            attribs = re.sub(
                                full_feature_id_string,
                                feature_id_label + ' "' + new_feature_id + '"',
                                attribs,
                            )
                            feature_counter += 1
                            values[8] = attribs

                gtf_out.write("\t".join(values))
                line = gtf_in.readline()
            else:
                logging.info(
                    "Feature type not recognised, will skip. Feature type: " + values[2]
                )
                line = gtf_in.readline()
        gtf_in.close()
    gtf_out.close()


def run_red(red_path, main_output_dir, genome_file):

    if not red_path:
        red_path = "Red"

    check_exe(red_path)
    red_dir = create_dir(main_output_dir, "red_output")
    red_mask_dir = create_dir(red_dir, "mask_output")
    red_repeat_dir = create_dir(red_dir, "repeat_output")
    red_genome_dir = create_dir(red_dir, "genome_dir")

    sym_link_genome_cmd = "ln -s " + genome_file + " " + red_genome_dir

    genome_file_name = os.path.basename(genome_file)
    red_genome_file = os.path.join(red_genome_dir, genome_file_name)
    masked_genome_file = os.path.join(
        red_mask_dir, os.path.splitext(genome_file_name)[0] + ".msk"
    )
    repeat_coords_file = os.path.join(
        red_repeat_dir, os.path.splitext(genome_file_name)[0] + ".rpt"
    )
    gtf_output_file_path = os.path.join(red_dir, "annotation.gtf")

    if os.path.exists(masked_genome_file):
        logging.warning(
            "Masked Genome file already found on the path to the Red mask output dir. Will not create a new file"
        )
        create_red_gtf(repeat_coords_file, gtf_output_file_path)
        return masked_genome_file

    if os.path.exists(red_genome_file):
        logging.warning(
            "Unmasked genome file already found on the path to the Red genome dir, will not create a sym link"
        )

    else:
        logging.info(
            "Preparing to sym link the genome file to the Red genome dir. Cmd\n%s"
            % sym_link_genome_cmd
        )
        subprocess.run(["ln", "-s", genome_file, red_genome_dir])

    if not os.path.exists(os.path.join(red_genome_dir, genome_file_name)):
        logging.error(
            "Could not find the genome file in the Red genome dir or sym link to the original file. Path expected:\n%s"
            % red_genome_file
        )

    logging.info("Running Red, this may take some time depending on the genome size")
    subprocess.run(
        [red_path, "-gnm", red_genome_dir, "-msk", red_mask_dir, "-rpt", red_repeat_dir]
    )

    logging.info("Completed running Red")

    create_red_gtf(repeat_coords_file, gtf_output_file_path)

    return masked_genome_file


def create_red_gtf(repeat_coords_file, gtf_output_file_path):

    red_in = open(repeat_coords_file, "r")
    red_out = open(gtf_output_file_path, "w+")
    line = red_in.readline()
    repeat_id = 1
    while line:
        result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
        if result_match:
            region_name = result_match.group(1)
            # Note that Red is 0-based, so add 1
            start = int(result_match.group(2)) + 1
            end = int(result_match.group(3)) + 1
            gtf_line = (
                region_name
                + "\tRed\trepeat\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t+\t.\t"
                + 'repeat_id "'
                + str(repeat_id)
                + '";\n'
            )
            red_out.write(gtf_line)
            repeat_id += 1
        line = red_in.readline()
    red_in.close()
    red_out.close()


def run_genblast_align(
    genblast_path,
    convert2blastmask_path,
    makeblastdb_path,
    genblast_dir,
    protein_file,
    masked_genome_file,
    max_intron_length,
    num_threads,
):

    genblast_timeout_secs = 10800

    if not genblast_path:
        genblast_path = "genblast"

    check_exe(genblast_path)

    if not convert2blastmask_path:
        convert2blastmask_path = "convert2blastmask"

    check_exe(convert2blastmask_path)

    if not makeblastdb_path:
        makeblastdb_path = "makeblastdb"

    check_exe(makeblastdb_path)

    create_dir(genblast_dir, None)

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(genblast_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Genblast gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    genblast_output_file = os.path.join(genblast_dir, "genblast")

    asnb_file = masked_genome_file + ".asnb"
    logging.info("ASNB file: %s" % asnb_file)

    if not os.path.exists("alignscore.txt"):
        subprocess.run(
            [
                "cp",
                os.environ["ENSCODE"]
                + "/ensembl-analysis/scripts/genebuild/gbiab/support_files/alignscore.txt",
                "./",
            ]
        )

    if not os.path.exists(masked_genome_file):
        raise IOError("Masked genome file does not exist: %s" % masked_genome_file)

    if not os.path.exists(protein_file):
        raise IOError("Protein file does not exist: %s" % protein_file)

    if not os.path.exists(asnb_file):
        run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file)
    else:
        logging.info("Found an existing asnb, so will skip convert2blastmask")

    if not os.path.exists(asnb_file):
        raise IOError("asnb file does not exist: %s" % asnb_file)

    run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file)

    batched_protein_files = split_protein_file(protein_file, genblast_dir, 20)

    pool = multiprocessing.Pool(int(num_threads))
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_genblast,
            args=(
                batched_protein_file,
                masked_genome_file,
                genblast_path,
                genblast_timeout_secs,
                max_intron_length,
            ),
        )
    pool.close()
    pool.join()

    logging.info("Completed running GenBlast")
    logging.info("Combining output into single GTF")
    generate_genblast_gtf(genblast_dir)


def multiprocess_genblast(
    batched_protein_file,
    masked_genome_file,
    genblast_path,
    genblast_timeout_secs,
    max_intron_length,
):

    batch_num = os.path.splitext(batched_protein_file)[0]
    batch_dir = os.path.dirname(batched_protein_file)
    logging.info("Running GenBlast on " + batched_protein_file + ":")

    genblast_cmd = [
        genblast_path,
        "-p",
        "genblastg",
        "-q",
        batched_protein_file,
        "-t",
        masked_genome_file,
        "-g",
        "T",
        "-pid",
        "-r",
        "1",
        "-P",
        "blast",
        "-gff",
        "-e",
        "1e-1",
        "-c",
        "0.8",
        "-W",
        "3",
        "-softmask",
        "-scodon",
        "50",
        "-i",
        "30",
        "-x",
        "10",
        "-n",
        "30",
        "-d",
        str(max_intron_length),
        "-o",
        batched_protein_file,
    ]

    logging.info(" ".join(genblast_cmd))
    try:
        subprocess.run(genblast_cmd, timeout=genblast_timeout_secs)
    except subprocess.TimeoutExpired:
        logging.error("Timeout reached for file:\n" + batched_protein_file)
        subprocess.run(["touch", (batched_protein_file + ".except")])

    files_to_delete = glob.glob(batched_protein_file + "*msk.blast*")
    files_to_delete.append(batched_protein_file)
    for file_to_delete in files_to_delete:
        subprocess.run(["rm", file_to_delete])


def generate_genblast_gtf(genblast_dir):
    file_out_name = os.path.join(genblast_dir, "annotation.gtf")
    file_out = open(file_out_name, "w+")
    genblast_extension = "_1.1c_2.3_s1_0_16_1"
    for root, dirs, files in os.walk(genblast_dir):
        for genblast_file in files:
            genblast_file = os.path.join(root, genblast_file)
            if genblast_file.endswith(".gff"):
                gtf_string = convert_gff_to_gtf(genblast_file)
                file_out.write(gtf_string)
            elif (
                genblast_file.endswith(".fa.blast")
                or genblast_file.endswith(".fa.blast.report")
                or genblast_file.endswith(genblast_extension)
            ):
                subprocess.run(["rm", genblast_file])
    file_out.close()


def convert_gff_to_gtf(gff_file):
    gtf_string = ""
    file_in = open(gff_file)
    line = file_in.readline()
    while line:
        #    match = re.search(r"genBlastG",line)
        #    if match:
        results = line.split()
        if not len(results) == 9:
            line = file_in.readline()
            continue
        if results[2] == "coding_exon":
            results[2] = "exon"
        attributes = set_attributes(results[8], results[2])
        results[8] = attributes
        converted_line = "\t".join(results)
        gtf_string += converted_line + "\n"
        line = file_in.readline()
    file_in.close()
    return gtf_string


def set_attributes(attributes, feature_type):

    converted_attributes = ""
    split_attributes = attributes.split(";")
    if feature_type == "transcript":
        match = re.search(r"Name\=(.+)$", split_attributes[1])
        name = match.group(1)
        converted_attributes = 'gene_id "' + name + '"; transcript_id "' + name + '";'
    elif feature_type == "exon":
        match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-", attributes)
        exon_rank = match.group(1)
        name = match.group(2)
        converted_attributes = (
            'gene_id "'
            + name
            + '"; transcript_id "'
            + name
            + '"; exon_number "'
            + exon_rank
            + '";'
        )

    return converted_attributes


# Example genBlast output
# 1       genBlastG       transcript      131128674       131137049       252.729 -       .       ID=259447-R1-1-A1;Name=259447;PID=84.65;Coverage=94.22;Note=PID:84.65-Cover:94.22
# 1       genBlastG       coding_exon     131137031       131137049       .       -       .       ID=259447-R1-1-A1-E1;Parent=259447-R1-1-A1
# 1       genBlastG       coding_exon     131136260       131136333       .       -       .       ID=259447-R1-1-A1-E2;Parent=259447-R1-1-A1
# 1       genBlastG       coding_exon     131128674       131130245       .       -       .       ID=259447-R1-1-A1-E3;Parent=259447-R1-1-A1
##sequence-region       1_group1        1       4534
# 1       genBlastG       transcript      161503457       161503804       30.94   +       .       ID=259453-R1-1-A1;Name=259453;PID=39.46;Coverage=64.97;Note=PID:39.46-Cover:64.97
# 1       genBlastG       coding_exon     161503457       161503804       .       +       .       ID=259453-R1-1-A1-E1;Parent=259453-R1-1-A1
##sequence-region       5_group1        1       4684
# 5       genBlastG       transcript      69461063        69461741        86.16   +       .       ID=259454-R1-1-A1;Name=259454;PID=82.02;Coverage=91.67;Note=PID:82.02-Cover:91.67
# 5       genBlastG       coding_exon     69461063        69461081        .       +       .       ID=259454-R1-1-A1-E1;Parent=259454-R1-1-A1
# 5       genBlastG       coding_exon     69461131        69461741        .       +       .       ID=259454-R1-1-A1-E2;Parent=259454-R1-1-A1


def split_protein_file(protein_file, protein_output_dir, batch_size):
    if batch_size is None:
        batch_size = 20

    batched_protein_files = []

    for i in range(0, 10):
        create_dir(protein_output_dir, ("bin_" + str(i)))

    file_in = open(protein_file)
    line = file_in.readline()
    seq_count = 0
    batch_count = 0
    current_record = ""
    initial_seq = 1
    while line:
        num_dir = random.randint(0, 9)
        match = re.search(r">(.+)$", line)
        if match and not initial_seq and seq_count % batch_size == 0:
            file_out_name = os.path.join(
                protein_output_dir,
                ("bin_" + str(random.randint(0, 9))),
                (str(batch_count) + ".fa"),
            )
            file_out = open(file_out_name, "w+")
            file_out.write(current_record)
            file_out.close()
            batch_count += 1
            seq_count += 1
            current_record = line
            batched_protein_files.append(file_out_name)
        elif match:
            current_record += line
            initial_seq = 0
            seq_count += 1
        else:
            current_record += line
        line = file_in.readline()
    file_in.close()

    if current_record:
        file_out_name = os.path.join(
            protein_output_dir,
            ("bin_" + str(random.randint(0, 9))),
            (str(batch_count) + ".fa"),
        )
        file_out = open(file_out_name, "w+")
        file_out.write(current_record)
        file_out.close()
        batched_protein_files.append(file_out_name)

    return batched_protein_files


def run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file):

    asnb_file = masked_genome_file + ".asnb"
    logging.info("Running convert2blastmask prior to GenBlast:")
    cmd = [
        convert2blastmask_path,
        "-in",
        masked_genome_file,
        "-parse_seqids",
        "-masking_algorithm",
        "other",
        "-masking_options",
        '"REpeatDetector, default"',
        "-outfmt",
        "maskinfo_asn1_bin",
        "-out",
        asnb_file,
    ]
    logging.info(" ".join(cmd))
    subprocess.run(cmd)
    logging.info("Completed running convert2blastmask")


def run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file):

    logging.info("Running makeblastdb prior to GenBlast")
    subprocess.run(
        [
            makeblastdb_path,
            "-in",
            masked_genome_file,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-mask_data",
            asnb_file,
            "-max_file_sz",
            "10000000000",
        ]
    )
    logging.info("Completed running makeblastdb")


def run_trimming(
    main_output_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads
):

    trim_galore_path = "trim_galore"
    check_exe(trim_galore_path)

    trim_dir = create_dir(main_output_dir, "trim_galore_output")

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    for file_type in file_types:
        fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    fastq_file_list = create_paired_paths(fastq_file_list)

    for fastq_file_path in fastq_file_list:
        logging.info("fastaq file path" + fastq_file_path)

    generic_trim_galore_cmd = [
        trim_galore_path,
        "--illumina",
        "--quality",
        "20",
        "--length",
        "50",
        "--output_dir",
        trim_dir,
    ]

    pool = multiprocessing.Pool(int(num_threads))
    for fastq_files in fastq_file_list:
        pool.apply_async(
            multiprocess_trim_galore,
            args=(
                generic_trim_galore_cmd,
                fastq_files,
                trim_dir,
            ),
        )
        if delete_pre_trim_fastq:
            for file_path in fastq_files:
                logging.info("Removing original fastq file post trimming:\n" + file_path)
                subprocess.run(["rm", file_path])

    pool.close()
    pool.join()

    trimmed_fastq_list = glob.glob(os.path.join(trim_dir, "*.fq.gz"))
    for trimmed_fastq_path in trimmed_fastq_list:
        logging.info("Trimmed file path:\n" + trimmed_fastq_path)
        sub_patterns = r"|".join(("_val_1.fq", "_val_2.fq", "_trimmed.fq"))
        updated_file_path = re.sub(sub_patterns, ".fq", trimmed_fastq_path)
        logging.info("Updated file path:\n" + updated_file_path)
        subprocess.run(["mv", trimmed_fastq_path, updated_file_path])

        files_to_delete_list = []
        for file_type in file_types:
            files_to_delete_list.extend(
                glob.glob(os.path.join(short_read_fastq_dir, file_type))
            )


def multiprocess_trim_galore(generic_trim_galore_cmd, fastq_files, trim_dir):

    fastq_file = fastq_files[0]
    fastq_file_pair = None
    if len(fastq_files) == 2:
        fastq_file_pair = fastq_files[1]

    trim_galore_cmd = generic_trim_galore_cmd
    if fastq_file_pair:
        trim_galore_cmd.append("--paired")

    trim_galore_cmd.append(fastq_file)

    if fastq_file_pair:
        trim_galore_cmd.append(fastq_file_pair)

    logging.info("Running Trim Galore with the following command:")
    logging.info(" ".join(trim_galore_cmd))
    subprocess.run(trim_galore_cmd)


def run_star_align(
    star_path,
    trim_fastq,
    subsample_script_path,
    main_output_dir,
    short_read_fastq_dir,
    genome_file,
    max_reads_per_sample,
    max_total_reads,
    max_short_read_intron_length,
    num_threads,
):
    # !!! Need to add in samtools path above instead of just using 'samtools' in command

    if not star_path:
        star_path = "STAR"

    check_exe(star_path)

    # If trimming has been enabled then switch the path for short_read_fastq_dir from the original location to the trimmed fastq dir
    if trim_fastq:
        short_read_fastq_dir = os.path.join(main_output_dir, "trim_galore_output")

    #  if not os.path.exists(subsample_script_path):
    subsample_script_path = "subsample_fastq.py"

    star_dir = create_dir(main_output_dir, "star_output")

    logging.info("Skip analysis if the final log file already exists")
    log_final_out = Path(os.path.join(star_dir, "Log.final.out"))
    log_out = Path(os.path.join(star_dir, "Log.out"))
    log_progress_out = Path(os.path.join(star_dir, "Log.progress.out"))
    if log_final_out.is_file() and log_out.is_file() and log_progress_out.is_file():
        logging.info("Star gtf file exists")
        return
    else:
        logging.info("No log files, go on with the analysis")

    star_tmp_dir = os.path.join(star_dir, "tmp")
    if os.path.exists(star_tmp_dir):
        subprocess.run(["rm", "-rf", star_tmp_dir])

    star_index_file = os.path.join(star_dir, "SAindex")

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    for file_type in file_types:
        fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    # This works out if the files are paired or not
    fastq_file_list = create_paired_paths(fastq_file_list)

    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        pool = multiprocessing.Pool(int(num_threads))
        for fastq_files in fastq_file_list:
            fastq_file = fastq_files[0]
            fastq_file_pair = ""
            if len(fastq_files) == 2:
                fastq_file_pair = fastq_files[1]

            if (
                fastq_file_pair
                and os.path.exists(fastq_file + ".sub")
                and os.path.exists(fastq_file_pair + ".sub")
            ):
                logging.info(
                    "Found an existing .sub files on the fastq path for both members of the pair, will use those instead of subsampling again. Files:"
                )
                logging.info(fastq_file + ".sub")
                logging.info(fastq_file_pair + ".sub")
            elif fastq_file_pair:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )
            elif os.path.exists(fastq_file + ".sub"):
                logging.info(
                    "Found an existing .sub file on the fastq path, will use that instead. File:"
                )
                logging.info(fastq_file + ".sub")
            else:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )

        pool.close()
        pool.join()

    fastq_file_list = check_for_fastq_subsamples(fastq_file_list)

    if not fastq_file_list:
        raise IndexError(
            "The list of fastq files is empty. Fastq dir:\n%s" % short_read_fastq_dir
        )

    if not os.path.exists(star_index_file):
        logging.info("Did not find an index file for Star. Will create now")
        seq_region_lengths = get_seq_region_lengths(genome_file, 0)
        genome_size = sum(seq_region_lengths.values())
        index_bases = min(14, math.floor((math.log(genome_size, 2) / 2) - 1))
        subprocess.run(
            [
                star_path,
                "--runThreadN",
                str(num_threads),
                "--runMode",
                "genomeGenerate",
                "--outFileNamePrefix",
                (star_dir + "/"),
                "--genomeDir",
                star_dir,
                "--genomeSAindexNbases",
                str(index_bases),
                "--genomeFastaFiles",
                genome_file,
            ]
        )

    if not star_index_file:
        raise IOError(
            "The index file does not exist. Expected path:\n%s" % star_index_file
        )

    logging.info("Running Star on the files in the fastq dir")
    for fastq_file_path in fastq_file_list:
        logging.info(fastq_file_path)
        fastq_file_name = os.path.basename(fastq_file_path)
        check_compression = re.search(r".gz$", fastq_file_name)

        # If there's a tmp dir already, the most likely cause is that STAR failed on the previous input file(s)
        # In this case STAR would effectively break for all the rest of the files, as it won't run if the tmp
        # dir exists. So clean it up and put out a warning. Another approach would be to just name the tmp dir
        # uniquely. Also there should be code checking the return on STAR anyway. So later when this is being
        # cleaned up to have proper tests, need to decide on the best implementation
        if os.path.exists(star_tmp_dir):
            logging.error(
                "Found an existing tmp dir, implies potential failure on previous file. Removing tmp dir"
            )
            try:
                shutil.rmtree(star_tmp_dir)
            except OSError as e:
                logging.error("Error: %s - %s." % (e.filename, e.strerror))

        sam_file_path = os.path.join(star_dir, (fastq_file_name + ".sam"))
        junctions_file_path = os.path.join(star_dir, (fastq_file_name + ".sj.tab"))
        sam_file_name = os.path.basename(sam_file_path)
        sam_temp_file_path = os.path.join(star_dir, (sam_file_name + ".tmp"))
        bam_sort_file_path = os.path.join(star_dir, re.sub(".sam", ".bam", sam_file_name))

        if (
            os.path.isfile(bam_sort_file_path)
            and not os.stat(bam_sort_file_path).st_size == 0
        ):
            logging.info(
                "Found an existing bam file for the fastq file, presuming the file has been processed, will skip"
            )
            continue

        logging.info("Processing %s" % fastq_file_path)
        #    star_command = [star_path,'--outFilterIntronMotifs','RemoveNoncanonicalUnannotated','--outSAMstrandField','intronMotif','--runThreadN',str(num_threads),'--twopassMode','Basic','--runMode','alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file_path,'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir,'--outSAMtype','SAM','--alignIntronMax',str(max_intron_length),'--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000','50000','50000','50000','50000','50000','100000']

        star_command = [
            star_path,
            "--outFilterIntronMotifs",
            "RemoveNoncanonicalUnannotated",
            "--outSAMstrandField",
            "intronMotif",
            "--runThreadN",
            str(num_threads),
            "--twopassMode",
            "Basic",
            "--runMode",
            "alignReads",
            "--genomeDir",
            star_dir,
            "--readFilesIn",
            fastq_file_path,
            "--outFileNamePrefix",
            (star_dir + "/"),
            "--outTmpDir",
            star_tmp_dir,
            "--outSAMtype",
            "SAM",
            "--alignIntronMax",
            str(max_intron_length),
        ]

        if check_compression:
            star_command.append("--readFilesCommand")
            star_command.append("gunzip")
            star_command.append("-c")

        subprocess.run(star_command)
        subprocess.run(["mv", os.path.join(star_dir, "Aligned.out.sam"), sam_file_path])
        subprocess.run(["mv", os.path.join(star_dir, "SJ.out.tab"), junctions_file_path])

        logging.info("Converting samfile into sorted bam file. Bam file:")
        logging.info(bam_sort_file_path)
        subprocess.run(
            [
                "samtools",
                "sort",
                "-@",
                str(num_threads),
                "-T",
                sam_temp_file_path,
                "-o",
                bam_sort_file_path,
                sam_file_path,
            ]
        )

        logging.info("Removing sam file")
        subprocess.run(["rm", sam_file_path])

    logging.info("Completed running STAR")


def run_subsample_script(fastq_file, fastq_file_pair, subsample_script_path):

    if fastq_file_pair:
        subprocess.run(
            [
                "python3",
                subsample_script_path,
                "--fastq_file",
                fastq_file,
                "--fastq_file_pair",
                fastq_file_pair,
            ]
        )
    else:
        subprocess.run(["python3", subsample_script_path, "--fastq_file", fastq_file])


def check_for_fastq_subsamples(fastq_file_list):
    # This should probably removed at some point as it is needlessly complicated
    # Would be better to just build into the previous step
    # Mainly just about making sure that if you have subsamples they're used and if you have pairs they're paired
    for idx, fastq_files in enumerate(fastq_file_list):
        fastq_file = fastq_files[0]
        subsample_file = fastq_file + ".sub"

        fastq_file_pair = ""
        subsample_file_pair = ""
        if len(fastq_files) == 2:
            fastq_file_pair = fastq_files[1]
            subsample_file_pair = fastq_file_pair + ".sub"

        # This bit will replace the list entry with a string, don't need a list after this function for each pair/file
        if os.path.exists(subsample_file):
            logging.info(
                "Found a subsampled file extension, will use that instead of the original file. Path:"
            )
            logging.info(subsample_file)
            fastq_file_list[idx] = subsample_file
        else:
            fastq_file_list[idx] = fastq_file

        # This bit just concats the paired file (or subsampled paired file) if it exists
        if os.path.exists(subsample_file_pair):
            logging.info(
                "Found a subsampled paired file extension, will use that instead of the original file. Path:"
            )
            logging.info(subsample_file_pair)
            fastq_file_list[idx] = subsample_file + "," + subsample_file_pair
        elif fastq_file_pair:
            fastq_file_list[idx] = fastq_file + "," + fastq_file_pair

        logging.info("Entry at current index:")
        logging.info(fastq_file_list[idx])

    return fastq_file_list


def run_minimap2_align(
    minimap2_path,
    paftools_path,
    main_output_dir,
    long_read_fastq_dir,
    genome_file,
    max_intron_length,
    num_threads,
):

    if not minimap2_path:
        minimap2_path = "minimap2"

    check_exe(minimap2_path)

    if not paftools_path:
        paftools_path = "paftools.js"

    check_exe(paftools_path)

    minimap2_output_dir = create_dir(main_output_dir, "minimap2_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(minimap2_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Minimap2 gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    genome_file_name = os.path.basename(genome_file)
    genome_file_index = genome_file_name + ".mmi"
    minimap2_index_file = os.path.join(minimap2_output_dir, genome_file_index)
    minimap2_hints_file = os.path.join(minimap2_output_dir, "minimap2_hints.gff")

    fastq_file_list = []
    for fastq_file in glob.glob(long_read_fastq_dir + "/*.fastq"):
        fastq_file_list.append(fastq_file)

    for fastq_file in glob.glob(long_read_fastq_dir + "/*.fq"):
        fastq_file_list.append(fastq_file)

    if not fastq_file_list:
        # NOTE: should update this to have a param that says it's okay to be empty if there's an override param
        # This is because it's important that an external user would have an exception if they provided a long read dir
        # but there was nothing in it. Whereas in the Ensembl pipeline there might not be long read data, but we don't
        # know this when the command line is being constructed. This is also true for the short read data. An alternative
        # would be to put in another analysis into the Ensembl pipeline to construct this commandline after the transcriptomic
        # data has been searched for
        return
    #    raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % long_read_fastq_dir)

    if not os.path.exists(minimap2_index_file):
        logging.info("Did not find an index file for minimap2. Will create now")
        subprocess.run(
            [
                minimap2_path,
                "-t",
                str(num_threads),
                "-d",
                os.path.join(minimap2_index_file),
                genome_file,
            ]
        )

    if not minimap2_index_file:
        raise IOError(
            "The minimap2 index file does not exist. Expected path:\n%s"
            % minimap2_index_file
        )

    logging.info("Running minimap2 on the files in the long read fastq dir")
    for fastq_file_path in fastq_file_list:
        fastq_file_name = os.path.basename(fastq_file_path)
        sam_file = os.path.join(minimap2_output_dir, (fastq_file_name + ".sam"))
        bed_file = os.path.join(minimap2_output_dir, (fastq_file_name + ".bed"))
        bed_file_out = open(bed_file, "w+")
        logging.info("Processing %s" % fastq_file)
        subprocess.run(
            [
                minimap2_path,
                "-G",
                str(max_intron_length),
                "-t",
                str(num_threads),
                "--cs",
                "--secondary=no",
                "-ax",
                "splice",
                "-u",
                "b",
                minimap2_index_file,
                fastq_file_path,
                "-o",
                sam_file,
            ]
        )
        logging.info("Creating bed file from SAM")
        subprocess.run([paftools_path, "splice2bed", sam_file], stdout=bed_file_out)
        bed_file_out.close()

    bed_to_gtf(minimap2_output_dir)

    logging.info("Completed running minimap2")


def bed_to_gtf(minimap2_output_dir):

    gtf_file_path = os.path.join(minimap2_output_dir, "annotation.gtf")
    gtf_out = open(gtf_file_path, "w+")
    exons_dict = {}
    gene_id = 1
    for bed_file in glob.glob(minimap2_output_dir + "/*.bed"):
        logging.info("Converting bed to GTF:")
        logging.info(bed_file)
        bed_in = open(bed_file)
        bed_lines = bed_in.readlines()
        for line in bed_lines:
            line = line.rstrip()
            elements = line.split("\t")
            seq_region_name = elements[0]
            offset = int(elements[1])
            hit_name = elements[3]
            strand = elements[5]
            block_sizes = elements[10].split(",")
            block_sizes = list(filter(None, block_sizes))
            block_starts = elements[11].split(",")
            block_starts = list(filter(None, block_starts))
            exons = bed_to_exons(block_sizes, block_starts, offset)
            transcript_line = [
                seq_region_name,
                "minimap",
                "transcript",
                0,
                0,
                ".",
                strand,
                ".",
                'gene_id "minimap_'
                + str(gene_id)
                + '"; '
                + 'transcript_id "minimap_'
                + str(gene_id)
                + '"',
            ]
            transcript_start = None
            transcript_end = None
            exon_records = []
            for i, exon_coords in enumerate(exons):
                if transcript_start is None or exon_coords[0] < transcript_start:
                    transcript_start = exon_coords[0]

                if transcript_end is None or exon_coords[1] > transcript_end:
                    transcript_end = exon_coords[1]

                exon_line = [
                    seq_region_name,
                    "minimap",
                    "exon",
                    str(exon_coords[0]),
                    str(exon_coords[1]),
                    ".",
                    strand,
                    ".",
                    'gene_id "minimap_'
                    + str(gene_id)
                    + '"; '
                    + 'transcript_id "minimap_'
                    + str(gene_id)
                    + '"; exon_number "'
                    + str(i + 1)
                    + '";',
                ]

                exon_records.append(exon_line)

            transcript_line[3] = str(transcript_start)
            transcript_line[4] = str(transcript_end)

            gtf_out.write("\t".join(transcript_line) + "\n")
            for exon_line in exon_records:
                gtf_out.write("\t".join(exon_line) + "\n")

            gene_id += 1

    gtf_out.close()


def bed_to_gff(input_dir, hints_file):

    gff_out = open(hints_file, "w+")
    exons_dict = {}
    for bed_file in glob.glob(input_dir + "/*.bed"):
        logging.info("Processing file for hints:")
        logging.info(bed_file)
        bed_in = open(bed_file)
        bed_lines = bed_in.readlines()
        for line in bed_lines:
            line = line.rstrip()
            elements = line.split("\t")
            seq_region_name = elements[0]
            offset = int(elements[1])
            hit_name = elements[3]
            strand = elements[5]
            block_sizes = elements[10].split(",")
            block_sizes = list(filter(None, block_sizes))
            block_starts = elements[11].split(",")
            block_starts = list(filter(None, block_starts))
            exons = bed_to_exons(block_sizes, block_starts, offset)
            for i, element in enumerate(exons):
                exon_coords = exons[i]
                exon_key = (
                    seq_region_name
                    + ":"
                    + exon_coords[0]
                    + ":"
                    + exon_coords[1]
                    + ":"
                    + strand
                )
                if exon_key in exons_dict:
                    exons_dict[exon_key][5] += 1
                else:
                    gff_list = [
                        seq_region_name,
                        "CDNA",
                        "exon",
                        exon_coords[0],
                        exon_coords[1],
                        1.0,
                        strand,
                        ".",
                    ]
                    exons_dict[exon_key] = gff_list

    for exon_key, gff_list in exons_dict.items():
        gff_list[5] = str(gff_list[5])
        gff_line = "\t".join(gff_list) + "\tsrc=W;mul=" + gff_list[5] + ";\n"
        gff_out.write(gff_line)

    gff_out.close()

    sorted_hints_out = open((hints_file + ".srt"), "w+")
    subprocess.run(
        ["sort", "-k1,1", "-k7,7", "-k4,4", "-k5,5", hints_file], stdout=sorted_hints_out
    )
    sorted_hints_out.close()


def bed_to_exons(block_sizes, block_starts, offset):
    exons = []
    for i, element in enumerate(block_sizes):
        block_start = offset + int(block_starts[i]) + 1
        block_end = block_start + int(block_sizes[i]) - 1

        if block_end < block_start:
            logging.warning("Warning: block end is less than block start, skipping exon")
            continue

        exon_coords = [str(block_start), str(block_end)]
        exons.append(exon_coords)

    return exons


def check_transcriptomic_output(main_output_dir):

    # This will check across the various transcriptomic dirs and check there's actually some data
    transcriptomic_dirs = ["scallop_output", "stringtie_output", "minimap2_output"]
    total_lines = 0
    min_lines = 100000
    for transcriptomic_dir in transcriptomic_dirs:
        full_file_path = os.path.join(
            main_output_dir, transcriptomic_dir, "annotation.gtf"
        )
        if not os.path.exists(full_file_path):
            logging.warning(
                "Warning, no annotation.gtf found for "
                + transcriptomic_dir
                + ". This might be fine, e.g. no long read data were provided"
            )
            continue
        num_lines = sum(1 for line in open(full_file_path))
        total_lines = total_lines + num_lines
        logging.info(
            "For "
            + transcriptomic_dir
            + " found a total of "
            + str(num_lines)
            + " in the annotation.gtf file"
        )
    if total_lines == 0:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, but the transcriptomic annotation files are empty"
        )
    elif total_lines <= min_lines:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, but the total number of lines in the output files were less than the min expected value"
            + "\n"
            "Found: " + str(total_lines) + "\n"
            "Min allowed: " + str(min_lines)
        )

    else:
        logging.info(
            "Found "
            + str(total_lines)
            + " total lines across the transcriptomic files. Checks passed"
        )


# start gene g1
# 1       AUGUSTUS        gene    1       33908   1       +       .       g1
# 1       AUGUSTUS        transcript      1       33908   .       +       .       g1.t1
# 1       AUGUSTUS        CDS     3291    3585    .       +       2       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    3291    3585    .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     11377   11510   .       +       1       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    11377   11510   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     30726   30871   .       +       2       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    30726   30871   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     32975   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    32975   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        stop_codon      33500   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        tts     33908   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# protein sequence = [GGRGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEKKEKKEEEEEEEEEEEEEEEK
# EEGEGMEGRDIMGIRGKQKQPRKQPELSSSFPTFLLTQKTPYCPESIKLLLILRNNIQTIVFYKGPKGVINDWRKFKLESEDGDSIPPSKKEILRQMS
# SPQSRDDSKERMSRKMSIQEYELIHQDKEDESCLRKYRRQCMQDMHQKLSFGPRYGFVYELETGEQFLETIEKEQKVTTIVVNIYEDGVRGCDALNSS
# LACLAVEYPMVKFCKIKASNTGAGDRFSTDVLPTLLVYKGGELISNFISVAEQFAEEFFAVDVESFLNEYGLLPEREIHDLEQTNMEDEDIE]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 0
# CDS exons: 0/4
# CDS introns: 0/4
# 5'UTR exons and introns: 0/0
# 3'UTR exons and introns: 0/1
# hint groups fully obeyed: 0
# incompatible hint groups: 0
# end gene g1
###


def augustus_output_to_gtf(augustus_output_dir, augustus_genome_dir):

    gtf_file_path = os.path.join(augustus_output_dir, "annotation.gtf")
    gtf_out = open(gtf_file_path, "w+")
    record_count = 1
    for gff_file_path in glob.glob(augustus_genome_dir + "/*.aug"):
        gff_file_name = os.path.basename(gff_file_path)
        match = re.search(r"\.rs(\d+)\.re(\d+)\.", gff_file_name)
        start_offset = int(match.group(1))

        exon_number = 1
        current_exon_hints_total = 0
        current_exon_hints_match = 0
        current_intron_hints_total = 0
        current_intron_hints_match = 0
        current_record = []
        gff_in = open(gff_file_path, "r")
        line = gff_in.readline()
        while line:
            match = re.search(r"# CDS exons\: (\d+)\/(\d+)", line)
            if match:
                current_exon_hints_match = match.group(1)
                current_exon_hints_total = match.group(2)

            match = re.search(r"# CDS introns\: (\d+)\/(\d+)", line)
            if match:
                current_introns_hints_match = match.group(1)
                current_introns_hints_total = match.group(2)

            if re.search(r"# end gene", line):
                for output_line in current_record:
                    gtf_out.write(output_line)
                current_record = []
                current_exon_hints_total = 0
                current_exon_hints_match = 0
                current_intron_hints_total = 0
                current_intron_hints_match = 0
                record_count += 1
                exon_number = 1

            values = line.split("\t")
            if (
                len(values) == 9
                and values[1] == "AUGUSTUS"
                and (values[2] == "transcript" or values[2] == "exon")
            ):
                values[3] = str(int(values[3]) + (start_offset - 1))
                values[4] = str(int(values[4]) + (start_offset - 1))
                values[8] = (
                    'gene_id "aug'
                    + str(record_count)
                    + '"; transcript_id "aug'
                    + str(record_count)
                    + '";'
                )
                if values[2] == "exon":
                    values[8] = values[8] + ' exon_number "' + str(exon_number) + '";'
                    exon_number += 1
                values[8] = values[8] + "\n"
                current_record.append("\t".join(values))

            line = gff_in.readline()
        gff_in.close()
    gtf_out.close()


def run_augustus_predict(augustus_path, main_output_dir, masked_genome_file, num_threads):

    min_seq_length = 1000

    if not augustus_path:
        augustus_path = "/hps/nobackup2/production/ensembl/jma/src/Augustus/bin/augustus"

    bam2hints_path = "/homes/fergal/bin/bam2hints"
    bam2wig_path = "/homes/fergal/bin/bam2wig"
    wig2hints_path = "/homes/fergal/bin/wig2hints"
    check_exe(augustus_path)

    # Run bam2hints, bam2wig, wig2hints, then combine the hints into a single file
    # Multiprocess with all three steps in the MP as that would be fastest

    augustus_dir = create_dir(main_output_dir, "augustus_output")
    augustus_hints_dir = create_dir(augustus_dir, "hints")
    augustus_genome_dir = create_dir(augustus_dir, "genome_dir")
    augustus_evidence_dir = create_dir(augustus_dir, "evidence")
    augustus_hints_file = os.path.join(augustus_evidence_dir, "augustus_hints.gff")
    star_dir = os.path.join(main_output_dir, "star_output")
    minimap2_output_dir = os.path.join(main_output_dir, "minimap2_output")

    if os.path.exists(star_dir):
        logging.info("Found a Star output dir, generating hints from any .sj.tab files")
        generate_hints(
            bam2hints_path,
            bam2wig_path,
            wig2hints_path,
            augustus_hints_dir,
            star_dir,
            num_threads,
        )
        hints_out = open(augustus_hints_file, "w+")
        for gff_file in glob.glob(augustus_hints_dir + "/*.bam.hints.gff"):
            gff_in = open(gff_file, "r")
            line = gff_in.readline()
            while line:
                hints_out.write(line)
                line = gff_in.readline()
            gff_in.close()
        hints_out.close()

    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 100000, 5000)

    generic_augustus_cmd = [
        augustus_path,
        "--species=human",
        "--UTR=on",
        (
            "--extrinsicCfgFile="
            + "/hps/nobackup2/production/ensembl/jma/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg"
        ),
    ]

    pool = multiprocessing.Pool(int(num_threads))
    tasks = []

    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_augustus_id,
            args=(
                generic_augustus_cmd,
                slice_id,
                masked_genome_file,
                augustus_hints_file,
                augustus_genome_dir,
            ),
        )

    pool.close()
    pool.join()

    augustus_output_to_gtf(augustus_dir, augustus_genome_dir)


def create_slice_ids(seq_region_lengths, slice_size, overlap, min_length):
    if not slice_size:
        slice_size = 1000000

    if not overlap:
        overlap = 0

    if not min_length:
        min_length = 0

    slice_ids = []

    for region in seq_region_lengths:
        region_length = int(seq_region_lengths[region])
        if region_length < min_length:
            continue

        if region_length <= slice_size:
            slice_ids.append([region, 1, region_length])
            continue

        start = 1
        end = start + slice_size - 1
        while end < region_length:
            start = start - overlap
            if start < 1:
                start = 1

            end = start + slice_size - 1
            if end > region_length:
                end = region_length
            if (end - start + 1) >= min_length:
                slice_ids.append([region, start, end])
            start = end + 1

    return slice_ids


def generate_hints(
    bam2hints_path,
    bam2wig_path,
    wig2hints_path,
    augustus_hints_dir,
    star_dir,
    num_threads,
):

    pool = multiprocessing.Pool(int(num_threads))
    for bam_file in glob.glob(star_dir + "/*.bam"):
        pool.apply_async(
            multiprocess_augustus_hints,
            args=(
                bam2hints_path,
                bam2wig_path,
                wig2hints_path,
                bam_file,
                augustus_hints_dir,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_augustus_hints(
    bam2hints_path, bam2wig_path, wig2hints_path, bam_file, augustus_hints_dir
):
    bam_file_name = os.path.basename(bam_file)
    logging.info("Processing " + bam_file_name + " for Augustus hints")

    bam2hints_file_name = bam_file_name + ".hints.gff"
    bam2hints_file_path = os.path.join(augustus_hints_dir, bam2hints_file_name)
    bam2hints_cmd = [
        bam2hints_path,
        ("--in=" + bam_file),
        ("--out=" + bam2hints_file_path),
        "--maxintronlen=100000",
    ]
    logging.info("bam2hints command:\n" + " ".join(bam2hints_cmd))
    subprocess.run(bam2hints_cmd)


#  bam2wig_cmd = [bam2wig_path,'-D',augustus_hints_dir,bam_file]
#  print("bam2wig command:\n" + ' '.join(bam2wig_cmd))
#  subprocess.run(bam2wig_cmd)


# wig2hints is odd in that it runs directly off STDIN and then just prints to STDOUT,
# so the code below is implemented in steps as it's not good practice to use pipes and
# redirects in a subprocess command
#  wig_file_name = re.sub('.bam','.wig',bam_file_name)
#  wig_file_path = os.path.join(augustus_hints_dir,wig_file_name)
#  wig_hints_file_name = (wig_file_name + '.hints.gff')
#  wig_hints_file_path =  os.path.join(augustus_hints_dir,wig_hints_file_name)
#  print("Writing wig file info to hints file:\n" + wig_hints_file_name)
#  wig2hints_out = open(wig_hints_file_path,'w+')
#  wigcat = subprocess.Popen(('cat',wig_file_path), stdout=subprocess.PIPE)
#  subprocess.run(wig2hints_path, stdin=wigcat.stdout, stdout=wig2hints_out)
#  wig2hints_out.close()


def multiprocess_augustus_id(cmd, slice_id, genome_file, hints_file, output_dir):

    region = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]
    seq = get_sequence(region, start, end, 1, genome_file, output_dir)

    region_fasta_file_name = region + ".rs" + str(start) + ".re" + str(end) + ".fa"
    region_fasta_file_path = os.path.join(output_dir, region_fasta_file_name)
    region_augustus_file_path = os.path.join(
        output_dir, (region_fasta_file_name + ".aug")
    )

    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region + "\n" + seq + "\n")
    region_fasta_out.close()

    region_hints_file = create_slice_hints_file(
        region, start, end, hints_file, region_fasta_file_path
    )

    aug_out = open(region_augustus_file_path, "w+")

    augustus_forward = cmd.copy()
    augustus_forward.append(("--hintsfile=" + region_hints_file))
    augustus_forward.append("--strand=forward")
    augustus_forward.append(region_fasta_file_path)
    subprocess.run(augustus_forward, stdout=aug_out)

    augustus_backward = cmd.copy()
    augustus_backward.append(("--hintsfile=" + region_hints_file))
    augustus_backward.append("--strand=backward")
    augustus_backward.append(region_fasta_file_path)
    subprocess.run(augustus_backward, stdout=aug_out)

    aug_out.close()


def create_slice_hints_file(region, start, end, hints_file, region_fasta_file_path):

    # Note this is trying to be memory and file efficient at the cost of speed
    # So files are only created as needed and the hints are being read line by line as written as needed
    # This comes with the downside of being slow, but it's only a very small amount of time relative
    # to how slow the step is in total. Given that this step in general eats up a low of memory, saving as much
    # as possible here is not a bad thing even if it's adding in an overhead by continuously reading the hints file

    region_hints_file_path = region_fasta_file_path + ".gff"
    hints_in = open(hints_file)
    hints_out = open(region_hints_file_path, "w+")
    hint_line = hints_in.readline()
    while hint_line:
        hint_line_values = hint_line.split("\t")
        if not len(hint_line_values) == 9:
            hint_line = hints_in.readline()
            continue

        hint_region = hint_line_values[0]
        hint_region_start = int(hint_line_values[3])
        hint_region_end = int(hint_line_values[4])

        if (
            hint_region == region
            and hint_region_start >= start
            and hint_region_end <= end
        ):
            hint_line_values[3] = str(int(hint_line_values[3]) - (start - 1))
            hint_line_values[4] = str(int(hint_line_values[4]) - (start - 1))
            hints_out.write("\t".join(hint_line_values))

        hint_line = hints_in.readline()
    hints_in.close()
    hints_out.close()

    return region_hints_file_path


def run_stringtie_assemble(
    stringtie_path, samtools_path, main_output_dir, genome_file, num_threads
):

    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    check_exe(stringtie_path)

    if not samtools_path:
        samtools_path = shutil.which("samtools")
    check_exe(samtools_path)

    stringtie_dir = create_dir(main_output_dir, "stringtie_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(stringtie_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Stringtie gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    stringtie_merge_input_file = os.path.join(stringtie_dir, "stringtie_assemblies.txt")
    stringtie_merge_output_file = os.path.join(stringtie_dir, "annotation.gtf")
    star_dir = os.path.join(main_output_dir, "star_output")

    if os.path.exists(star_dir):
        logging.info("Found a Star output dir, will load sam file")

    sorted_bam_files = []
    for bam_file in glob.glob(star_dir + "/*.bam"):
        sorted_bam_files.append(bam_file)

    if not sorted_bam_files:
        raise IndexError(
            "The list of sorted bam files is empty, expected them in Star output dir. Star dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was fast enough in serial. But consider multiprocessing if
    # the mem usage is low
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = os.path.basename(sorted_bam_file)
        transcript_file_name = re.sub(".bam", ".stringtie.gtf", sorted_bam_file_name)
        transcript_file_path = os.path.join(stringtie_dir, transcript_file_name)

        if os.path.exists(transcript_file_path):
            logging.info(
                "Found an existing stringtie gtf file, will not overwrite. File found:"
            )
            logging.info(transcript_file_path)
        else:
            logging.info("Running Stringtie on: " + sorted_bam_file_name)
            logging.info("Writing output to: " + transcript_file_path)
            subprocess.run(
                [
                    stringtie_path,
                    sorted_bam_file,
                    "-o",
                    transcript_file_path,
                    "-p",
                    str(num_threads),
                    "-t",
                    "-a",
                    "15",
                ]
            )

    # Now need to merge
    logging.info("Creating Stringtie merge input file: " + stringtie_merge_input_file)

    # Now need to merge
    logging.info("Creating Stringtie merge input file: " + stringtie_merge_input_file)
    gtf_list_out = open(stringtie_merge_input_file, "w+")
    for gtf_file in glob.glob(stringtie_dir + "/*.stringtie.gtf"):
        transcript_count = check_gtf_content(gtf_file, "transcript")
        if transcript_count > 0:
            gtf_list_out.write(gtf_file + "\n")
        else:
            logging.warning(
                "Warning, skipping file with no transcripts. Path:\n" + gtf_file
            )
    gtf_list_out.close()

    logging.info("Merging Stringtie results. Writing to the following file:")
    logging.info(stringtie_merge_output_file)
    subprocess.run(
        [
            stringtie_path,
            "--merge",
            "-o",
            stringtie_merge_output_file,
            stringtie_merge_input_file,
        ]
    )


def run_scallop_assemble(scallop_path, stringtie_path, main_output_dir):

    if not scallop_path:
        scallop_path = shutil.which("scallop")
    check_exe(scallop_path)

    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    check_exe(stringtie_path)

    scallop_dir = create_dir(main_output_dir, "scallop_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(scallop_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Scallop gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")

    stringtie_merge_input_file = os.path.join(scallop_dir, "scallop_assemblies.txt")
    stringtie_merge_output_file = os.path.join(scallop_dir, "annotation.gtf")
    star_dir = os.path.join(main_output_dir, "star_output")
    memory_limit = 40 * 1024**3

    if os.path.exists(star_dir):
        logging.info("Found a Star output dir, will load bam files")

    sorted_bam_files = []
    for bam_file in glob.glob(star_dir + "/*.bam"):
        sorted_bam_files.append(bam_file)

    if not sorted_bam_files:
        raise IndexError(
            "The list of sorted bam files is empty, expected them in Star output dir. Star dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was fast enough in serial. But consider multiprocessing if
    # the mem usage is low
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = os.path.basename(sorted_bam_file)
        transcript_file_name = re.sub(".bam", ".scallop.gtf", sorted_bam_file_name)
        transcript_file_path = os.path.join(scallop_dir, transcript_file_name)

        if os.path.exists(transcript_file_path):
            logging.info(
                "Found an existing scallop gtf file, will not overwrite. File found:"
            )
            logging.info(transcript_file_path)
        else:
            logging.info("Running Scallop on: " + sorted_bam_file_name)
            logging.info("Writing output to: " + transcript_file_path)
            scallop_cmd = [
                scallop_path,
                "-i",
                sorted_bam_file,
                "-o",
                transcript_file_path,
                "--min_flank_length",
                "10",
            ]
            if memory_limit is not None:
                scallop_cmd = prlimit_command(scallop_cmd, memory_limit)

            return_value = None
            try:
                return_value = subprocess.check_output(scallop_cmd)
            except subprocess.CalledProcessError as ex:
                logging.error("Issue processing the following region with scallop")
                logging.error("Return value: " + str(return_value))

    #      subprocess.run([scallop_path,'-i',sorted_bam_file,'-o',transcript_file_path,'--min_flank_length','10'])

    # Now need to merge
    logging.info("Creating Stringtie merge input file: " + stringtie_merge_input_file)

    gtf_list_out = open(stringtie_merge_input_file, "w+")
    for gtf_file in glob.glob(scallop_dir + "/*.scallop.gtf"):
        transcript_count = check_gtf_content(gtf_file, "transcript")
        if transcript_count > 0:
            gtf_list_out.write(gtf_file + "\n")
        else:
            logging.warning(
                "Warning, skipping file with no transcripts. Path:\n" + gtf_file
            )
    gtf_list_out.close()

    logging.info("Merging Scallop results. Writing to the following file:")
    logging.info(stringtie_merge_output_file)
    subprocess.run(
        [
            stringtie_path,
            "--merge",
            "-o",
            stringtie_merge_output_file,
            stringtie_merge_input_file,
        ]
    )


def check_gtf_content(gtf_file, content_obj):
    logging.info("check gtf transcript function")
    # This just checks how many transcript lines are in a GTF
    transcript_count = 0
    gtf_in = open(gtf_file)
    line = gtf_in.readline()
    while line:
        eles = line.split("\t")
        if not len(eles) == 9:
            line = gtf_in.readline()
            continue
        if eles[2] == content_obj:
            transcript_count += 1
        line = gtf_in.readline()
    gtf_in.close()
    logging.info(transcript_count)
    return transcript_count


def splice_junction_to_gff(input_dir, hints_file):

    sjf_out = open(hints_file, "w+")

    for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
        sjf_in = open(sj_tab_file)
        sjf_lines = sjf_in.readlines()
        for line in sjf_lines:
            elements = line.split("\t")
            strand = "+"
            # If the strand is undefined then skip, Augustus expects a strand
            if elements[3] == "0":
                continue
            elif elements[3] == "2":
                strand = "-"

            junction_length = int(elements[2]) - int(elements[1]) + 1
            if junction_length < 100:
                continue

            if not elements[4] and elements[7] < 10:
                continue

            # For the moment treat multimapping and single mapping things as a combined score
            score = float(elements[6]) + float(elements[7])
            score = str(score)
            output_line = [
                elements[0],
                "RNASEQ",
                "intron",
                elements[1],
                elements[2],
                score,
                strand,
                ".",
                ("src=W;mul=" + score + ";"),
            ]
            sjf_out.write("\t".join(output_line) + "\n")

    sjf_out.close()


def model_builder(work_dir):

    star_output_dir = os.path.join(work_dir, "star_output")

    all_junctions_file = os.path.join(star_output_dir, "all_junctions.sj")
    sjf_out = open(all_junctions_file, "w+")

    for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
        sjf_in = open(sj_tab_file)
        sjf_lines = sjf_in.readlines()
        for line in sjf_lines:
            elements = line.split("\t")
            strand = "+"

            #    my $slice_name = $eles[0];
            #    my $start = $eles[1];
            #    my $end = $eles[2];
            #    my $strand = $eles[3];

            # If the strand is undefined then skip, Augustus expects a strand
            if elements[3] == "0":
                continue
            elif elements[3] == "2":
                strand = "-"

            junction_length = int(elements[2]) - int(elements[1]) + 1
            if junction_length < 100:
                continue

            if not elements[4] and elements[7] < 10:
                continue

            # For the moment treat multimapping and single mapping things as a combined score
            score = float(elements[6]) + float(elements[7])
            score = str(score)
            output_line = [
                elements[0],
                "RNASEQ",
                "intron",
                elements[1],
                elements[2],
                score,
                strand,
                ".",
                ("src=W;mul=" + score + ";"),
            ]
            sjf_out.write("\t".join(output_line) + "\n")

    sjf_out.close()


def split_genome(genome_file, target_dir, min_seq_length):
    # This is the lazy initial way of just splitting into a dir of files based on the toplevel sequence with a min sequence length filter
    # There are a couple of obvious improvements:
    # 1) Instead of making files for all seqs, just process N seqs parallel, where N = num_threads. Then you could clean up the seq file
    #    after each seq finishes, thus avoiding potentially having thousands of file in a dir
    # 2) Split the seq into even slices and process these in parallel (which the same cleanup as in 1). For sequences smaller than the
    #    target slice size, bundle them up together into a single file. Vastly more complex, partially implemented in the splice_genome method
    #    Allows for more consistency with parallelisation (since there should be no large outliers). But require a mapping strategy for the
    #    coords and sequence names and all the hints will need to be adjusted
    current_header = ""
    current_seq = ""

    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">(.+)$", line)
        if match and current_header:
            if len(current_seq) > min_seq_length:
                file_out_name = os.path.join(target_dir, (current_header + ".split.fa"))
                if not os.path.exists(file_out_name):
                    file_out = open(file_out_name, "w+")
                    file_out.write(">" + current_header + "\n" + current_seq + "\n")
                    file_out.close()

                else:
                    print(
                        "Found an existing split file, so will not overwrite. File found:"
                    )
                    print(file_out_name)

            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

        line = file_in.readline()

    if len(current_seq) > min_seq_length:
        file_out_name = os.path.join(target_dir, (current_header + ".split.fa"))
        if not os.path.exists(file_out_name):
            file_out = open(file_out_name, "w+")
            file_out.write(">" + current_header + "\n" + current_seq + "\n")
            file_out.close()

        else:
            logging.info(
                "Found an existing split file, so will not overwrite. File found:"
            )
            logging.info(file_out_name)

    file_in.close()


def get_seq_region_lengths(genome_file, min_seq_length):
    current_header = ""
    current_seq = ""

    seq_regions = {}
    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">(.+)$", line)
        if match and current_header:
            if len(current_seq) > min_seq_length:
                seq_regions[current_header] = len(current_seq)

            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

        line = file_in.readline()

    if len(current_seq) > min_seq_length:
        seq_regions[current_header] = len(current_seq)

    return seq_regions


def run_finalise_geneset(
    main_script_dir,
    main_output_dir,
    genome_file,
    seq_region_names,
    validation_type,
    diamond_validation_db,
    num_threads,
):

    if validation_type is None:
        logging.info("Setting validation type to relaxed")
    else:
        logging.info("Setting validation type to " + validation_type)

    final_annotation_dir = create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = create_dir(final_annotation_dir, "final_region_gtfs")
    utr_region_annotation_dir = create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = get_seq_region_lengths(genome_file, 0)

    # This used to be a list of output dirs and a loop which was neat, I'm coverting to a list of conditions as
    # it's more straightforward with the renaming and having to merge scallop and stringtie
    protein_annotation_raw = os.path.join(
        main_output_dir, "genblast_output", "annotation.gtf"
    )
    minimap2_annotation_raw = os.path.join(
        main_output_dir, "minimap2_output", "annotation.gtf"
    )
    stringtie_annotation_raw = os.path.join(
        main_output_dir, "stringtie_output", "annotation.gtf"
    )
    scallop_annotation_raw = os.path.join(
        main_output_dir, "scallop_output", "annotation.gtf"
    )
    busco_annotation_raw = os.path.join(main_output_dir, "busco_output", "annotation.gtf")

    transcript_selector_script = os.path.join(
        main_script_dir, "support_scripts_perl", "select_best_transcripts.pl"
    )
    finalise_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "finalise_geneset.pl"
    )
    clean_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_geneset.pl"
    )
    clean_utrs_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_utrs_and_lncRNAs.pl"
    )
    gtf_to_seq_script = os.path.join(
        main_script_dir, "support_scripts_perl", "gtf_to_seq.pl"
    )

    transcriptomic_annotation_raw = os.path.join(
        final_annotation_dir, "transcriptomic_raw.gtf"
    )
    file_out = open(transcriptomic_annotation_raw, "w+")
    for transcriptomic_file in [
        minimap2_annotation_raw,
        scallop_annotation_raw,
        stringtie_annotation_raw,
    ]:

        if not os.path.exists(transcriptomic_file):
            logging.info(
                "No annotation.gtf file found in " + transcriptomic_file + ", skipping"
            )
            continue

        file_in = open(transcriptomic_file)
        line = file_in.readline()
        while line:
            print(line.rstrip(), file=file_out)
            line = file_in.readline()
        file_in.close()
    file_out.close()

    # Copy the raw files into the annotation dir, this is not needed as such, but collecting them in one place and relabelling is
    # helpful for a user
    if os.path.exists(busco_annotation_raw):
        subprocess.run(
            [
                "cp",
                busco_annotation_raw,
                os.path.join(final_annotation_dir, "busco_raw.gtf"),
            ]
        )

    if os.path.exists(protein_annotation_raw):
        subprocess.run(
            [
                "cp",
                protein_annotation_raw,
                os.path.join(final_annotation_dir, "protein_raw.gtf"),
            ]
        )

    gtf_files = ["transcriptomic_raw.gtf", "protein_raw.gtf", "busco_raw.gtf"]
    generic_select_cmd = ["perl", transcript_selector_script, "-genome_file", genome_file]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on whether the seqs are from transcriptomic data or not
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        transcriptomic_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".trans.gtf")
        )
        busco_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".busco.gtf")
        )
        protein_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".protein.gtf")
        )

        if os.path.exists(transcriptomic_annotation_raw):
            logging.info("Finalising transcriptomic data for: " + seq_region_name)
            transcriptomic_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", transcriptomic_annotation_raw
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    transcriptomic_annotation_raw,
                    "-output_gtf_file",
                    transcriptomic_region_gtf_path,
                    "-cds_search",
                    "-final_biotype",
                    "transcriptomic",
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

        if os.path.exists(busco_annotation_raw):
            logging.info("Finalising BUSCO data for: " + seq_region_name)
            busco_annotation_select = re.sub("_raw.gtf", "_sel.gtf", busco_annotation_raw)
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    busco_annotation_raw,
                    "-output_gtf_file",
                    busco_region_gtf_path,
                    "-all_cds_exons",
                    "-final_biotype",
                    "busco",
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

        if os.path.exists(protein_annotation_raw):
            logging.info("Finalising protein data for: " + seq_region_name)
            protein_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", protein_annotation_raw
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    protein_annotation_raw,
                    "-output_gtf_file",
                    protein_region_gtf_path,
                    "-clean_transcripts",
                    "-all_cds_exons",
                    "-final_biotype",
                    "protein",
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    # At this point we will have the region files for all the,
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".trans.gtf", "transcriptomic"
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".busco.gtf", "busco"
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".protein.gtf", "protein"
    )

    # Create a single GTF file with all the selected transcripts now that they have proper ids
    fully_merged_gtf_path = os.path.join(
        final_annotation_dir, "all_selected_transcripts.gtf"
    )
    fully_merged_gtf_out = open(fully_merged_gtf_path, "w+")

    merge_gtf_cmd = ["cat"]
    merge_gtf_cmd.extend(glob.glob(final_annotation_dir + "/*_sel.gtf"))
    subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out)
    fully_merged_gtf_out.close()

    # Now collapse the gene set
    generic_finalise_cmd = ["perl", finalise_geneset_script, "-genome_file", genome_file]

    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        final_region_gtf_path = os.path.join(
            final_region_annotation_dir, (region_details + ".final.gtf")
        )

        cmd = generic_finalise_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                fully_merged_gtf_path,
                "-output_gtf_file",
                final_region_gtf_path,
            ]
        )
        pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir, final_region_annotation_dir, ".final.gtf", "prevalidation"
    )
    merged_gtf_file = os.path.join(final_annotation_dir, ("prevalidation_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, ("prevalidation_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, ("prevalidation_sel.prot.fa")
    )
    updated_gtf_lines = validate_coding_transcripts(
        merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,
    )
    postvalidation_gtf_file = os.path.join(final_annotation_dir, ("postvalidation.gtf"))
    file_out = open(postvalidation_gtf_file, "w+")
    for line in updated_gtf_lines:
        file_out.write(line)
    file_out.close()

    cleaned_initial_gtf_file = os.path.join(final_annotation_dir, ("cleaned_pre_utr.gtf"))
    cleaned_utr_gtf_file = os.path.join(final_annotation_dir, ("annotation.gtf"))

    logging.info("Cleaning initial set")
    cleaning_cmd = [
        "perl",
        clean_geneset_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        postvalidation_gtf_file,
        "-output_gtf_file",
        cleaned_initial_gtf_file,
    ]
    logging.info(" ".join(cleaning_cmd))
    subprocess.run(cleaning_cmd)

    # Clean UTRs
    generic_clean_utrs_cmd = [
        "perl",
        clean_utrs_script,
        "-genome_file",
        genome_file,
        "-input_gtf_file",
        cleaned_initial_gtf_file,
    ]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        utr_region_gtf_path = os.path.join(
            utr_region_annotation_dir, (region_details + ".utr.gtf")
        )

        cmd = generic_clean_utrs_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                cleaned_initial_gtf_file,
                "-output_gtf_file",
                utr_region_gtf_path,
            ]
        )
        pool.apply_async(multiprocess_generic, args=(cmd,))
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir, utr_region_annotation_dir, ".utr.gtf", "annotation"
    )
    subprocess.run(
        [
            "mv",
            os.path.join(final_annotation_dir, "annotation_sel.gtf"),
            cleaned_utr_gtf_file,
        ]
    )

    logging.info("Dumping transcript and translation sequences")
    dumping_cmd = [
        "perl",
        gtf_to_seq_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        cleaned_utr_gtf_file,
    ]
    logging.info(" ".join(dumping_cmd))
    subprocess.run(dumping_cmd)

    logging.info("Finished creating gene set")


def validate_coding_transcripts(
    cdna_file,
    amino_acid_file,
    validation_dir,
    validation_type,
    diamond_validation_db,
    gtf_file,
    num_threads,
):

    logging.info("Running CDS validation with RNAsamba and CPC2")
    rnasamba_weights = "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"
    rnasamba_output_path = os.path.join(validation_dir, "rnasamba.tsv.txt")
    cpc2_output_path = os.path.join(validation_dir, "cpc2.tsv")
    rnasamba_volume = validation_dir + "/:/app:rw"
    rnasamba_cmd = [
        "singularity",
        "exec",
        "--bind",
        rnasamba_volume,
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif",
        "rnasamba",
        "classify",
        rnasamba_output_path,
        cdna_file,
        rnasamba_weights,
    ]
    logging.info(" ".join(rnasamba_cmd))
    subprocess.run(rnasamba_cmd)
    cpc2_volume = validation_dir + "/:/app:rw"
    cpc2_cmd = [
        "singularity",
        "exec",
        "--bind",
        cpc2_volume,
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif",
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        cdna_file,
        "--ORF",
        "-o",
        cpc2_output_path,
    ]
    logging.info(" ".join(cpc2_cmd))
    subprocess.run(cpc2_cmd)
    cpc2_output_path = cpc2_output_path + ".txt"

    check_file(rnasamba_output_path)
    check_file(cpc2_output_path)

    logging.info("diamond validation")
    diamond_results = None
    if diamond_validation_db is not None:
        diamond_output_dir = create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db, amino_acid_file, diamond_output_dir, num_threads
        )
        diamond_results = read_diamond_results(diamond_output_dir)

    logging.info("read results")
    rnasamba_results = read_rnasamba_results(rnasamba_output_path)
    cpc2_results = read_cpc2_results(cpc2_output_path)
    combined_results = combine_results(rnasamba_results, cpc2_results, diamond_results)
    logging.info("read gtf genes")
    parsed_gtf_genes = read_gtf_genes(gtf_file)
    updated_gtf_lines = update_gtf_genes(
        parsed_gtf_genes, combined_results, validation_type
    )

    return updated_gtf_lines


def diamond_validation(
    diamond_validation_db, amino_acid_file, diamond_output_dir, num_threads
):

    batched_protein_files = split_protein_file(amino_acid_file, diamond_output_dir, 100)

    pool = multiprocessing.Pool(int(num_threads))
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_diamond,
            args=(
                batched_protein_file,
                diamond_output_dir,
                diamond_validation_db,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_diamond(
    batched_protein_file,
    diamond_output_dir,
    diamond_validation_db,
):

    batch_num = os.path.splitext(batched_protein_file)[0]
    batch_dir = os.path.dirname(batched_protein_file)
    diamond_output_file = batched_protein_file + ".dmdout"
    logging.info("Running diamond on " + batched_protein_file + ":")

    diamond_cmd = [
        "diamond",
        "blastp",
        "--query",
        batched_protein_file,
        "--db",
        diamond_validation_db,
        "--out",
        diamond_output_file,
    ]

    logging.info(" ".join(diamond_cmd))
    subprocess.run(diamond_cmd)
    subprocess.run(["mv", diamond_output_file, diamond_output_dir])


def update_gtf_genes(parsed_gtf_genes, combined_results, validation_type):

    output_lines = []

    for gene_id in parsed_gtf_genes.keys():
        transcript_ids = parsed_gtf_genes[gene_id].keys()
        for transcript_id in transcript_ids:
            transcript_line = parsed_gtf_genes[gene_id][transcript_id]["transcript"]
            single_cds_exon_transcript = 0
            translation_match = re.search(
                r'; translation_coords "([^"]+)";', transcript_line
            )
            if translation_match:
                translation_coords = translation_match.group(1)
                translation_coords_list = translation_coords.split(":")
                # If the start exon coords of both exons are the same, then it's the same exon and thus a single exon cds
                if translation_coords_list[0] == translation_coords_list[3]:
                    single_cds_exon_transcript = 1

            exon_lines = parsed_gtf_genes[gene_id][transcript_id]["exons"]
            validation_results = combined_results[transcript_id]
            rnasamba_coding_probability = float(validation_results[0])
            rnasamba_coding_potential = validation_results[1]
            cpc2_coding_probability = float(validation_results[2])
            cpc2_coding_potential = validation_results[3]
            transcript_length = int(validation_results[4])
            peptide_length = int(validation_results[5])
            diamond_e_value = None
            if len(validation_results) == 7:
                diamond_e_value = validation_results[6]

            avg_coding_probability = (
                rnasamba_coding_probability + cpc2_coding_probability
            ) / 2
            max_coding_probability = max(
                rnasamba_coding_probability, cpc2_coding_probability
            )

            match = re.search(r'; biotype "([^"]+)";', transcript_line)
            biotype = match.group(1)
            if biotype == "busco" or biotype == "protein":
                transcript_line = re.sub(
                    '; biotype "' + biotype + '";',
                    '; biotype "protein_coding";',
                    transcript_line,
                )
                output_lines.append(transcript_line)
                output_lines.extend(exon_lines)
                continue

            min_single_exon_pep_length = 100
            min_multi_exon_pep_length = 75
            min_single_source_probability = 0.8
            min_single_exon_probability = 0.9

            # Note that the below looks at validating things under different levels of strictness
            # There are a few different continue statements, where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            elif single_cds_exon_transcript == 1 and validation_type == "moderate":
                if (
                    diamond_e_value is not None
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        and cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and avg_coding_probability >= min_single_exon_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_multi_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "lncRNA";',
                        transcript_line,
                    )
                    transcript_line = re.sub(
                        ' translation_coords "[^"]+";', "", transcript_line
                    )
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines


def read_rnasamba_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^sequence_name", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 3:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        coding_probability = eles[1]
        coding_potential = eles[2]
        results.append([transcript_id, coding_probability, coding_potential])
        line = file_in.readline()
    file_in.close()

    return results


def read_cpc2_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^#ID", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 9:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        transcript_length = eles[1]
        peptide_length = eles[2]
        coding_probability = eles[7]
        coding_potential = eles[8]
        results.append(
            [
                transcript_id,
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )
        line = file_in.readline()
    file_in.close()

    return results


def read_diamond_results(diamond_output_dir):

    results = []
    diamond_files = glob.glob(diamond_output_dir + "/*.dmdout")
    for file_path in diamond_files:
        file_in = open(file_path)
        line = file_in.readline()
        while line:
            line = line.rstrip()

            eles = line.split("\t")
            if not len(eles) == 12:
                line = file_in.readline()
                continue

            transcript_id = eles[0]
            e_value = eles[10]
            results.append([transcript_id, e_value])
            line = file_in.readline()
    file_in.close()

    return results


def combine_results(rnasamba_results, cpc2_results, diamond_results):

    transcript_ids = {}

    for result in rnasamba_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]

        if transcript_id not in transcript_ids:
            transcript_ids[transcript_id] = [coding_probability, coding_potential]

    for result in cpc2_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]
        transcript_length = result[3]
        peptide_length = result[4]
        transcript_ids[transcript_id].extend(
            [coding_probability, coding_potential, transcript_length, peptide_length]
        )

    if diamond_results is not None:
        for result in diamond_results:
            transcript_id = result[0]
            e_value = result[1]
            # There seems to be an issue where there are a small number of sequences that don't make it into the cpc2/rnasamba output
            # Should code in a system for this, but it would be good to understand why it happens to begin with. Seems to be the same
            # number of missing seqs in both, so maybe a shared cut-off
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].extend([e_value])

    return transcript_ids


def read_gtf_genes(gtf_file):

    gtf_genes = {}
    gtf_in = open(gtf_file)
    line = gtf_in.readline()
    while line:
        eles = line.split("\t")
        if not len(eles) == 9:
            line = gtf_in.readline()
            continue

        match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

        if not match:
            line = gtf_in.readline()
            continue

        gene_id = match.group(1)
        transcript_id = match.group(2)
        feature_type = eles[2]
        if gene_id not in gtf_genes:
            gtf_genes[gene_id] = {}
        if feature_type == "transcript":
            gtf_genes[gene_id][transcript_id] = {}
            gtf_genes[gene_id][transcript_id]["transcript"] = line
            gtf_genes[gene_id][transcript_id]["exons"] = []
        elif feature_type == "exon":
            gtf_genes[gene_id][transcript_id]["exons"].append(line)
        line = gtf_in.readline()
    gtf_in.close()

    return gtf_genes


def merge_finalise_output_files(
    final_annotation_dir, region_annotation_dir, extension, id_label
):

    gtf_files = glob.glob(region_annotation_dir + "/*" + extension)

    merged_gtf_file = os.path.join(final_annotation_dir, (id_label + "_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, (id_label + "_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, (id_label + "_sel.prot.fa")
    )

    # The below is not great, it's a bit messy because there might be some cases where there aren't
    # translations. So it's not as straightforward as reading the records across all three files
    # in parallel. The solution is to just load the seqs into memory and index them on the current
    # header, which should correspond to a transcript/gene id in the GTF. When writing the results
    # into the single merged files the ids will be updated to be unique and consistent across the
    # three file types

    gene_id_counter = 0
    transcript_id_counter = 0
    gtf_out = open(merged_gtf_file, "w+")
    cdna_out = open(merged_cdna_file, "w+")
    amino_acid_out = open(merged_amino_acid_file, "w+")
    for gtf_file in gtf_files:
        logging.info("GTF file: " + gtf_file)
        cdna_seq_index = {}
        amino_acid_seq_index = {}
        cdna_file = gtf_file + ".cdna"
        amino_acid_file = gtf_file + ".prot"
        cdna_in = open(cdna_file)
        amino_acid_in = open(amino_acid_file)
        cdna_seq_index = fasta_to_dict(cdna_in.readlines())
        amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())
        cdna_in.close()
        amino_acid_in.close()

        current_gene_id = ""
        gtf_in = open(gtf_file)
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)
            if match and eles[2] == "transcript":
                transcript_id_counter += 1

            gene_id = match.group(1)
            transcript_id = match.group(2)

            if not current_gene_id:
                gene_id_counter += 1
                current_gene_id = gene_id

            if not gene_id == current_gene_id:
                gene_id_counter += 1
                current_gene_id = gene_id

            new_gene_id = id_label + "_" + str(gene_id_counter)
            new_transcript_id = id_label + "_" + str(transcript_id_counter)
            line = re.sub(
                'gene_id "' + gene_id + '"', ('gene_id "' + new_gene_id + '"'), line
            )
            line = re.sub(
                'transcript_id "' + transcript_id + '"',
                ('transcript_id "' + new_transcript_id + '"'),
                line,
            )
            gtf_out.write(line)
            line = gtf_in.readline()

            if eles[2] == "transcript":
                new_header = ">" + new_transcript_id + "\n"
                cdna_out.write(new_header + cdna_seq_index[transcript_id])

                if transcript_id in amino_acid_seq_index:
                    amino_acid_out.write(new_header + amino_acid_seq_index[transcript_id])

    gtf_out.close()
    cdna_out.close()
    amino_acid_out.close()


def fasta_to_dict(fasta_list):

    index = {}
    it = iter(fasta_list)
    for header in it:
        match = re.search(r">(.+)\n$", header)
        header = match.group(1)
        seq = next(it)
        index[header] = seq
    return index


# def merge_gtf_files(file_paths,id_label):

#  gtf_file_path = os.path.join(output_dir,'annotation.gtf')
#  gtf_out = open(gtf_file_path,'w+')
#  for gtf_file_path in gtf_files:
#    gtf_file_name = os.path.basename(gtf_file_path)
#    match = re.search(r'\.rs(\d+)\.re(\d+)\.',gtf_file_name)
#    start_offset = int(match.group(1))
#    gtf_in = open(gtf_file_path,'r')
#    line = gtf_in.readline()
#    while line:
#      values = line.split("\t")
#      if len(values) == 9 and (values[2] in feature_types):
#        values[3] = str(int(values[3]) + (start_offset - 1))
#        values[4] = str(int(values[4]) + (start_offset - 1))
#        gtf_out.write("\t".join(values))
#        line = gtf_in.readline()
#    gtf_in.close()
#  gtf_out.close()


def multiprocess_generic(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd)


def multiprocess_finalise_geneset(cmd):

    print(" ".join(cmd))
    subprocess.run(cmd)


def get_sequence(seq_region, start, end, strand, fasta_file, output_dir):
    start = int(start)
    end = int(end)
    strand = int(strand)
    start -= 1
    bedtools_path = "bedtools"

    # This creates a tempfile and writes the bed info to it based on whatever information
    # has been passed in about the sequence. Then runs bedtools getfasta. The fasta file
    # should have a faidx. This can be created with the create_faidx static method prior
    # to fetching sequence
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as bed_temp_file:
        bed_temp_file.writelines(seq_region + "\t" + str(start) + "\t" + str(end))
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
            if strand == 1:
                sequence = line.rstrip()
            else:
                sequence = reverse_complement(line.rstrip())

    os.remove(bed_temp_file.name)
    return sequence


def reverse_complement(sequence):
    rev_matrix = str.maketrans("atgcATGC", "tacgTACG")
    return sequence.translate(rev_matrix)[::-1]


def seq_region_names(genome_file):
    region_list = []

    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">([^\s]+)", line)
        if match:
            region_name = match.group(1)
            if region_name == "MT":
                logging.info("Skipping region named MT")
                line = file_in.readline()
                continue
            else:
                region_list.append(match.group(1))
        line = file_in.readline()

    return region_list


def slice_genome(genome_file, target_dir, target_slice_size):
    # The below is sort of tested
    # Without the
    target_seq_length = 50000000
    min_seq_length = 1000
    current_header = ""
    current_seq = ""
    seq_dict = {}
    for line in seq:
        match = re.search(r">(.+)$", line)
        if match and current_header:
            seq_dict[current_header] = current_seq
            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

    seq_dict[current_header] = current_seq

    seq_buffer = 0
    file_number = 0
    file_name = "genome_file_" + str(file_number)

    for header in seq_dict:
        seq_iterator = 0
        seq = seq_dict[header]

        while len(seq) > target_seq_length:
            file_out = open(os.path.join(target_dir, file_name), "w+")
            subseq = seq[0:target_seq_length]
            file_out.write(
                ">" + header + "_sli" + str(seq_iterator) + "\n" + subseq + "\n"
            )
            file_out.close()
            seq = seq[target_seq_length:]
            seq_iterator += 1
            file_number += 1
            file_name = "genome_file_" + str(file_number)

        if len(seq) >= min_seq_length:
            file_name = "genome_file_" + str(file_number)
            file_out = open(os.path.join(file_name), "w+")
            file_out.write(">" + header + "_sli" + str(seq_iterator) + "\n" + seq + "\n")
            file_out.close()
            file_number += 1
            file_name = "genome_file_" + str(file_number)


def create_paired_paths(fastq_file_paths):
    path_dict = {}
    final_list = []

    for path in fastq_file_paths:
        match = re.search(r"(.+)_\d+\.(fastq|fq)", path)
        if not match:
            logging.error(
                "Could not find _1 or _2 at the end of the prefix for file. Assuming file is not paired:"
            )
            logging.error(path)
            final_list.append([path])
            continue

        prefix = match.group(1)
        if prefix in path_dict:
            #      path_dict[prefix] = path_dict[prefix] + ',' + path
            path_dict[prefix].append(path)
        else:
            path_dict[prefix] = [path]

    for pair in path_dict:
        final_list.append(path_dict[pair])

    return final_list


def check_exe(exe_path):

    if not shutil.which(exe_path):
        raise OSError("Exe does not exist. Path checked: %s" % exe_path)


def check_file(file_path):

    if not os.path.exists(file_path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)


def coallate_results(main_output_dir):

    results_dir = create_dir(main_output_dir, "results")
    output_dirs = [
        "augustus_output",
        "cpg_output",
        "dust_output",
        "eponine_output",
        "red_output",
        "rfam_output",
        "trf_output",
        "trnascan_output",
    ]
    for output_dir in output_dirs:
        match = re.search(r"(.+)_output", output_dir)
        result_type = match.group(1)
        results_path = os.path.join(main_output_dir, output_dir, "annotation.gtf")
        copy_path = os.path.join(results_dir, (result_type + ".gtf"))
        if os.path.exists(results_path):
            cpy_cmd = ["cp", results_path, copy_path]
            subprocess.run(cpy_cmd)


def _configure_logging(work_dir):
    filename = os.path.join(work_dir, "AnnoLogFile.log")
    print("FILENAME", filename)
    # Remove all handlers associated with the root logger object.
    # for handler in logging.root.handlers[:]:
    #    logging.root.removeHandler(handler)
    logging.basicConfig(
        filename=filename,
        filemode="w",
        format="%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        level=logging.DEBUG,
    )
    logging.warning("This will get logged to a file")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        help="Path where the output and temp files will write to. Uses current dir by default",
        required=False,
    )
    parser.add_argument(
        "--genome_file", help="Path to the fasta genome file", required=True
    )
    parser.add_argument(
        "--num_threads", type=int, help="Number of threads to use", required=False
    )
    parser.add_argument(
        "--run_masking",
        help="Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome",
        required=False,
    )
    parser.add_argument(
        "--red_path",
        help="Path to Red executable. See http://toolsmith.ens.utulsa.edu",
        required=False,
    )
    parser.add_argument(
        "--genblast_path",
        help="Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html",
        required=False,
    )
    parser.add_argument(
        "--convert2blastmask_path",
        help="Path to convert2blastmask executable",
        required=False,
    )
    parser.add_argument(
        "--makeblastdb_path", help="Path to makeblastdb executable", required=False
    )
    parser.add_argument(
        "--run_genblast", help="Run GenBlast to align protein sequences", required=False
    )
    parser.add_argument(
        "--run_busco",
        help="Run GenBlast to align BUSCO protein sequences",
        required=False,
    )
    parser.add_argument(
        "--protein_file",
        help="Path to a fasta file with protein sequences",
        required=False,
    )
    parser.add_argument(
        "--busco_protein_file",
        help="Path to a fasta file with BUSCO protein sequences",
        required=False,
    )
    parser.add_argument(
        "--rfam_accessions_file",
        help="Path to a file with Rfam CM accessions, one accession per line, to use with cmsearch",
        required=False,
    )
    parser.add_argument(
        "--run_star", help="Run Star for short read alignment", required=False
    )
    parser.add_argument(
        "--star_path", help="Path to Star for short read alignment", required=False
    )
    parser.add_argument(
        "--max_reads_per_sample",
        nargs="?",
        const=0,
        type=int,
        help="The maximum number of reads to use per sample. Default=0 (unlimited)",
        required=False,
    )
    parser.add_argument(
        "--max_total_reads",
        nargs="?",
        const=0,
        type=int,
        help="The maximum total number of reads. Default=0 (unlimited)",
        required=False,
    )
    parser.add_argument(
        "--short_read_fastq_dir",
        help="Path to short read fastq dir for running with Star",
        required=False,
    )
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        const=100000,
        type=int,
        help="The maximum intron size for alignments. Default=100000",
        required=False,
    )
    parser.add_argument(
        "--run_minimap2", help="Run minimap2 for long read alignment", required=False
    )
    parser.add_argument(
        "--minimap2_path", help="Path to minimap2 for long read alignment", required=False
    )
    parser.add_argument(
        "--paftools_path",
        help="Path to paftools for SAM to BED conversion",
        required=False,
    )
    parser.add_argument(
        "--long_read_fastq_dir",
        help="Path to long read fastq dir for running with minimap2",
        required=False,
    )
    parser.add_argument(
        "--run_augustus",
        help="Run Augustus with hints for gene/transcript prediction",
        required=False,
    )
    parser.add_argument("--augustus_path", help="Path to Augustus", required=False)
    parser.add_argument(
        "--run_stringtie",
        help="Run Stringtie on the results from the STAR alignments",
        required=False,
    )
    parser.add_argument(
        "--run_scallop",
        help="Run Scallop on the results from the STAR alignments",
        required=False,
    )
    parser.add_argument("--stringtie_path", help="Path to Stringtie", required=False)
    parser.add_argument("--scallop_path", help="Path to Scallop", required=False)
    parser.add_argument(
        "--subsample_script_path", help="Path to gbiab subsampling script", required=False
    )
    parser.add_argument(
        "--samtools_path", help="Path to subsampling script", required=False
    )
    parser.add_argument(
        "--finalise_geneset",
        help="Used to finalise the gene set from the various GTF files generated",
        required=False,
    )
    parser.add_argument(
        "--db_details",
        help="A comma separated stinrg of dbname,host,port,user,pass",
        required=False,
    )
    parser.add_argument(
        "--run_cmsearch",
        help="Search for sncRNA structures using Rfam and cmsearch",
        required=False,
    )
    parser.add_argument(
        "--run_trf", help="Run TRF to find tandem repeats", required=False
    )
    parser.add_argument("--trf_path", help="Path to TRF", required=False)
    parser.add_argument(
        "--run_dust", help="Run Dust to find low complexity regions", required=False
    )
    parser.add_argument("--dust_path", help="Path to Dust", required=False)
    parser.add_argument(
        "--run_repeatmasker",
        help="Run RepeatMasker to find repeat regions",
        required=False,
    )
    parser.add_argument(
        "--repeatmasker_path", help="Path to RepeatMasker", required=False
    )
    parser.add_argument(
        "--run_trnascan", help="Run tRNAscan-SE to find tRNAs", required=False
    )
    parser.add_argument("--trnascan_path", help="Path to tRNAscan-SE", required=False)
    parser.add_argument(
        "--trnascan_filter_path",
        help="Path to tRNAscan-SE high confidence filter",
        required=False,
    )
    parser.add_argument(
        "--run_cpg", help="Run cpg_lh to find CpG islands", required=False
    )
    parser.add_argument("--cpg_path", help="Path to cpg_lh", required=False)
    parser.add_argument(
        "--run_eponine",
        help="Run Eponine to find transcription start sites",
        required=False,
    )
    parser.add_argument("--eponine_path", help="Path to Eponine jar file", required=False)
    parser.add_argument(
        "--java_path", help="Path to Java for use with Eponine", required=False
    )
    parser.add_argument(
        "--run_full_annotation",
        help="Run a full annotation, will automatically check for input data and run tools based on that",
        required=False,
    )
    parser.add_argument("--run_repeats", help="Run Red, Dust, TRF", required=False)
    parser.add_argument("--run_simple_features", help="Run CpG, Eponine", required=False)
    parser.add_argument("--run_sncrnas", help="Run Rfam, tRNAscan-SE", required=False)
    parser.add_argument(
        "--run_transcriptomic",
        help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir and/or long_read_fastq_dir are provided)",
        required=False,
    )
    parser.add_argument(
        "--run_proteins",
        help="Run GenBlast if protein_file and/or busco_protein_file",
        required=False,
    )
    parser.add_argument(
        "--diamond_validation_db",
        help="Use a Diamond db with blastp mode to help validate cds sequences",
        required=False,
    )
    parser.add_argument(
        "--validation_type",
        help='The strength of evidence needed to validate and ORF as protein coding, can be "relaxed" or "moderate"',
        required=False,
    )
    parser.add_argument(
        "--load_to_ensembl_db",
        help="Load results to an Ensembl db, must also provide the db_details flag",
        required=False,
    )
    parser.add_argument(
        "--trim_fastq", help="Trim the short read files using Trim Galore", required=False
    )
    parser.add_argument(
        "--delete_pre_trim_fastq",
        help="Delete the original fastq files after trimming",
        required=False,
    )
    parser.add_argument(
        "--repeatmasker_library", help="Specify library for repeatmasker ", required=False
    )
    parser.add_argument(
        "--repeatmasker_species",
        help="Specify species for repeatmasker (default homo) ",
        required=False,
    )
    args = parser.parse_args()

    work_dir = args.output_dir
    genome_file = args.genome_file
    num_threads = args.num_threads
    masked_genome_file = genome_file  # This will be updated later if Red is run
    run_masking = args.run_masking
    red_path = args.red_path
    genblast_path = args.genblast_path
    convert2blastmask_path = args.convert2blastmask_path
    makeblastdb_path = args.makeblastdb_path
    run_genblast = args.run_genblast
    run_busco = args.run_busco
    protein_file = args.protein_file
    busco_protein_file = args.busco_protein_file
    rfam_accessions_file = args.rfam_accessions_file
    run_star = args.run_star
    star_path = args.star_path
    short_read_fastq_dir = args.short_read_fastq_dir
    max_intron_length = args.max_intron_length
    max_reads_per_sample = args.max_reads_per_sample
    max_total_reads = args.max_total_reads
    run_minimap2 = args.run_minimap2
    minimap2_path = args.minimap2_path
    paftools_path = args.paftools_path
    long_read_fastq_dir = args.long_read_fastq_dir
    run_augustus = args.run_augustus
    augustus_path = args.augustus_path
    run_stringtie = args.run_stringtie
    run_scallop = args.run_scallop
    stringtie_path = args.stringtie_path
    scallop_path = args.scallop_path
    subsample_script_path = args.subsample_script_path
    samtools_path = args.samtools_path
    finalise_geneset = args.finalise_geneset
    db_details = args.db_details
    run_cmsearch = args.run_cmsearch
    run_trf = args.run_trf
    trf_path = args.trf_path
    run_dust = args.run_dust
    dust_path = args.dust_path
    run_trnascan = args.run_trnascan
    trnascan_path = args.trnascan_path
    trnascan_filter_path = args.trnascan_filter_path
    run_cpg = args.run_cpg
    cpg_path = args.cpg_path
    run_eponine = args.run_eponine
    eponine_path = args.eponine_path
    java_path = args.java_path
    run_repeatmasker = args.run_repeatmasker
    repeatmasker_path = args.repeatmasker_path
    run_full_annotation = args.run_full_annotation
    run_repeats = args.run_repeats
    run_simple_features = args.run_simple_features
    run_sncrnas = args.run_sncrnas
    run_transcriptomic = args.run_transcriptomic
    run_proteins = args.run_proteins
    diamond_validation_db = args.diamond_validation_db
    validation_type = args.validation_type
    load_to_ensembl_db = args.load_to_ensembl_db
    trim_fastq = args.trim_fastq
    delete_pre_trim_fastq = args.delete_pre_trim_fastq
    library = args.repeatmasker_library
    species = args.repeatmasker_species

    main_script_dir = os.path.dirname(os.path.realpath(__file__))
    # work_dir=glob.glob(work_dir)
    if not os.path.exists(genome_file):
        raise IOError("File does not exist: %s" % genome_file)

    if not work_dir:
        work_dir = os.getcwd()
        # work_dir=glob.glob(work_dir)

    _configure_logging(work_dir)
    logging.info("Work dir is: %s" % work_dir)

    if not os.path.exists(work_dir):
        logging.info("Work dir does not exist, will create")
        create_dir(work_dir, None)

    if not num_threads:
        logging.info("No thread count specified, so defaulting to 1. This might be slow")
        num_threads = 1

    # If the run_full_annotation flag is set then we want to set a standardised set of analyses
    if run_full_annotation:
        run_repeats = 1
        run_simple_features = 1
        run_sncrnas = 1
        run_transcriptomic = 1
        run_proteins = 1
        finalise_geneset = 1

    # These are subsets of the analyses that can be run, group by type
    if run_repeats:
        run_masking = 1
        run_dust = 1
        run_trf = 1

    if run_simple_features:
        run_cpg = 1
        run_eponine = 1

    if run_sncrnas:
        if rfam_accessions_file:
            run_cmsearch = 1
        run_trnascan = 1

    if run_transcriptomic:
        if short_read_fastq_dir:
            run_star = 1
            run_scallop = 1
            run_stringtie = 1
        if long_read_fastq_dir:
            run_minimap2 = 1

    if run_proteins:
        if protein_file:
            run_genblast = 1
        if busco_protein_file:
            run_busco = 1

    # Collect a list of seq region names, most useful for multiprocessing regions
    seq_region_names = seq_region_names(genome_file)
    for i in seq_region_names:
        logging.info(i)

    #################################
    # Repeat analyses
    #################################
    if run_masking:
        logging.info("Running masking via Red")
        masked_genome_file = run_red(red_path, work_dir, genome_file)
        logging.info("Masked genome file: " + masked_genome_file)

    else:
        logging.info("Not running masking, presuming the genome file is softmasked")

    if run_dust:
        logging.info("Annotating low complexity regions")
        run_dust_regions(genome_file, dust_path, work_dir, num_threads)

    if run_trf:
        logging.info("Annotating tandem repeats")
        run_trf_repeats(genome_file, trf_path, work_dir, num_threads)

    if run_repeatmasker:
        logging.info("Annotating repeats with RepeatMasker")
        run_repeatmasker_regions(
            genome_file, repeatmasker_path, library, species, work_dir, num_threads
        )

    #################################
    # Simple feature analyses
    #################################
    if run_cpg:
        logging.info("Annotating CpG islands")
        run_cpg_regions(genome_file, cpg_path, work_dir, num_threads)

    if run_eponine:
        logging.info("Running Eponine to find transcription start sites")
        run_eponine_regions(genome_file, java_path, eponine_path, work_dir, num_threads)

    #################################
    # sncRNA analyses
    #################################
    # Search Rfam with cmsearch
    if run_cmsearch:
        logging.info("Annotating sncRNAs")
        run_cmsearch_regions(
            genome_file, None, None, None, rfam_accessions_file, work_dir, num_threads
        )

    if run_trnascan:
        logging.info("Annotating tRNAs")
        run_trnascan_regions(
            genome_file, trnascan_path, trnascan_filter_path, work_dir, num_threads
        )

    #################################
    # Transcriptomic analyses
    #################################
    if trim_fastq:
        run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    # Run STAR
    if run_star:
        logging.info("Running Star")
        run_star_align(
            star_path,
            trim_fastq,
            subsample_script_path,
            work_dir,
            short_read_fastq_dir,
            genome_file,
            max_reads_per_sample,
            max_total_reads,
            max_intron_length,
            num_threads,
        )

    # Run Scallop
    if run_scallop:
        logging.info("Running Scallop")
        run_scallop_assemble(scallop_path, stringtie_path, work_dir)

    # Run Stringtie
    if run_stringtie:
        logging.info("Running Stringtie")
        run_stringtie_assemble(
            stringtie_path, samtools_path, work_dir, genome_file, num_threads
        )

    # Run minimap2
    if run_minimap2:
        logging.info("Running minimap2")
        run_minimap2_align(
            minimap2_path,
            paftools_path,
            work_dir,
            long_read_fastq_dir,
            genome_file,
            max_intron_length,
            num_threads,
        )

    if run_transcriptomic:
        check_transcriptomic_output(work_dir)

    #################################
    # Protein analyses
    #################################
    # Run GenBlast
    if run_genblast:
        logging.info("Running GenBlast")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            os.path.join(work_dir, "genblast_output"),
            protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
        )

    # Run GenBlast on BUSCO set, gives higher priority when creating the final genes in cases where transcriptomic data are missing or fragmented
    if run_busco:
        logging.info("Running GenBlast of BUSCO proteins")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            os.path.join(work_dir, "busco_output"),
            busco_protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
        )

    #################################
    # Finalisation analyses
    #################################
    # Do some magic
    if finalise_geneset:
        logging.info("Finalise geneset")
        run_finalise_geneset(
            main_script_dir,
            work_dir,
            genome_file,
            seq_region_names,
            validation_type,
            diamond_validation_db,
            num_threads,
        )

    #################################
    # Other analyses
    #################################
    # Run Augustus
    if run_augustus:
        logging.info("Running Augustus")
        run_augustus_predict(augustus_path, work_dir, masked_genome_file, num_threads)

    if load_to_ensembl_db:
        load_results_to_ensembl_db(
            main_script_dir,
            load_to_ensembl_db,
            genome_file,
            work_dir,
            db_details,
            num_threads,
        )

#  coallate_results(work_dir)

#  run_find_orfs(genome_file,work_dir)
