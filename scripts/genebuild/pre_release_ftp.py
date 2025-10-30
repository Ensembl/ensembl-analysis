# Copyright [2025] EMBL-European Bioinformatics Institute
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


#!/usr/bin/env python3
"""
Pre-release FTP processing script for Ensembl genebuild pipeline.

This script gathers GTF, GFF3, and reheadered toplevel FASTA files from a 
pre-release directory structure and organizes them into an FTP-ready directory
structure with standardized naming and MD5 checksums.
"""

import argparse
import os
import sys
import glob
import logging
import shutil
import hashlib
import re
from typing import Dict, List, Optional, Tuple


def setup_logging() -> logging.Logger:
    """
    Set up logging configuration with timestamps and appropriate formatting.

    Returns:
        logging.Logger: Configured logger instance
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def parse_gca_id(gca_string: str) -> Tuple[str, str]:
    """
    Parse GCA ID to extract accession number and version information.

    Args:
        gca_string: GCA identifier in format GCA_030222105.1

    Returns:
        Tuple[str, str]: GCA number and version (e.g., ('030222105', '1'))

    Raises:
        ValueError: If GCA format is invalid
    """
    # Extract GCA number and version from format like GCA_030222105.1
    match = re.match(r"GCA_(\d+)\.(\d+)", gca_string)
    if not match:
        raise ValueError(
            f"Invalid GCA format: {gca_string}. Expected format: GCA_NNNNNNNN.N"
        )

    gca_number = match.group(1)
    version = match.group(2)

    return gca_number, version


def find_reheadered_fasta(output_path: str) -> Tuple[str, str]:
    """
    Find the reheadered toplevel FASTA file and its corresponding FAI index file.

    Args:
        output_path: Directory path to search for FASTA files

    Returns:
        tuple[str, str]: Full paths to (FASTA file, FAI file)

    Raises:
        FileNotFoundError: If reheadered FASTA file or its FAI index is not found
    """
    # Find FASTA file
    fasta_pattern = os.path.join(output_path, "*_softmasked_toplevel.fa.gz")
    fasta_matches = glob.glob(fasta_pattern)
    if not fasta_matches:
        raise FileNotFoundError(
            f"No reheadered toplevel FASTA file found in {output_path}"
        )
    
    return fasta_matches[0]


def find_2bit(output_path: str) -> Optional[str]:
    """
    Find the 2bit file in the output directory.

    Args:
        output_path: Path to the base output directory
    Returns:
        Optional[str]: Path to the 2bit file or None if not found
    """
    two_bit_pattern = os.path.join(output_path, "*_softmasked_toplevel.2bit")
    two_bit_matches = glob.glob(two_bit_pattern)
    if two_bit_matches:
        return two_bit_matches[0]  # Return first match
    return None

def find_annotation_files(output_path: str) -> Dict[str, Optional[str]]:
    """
    Find GFF3 and GTF annotation files in the clade (vertabrates/ etc) directory structure.
    Prioritizes main annotation files over specialized variants like abinitio.

    Args:
        output_path: Path to the base output directory

    Returns:
        Dict[str, Optional[str]]: Dictionary with 'gff3' and 'gtf' keys mapping
                                 to file paths or None if not found
    """
    annotation_files: Dict[str, Optional[str]] = {"gff3": None, "gtf": None}

    # Look for GFF3 files
    gff3_pattern = os.path.join(output_path, '*', "gff3", "**", "*.gff3.gz")
    gff3_matches = glob.glob(gff3_pattern, recursive=True)
    if gff3_matches:
        # Prioritize main GFF3 files over specialized variants
        main_gff3 = _select_main_annotation_file(gff3_matches, "gff3")
        annotation_files["gff3"] = main_gff3

    # Look for GTF files
    gtf_pattern = os.path.join(output_path, '*', "gtf", "**", "*.gtf.gz")
    gtf_matches = glob.glob(gtf_pattern, recursive=True)
    if gtf_matches:
        # Prioritize main GTF files over specialized variants
        main_gtf = _select_main_annotation_file(gtf_matches, "gtf")
        annotation_files["gtf"] = main_gtf

    return annotation_files


def _select_main_annotation_file(file_list: List[str], file_type: str) -> str:
    """
    Select the main annotation file from a list of candidates.
    
    Prioritizes files in this order:
    1. Files ending with .{file_type}.gz (main annotation)
    2. Files ending with .chr.{file_type}.gz (chromosome-level)
    3. Any other files (avoiding abinitio, which are typically minimal)
    
    Args:
        file_list: List of file paths to choose from
        file_type: Type of file ('gtf' or 'gff3')
    
    Returns:
        str: Path to the selected main annotation file
    """
    if not file_list:
        raise ValueError("Empty file list provided")
    
    # Sort files by priority
    main_files = []
    chr_files = []
    other_files = []
    abinitio_files = []
    
    for filepath in file_list:
        filename = os.path.basename(filepath)
        
        if 'abinitio' in filename.lower():
            abinitio_files.append(filepath)
        elif filename.endswith(f'.{file_type}.gz'):
            main_files.append(filepath)
        elif filename.endswith(f'.chr.{file_type}.gz'):
            chr_files.append(filepath)
        else:
            other_files.append(filepath)
    
    # Return in order of preference
    if main_files:
        return main_files[0]
    elif chr_files:
        return chr_files[0]
    elif other_files:
        return other_files[0]
    else:
        # Only abinitio files found - return first one but this might not be ideal
        return abinitio_files[0]

def is_compressed_file(filepath: str) -> bool:
    """
    Check if a file is actually compressed by examining its magic bytes.
    
    Checks for gzip magic bytes (0x1f, 0x8b) at the beginning of the file.
    
    Args:
        filepath: Path to the file to check
        
    Returns:
        bool: True if file is actually gzip compressed
    """
    try:
        with open(filepath, 'rb') as f:
            magic_bytes = f.read(2)
            # Check for gzip magic bytes (0x1f, 0x8b)
            return len(magic_bytes) == 2 and magic_bytes[0] == 0x1f and magic_bytes[1] == 0x8b
    except (IOError, OSError):
        # If we can't read the file, assume it's not compressed
        return False    

def extract_species_from_filename(filename: str) -> str:
    """
    Extract species name from Ensembl filename format.

    Parses filenames like:
    - Bos_taurus_gca002263795v4.ARS-UCD2.0.114.primary_assembly.1.gff3.gz
    - macaca_mulatta_gca030222105v1.gtf.gz

    Args:
        filename: Full filename or path to parse

    Returns:
        str: Species name (e.g., 'Bos_taurus', 'macaca_mulatta')
    """
    basename = os.path.basename(filename)
    # Split on first occurrence of version pattern (gca followed by numbers)
    match = re.match(r"([a-zA-Z_]+)_gca\d+", basename, re.IGNORECASE)
    if match:
        return match.group(1)

    # Fallback: take everything before the first dot
    return basename.split(".")[0]


def generate_new_filename(
    original_file: str, species_name: str, gca_number: str, version: str, file_type: str
) -> str:
    """
    Generate standardized filename for FTP release.

    Creates filenames in the format: species_gca{number}v{version}.{type}.gz

    Args:
        original_file: Original file path (used for reference)
        species_name: Species name (e.g., 'Macaca_mulatta')
        gca_number: GCA accession number (e.g., '030222105')
        version: Version number (e.g., '1')
        file_type: File type (e.g., 'gff3', 'gtf', 'dna.toplevel.fa')

    Returns:
        str: Standardized filename
    """
    species_lower = species_name.lower()
    new_name = f"{species_lower}_gca{gca_number}v{version}.{file_type}"

    if is_compressed_file(original_file) and not new_name.endswith('.gz'):
        new_name += '.gz'
    return new_name


def calculate_md5(filepath: str) -> str:
    """
    Calculate MD5 checksum of a file.

    Args:
        filepath: Path to the file to checksum

    Returns:
        str: Hexadecimal MD5 hash

    Raises:
        IOError: If file cannot be read
    """
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def format_species_name(species_name: str) -> str:
    """
    Reformat species name into typical binomial nomenclature format.

    Converts species names to the standard scientific naming convention where
    the genus is capitalized, subsequent parts are lowercase, and all parts
    are separated by underscores. Replaces spaces, hyphens, and periods with
    underscores.

    Args:
        species_name: Species name to format (e.g., 'homo sapiens', 'canis_lupus')

    Returns:
        str: Formatted species name (e.g., 'Homo_sapiens', 'Canis_lupus')
    """
    parts = species_name.replace("_", " ").split()
    if len(parts) >= 2:
        species_title = parts[0].capitalize() + "_" + parts[1].lower()
        # Add any additional parts if they exist
        if len(parts) > 2:
            species_title += "_" + "_".join(part.lower() for part in parts[2:])
    else:
        species_title = species_name.capitalize()
    return species_title.replace(" ", "_").replace("-", "_").replace(".", "_")


def create_ftp_directory_structure(
    output_path: str,
    species_name: str,
    gca_string: str,
    annotation_files: Dict[str, Optional[str]],
    fasta_file: str,
    two_bit_file: str,
    logger: logging.Logger,
) -> Tuple[str, List[str]]:
    """
    Create FTP directory structure and copy files with standardized naming.

    Creates directory structure: {output_path}/ftp_release/{Species_Name}/{GCA_ID}/
    and copies files with standardized names.

    Args:
        output_path: Base output directory
        species_name: Species name for directory creation
        gca_string: Full GCA identifier (e.g., 'GCA_030222105.1')
        annotation_files: Dictionary of annotation file paths
        fasta_file: Path to reheadered FASTA file
        two_bit_file: Path to 2bit file
        logger: Logger instance for status messages

    Returns:
        Tuple[str, List[str]]: FTP directory path and list of copied file paths

    Raises:
        OSError: If directory creation or file copying fails
    """
    # Parse GCA information
    gca_number, version = parse_gca_id(gca_string)

    # Create FTP directory structure: species_name/GCA_ID/
    species_title = format_species_name(species_name)

    ftp_base = os.path.join(output_path, "ftp_release")
    ftp_species_dir = os.path.join(ftp_base, species_title)
    ftp_final_dir = os.path.join(ftp_species_dir, gca_string)

    # Create directories
    os.makedirs(ftp_final_dir, exist_ok=True)
    logger.info(f"Created FTP directory: {ftp_final_dir}")

    copied_files: List[str] = []

    # Copy and rename annotation files
    for file_type, source_file in annotation_files.items():
        if source_file and os.path.exists(source_file):
            new_filename = generate_new_filename(
                source_file, species_name, gca_number, version, file_type
            )
            dest_file = os.path.join(ftp_final_dir, new_filename)

            logger.info(
                f"Copying {file_type}: {os.path.basename(source_file)} -> {new_filename}"
            )
            shutil.copy2(source_file, dest_file)
            copied_files.append(dest_file)
        else:
            logger.warning(f"No {file_type} file found to copy")

    # Copy and rename FASTA file
    if fasta_file and os.path.exists(fasta_file):
        new_fasta_name = generate_new_filename(
            fasta_file, species_name, gca_number, version, "dna.softmasked.fa"
        )
        dest_fasta = os.path.join(ftp_final_dir, new_fasta_name)

        logger.info(
            f"Copying FASTA: {os.path.basename(fasta_file)} -> {new_fasta_name}"
        )
        shutil.copy2(fasta_file, dest_fasta)
        copied_files.append(dest_fasta)
    else:
        logger.warning("No reheadered FASTA file found to copy")

    # Copy 2bit file if it exists
    if two_bit_file and os.path.exists(two_bit_file):
        new_2bit_name = generate_new_filename(
            two_bit_file, species_name, gca_number, version, "dna.toplevel.2bit"
        )
        dest_2bit = os.path.join(ftp_final_dir, new_2bit_name)

        logger.info(f"Copying 2bit file: {os.path.basename(two_bit_file)} -> {new_2bit_name}")
        shutil.copy2(two_bit_file, dest_2bit)
        copied_files.append(dest_2bit)
    else:
        logger.warning("No 2bit file found to copy")
    return ftp_final_dir, copied_files


def generate_md5_checksums(
    directory: str, files: List[str], logger: logging.Logger
) -> str:
    """
    Generate MD5 checksums for files and write to CHECKSUMS file.

    Creates a CHECKSUMS file in the specified directory with MD5 hashes
    for all provided files in the format: {hash}  {filename}

    Args:
        directory: Directory where CHECKSUMS file will be created
        files: List of file paths to checksum
        logger: Logger instance for status messages

    Returns:
        str: Path to the created CHECKSUMS file

    Raises:
        IOError: If CHECKSUMS file cannot be written
    """
    md5_file = os.path.join(directory, "CHECKSUMS")

    logger.info(f"Generating MD5 checksums in: {md5_file}")

    with open(md5_file, "w") as f:
        for filepath in sorted(files):
            if os.path.exists(filepath):
                filename = os.path.basename(filepath)
                md5_hash = calculate_md5(filepath)
                f.write(f"{md5_hash}  {filename}\n")
                logger.info(f"  {filename}: {md5_hash}")

    return md5_file


def main() -> None:
    """
    Main function to process pre-release FTP data.

    Parses command line arguments, validates inputs, and orchestrates the
    file gathering and FTP structure creation process.

    Raises:
        SystemExit: On validation errors or processing failures
    """
    parser = argparse.ArgumentParser(
        description="Gather GTF, GFF3, and FASTA files into FTP directory structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -p /path/to/GCA_030222105.1/ -g GCA_030222105.1 -s Macaca_mulatta
  %(prog)s -p /path/to/data/ -g GCA_002263795.4
        """,
    )
    parser.add_argument(
        "-p",
        "--path",
        required=True,
        help="Output path containing the genome data directory structure",
    )
    parser.add_argument(
        "-g",
        "--gca",
        required=True,
        help="GCA identifier in format GCA_NNNNNNNN.N (e.g., GCA_030222105.1)",
    )
    parser.add_argument(
        "-s",
        "--species",
        help="Species name (e.g., Macaca_mulatta). If not provided, will extract from filenames",
    )

    args = parser.parse_args()
    logger = setup_logging()

    output_path: str = args.path
    gca_string: str = args.gca
    species_name: Optional[str] = args.species

    # Validate input path
    if not os.path.exists(output_path):
        logger.error(f"Output path does not exist: {output_path}")
        sys.exit(1)

    if not os.path.isdir(output_path):
        logger.error(f"Output path is not a directory: {output_path}")
        sys.exit(1)

    logger.info(f"Processing pre-release FTP data in: {output_path}")
    logger.info(f"GCA identifier: {gca_string}")

    try:
        # Find reheadered FASTA file
        try:
            fasta_file = find_reheadered_fasta(output_path)
            logger.info(f"Found reheadered FASTA file: {os.path.basename(fasta_file)}")

        except FileNotFoundError as e:
            logger.error(str(e))
            sys.exit(1)

        annotation_files = find_annotation_files(output_path)
        logger.info(f"Found annotation files: {annotation_files}")

        # Find 2bit file if it exists
        two_bit_file = find_2bit(output_path)
        if two_bit_file:
            logger.info(f"Found 2bit file: {os.path.basename(two_bit_file)}")
        else:
            logger.warning("No 2bit file found in the output directory")

        # Extract species name if not provided
        if not species_name:
            # Try to extract from GFF3 file first, then GTF
            for file_type, filepath in annotation_files.items():
                if filepath:
                    species_name = extract_species_from_filename(filepath)
                    logger.info(
                        f"Extracted species name from {file_type} file: {species_name}"
                    )
                    break

            if not species_name:
                logger.error(
                    "Could not determine species name. Please provide with -s option"
                )
                sys.exit(1)

        # Create FTP directory structure and copy files
        ftp_dir, copied_files = create_ftp_directory_structure(
            output_path, species_name, gca_string, annotation_files, fasta_file, two_bit_file, logger
        )

        # Generate MD5 checksums
        if copied_files:
            checksum_file = generate_md5_checksums(ftp_dir, copied_files, logger)
            logger.info(f"MD5 checksums written to: {checksum_file}")

        logger.info(f"FTP structure created successfully in: {ftp_dir}")
        logger.info("Pre-release FTP processing completed successfully")

    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
