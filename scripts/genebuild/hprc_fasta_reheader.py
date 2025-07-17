#!/usr/bin/env python3

import argparse
import re
import hashlib
import sys
from pathlib import Path


def calculate_sequence_md5s(fasta_file):
    """Calculate MD5 hash for each sequence in a FASTA file."""
    sequence_md5s = {}
    current_name = None
    current_sequence = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Process previous sequence if it exists
                if current_name is not None and current_sequence:
                    seq_string = "".join(current_sequence)
                    sequence_md5s[current_name] = hashlib.md5(
                        seq_string.encode()
                    ).hexdigest()

                # Start new sequence
                current_name = line[1:]  # Remove '>' character
                current_sequence = []
            else:
                # Add to current sequence
                current_sequence.append(line)

        # Process the last sequence
        if current_name is not None and current_sequence:
            seq_string = "".join(current_sequence)
            sequence_md5s[current_name] = hashlib.md5(seq_string.encode()).hexdigest()

    return sequence_md5s


def load_verification_data(tsv_file):
    """Load name and MD5 pairs from TSV file."""
    verification_data = {}
    with open(tsv_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                print(
                    f"Warning: Line {line_num} in {tsv_file} doesn't have exactly 2 columns: {line}",
                    file=sys.stderr,
                )
                continue
            name, md5_hash = parts
            verification_data[name] = md5_hash.lower()
    return verification_data


def process_fasta(input_file, output_file):
    """Process FASTA file according to the Perl one-liner logic."""
    sequence_names = []

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.rstrip("\n\r")

            if line.startswith(">"):
                chromosome_match = re.search(r"chromosome ([^,\s]+)", line)
                if chromosome_match:
                    new_header = f">{chromosome_match.group(1)}"
                    sequence_names.append(chromosome_match.group(1))
                else:
                    # Extract contig pattern like ABC123.1
                    accession_match = re.match(r"^>([A-Z]+\d+\.\d+)", line)
                    if accession_match:
                        new_header = f">{accession_match.group(1)}"
                        sequence_names.append(accession_match.group(1))
                    else:
                        clean_header = re.sub(r"^>([^\s]+).*$", r">\1", line)
                        new_header = clean_header
                        name = (
                            clean_header[1:]
                            if clean_header.startswith(">")
                            else clean_header
                        )
                        sequence_names.append(name)

                outfile.write(new_header + "\n")
            else:
                outfile.write(line.upper() + "\n")

    return sequence_names


def verify_output(output_file, verification_file, sequence_names):
    """Verify the output file against the verification TSV."""
    if not verification_file:
        print("No verification file provided. Skipping verification.")
        return True

    print(f"Loading verification data from {verification_file}...")
    verification_data = load_verification_data(verification_file)

    print(f"Calculating MD5 of output file {output_file}...")
    actual_md5s = calculate_sequence_md5s(output_file)

    all_passed = True
    matches_found = 0

    for seq_name, actual_md5 in actual_md5s.items():
        if seq_name in verification_data:
            matches_found += 1
            expected_md5 = verification_data[seq_name]
            print(f"\nVerification for sequence '{seq_name}':")
            print(f"  Expected MD5: {expected_md5}")
            print(f"  Actual MD5:   {actual_md5}")

            if actual_md5.lower() == expected_md5.lower():
                print(f"  PASSED")
            else:
                print(f"  FAILED")
                all_passed = False
        else:
            print(f"\nWarning: Sequence '{seq_name}' not found in verification file")

    for verify_name in verification_data:
        if verify_name not in actual_md5s:
            print(
                f"\nWarning: Expected sequence '{verify_name}' not found in output file"
            )
            all_passed = False

    if matches_found == 0:
        print(f"\nError: No sequences matched between output and verification files")
        print(f"Sequences in output: {list(actual_md5s.keys())}")
        print(f"Sequences in verification: {list(verification_data.keys())}")
        return False

    print(f"\nSummary: {matches_found} sequences verified")
    return all_passed


def main():
    parser = argparse.ArgumentParser(
        description="Process FASTA file by cleaning headers and converting sequences to uppercase",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.fna output.fa
  %(prog)s input.fna output.fa --verify checksums.tsv

The verification TSV file should have format: name<TAB>MD5
        """,
    )

    parser.add_argument("input_file", help="Input FASTA file to process")
    parser.add_argument(
        "output",
        nargs="?",
        help="Output FASTA file (default: input_file with _reheadered suffix)",
    )
    parser.add_argument(
        "-v",
        "--verify",
        help="TSV file with name and MD5 for verification (format: name<TAB>MD5)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually processing",
    )

    args = parser.parse_args()

    if args.output:
        output_file = args.output
    else:
        input_path = Path(args.input_file)
        output_file = (
            input_path.parent / f"{input_path.stem}_reheadered{input_path.suffix}"
        )

    if not Path(args.input_file).exists():
        raise Exception(f"Error: Input file '{args.input_file}' does not exist.")

    if args.verify and not Path(args.verify).exists():
        raise Exception(f"Error: Verification file '{args.verify}' does not exist.")

    if args.dry_run:
        print(f"Would process: {args.input_file} -> {output_file}")
        if args.verify:
            print(f"Would verify against: {args.verify}")
        sys.exit(0)

    try:
        print(f"Processing {args.input_file} -> {output_file}")
        sequence_names = process_fasta(args.input_file, output_file)
        print(f"Successfully processed {len(sequence_names)} sequences")

        if args.verify:
            success = verify_output(output_file, args.verify, sequence_names)
            if not success:
                print("Verification failed!", file=sys.stderr)
                sys.exit(1)
            else:
                print("Verification passed!")

        print("Done!")

    except Exception as e:
        raise Exception(f"Error processing file: {e}") from e


if __name__ == "__main__":
    main()
