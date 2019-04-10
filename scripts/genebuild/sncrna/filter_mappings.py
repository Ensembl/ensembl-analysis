import pandas as pd
import argparse
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s.%(msecs)03d %(name)-4s [%(levelname)-4s] %(message)s',
                    datefmt='%d-%m-%Y %H:%M:%S')


def collapse_mappings(alignments, working_dir):
    mappings = {}
    alignments.sort_values(3, inplace=True)

    for idx, sam in alignments.iterrows():
        target = sam[2]
        starts = mappings[target] if target in mappings else []

        overlapped = False
        sam_range = range(sam[3], sam[3] + len(sam[9]))
        for interval in starts:
            overlapped = bool(set(sam_range) & set(interval))  # check overlapped coords
            sam_range = interval if overlapped and len(interval) > len(
                sam_range) else sam_range  # update interval to choose longest mapped

        if not overlapped:
            starts.append(sam_range)

        mappings[target] = starts

    logging.info("Finished extracting mapped mature products - Number of stem-loop structures: {}"
                 .format(len(mappings)))

    # write_output
    with open(working_dir + "/mirna_coords.txt", 'w') as f:
        for target in mappings:
            target_size = abs(int(target.split(":")[2].split("-")[1]) - int(target.split(":")[2].split("-")[0]))
            for interval in mappings[target]:
                arm = 5 if (target_size - min(interval)) > min(interval) else 3
                f.write("{T}\t{S}\t{E}\t{A}\n".format(T=target, S=min(interval), E=max(interval), A=arm))


if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    parser.add_argument('-s', '--sam_alignment', action="store", dest="alignment",
                        required=True, help="Mapped mature products in SAM format")
    parser.add_argument('-w', '--working_dir', action="store", dest="output_dir",
                        required=True, help="Output directory")

    args = parser.parse_args()
    alignments = pd.read_csv(args.alignment, sep="\t", header=None)

    logging.info("Finished reading SAM file - Number of alignments: {}"
                 .format(len(alignments)))

    collapse_mappings(alignments, args.output_dir)


