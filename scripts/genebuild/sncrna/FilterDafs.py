#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

import pandas as pd
import argparse
import logging
from sklearn.externals import joblib

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s.%(msecs)03d %(name)-4s [%(levelname)-4s] %(message)s',
                    datefmt='%d-%m-%Y %H:%M:%S')


def get_input_files(working_dir):
    dafs_path = working_dir + "/annotated_dafs.tsv"
    mfe_path = working_dir + "/rna_fold_results.txt"

    dafs = pd.read_csv(dafs_path, sep="\t", index_col=3, names=list("ABCDEFGHIJKL"))
    dafs.sort_index(inplace=True)

    mfe = pd.read_csv(mfe_path, sep="\t", index_col=3, names=list("ABCDEFG"))

    temp_df = pd.merge(dafs, pd.DataFrame(mfe['E']), left_index=True, right_index=True, how="inner")
    temp_df.rename(columns={'E_x': 'old_score', 'G': 'mir_id', 'H': 'blast_evalue', 'I': 'pid', 'J': 'cigar',
                            'K': 'repeat_coverage',
                            'L': 'gc_perc', 'E_y': 'mfe'}, inplace=True)

    temp_df['size'] = temp_df.apply(lambda x: x['C'] - x['B'], axis = 1)
    temp_df['norm_mfe'] = temp_df.apply(lambda x: x['mfe'] / x['size'], axis = 1)

    df = temp_df[['size', 'blast_evalue', 'pid', 'repeat_coverage', 'gc_perc', 'mfe', 'norm_mfe']]

    logging.info("Input files fetched successfully; # of predictions: {}".format(len(df)))

    return df


def get_predictions(scaler_path, model_path, df):
    scaler = joblib.load(scaler_path)
    model = joblib.load(model_path)
    scaled = scaler.transform(df)
    predictions = model.predict(scaled)

    logging.info("Classifying stem-loops with RFC model: {}".format(model_path))

    return predictions


def filter_predictions(working_dir, df, predictions):
    df['predictions'] = predictions
    removed_dafs = df.query('predictions == 0')
    removed_dafs['coords'] = removed_dafs.index.astype(str)

    putative_stem_loops_df = pd.read_csv(working_dir + "/identified_mirnas.bed",
                                         sep="\t", index_col=3, names=list("ABCDEFG"))
    putative_stem_loops_df['coords'] = putative_stem_loops_df.index.astype(str)

    merged = pd.merge(putative_stem_loops_df, removed_dafs,
                      left_on='coords', right_on='coords', how="inner")
    merged.drop_duplicates(subset='coords', inplace=True)
    merged['G'].to_csv(working_dir + "/mirnas_to_delete.txt", header=False, index=False)
    removed_dafs.to_csv(working_dir + "/excluded_dafs.tsv", sep="\t")
    df.to_csv(working_dir + "/labelled_dafs.tsv", sep="\t")

    logging.info("Flagged {0} stem-loops for removal, see {1}/mirnas_to_delete.txt ".
                 format(len(removed_dafs), working_dir))


if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    parser.add_argument('-s', '--scaler_path', action="store", dest="scaler_path",
                        required=True, help="Scaler pickle dump generated during pre-processing")
    parser.add_argument('-m', '--model_path', action="store", dest="model_path",
                        required=True, help="RFC model pickle path")
    parser.add_argument('-w', '--working_dir', action="store", dest="output_dir",
                        required=True, help="Output directory")

    args = parser.parse_args()
    df = get_input_files(args.output_dir)
    predictions = get_predictions(args.scaler_path, args.model_path, df)
    filter_predictions(args.output_dir, df, predictions)

