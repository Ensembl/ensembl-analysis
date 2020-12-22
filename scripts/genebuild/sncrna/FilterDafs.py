# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
import numpy as np
import sys, os

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


model_path = sys.argv[1]
scaler_path = sys.argv[2]
working_dir = sys.argv[3]

dafs_path = sys.argv[4]
mfe_path = sys.argv[5]
putative_stem_loops = sys.argv[6]
output_file = sys.argv[7]

scaler = joblib.load(scaler_path)
model = joblib.load(model_path)

dafs = pd.read_csv(dafs_path, sep = "\t", index_col=3, names = list("ABCDEFGHIJKL"), low_memory = False)
dafs.sort_index(inplace=True)

mfe = pd.read_csv(mfe_path, sep="\t", index_col=3, names=list("ABCDEFG"), low_memory = False)

mfe.sort_index(inplace=True)
temp_df = pd.merge(dafs, pd.DataFrame(mfe['E']), left_index=True, right_index=True, how="inner")
temp_df.rename(columns={'E_x':'old_score', 'G':'mir_id', 'H':'blast_evalue', 'I':'pid', 'J':'cigar', 'K':'repeat_coverage', 
                   'L':'gc_perc', 'E_y':'mfe'}, inplace=True)

temp_df['size'] = temp_df.apply(lambda x: x['C'] - x['B'], axis = 1)
temp_df['norm_mfe'] = temp_df.apply(lambda x: float(x['mfe']) / x['size'], axis = 1)

df = temp_df[['size', 'blast_evalue', 'pid', 'repeat_coverage', 'gc_perc', 'mfe', 'norm_mfe']]
scaled = scaler.transform(df)
preds = model.predict(scaled)

df.insert(len(df.columns), 'predictions', preds)

removed_dafs = df.query('predictions == 0')

removed_dafs.to_csv(working_dir + "/excluded_dafs.tsv", sep="\t")
df.to_csv(working_dir + "/labelled_dafs.tsv", sep="\t")

########################### db ids
putative_stem_loops_df = pd.read_csv(putative_stem_loops, sep="\t", index_col = 3, names=list("ABCDEFG"), low_memory = False)
removed_dafs.insert(len(removed_dafs.columns), 'coords', removed_dafs.index.astype(str))
putative_stem_loops_df['coords'] = putative_stem_loops_df.index.astype(str)

merged = pd.merge(putative_stem_loops_df, removed_dafs, left_on = 'coords', right_on = 'coords', how="inner")
merged.drop_duplicates(subset='coords',inplace=True)
merged['G'].to_csv(output_file, header=False, index=False)


