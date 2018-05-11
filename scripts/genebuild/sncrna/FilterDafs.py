import pandas as pd
import numpy as np
import sys, os

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


model_path = sys.argv[1]
scaler_path = sys.argv[2]
working_dir = sys.argv[3]

dafs_path = working_dir + "/annotated_dafs.tsv"
mfe_path = working_dir + "/rna_fold_results.txt"
putative_stem_loops = working_dir + "/identified_mirnas.bed"

scaler = joblib.load(scaler_path)
# scaler = StandardScaler()
model = joblib.load(model_path)

dafs = pd.read_csv(dafs_path, sep = "\t", index_col=3, names = list("ABCDEFGHIJKL"))
dafs.sort_index(inplace=True)

mfe = pd.read_csv(mfe_path, sep="\t", index_col=3, names=list("ABCDEFG"))

print(dafs.head())
print(mfe.head())
mfe.sort_index(inplace=True)
temp_df = pd.merge(dafs, pd.DataFrame(mfe['E']), left_index=True, right_index=True, how="inner")
print(temp_df.head())
temp_df.rename(columns={'E_x':'old_score', 'G':'mir_id', 'H':'blast_evalue', 'I':'pid', 'J':'cigar', 'K':'repeat_coverage', 
                   'L':'gc_perc', 'E_y':'mfe'}, inplace=True)

temp_df['size'] = temp_df.apply(lambda x: x['C'] - x['B'], axis = 1)
temp_df['norm_mfe'] = temp_df.apply(lambda x: x['mfe'] / x['size'], axis = 1)

df = temp_df[['size', 'blast_evalue', 'pid', 'repeat_coverage', 'gc_perc', 'mfe', 'norm_mfe']]
# scaled = scaler.fit_transform(df)
scaled = scaler.transform(df)
preds = model.predict(scaled)

df['predictions'] = preds

removed_dafs = df.query('predictions == 0')

removed_dafs.to_csv(working_dir + "/excluded_dafs.tsv", sep="\t")
df.to_csv(working_dir + "/labelled_dafs.tsv", sep="\t")

########################### db ids
putative_stem_loops_df = pd.read_csv(putative_stem_loops, sep="\t", index_col = 3, names=list("ABCDEFG"))
merged = pd.merge(putative_stem_loops_df, removed_dafs, left_index = True, right_index = True, how = "inner")

merged.drop_duplicates(inplace=True)
merged['G'].to_csv(working_dir + "mirnas_to_delete.txt", header=False, index=False)


