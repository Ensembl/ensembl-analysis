#!/usr/bin/env perl

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
import numpy as np
import logging
import sys
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.externals import joblib


def pre_process(_df):
	scaler = StandardScaler()
	tdf = _df.drop('in_mirbase', inplace=False, axis = 1)
	scaled_df = pd.DataFrame(scaler.fit_transform(tdf))
	scaled_df.columns = _df.columns[:-1]
	scaled_df['in_mirbase'] = list(_df.in_mirbase)
	
	logging.info("Finished pre-processing data, splitting data, train/test split : 75/25 %")
	return scaler, scaled_df


def build_rfc_model(X, y, _cv=10, _njobs=10):
	rfc = RandomForestClassifier(n_estimators = 100)
	rfc.fit(X, y)
	
	cv_scores = cross_val_score(rfc, X, y, scoring = 'f1_weighted', cv = _cv, n_jobs = _njobs)
	logging.info("Finished building RFC model:\n{M}\n\n>>> Accuracy on train data: {S}".format(M=rfc, S=cv_scores.mean()))

	return rfc


def train_rfc_model(X, y, _cv=5, _njobs=5):
	rfc = RandomForestClassifier()
	params = {'n_estimators': [10, 50, 100, 250, 500, 1000], 'criterion':['gini', 'entropy'], 'max_features': ['sqrt','log2', None], 'min_samples_split':[2, 5, 10, 20, 50]}

	try:
		model = GridSearchCV(rfc, params, cv = _cv, scoring = 'f1_weighted', error_score = 0, n_jobs = _njobs)
		model.fit(X, y)
		logging.info("Finished tuning parameters for RFC  model:\n{M}\n Accuracy on train data: {S}".format(M = model, S = model.best_score_))
		return model.best_estimator_		
	except:
		logging.error("Error tuning parameters for RFC model, will try building with default parameters")
		build_rfc_model(X, y)


def evaluate_model(model, X_test, y_test, column_names, wd):
	predictions = model.predict(X_test)
	importances = model.feature_importances_
	indices = np.argsort(importances)[::-1]
	
	with open(wd + "/classification_model_report.log", 'w') as f:
		f.write("[ Classification report ]\n")
		f.write(str(classification_report(y_test, predictions))) 
		f.write("\n[ Confusion Matrix ]\n")
		f.write(str(confusion_matrix(y_test, predictions))) 
		f.write("\n[ Feature importance ]\n")
		for i in range(X_test.shape[1]):
			f.write("%d. %s (%f)" % (i + 1, column_names[indices[i]], importances[indices[i]]))
			f.write("\n")
			

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--dafs_path', action="store", dest="dafs_path", required=True,
						help="Path to MiRNA DAFs (BLAST alignments) in TSV format required")

	parser.add_argument('-m', '--mfe_path', action="store", dest="mfe_path", required=True,
						help="Path to MiRNA MFE (metrics from RNAfold) in TSV format required")

	parser.add_argument('-p', '--tune_params', action="store_true", dest="tune_params",
						help="Tune parameters?")

	parser.add_argument('-s', '--species_name', action="store", dest="species_name", required=True,
						help="Species name required")

	parser.add_argument('-w', '--working_dir', action="store", dest="working_dir", required=True,
						help="Output directory required")

	parser.set_defaults(tune_params=False)

	args = parser.parse_args()

	dafs = pd.read_csv(args.dafs_path, sep="\t", index_col=3, names = list('ABCDEFGHIJKLMN')); print(dafs.head()); dafs.index = dafs.index.str.replace('chr', '')
	mfe = pd.read_csv(args.mfe_path, sep="\t", index_col=3, names=list('ABCDEFG')); print(mfe.head());


	# tune_params_flag = True if int(sys.argv[3]) > 0 else False

	species_name = args.species_name

	working_dir = args.working_dir

	# ============================================================================= [ Prepare data ]

	tdf = pd.merge(dafs, pd.DataFrame(mfe['E']), left_index=True, right_index=True, how="inner")

	tdf.rename(columns={'E_x':'old_score', 'G':'mir_id', 'H':'blast_evalue', 'I':'pid', 'J':'cigar',
						'L':'repeat_coverage', 'M':'gc_perc', 'N':'gb_reported', 'E_y':'mfe'}, inplace=True)

	tdf['in_mirbase'] = tdf.apply(lambda x: 1 if x['K'] > 0 else 0 ,axis = 1)
	tdf['size'] = tdf.apply(lambda x: x['C'] - x['B'], axis = 1)
	tdf['norm_mfe'] = tdf.apply(lambda x: x['mfe'] / float(x['size']), axis = 1)

	df = tdf[['size', 'blast_evalue', 'pid', 'repeat_coverage', 'gc_perc', 'mfe', 'norm_mfe', 'in_mirbase']]

	scaler, scaled_df = pre_process(df)

	# ============================================================================ [ Build estimators ]
	X = scaled_df.iloc[:,:-1]
	y = scaled_df['in_mirbase']

	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)

	model = train_rfc_model(X_train, y_train) if args.tune_params else build_rfc_model(X_train, y_train)

	evaluate_model(model, X_test, y_test, scaled_df.columns, working_dir)

	# ============================================================================= [ Save models ]
	joblib.dump(scaler, '{0}/filter_dafs_scaler_{1}.pkl'.format(working_dir, species_name))
	joblib.dump(model, '{0}/filter_dafs_rfc_model_{1}.pkl'.format(working_dir, species_name))


	
