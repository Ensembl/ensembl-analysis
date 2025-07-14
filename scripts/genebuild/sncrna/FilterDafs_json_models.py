# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
import json

class StandardScalerFromJSON:
    """Version-independent StandardScaler from JSON parameters"""
    
    def __init__(self, params_file):
        with open(params_file, 'r') as f:
            params = json.load(f)
        
        self.mean_ = np.array(params['mean_'])
        self.scale_ = np.array(params['scale_'])
        self.with_mean = params['with_mean']
        self.with_std = params['with_std']
    
    def transform(self, X):
        X = np.array(X, dtype=np.float64)
        if self.with_mean:
            X = X - self.mean_
        if self.with_std:
            X = X / self.scale_
        return X

class DecisionTreeFromJSON:
    """Version-independent DecisionTree from JSON parameters"""
    
    def __init__(self, tree_params):
        self.children_left = np.array(tree_params['children_left'])
        self.children_right = np.array(tree_params['children_right'])
        self.feature = np.array(tree_params['feature'])
        self.threshold = np.array(tree_params['threshold'])
        self.value = np.array(tree_params['value'])
        self.n_classes_ = tree_params['n_classes_']
    
    def predict(self, X):
        n_samples = X.shape[0]
        predictions = np.zeros(n_samples, dtype=np.int64)
        
        for i in range(n_samples):
            node = 0
            while self.children_left[node] != self.children_right[node]:
                if X[i, self.feature[node]] <= self.threshold[node]:
                    node = self.children_left[node]
                else:
                    node = self.children_right[node]
            
            leaf_values = self.value[node][0]
            predictions[i] = np.argmax(leaf_values)
        
        return predictions

class RandomForestFromJSON:
    """Version-independent RandomForest from JSON parameters"""
    
    def __init__(self, params_file):
        with open(params_file, 'r') as f:
            params = json.load(f)
        
        self.n_estimators = params['n_estimators']
        self.classes_ = np.array(params['classes_'])
        self.n_classes_ = params['n_classes_']
        
        self.trees = []
        for tree_params in params['trees']:
            self.trees.append(DecisionTreeFromJSON(tree_params))
    
    def predict(self, X):
        X = np.array(X, dtype=np.float64)
        n_samples = X.shape[0]
        
        tree_predictions = np.zeros((n_samples, self.n_estimators), dtype=np.int64)
        for i, tree in enumerate(self.trees):
            tree_predictions[:, i] = tree.predict(X)
        
        predictions = np.zeros(n_samples, dtype=np.int64)
        for i in range(n_samples):
            unique, counts = np.unique(tree_predictions[i], return_counts=True)
            predictions[i] = unique[np.argmax(counts)]
        
        return predictions

# Main script starts here
if __name__ == "__main__":
    model_path = sys.argv[1]
    scaler_path = sys.argv[2]
    working_dir = sys.argv[3]
    dafs_path = sys.argv[4]
    mfe_path = sys.argv[5]
    putative_stem_loops = sys.argv[6]
    output_file = sys.argv[7]
    
    # Load models from JSON
    scaler = StandardScalerFromJSON(scaler_path.replace('.pkl', '_params.json'))
    model = RandomForestFromJSON(model_path.replace('.pkl', '_params.json'))
    
    # Load and process data
    dafs = pd.read_csv(dafs_path, sep="\t", index_col=3, names=list("ABCDEFGHIJKL"), low_memory=False)
    dafs.sort_index(inplace=True)
    mfe = pd.read_csv(mfe_path, sep="\t", index_col=3, names=list("ABCDEFG"), low_memory=False)
    mfe.sort_index(inplace=True)
    
    temp_df = pd.merge(dafs, pd.DataFrame(mfe['E']), left_index=True, right_index=True, how="inner")
    temp_df.rename(columns={'E_x':'old_score', 'G':'mir_id', 'H':'blast_evalue', 'I':'pid', 'J':'cigar', 'K':'repeat_coverage', 
                       'L':'gc_perc', 'E_y':'mfe'}, inplace=True)
    temp_df['size'] = temp_df.apply(lambda x: x['C'] - x['B'], axis=1)
    temp_df['norm_mfe'] = temp_df.apply(lambda x: float(x['mfe']) / x['size'], axis=1)
    
    df = temp_df[['size', 'blast_evalue', 'pid', 'repeat_coverage', 'gc_perc', 'mfe', 'norm_mfe']]
    
    # Apply models
    scaled = scaler.transform(df)
    preds = model.predict(scaled)
    
    df.insert(len(df.columns), 'predictions', preds)
    removed_dafs = df.query('predictions == 0')
    removed_dafs.to_csv(working_dir + "/excluded_dafs.tsv", sep="\t")
    df.to_csv(working_dir + "/labelled_dafs.tsv", sep="\t")
    
    # Database IDs processing
    putative_stem_loops_df = pd.read_csv(putative_stem_loops, sep="\t", index_col=3, names=list("ABCDEFG"), low_memory=False)
    removed_dafs.insert(len(removed_dafs.columns), 'coords', removed_dafs.index.astype(str))
    putative_stem_loops_df['coords'] = putative_stem_loops_df.index.astype(str)
    merged = pd.merge(putative_stem_loops_df, removed_dafs, left_on='coords', right_on='coords', how="inner")
    merged.drop_duplicates(subset='coords', inplace=True)
    merged['G'].to_csv(output_file, header=False, index=False)
