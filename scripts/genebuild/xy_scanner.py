# Copyright [2022] EMBL-European Bioinformatics Institute
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

import csv
import argparse
import os
import tempfile
import subprocess
import re


def run_minimap(genome_index,marker_file,output_dir):

  temp_name = next(tempfile._get_candidate_names())
  temp_path = os.path.join(output_dir,(temp_name + '.mm'))
  subprocess.run(['minimap2','--secondary=no',genome_index,marker_file,'-o',temp_path])
  chromosome_present = process_results(temp_path,marker_file)
  return chromosome_present

  
def process_results(temp_path,marker_file):

  presence_threshold = 0.75
  results_dict = {}
  marker_in = open(marker_file,'r')
  line = marker_in.readline()
  while line:
    match = re.search(r'^>(.+)$',line)
    if match:
      label = match.group(1)
      results_dict[label] = 0
    line = marker_in.readline()

  with open(temp_path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter="\t")
    for row in csv_reader:
      label = row[0]
      query_length = row[1]
      matching_bases = row[9]
      percent_matching = float(matching_bases) / float(query_length)
      if results_dict[label] < percent_matching:
        results_dict[label] = percent_matching
  csv_file.close()

  for label in results_dict:
    percent_matching = results_dict[label]
    print("Result for " + label + ": " + str(percent_matching))

  num_markers = len(results_dict.keys())
  average_matching_bases = sum(results_dict.values()) / num_markers
  average_matching_bases = "{:.2f}".format(average_matching_bases)
  print("Average matching bases: " + average_matching_bases)  

  if float(average_matching_bases) >= presence_threshold:
    return 1
  else:
    return 0


if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--genome_index', help='Path to the genome index file (mmi)', required=True)
  parser.add_argument('--x_marker_file', help='Number of threads to use', required=True)
  parser.add_argument('--y_marker_file', help='Number of threads to use', required=True)
  parser.add_argument('--output_dir', help='Path to output dir to write minimap output to',required=False)
  parser.add_argument('--output_file_name', help='File name for the result, will be written to output dir',required=False)
  args = parser.parse_args()
  genome_index = args.genome_index
  x_marker_file = args.x_marker_file
  y_marker_file = args.y_marker_file
  output_dir = args.output_dir
  output_file_name = args.output_file_name

  if output_file_name is None:
    output_file_name = 'xy_scanner.out'

  x_present = run_minimap(genome_index,x_marker_file,output_dir)
  y_present = run_minimap(genome_index,y_marker_file,output_dir)

  result_out = open(os.path.join(output_dir,output_file_name),'w+')
  if(x_present and y_present):
    print("XY present")
    result_out.write("XY")
  elif(x_present):
    print("X present")
    result_out.write("X")
  elif(y_present):
    print("Y present")
    result_out.write("Y")
  else:
    print("None present")
    result_out.write("None")
  result_out.close()
