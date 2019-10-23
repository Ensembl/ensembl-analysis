import argparse
import os
import re
import random
import multiprocessing
import gzip

def subsample(fastq_files,output_files,subsample_read_limit,num_threads,compressed):

  fastq_file = fastq_files[0]
  fastq_file_pair = fastq_files[1]
  output_file = output_files[0]
  output_file_pair = output_files[1]

  check_compression = re.search(r'\.gz$',fastq_file)
  if check_compression:
    print("Found a .gz extension, so assuming compression")
    compressed = 1

  # Count the file to begin with
  if compressed:
    num_lines = sum(1 for line in gzip.open(fastq_file))
  else:
    num_lines = sum(1 for line in open(fastq_file))

  range_limit = int(num_lines/4)


  if range_limit <= subsample_read_limit:
    print("Number of reads (" + str(range_limit) + ") is less than the max allowed read count (" + str(subsample_read_limit) + "), no need to subsample")
    return

  random_indices = {}

  rand_list = random.sample(range(0,range_limit - 1), subsample_read_limit)
  rand_count = 0
  for idx,item in enumerate(rand_list):
    random_indices[rand_list[idx] * 4] = 1

  # Note that because of the checks in main this will only be true if there if a paired file that exists
  if(num_threads == 2):
    print("Processing paired files in parallel")
    pool = multiprocessing.Pool(int(num_threads))
    pool.apply_async(print_subsample, args=(fastq_file,output_file,random_indices,compressed,))
    pool.apply_async(print_subsample, args=(fastq_file_pair,output_file_pair,random_indices,compressed,))
    pool.close()
    pool.join()

  else:
    print_subsample(fastq_file,output_file,random_indices,compressed)
    if fastq_file_pair:
      print_subsample(fastq_file_pair,output_file_pair,random_indices,compressed)  


def print_subsample(fastq_file,output_file,random_indices,compressed):

  line_index = 0

  if compressed:
    file_in = gzip.open(fastq_file,'rt')
  else:
    file_in = open(fastq_file)

  file_out = open(output_file,"w+")

  l1 = file_in.readline()
  l2 = file_in.readline()
  l3 = file_in.readline()
  l4 = file_in.readline()

  if(line_index in random_indices):
    file_out.write(l1)
    file_out.write(l2)
    file_out.write(l3)
    file_out.write(l4)
  line_index += 4

  while(l4):
    if(line_index in random_indices):
      file_out.write(l1)
      file_out.write(l2)
      file_out.write(l3)
      file_out.write(l4)

    l1 = file_in.readline()
    l2 = file_in.readline()
    l3 = file_in.readline()
    l4 = file_in.readline()
    line_index += 4

  file_in.close()
  file_out.close()


if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--fastq_file', help='Path to the fastq file', required=True)
  parser.add_argument('--fastq_file_pair', help='Path to the paired file if it exists', required=False)
  parser.add_argument('--output_file', help='Designated output file. Defaults to appending ".sub" to input file', required=False)
  parser.add_argument('--output_file_pair', help='Designated output file for the paired file if it exists. Defaults to appending ".sub" to input file', required=False)
  parser.add_argument('--subsample_read_limit', type=int, help='Maximum number of reads to subsample', required=False)
  parser.add_argument('--num_threads', type=int, help='Can use two threads if a paired file is provided. Default is 1', required=False)
  parser.add_argument('--compressed', type=int, help='Use if the files are compressed with gzip. Defaults on 0', required=False)

  args = parser.parse_args()

  fastq_file = args.fastq_file
  fastq_file_pair = args.fastq_file_pair
  output_file = args.output_file
  output_file_pair = args.output_file_pair
  subsample_read_limit = args.subsample_read_limit
  num_threads = args.num_threads
  compressed = args.compressed

  if not os.path.exists(fastq_file):
    raise OSError('Fastq file does not exist. Path checked: %s' % fastq_file)

  if fastq_file_pair and not os.path.exists(fastq_file_pair):
    raise OSError('Paired fastq file does not exist. Path checked: %s' % fastq_file_path)

  if not output_file:
    output_file = fastq_file + '.sub'
    print("No output file designated. Will write to:")
    print(output_file)

  if fastq_file_pair and not output_file_pair:
    output_file_pair = fastq_file_pair + '.sub'
    print("No output file for the paired file designated. Will write to:")
    print(output_file_pair)

  if not subsample_read_limit:
    subsample_read_limit = 100000000
    print("subsample_read_limit not set, defaulting to",str(subsample_read_limit))

  if not num_threads:
    num_threads = 1
  elif num_threads > 2 and fastq_file_pair:
    num_threads = 2
    print("Maximum number of usable threads is 2, so setting num_threads to 2")
  elif not fastq_file_pair and num_threads > 1:
    num_threads = 1
    print("No paired file provided therefore maximum number of usable threads is 1, so setting num_threads to 1")

  fastq_files = []
  output_files = []  
  fastq_files.append(fastq_file)
  fastq_files.append(fastq_file_pair)
  output_files.append(output_file)
  output_files.append(output_file_pair)

  subsample(fastq_files,output_files,subsample_read_limit,num_threads,compressed)
