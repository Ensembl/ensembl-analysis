# Copyright [2019] EMBL-European Bioinformatics Institute
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

import logging
import time
import sys
sys.path.append('/hps/nobackup2/production/ensembl/fergal/coding/gbiab/python_threading/')

import argparse
import os
import shutil
import subprocess
import glob
import re
import multiprocessing
import gbiab_multiprocess



def create_dir(main_output_dir,dir_name):

  target_dir = os.path.join(main_output_dir,dir_name)

  if os.path.exists(target_dir):
    print ("Directory already exists, will not create again")
    return target_dir

  print ("Attempting to create target dir: %s" % target_dir)

  try:
    os.mkdir(target_dir)

  except OSError:
    print ("Creation of the dir failed, path used: %s" % target_dir)
  else:
    print ("Successfully created the dir on the following path: %s" % target_dir)

  return target_dir


def run_red(red_path,main_output_dir,genome_file):

  if not red_path:
    red_path = 'Red'

  check_exe(red_path)
  red_dir = create_dir(main_output_dir,'red_output')
  red_mask_dir = create_dir(red_dir,'mask_output')
  red_repeat_dir = create_dir(red_dir,'repeat_output')
  red_genome_dir = create_dir(red_dir,'genome_dir')

  sym_link_genome_cmd = 'ln -s ' + genome_file + ' ' + red_genome_dir

  genome_file_name = os.path.basename(genome_file)
  red_genome_file = os.path.join(red_genome_dir,genome_file_name)

  if os.path.exists(red_genome_file):
    print ('Genome file already found on the path to the Red genome dir, will not create a sym link')

  else:
    print ('Preparing to sym link the genome file to the Red genome dir. Cmd\n%s' % sym_link_genome_cmd)
    subprocess.run(['ln','-s',genome_file,red_genome_dir])

  if not os.path.exists(os.path.join(red_genome_dir,genome_file_name)):
    print ('Could not find the genome file in the Red genome dir or sym link to the original file. Path expected:\n%s' % red_genome_file)

  print ('Running Red, this may take some time depending on the genome size')
  subprocess.run([red_path,'-gnm',red_genome_dir,'-msk',red_mask_dir,'-rpt',red_repeat_dir])

  print ('Completed running Red')

  return red_genome_file


def run_genblast(genblast_path,convert2blastmask_path,makeblastdb_path,main_output_dir,protein_file,masked_genome_file):

  if not genblast_path:
    genblast_path = 'genblast'

  check_exe(genblast_path)

  if not convert2blastmask_path:
    convert2blastmask_path = 'convert2blastmask'

  check_exe(convert2blastmask_path)

  if not makeblastdb_path:
    makeblastdb_path = 'makeblastdb'

  check_exe(makeblastdb_path)

  genblast_dir = create_dir(main_output_dir,'genblast_output')

  genblast_output_file = os.path.join(genblast_dir,'genblast')

  asnb_file = masked_genome_file + '.asnb'
  print ("ASNB file: %s" % asnb_file)

  if not os.path.exists(masked_genome_file):
    raise IOError('Masked genome file does not exist: %s' % masked_genome_file)

  if not os.path.exists(protein_file):
    raise IOError('Protein file does not exist: %s' % protein_file)

  if not os.path.exists(asnb_file):
    run_convert2blastmask(convert2blastmask_path,masked_genome_file,asnb_file)
  else:
    print ('Found an existing asnb, so will skip convert2blastmask')

  if not os.path.exists(asnb_file):
    raise IOError('asnb file does not exist: %s' % asnb_file)

  run_makeblastdb(makeblastdb_path,masked_genome_file,asnb_file)

  print ('Running genblast, this may take some time depending on the genome size, number of proteins and level of softmasking')
  subprocess.run([genblast_path,'-p','genblastg','-q',protein_file,'-t',masked_genome_file,'-o',genblast_output_file,'-g','T','-pid','-r','5','-P','blast','-gff','-e','1e-1','-c','0.5','-W','3','-softmask','-scodon','50','-i','30','-x','10','-n','30','-d','200000'])

  print ('Completed running GenBlast')


def run_convert2blastmask(convert2blastmask_path,masked_genome_file,asnb_file):

  asnb_file = masked_genome_file + '.asnb'
  print ('Running convert2blastmask prior to GenBlast')
  subprocess.run([convert2blastmask_path,'-in',masked_genome_file,'-parse_seqids','-masking_algorithm','other','-masking_options','"REpeatDetector, default"','-outfmt','maskinfo_asn1_bin','-out',asnb_file])
  print ('Completed running convert2blastmask')


def run_makeblastdb(makeblastdb_path,masked_genome_file,asnb_file):

  print ('Running makeblastdb prior to GenBlast')
  subprocess.run([makeblastdb_path,'-in',masked_genome_file,'-dbtype','nucl','-parse_seqids','-mask_data',asnb_file,'-max_file_sz','10000000000'])
  print ('Completed running makeblastdb')


def run_star_align(star_path,main_output_dir,short_read_fastq_dir,genome_file,num_threads):

  if not star_path:
    star_path = 'star'

  check_exe(star_path)

  star_dir = create_dir(main_output_dir,'star_output')
  star_tmp_dir = os.path.join(star_dir,'tmp')
  if os.path.exists(star_tmp_dir):
    subprocess.run(['rm','-rf',star_tmp_dir])

  star_index_file = os.path.join(star_dir,'SAindex')

  fastq_file_list = []
  for fastq_file in glob.glob(short_read_fastq_dir + "/*.fastq"):
    fastq_file_list.append(fastq_file)

  for fastq_file in glob.glob(short_read_fastq_dir + "/*.fq"):
    fastq_file_list.append(fastq_file)

  if not fastq_file_list:
    raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % short_read_fastq_dir) 

  if not os.path.exists(star_index_file):
    print ('Did not find an index file for Star. Will create now')
    subprocess.run([star_path,'--runThreadN',num_threads,'--runMode','genomeGenerate','--outFileNamePrefix',(star_dir + '/'),'--genomeDir',star_dir,'--genomeFastaFiles',genome_file])

  if not star_index_file:
    raise IOError('The index file does not exist. Expected path:\n%s' % star_index_file)

  print ('Running Star on the files in the fastq dir')
  for fastq_file_path in fastq_file_list:
    fastq_file_name = os.path.basename(fastq_file_path)
    print ("Processing %s" % fastq_file)
    subprocess.run([star_path,'--runThreadN',num_threads,'--twopassMode','Basic','--runMode','alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file_path,'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir])
    subprocess.run(['mv',os.path.join(star_dir,'Aligned.out.sam'),os.path.join(star_dir,(fastq_file_name + '.sam'))])
    subprocess.run(['mv',os.path.join(star_dir,'SJ.out.tab'),os.path.join(star_dir,(fastq_file_name + '.sj.tab'))])

  print ('Completed running STAR')


def run_minimap2_align(minimap2_path,main_output_dir,long_read_fastq_dir,genome_file,num_threads):

  if not minimap2_path:
    minimap2_path = 'minimap2'

  check_exe(minimap2_path)

  minimap2_dir = create_dir(main_output_dir,'minimap2_output')

  genome_file_name = os.path.basename(genome_file)
  genome_file_index = (genome_file_name + '.mmi')
  minimap2_index_file = os.path.join(minimap2_dir,genome_file_index)

  fastq_file_list = []
  for fastq_file in glob.glob(long_read_fastq_dir + "/*.fastq"):
    fastq_file_list.append(fastq_file)

  for fastq_file in glob.glob(long_read_fastq_dir + "/*.fq"):
    fastq_file_list.append(fastq_file)

  if not fastq_file_list:
    raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % long_read_fastq_dir) 

  if not os.path.exists(minimap2_index_file):
    print ('Did not find an index file for minimap2. Will create now')
    subprocess.run([minimap2_path,'-t',num_threads,'-d',os.path.join(minimap2_index_file),genome_file])

  if not minimap2_index_file:
    raise IOError('The minimap2 index file does not exist. Expected path:\n%s' % minimap2_index_file)

  print ('Running minimap2 on the files in the long read fastq dir')
  for fastq_file_path in fastq_file_list:
    fastq_file_name = os.path.basename(fastq_file_path)
    print ("Processing %s" % fastq_file)
    subprocess.run([minimap2_path,'-t',num_threads,'--cs','-N','1','-ax','splice:hq','-u','b',minimap2_index_file,fastq_file_path,'-o',os.path.join(minimap2_dir,(fastq_file_name + '.sam'))])

  print ('Completed running minimap2')


def run_augustus_predict(augustus_path,main_output_dir,genome_file,num_threads):

  min_seq_length = 1000

  if not augustus_path:
    augustus_path = 'augustus'

  check_exe(augustus_path)

  augustus_dir = create_dir(main_output_dir,'augustus_output')
  augustus_genome_dir = create_dir(augustus_dir,'genome_dir')
  augustus_evidence_dir = create_dir(augustus_dir,'evidence')
  augustus_hints_file = os.path.join(augustus_evidence_dir,'augustus_hints.gff')
  star_dir = os.path.join(main_output_dir,'star_output')
  minimap2_dir = os.path.join(main_output_dir,'minimap2_output')
  
  if(os.path.exists(star_dir)):
    print("Found a Star output dir, generating hints from any .sj.tab files")
    splice_junction_to_gff(star_dir,augustus_hints_file)

  print("Splitting the genome into separate files for Augustus. Will ingore sequences of less than",min_seq_length,"in length")
#  split_genome(genome_file,augustus_genome_dir,min_seq_length)

  generic_augustus_cmd = [augustus_path,'--species=human',('--hintsfile=' + augustus_hints_file),'--UTR=on','--alternatives-from-evidence=true','--allow_hinted_splicesites=atac',('--extrinsicCfgFile=' + '/homes/thibaut/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg')]
  pool = multiprocessing.Pool(int(num_threads))
  tasks = []
  for seq_file in glob.glob(augustus_genome_dir + "/*.split.fa"):
    augustus_forward = generic_augustus_cmd.copy()
    augustus_forward.append('--strand=forward')
    augustus_forward.append(seq_file)
#    augustus_forward.append('>')
#    augustus_forward.append('> ' + seq_file + '.forward.gff')
#    print(augustus_forward)
    augustus_backward = generic_augustus_cmd.copy()
    augustus_backward.append('--strand=backward')
    augustus_backward.append(seq_file)
#    augustus_backward.append('>')
 #   augustus_backward.append(seq_file + '.backward.gff')

#    print("Augustus forward command:")
#    print(' '.join(augustus_forward))
#    print("Augustus backward command:")
#    print(' '.join(augustus_backward))

#    tasks.append(augustus_forward)
#    tasks.append(augustus_backward)
#  results = [pool.apply_async(test_mulit, t) for t in tasks]

#    pool.apply_async(test_mulit, args=(augustus_forward))
#    i = 10
    # NB!!!!! call should be args=(augustus_forward,), removing the ',' seems to lead to the method not getting called
    pool.apply_async(multiprocess_augustus, args=(augustus_forward,(seq_file + '.forward.gff'),))
    pool.apply_async(multiprocess_augustus, args=(augustus_backward,(seq_file + '.backward.gff'),))
#    results = pool.apply_async(test_mulit, args=(augustus_forward), callback=collect_result)
#    results = pool.apply_async(test_mulit, args=(augustus_backward), callback=collect_result)

  pool.close()
  pool.join()

#  augustus_cmd = [augustus_path,'--species=human',('--hintsfile=' + augustus_hints_file),'--UTR=on','--alternatives-from-evidence=true','--allow_hinted_splicesites=atac',('--extrinsicCfgFile=' + '/homes/thibaut/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg'),genome_file]
#  print("Running Augustus with the following command")
#  print(' '.join(augustus_cmd))
#  subprocess.run(augustus_cmd)

#results = []
#def collect_result(result):
#    global results
#    results.append(result)


def multiprocess_augustus(cmd,output_file):

#  format = "%(asctime)s: %(message)s"
#  logging.basicConfig(format=format, level=logging.INFO,
#                      datefmt="%H:%M:%S")

#  logging.info("Subprocess: starting", cmd)

#  file_out = open('/hps/nobackup2/production/ensembl/fergal/coding/gbiab/python_threading/testcheck.tct','w+')
#  file_out.write("TEST!!!!!")
#  file_out.close()

  file_out = open(output_file,'w+')
  print('Running Augustus with the following command:')
  print(' '.join(cmd))
  print('Output will be directed to:')
  print(output_file)
  subprocess.run(cmd, stdout=file_out)
  file_out.close()

#  logging.info("Subprocess: fishi", cmd)

#  return(cmd)


def splice_junction_to_gff(input_dir,augustus_hints_file):
  
  sjf_out = open(augustus_hints_file,"a+")
  
  for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
    sjf_in = open(sj_tab_file)
    sjf_lines = sjf_in.readlines()
    for line in sjf_lines:
      elements = line.split('\t')
      strand = '+'
      # If the strand is undefined then skip, Augustus expects a strand
      if elements[3] == '0':
        continue;
      elif elements[3] == '2':
        strand = '-'

      junction_length = int(elements[2]) - int(elements[1]) + 1
      if junction_length < 100:
        continue;

      if not elements[4] and elements[7] < 10:
        continue;
       
      # For the moment treat multimapping and single mapping things as a combined score
      score = float(elements[6]) + float(elements[7]);
      score = str(score)
      output_line = [elements[0],'RNASEQ','intron',elements[1],elements[2],score,strand,'.',('src=W;mul=' + score + ';')]
      sjf_out.write('\t'.join(output_line) + '\n')

  sjf_out.close()


def split_genome(genome_file,target_dir,min_seq_length):
  # This is the lazy initial way of just splitting into a dir of files based on the toplevel sequence with a min sequence length filter
  # There are a couple of obvious improvements:
  # 1) Instead of making files for all seqs, just process N seqs parallel, where N = num_threads. Then you could clean up the seq file
  #    after each seq finishes, thus avoiding potentially having thousands of file in a dir
  # 2) Split the seq into even slices and process these in parallel (which the same cleanup as in 1). For sequences smaller than the
  #    target slice size, bundle them up together into a single file. Vastly more complex, partially implemented in the splice_genome method
  #    Allows for more consistency with parallelisation (since there should be no large outliers). But require a mapping strategy for the
  #    coords and sequence names and all the hints will need to be adjusted
  current_header = ""
  current_seq = ""

  file_in = open(genome_file)
  line = file_in.readline()
  while line:
    match = re.search(r'>(.+)$',line)
    if match and current_header:
      if len(current_seq) > min_seq_length:
        file_out_name = os.path.join(target_dir,(current_header + '.split.fa'))
        if not os.path.exists(file_out_name):
          file_out = open(file_out_name,'w+')
          file_out.write(">" + current_header + "\n" + current_seq + "\n")
          file_out.close()

        else:
          print("Found an existing split file, so will not overwrite. File found:")
          print(file_out_name)
 
      current_seq = ""
      current_header = match.group(1)
    elif match:
      current_header = match.group(1)
    else:
      current_seq += line.rstrip()

    line = file_in.readline()

  if len(current_seq) > min_seq_length:
    file_out_name =os.path.join(target_dir,(current_header + '.split.fa'))
    if not os.path.exists(file_out_name):
      file_out = open(file_out_name,'w+')
      file_out.write(">" + current_header + "\n" + current_seq + "\n")
      file_out.close()

    else:
      print("Found an existing split file, so will not overwrite. File found:")
      print(file_out_name)

  file_in.close()


def slice_genome(genome_file,target_dir,target_slice_size):
  # The below is sort of tested
  # Without the 
  target_seq_length = 50000000
  min_seq_length = 1000
  current_header = ""
  current_seq = ""
  seq_dict = {}
  for line in seq:
    match = re.search(r'>(.+)$',line)
    if match and current_header:
      seq_dict[current_header] = current_seq
      current_seq = ""
      current_header = match.group(1)
    elif match:
      current_header = match.group(1)
    else:
      current_seq += line.rstrip()

  seq_dict[current_header] = current_seq

  seq_buffer = 0
  file_number = 0
  file_name = 'genome_file_' + str(file_number)

  for header in seq_dict:
    seq_iterator = 0
    seq = seq_dict[header]

    while len(seq) > target_seq_length:
      file_out = open(os.path.join(target_dir,file_name),"w+")
      subseq = seq[0:target_seq_length]
      file_out.write(">" + header + "_sli" + str(seq_iterator) + "\n" + subseq + "\n")
      file_out.close()
      seq = seq[target_seq_length:]
      seq_iterator += 1
      file_number += 1
      file_name = 'genome_file_' + str(file_number)

    if len(seq) >= min_seq_length:
      file_name = 'genome_file_' + str(file_number)
      file_out = open(os.path.join(file_name),"w+")
      file_out.write(">" + header + "_sli" + str(seq_iterator) + "\n" + seq + "\n")
      file_out.close()
      file_number += 1
      file_name = 'genome_file_' + str(file_number)
  
  

def check_exe(exe_path):

  if not shutil.which(exe_path):
    raise OSError('Exe does not exist. Path checked: %s' % exe_path)

  

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--output_dir', help='Path where the output and temp files will write to. Uses current dir by default', required=False)
  parser.add_argument('--genome_file', help='Path to the fasta genome file', required=True)
  parser.add_argument('--num_threads', help='Number of threads to use', required=False)
  parser.add_argument('--run_masking', help='Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome', required=False)
  parser.add_argument('--red_path', help='Path to Red executable. See http://toolsmith.ens.utulsa.edu', required=False)
  parser.add_argument('--genblast_path', help='Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html', required=False)
  parser.add_argument('--convert2blastmask_path', help='Path to convert2blastmask executable', required=False)
  parser.add_argument('--makeblastdb_path', help='Path to makeblastdb executable', required=False)
  parser.add_argument('--protein_file', help='Path to a fasta file with protein sequences', required=False)
  parser.add_argument('--run_star', help='Run Star for short read alignment', required=False)
  parser.add_argument('--star_path', help='Path to Star for short read alignment', required=False)
  parser.add_argument('--short_read_fastq_dir', help='Path to short read fastq dir for running with Star', required=False)
  parser.add_argument('--run_minimap2', help='Run minimap2 for long read alignment', required=False)
  parser.add_argument('--minimap2_path', help='Path to minimap2 for long read alignment', required=False)
  parser.add_argument('--long_read_fastq_dir', help='Path to long read fastq dir for running with minimap2', required=False)
  parser.add_argument('--run_augustus', help='Run Augustus with hints for gene/transcript prediction', required=False)
  parser.add_argument('--augustus_path', help='Path to Augustus', required=False)

  args = parser.parse_args()

  work_dir = args.output_dir
  genome_file = args.genome_file
  num_threads = args.num_threads
  masked_genome_file = genome_file # This will be updated later if Red is run
  run_masking = args.run_masking
  red_path = args.red_path
  genblast_path = args.genblast_path
  convert2blastmask_path = args.convert2blastmask_path
  makeblastdb_path = args.makeblastdb_path
  protein_file = args.protein_file
  run_star = args.run_star
  star_path = args.star_path
  short_read_fastq_dir = args.short_read_fastq_dir
  run_minimap2 = args.run_minimap2
  minimap2_path = args.minimap2_path
  long_read_fastq_dir = args.long_read_fastq_dir
  run_augustus = args.run_augustus
  augustus_path = args.augustus_path
  
  if not os.path.exists(genome_file):
    raise IOError('File does not exist: %s' % genome_file)

  if not work_dir:
    work_dir = os.getcwd()
  
  print ('Work dir is: %s' % work_dir)

  if not os.path.exists(work_dir):
    print ("Work dir does not exist, will create")
    create_dir(work_dir)

  if not num_threads:
    print ("No thread count specified, so defaulting to 1. This might be slow")
    num_threads = 1

  # Run masking
  if run_masking:
    print ("Running masking via Red")
    masked_genome_file = run_red(red_path,work_dir,genome_file)

  else:
    print ("Not running masking a presuming the genome file is softmasked")

  # Run GenBlast
  if protein_file:
    print ("Running GenBlast")
    run_genblast(genblast_path,convert2blastmask_path,makeblastdb_path,work_dir,protein_file,masked_genome_file)

  # Run STAR
  if run_star:
     print ("Running Star")
     run_star_align(star_path,work_dir,short_read_fastq_dir,genome_file,num_threads)

  # Run minimap2
  if run_minimap2:
     print ("Running minimap2")
     run_minimap2_align(minimap2_path,work_dir,long_read_fastq_dir,genome_file,num_threads)

  # Run Augustus
  if run_augustus:
     print ("Running Augustus")
     run_augustus_predict(augustus_path,work_dir,genome_file,num_threads)

  # Do some magic

