# Copyright [2019-2021] EMBL-European Bioinformatics Institute
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

import argparse
import os
import shutil
import subprocess
import glob
import re
import multiprocessing
import random
import tempfile
import io


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



def search_rfam(genome_file,cmsearch_path,rfam_cm_db_path,rfam_seeds_file_path,main_output_dir,clade,num_threads):

  if not cmsearch_path:
    cmsearch_path = 'cmsearch'

  if not cmsearch_path:
    cmsearch_path = 'cmsearch'

  check_exe(cmsearch_path)
  rfam_output_dir = create_dir(main_output_dir,'rfam_output')

  rfam_dbname = 'Rfam'
  rfam_user = 'rfamro'
  rfam_host = 'mysql-rfam-public.ebi.ac.uk'
  rfam_port = '4497'
#  rfam_accession_query_cmd = ["mysql -h", rfam_host,"-u",rfam_user,"-P",rfam_port,"-NB -e",rfam_dbname,"'select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE \'%" + clade + "%\' AND (f.type LIKE \'%snRNA%\' OR f.type LIKE \'%rRNA%\' OR LOWER(f.rfam_id) LIKE \'%rnase%\' OR LOWER(f.rfam_id) LIKE \'%vault%\' OR LOWER(f.rfam_id) LIKE \'%y_rna%\' OR f.rfam_id LIKE \'%Metazoa_SRP%\') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE \'%bacteria%\' AND rfam_id NOT LIKE \'%archaea%\' AND rfam_id NOT LIKE \'%microsporidia%\';'"]

  rfam_cm_db_path = '/hps/nobackup2/production/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.cm'
  rfam_seeds_file_path = '/hps/nobackup2/production/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed'
  rfam_accession_file = os.path.join(main_output_dir,'rfam_accessions.txt') #!!!!!
  rfam_selected_models_file = os.path.join(rfam_output_dir,'rfam_models.cm')
  with open(rfam_accession_file) as rfam_accessions_in:
    rfam_accessions = rfam_accessions_in.read().splitlines()

  with open(rfam_cm_db_path, 'r') as rfam_cm_in:
    rfam_data = rfam_cm_in.read()

  rfam_models = rfam_data.split("//\n")
  rfam_cm_out = open(rfam_selected_models_file, 'w+')

  for model in rfam_models:
    # The Rfam.cm file has INFERNAL and HMMR models
    # We just want the INFERNAL ones
    match = re.search(r"(RF\d+)",model)
    match_infernal = re.search(r"INFERNAL",model)
    if match: # and match_infernal:
      model_accession = match.group(1)
      if model_accession in rfam_accessions:
        rfam_cm_out.write(model + "//\n")
  rfam_cm_out.close()

  seed_descriptions = get_rfam_seed_descriptions(rfam_seeds_file_path)
  cv_models = extract_rfam_metrics(rfam_selected_models_file)

  print("Creating list of genomic slices")
  seq_region_lengths = get_seq_region_lengths(genome_file,5000)
  slice_ids = create_slice_ids(seq_region_lengths,1000000,0,5000)

  print("Running Rfam")
#  pool = multiprocessing.Pool(int(num_threads))
#  tasks = []
  generic_cmsearch_cmd = [cmsearch_path,'--rfam','--cpu','1','--nohmmonly','--cut_ga','--tblout']
  for slice_id in slice_ids:
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    print("Processing Rfam data using cmsearch against slice: " + region_name + ":" + str(start) + ":" + str(end))
    seq = get_sequence(region_name,start,end,1,genome_file,rfam_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(rfam_output_dir,region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path,'w+')
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")  
    region_fasta_out.close()

    region_tblout_file_name = slice_file_name + ".tblout"
    region_tblout_file_path = os.path.join(rfam_output_dir,region_tblout_file_name)

    region_results_file_name = slice_file_name + ".rfam.gtf"
    region_results_file_path = os.path.join(rfam_output_dir,region_results_file_name)

  
    cmsearch_cmd = generic_cmsearch_cmd.copy()
    cmsearch_cmd.append(region_tblout_file_path)
    cmsearch_cmd.append(rfam_selected_models_file)
    cmsearch_cmd.append(region_fasta_file_path)
    print(" ".join(cmsearch_cmd))
    subprocess.run(cmsearch_cmd)

#    rfam_results_out = open(region_results_file_path,'w+') !!!!!!!! remove?
#    rfam_results_out.close() !!!!!!!!! remove?

    initial_table_results = parse_rfam_tblout(region_tblout_file_path,region_name)
    unique_table_results = remove_rfam_overlap(initial_table_results)
    filtered_table_results = filter_rfam_results(unique_table_results,cv_models)
    create_rfam_gtf(filtered_table_results,cv_models,seed_descriptions,region_name,region_results_file_path)

#    pool.apply_async(multiprocess_cmsearch, args=(generic_cmsearch_cmd,slice_id,masked_genome_file,rfam_output_dir,rfam_selected_models_file,))

#  pool.close()
#  pool.join()


def multiprocess_cmsearch(generic_cmsearch_cmd,slice_id,masked_genome_file,rfam_output_dir,rfam_selected_models_file):

  region = slice_id[0]
  start = slice_id[1]
  end = slice_id[2]
  seq = get_sequence(region,start,end,1,genome_file,output_dir)

  slice_file_name = region + ".rs" + str(start) + ".re" + str(end)

  region_fasta_file_name = slice_file_name + ".fa"
  region_fasta_file_path = os.path.join(output_dir,region_fasta_file_name)
  region_fasta_out = open(region_fasta_file_path,'w+')
  region_fasta_out.write(">" + region + "\n" + seq + "\n")  
  region_fasta_out.close()

  region_tblout_file_name = slice_file_name + ".tblout"
  region_tblout_file_path = os.path.join(output_dir,region_tblout_file_name)

  region_results_file_name = slice_file_name + ".rfam.gtf"
  region_results_file_path = os.path.join(output_dir,region_results_file_name)

  rfam_results_out = open(region_results_file_path,'w+')
  
  cmsearch_cmd = cmd.copy()
  cmsearch_cmd.append(region_tblout_file_path)
  cmsearch_cmd.append(region_fasta_file_path)
  subprocess.run(cmsearch_cmd,stdout=rfam_results_out)

  rfam_results_out.close()


def get_rfam_seed_descriptions(rfam_seeds_path):

  descriptions = {}
  rfam_seeds = []

  # NOTE: for some reason the decoder breaks on the seeds file, so I have made this ignore errors
  with open(rfam_seeds_path,encoding='utf-8',errors='ignore') as rfam_seeds_in:
    rfam_seeds = rfam_seeds_in.read().splitlines()

  for seed in rfam_seeds:
    domain_match = re.search("^\#=GF AC   (.+)",seed)
    if domain_match:
      domain = domain_match.group(1)
      descriptions[domain] = {}
      continue

    description_match = re.search("^\#=GF DE   (.+)",seed)
    if description_match:
      description = description_match.group(1)
      descriptions[domain]['description'] = description
      continue

    name_match = re.search("^\#=GF ID   (.+)",seed)
    if name_match:
      name = name_match.group(1)
      descriptions[domain]['name'] = name
      continue

    type_match = re.search("^\#=GF TP   Gene; (.+)",seed)
    if type_match:
      rfam_type = type_match.group(1)
      descriptions[domain]['type'] = rfam_type
      continue

  return descriptions


def extract_rfam_metrics(rfam_selected_models_file):

  with open(rfam_selected_models_file, 'r') as rfam_cm_in:
    rfam_data = rfam_cm_in.read()

  rfam_models = rfam_data.split("//\n")
  parsed_cm_data = {}
  for model in rfam_models:
    temp = model.split("\n")
    model_name_match = re.search(r"NAME\s+(\S+)",model)
    match_infernal = re.search(r"INFERNAL",model)
    if model_name_match and match_infernal:
      model_name = model_name_match.group(1)
      parsed_cm_data[model_name] = {}
      for line in temp:
        name_match = re.search(r"^NAME\s+(\S+)",line)
        if name_match:
          parsed_cm_data[model_name]['-name'] = name_match.group(1)
          continue

        description_match = re.search(r"^DESC\s+(\S+)",line)
        if description_match:
          parsed_cm_data[model_name]['-description'] = description_match.group(1)
          continue

        length_match = re.search(r"^CLEN\s+(\d+)",line)
        if length_match:
          parsed_cm_data[model_name]['-length'] = length_match.group(1)
          continue

        max_length_match = re.search(r"^W\s+(\d+)",line)
        if max_length_match:
          parsed_cm_data[model_name]['-maxlength'] = max_length_match.group(1)
          continue

        accession_match = re.search(r"^ACC\s+(\S+)",line)
        if accession_match:
          parsed_cm_data[model_name]['-accession'] = accession_match.group(1)
          continue

        threshold_match = re.search(r"^GA\s+(\d+)",line)
        if threshold_match:
          parsed_cm_data[model_name]['-threshold'] = threshold_match.group(1)
          continue

  return parsed_cm_data


def parse_rfam_tblout(region_tblout_file_path,region_name):

  with open(region_tblout_file_path, 'r') as rfam_tbl_in:
    rfam_tbl_data = rfam_tbl_in.read()

  tbl_results = rfam_tbl_data.split("\n")

  all_parsed_results = []
  for result in tbl_results:
    parsed_tbl_data = {}
    if not re.match(r"^" + region_name,result):
      continue

    hit = result.split()
    accession = hit[3]
    target_name = hit[0]
    query_name = hit[2]
    hstart = hit[5];
    hend = hit[6];
    start = hit[7];
    end = hit[8];
    if hit[9] == "+":
      strand = 1
    else:
      strand = -1
    evalue = hit[15];
    score = hit[14];

    parsed_tbl_data['accession'] = accession
    parsed_tbl_data['start'] = start
    parsed_tbl_data['end'] = end
    parsed_tbl_data['strand'] = strand
    parsed_tbl_data['query_name'] = query_name
    parsed_tbl_data['score'] = score
    all_parsed_results.append(parsed_tbl_data)
  return all_parsed_results

# NOTE some of the code above and the code commented out here is to do with creating
# a DAF. As we don't have a python concept of this I'm leaving it out for the moment
# but the code below is a reference
#    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
#      -slice          => $self->queries,
#      -start          => $strand == 1 ? $start : $end,
#      -end            => $strand == 1 ? $end : $start,
#      -strand         => $strand,
#      -hstart         => $hstart,
#      -hend           => $hend,
#      -hstrand        => $strand,
#      -score          => $score,
#      -hseqname       => length($target_name) > 39 ? substr($target_name, 0, 39) : $target_name,,
#      -p_value  => $evalue,
#      -align_type => 'ensembl',
#      -cigar_string  => abs($hend - $hstart) . "M",
#   );


def remove_rfam_overlap(parsed_tbl_data):

  excluded_structures = {}
  chosen_structures = []
  for structure_x in parsed_tbl_data:
    chosen_structure = structure_x
    structure_x_start = int(structure_x['start'])
    structure_x_end = int(structure_x['end'])
    structure_x_score = float(structure_x['score'])
    structure_x_accession = structure_x['accession']
    structure_x_string = ":".join([str(structure_x_start),str(structure_x_end),str(structure_x_score),structure_x_accession])
    for structure_y in parsed_tbl_data:
      structure_y_start = int(structure_y['start'])
      structure_y_end = int(structure_y['end'])
      structure_y_score = float(structure_y['score'])
      structure_y_accession = structure_y['accession']
      structure_y_string = ":".join([str(structure_y_start),str(structure_y_end),str(structure_y_score),structure_y_accession])
      if structure_y_string in excluded_structures:
        continue

      if structure_x_start <= structure_y_end and structure_x_end >= structure_y_start:
        if structure_x_score < structure_y_score:
          chosen_structure = structure_y
          excluded_structures[structure_x_string] = 1
        else:
          excluded_structures[structure_y_string] = 1

    chosen_structures.append(chosen_structure)

  return chosen_structures


def filter_rfam_results(unfiltered_tbl_data,cv_models):

  filtered_results = []
  for structure in unfiltered_tbl_data:
    query = structure['query_name']
    if query in cv_models:
      threshold = cv_models[query]['-threshold']
      if query == 'LSU_rRNA_eukarya':
        threshold = 1700;
  
      elif query == 'LSU_rRNA_archaea':
        continue

      elif query == 'LSU_rRNA_bacteria':
        continue

      elif query == 'SSU_rRNA_eukarya':
        threshold = 1600;

      elif query == '5_8S_rRNA':
        threshold = 85;

      elif query == '5S_rRNA':
        threshold = 75;

      if threshold and float(structure['score']) >= float(threshold):
        filtered_results.append(structure)
  return filtered_results

# NOTE: The below are notes from the perl code (which has extra code) about possible improvements that are not implemented there
# Although not included in RefSeq filters, additional filters that consider sizes and score_to_size ratios can be applied
# in future work to further exclude FPs
#
# my $is_valid_size = $mapping_length > $min_length && $mapping_length < $max_length ? 1 : 0;
# my $score_size_ratio = $result->{'score'} / $mapping_length;


def create_rfam_gtf(filtered_results,cm_models,descriptions,region_name,region_results_file_path):

  if not filtered_results:
    return

  rfam_gtf_out = open(region_results_file_path,'w+')  
  gene_counter = 1
  for structure in filtered_results:
    query = structure['query_name']
    accession = structure['accession']
    if query in cm_models:
      model = cm_models[query]
      description = descriptions[accession]
      domain = structure['query_name']
      padding = model['-length']
      rfam_type = description['type']  
      strand = structure['strand']
      if strand == 1:
        start = structure['start']
        end = structure['end']
        strand = '+'
      else:
        start = structure['end']
        end = structure['start']
        score = structure['score']
        strand = '-'

      biotype = "misc_RNA";
      if re.match(r"^snRNA;",rfam_type):
        biotype = "snRNA"
      if re.match(r"^snRNA; snoRNA",rfam_type):
        biotype = "snoRNA"
      if re.match(r"^snRNA; snoRNA; scaRNA;",rfam_type):
        biotype = "scaRNA"
      if re.match(r"rRNA;",rfam_type):
        biotype = "rRNA"
      if re.match(r"antisense;",rfam_type):
        biotype = "antisense"
      if re.match(r"antitoxin;",rfam_type):
        biotype = "antitoxin"
      if re.match(r"ribozyme;",rfam_type):
        biotype = "ribozyme"
      if re.match(r"" + domain,rfam_type):
        biotype = domain
      if re.match(r"" + domain,rfam_type):
        biotype = "Vault_RNA"
      if re.match(r"" + domain,rfam_type):
        biotype = "Y_RNA"

      transcript_string = (region_name + "\tRfam\ttranscript\t" + str(start) + "\t" + str(end) + "\t" + strand + "\t.\t" + 'gene_id "' +
                           str(gene_counter) + '"; transcript_id "' + str(gene_counter) + '"; biotype "' + biotype + '";\n')
      exon_string = (region_name + "\tRfam\texon\t" + str(start) + "\t" + str(end) + "\t" + strand + "\t.\t" + 'gene_id "' +
                     str(gene_counter) + '"; transcript_id "' + str(gene_counter) + '"; exon_number "1"; biotype "' + biotype + '";\n')

      rfam_gtf_out.write(transcript_string)
      rfam_gtf_out.write(exon_string)
      gene_counter += 1

  rfam_gtf_out.close()      

#  return unless ($exon->start >= 1);
#  return unless ($exon->overlaps($daf));

#  my $RNAfold = Bio::EnsEMBL::Analysis::Runnable::RNAFold->new
#  (
#    -analysis  => $self->analysis,
#    -sequence  => $seq,
#  );
#  $RNAfold->run;



# blast_db_path = /hps/nobackup2/production/ensembl/genebuild/blastdb/ncrna/ncrna_2016_05/
# cm_search_path = /nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/bin/cmsearch
# rfam_cm = /hps/nobackup2/production/ensembl/genebuild/production/African_pygmy_mouse-Mus_minutoides-GCA_902729475.1/mus_minutoides/GCA_902729475.1/ncrna/Rfam.cm
# rfam_seeds = /hps/nobackup2/production/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed
#  my $options = " --rfam --cpu 4 --nohmmonly --cut_ga --tblout $tblout_path ";
#  $options =~ s/\-\-toponly//;
#  my $command = "$program $options $cm_path $filename > $results_path ";


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
  masked_genome_file = os.path.join(red_mask_dir,os.path.splitext(genome_file_name)[0] + ".msk")

  if os.path.exists(masked_genome_file):
    print ('Masked Genome file already found on the path to the Red mask output dir. Will not create a new file')
    return masked_genome_file

  if os.path.exists(red_genome_file):
    print ('Unmasked genome file already found on the path to the Red genome dir, will not create a sym link')

  else:
    print ('Preparing to sym link the genome file to the Red genome dir. Cmd\n%s' % sym_link_genome_cmd)
    subprocess.run(['ln','-s',genome_file,red_genome_dir])

  if not os.path.exists(os.path.join(red_genome_dir,genome_file_name)):
    print ('Could not find the genome file in the Red genome dir or sym link to the original file. Path expected:\n%s' % red_genome_file)

  print ('Running Red, this may take some time depending on the genome size')
  subprocess.run([red_path,'-gnm',red_genome_dir,'-msk',red_mask_dir,'-rpt',red_repeat_dir])

  print ('Completed running Red')

  return masked_genome_file


def run_genblast_align(genblast_path,convert2blastmask_path,makeblastdb_path,main_output_dir,protein_file,masked_genome_file,num_threads):

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

  if not os.path.exists('alignscore.txt'):
    subprocess.run(['cp','/homes/fergal/enscode/ensembl-analysis/scripts/genebuild/gbiab/support_files/alignscore.txt','./'])

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

  batched_protein_files = split_protein_file(protein_file,genblast_dir)

  pool = multiprocessing.Pool(int(num_threads))
  for batched_protein_file in batched_protein_files:
    pool.apply_async(multiprocess_genblast, args=(batched_protein_file,masked_genome_file,genblast_path,))

  pool.close()
  pool.join()

  print ('Completed running GenBlast')
  print ('Combining output into single GTF')
  generate_genblast_gtf(genblast_dir)


def multiprocess_genblast(batched_protein_file,masked_genome_file,genblast_path):

  batch_num = os.path.splitext(batched_protein_file)[0]
  batch_dir = os.path.dirname(batched_protein_file)
  print("Running GenBlast on " + batched_protein_file + ":")
  
  genblast_cmd = [genblast_path,'-p','genblastg','-q',batched_protein_file,'-t',masked_genome_file,'-g','T','-pid','-r','1','-P','blast','-gff','-e','1e-1','-c','0.5','-W','3','-softmask','-scodon','50','-i','30','-x','10','-n','30','-d','100000','-o',batched_protein_file]

  print(" ".join(genblast_cmd))
  subprocess.run(genblast_cmd)


def generate_genblast_gtf(genblast_dir):
  file_out_name = os.path.join(genblast_dir,"annotation.gtf")
  file_out = open(file_out_name,'w+')
  genblast_extension = '_1.1c_2.3_s1_0_16_1'
  for root, dirs, files in os.walk(genblast_dir):
    for genblast_file in files:
      genblast_file = os.path.join(root,genblast_file)
      if genblast_file.endswith(".gff"):
        gtf_string = convert_gff_to_gtf(genblast_file)
        file_out.write(gtf_string)      
      elif genblast_file.endswith(".fa.blast") or genblast_file.endswith(".fa.blast.report") or genblast_file.endswith(genblast_extension):
        subprocess.run(['rm',genblast_file])
  file_out.close()


def convert_gff_to_gtf(gff_file):
  gtf_string = ""
  file_in = open(gff_file)
  line = file_in.readline()
  while line:
    match = re.search(r"genBlastG",line)
    if match:
      results = line.split()
      if results[2] == "coding_exon":
        results[2] = "exon"
      attributes = set_attributes(results[8],results[2])
      results[8] = attributes
      converted_line = "\t".join(results)
      gtf_string += converted_line + "\n"
    line = file_in.readline()
  file_in.close()  
  return gtf_string


def set_attributes(attributes,feature_type):

  converted_attributes = ""
  split_attributes = attributes.split(";")
  if feature_type == "transcript":
    match = re.search(r"Name\=(.+)$",split_attributes[1])
    name = match.group(1)
    converted_attributes = 'gene_id "' + name + '"; transcript_id "' + name + '";'
  elif feature_type == "exon":
    match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-",attributes)
    exon_rank = match.group(1)
    name = match.group(2)
    converted_attributes = 'gene_id "' + name + '"; transcript_id "' + name + '"; exon_number "' + exon_rank + '";'

  return converted_attributes

# Example genBlast output
#1       genBlastG       transcript      131128674       131137049       252.729 -       .       ID=259447-R1-1-A1;Name=259447;PID=84.65;Coverage=94.22;Note=PID:84.65-Cover:94.22
#1       genBlastG       coding_exon     131137031       131137049       .       -       .       ID=259447-R1-1-A1-E1;Parent=259447-R1-1-A1
#1       genBlastG       coding_exon     131136260       131136333       .       -       .       ID=259447-R1-1-A1-E2;Parent=259447-R1-1-A1
#1       genBlastG       coding_exon     131128674       131130245       .       -       .       ID=259447-R1-1-A1-E3;Parent=259447-R1-1-A1
##sequence-region       1_group1        1       4534
#1       genBlastG       transcript      161503457       161503804       30.94   +       .       ID=259453-R1-1-A1;Name=259453;PID=39.46;Coverage=64.97;Note=PID:39.46-Cover:64.97
#1       genBlastG       coding_exon     161503457       161503804       .       +       .       ID=259453-R1-1-A1-E1;Parent=259453-R1-1-A1
##sequence-region       5_group1        1       4684
#5       genBlastG       transcript      69461063        69461741        86.16   +       .       ID=259454-R1-1-A1;Name=259454;PID=82.02;Coverage=91.67;Note=PID:82.02-Cover:91.67
#5       genBlastG       coding_exon     69461063        69461081        .       +       .       ID=259454-R1-1-A1-E1;Parent=259454-R1-1-A1
#5       genBlastG       coding_exon     69461131        69461741        .       +       .       ID=259454-R1-1-A1-E2;Parent=259454-R1-1-A1


def split_protein_file(protein_file,genblast_dir):
  batch_size = 20
  batched_protein_files = []

  for i in range(0,10):
    create_dir(genblast_dir,('bin_' + str(i)))

  file_in = open(protein_file)
  line = file_in.readline()
  seq_count = 0
  batch_count = 0
  current_record = ""
  initial_seq = 1
  while line:
    num_dir = random.randint(0,9)
    match = re.search(r'>(.+)$',line)
    if match and not initial_seq and seq_count % batch_size == 0:
      file_out_name = os.path.join(genblast_dir,('bin_' + str(random.randint(0,9))),(str(batch_count) + '.fa'))
      file_out = open(file_out_name,'w+')
      file_out.write(current_record)
      file_out.close()
      batch_count += 1
      seq_count += 1
      current_record = line
      batched_protein_files.append(file_out_name)
    elif match:
      current_record += line
      initial_seq = 0
      seq_count += 1
    else:
      current_record += line
    line = file_in.readline()
  file_in.close()

  if current_record:
    file_out_name = os.path.join(genblast_dir,('bin_' + str(random.randint(0,9))),(str(batch_count) + '.fa'))
    file_out = open(file_out_name,'w+')
    file_out.write(current_record)
    file_out.close()
    batched_protein_files.append(file_out_name)

  return batched_protein_files


def run_convert2blastmask(convert2blastmask_path,masked_genome_file,asnb_file):

  asnb_file = masked_genome_file + '.asnb'
  print ('Running convert2blastmask prior to GenBlast:')
  cmd = [convert2blastmask_path,'-in',masked_genome_file,'-parse_seqids','-masking_algorithm','other','-masking_options','"REpeatDetector, default"','-outfmt','maskinfo_asn1_bin','-out',asnb_file]
  print(' '.join(cmd))
  subprocess.run(cmd)
  print ('Completed running convert2blastmask')


def run_makeblastdb(makeblastdb_path,masked_genome_file,asnb_file):

  print ('Running makeblastdb prior to GenBlast')
  subprocess.run([makeblastdb_path,'-in',masked_genome_file,'-dbtype','nucl','-parse_seqids','-mask_data',asnb_file,'-max_file_sz','10000000000'])
  print ('Completed running makeblastdb')


def run_star_align(star_path,subsample_script_path,main_output_dir,short_read_fastq_dir,genome_file,max_reads_per_sample,max_total_reads,num_threads):
  # !!! Need to add in samtools path above instead of just using 'samtools' in command

  if not star_path:
    star_path = 'STAR'

  check_exe(star_path)

  if not os.path.exists(subsample_script_path):
    subsample_script_path = 'subsample_fastq.py'

  star_dir = create_dir(main_output_dir,'star_output')
  star_tmp_dir = os.path.join(star_dir,'tmp')
  if os.path.exists(star_tmp_dir):
    subprocess.run(['rm','-rf',star_tmp_dir])

  star_index_file = os.path.join(star_dir,'SAindex')

  fastq_file_list = []
  file_types = ('*.fastq','*.fq','*.fastq.gz','*.fq.gz')
  for file_type in file_types:
    fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir,file_type)))

  # This works out if the files are paired or not
  fastq_file_list = create_paired_paths(fastq_file_list)

  # Subsamples in parallel if there's a value set
  if max_reads_per_sample:
    pool = multiprocessing.Pool(int(num_threads))
    for fastq_files in fastq_file_list:
      fastq_file = fastq_files[0]
      fastq_file_pair = ''
      if(len(fastq_files) == 2):
        fastq_file_pair = fastq_files[1]

      if fastq_file_pair and os.path.exists(fastq_file + '.sub') and os.path.exists(fastq_file_pair + '.sub'):
        print("Found an existing .sub files on the fastq path for both members of the pair, will use those instead of subsampling again. Files:")
        print(fastq_file + '.sub')
        print(fastq_file_pair + '.sub')
      elif fastq_file_pair:
        pool.apply_async(run_subsample_script, args=(fastq_file,fastq_file_pair,subsample_script_path,))
      elif os.path.exists(fastq_file + '.sub'):
        print("Found an existing .sub file on the fastq path, will use that instead. File:")
        print(fastq_file + '.sub')
      else:
        pool.apply_async(run_subsample_script, args=(fastq_file,fastq_file_pair,subsample_script_path,))

    pool.close()
    pool.join()

  fastq_file_list = check_for_fastq_subsamples(fastq_file_list)

  if not fastq_file_list:
    raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % short_read_fastq_dir) 

  if not os.path.exists(star_index_file):
    print ('Did not find an index file for Star. Will create now')
    subprocess.run([star_path,'--runThreadN',str(num_threads),'--runMode','genomeGenerate','--outFileNamePrefix',(star_dir + '/'),'--genomeDir',star_dir,'--genomeFastaFiles',genome_file])

  if not star_index_file:
    raise IOError('The index file does not exist. Expected path:\n%s' % star_index_file)

  print ('Running Star on the files in the fastq dir')
  for fastq_file_path in fastq_file_list:
    print(fastq_file_path)
    fastq_file_name = os.path.basename(fastq_file_path)
    check_compression= re.search(r'.gz$',fastq_file_name)
    print ("Processing %s" % fastq_file_path)

    star_command = [star_path,'--outFilterIntronMotifs','RemoveNoncanonicalUnannotated','--outSAMstrandField','intronMotif','--runThreadN',str(num_threads),'--twopassMode','Basic','--runMode','alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file_path,'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir,'--outSAMtype','SAM','--alignIntronMax','100000','--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000','50000','50000','50000','50000','50000','100000']

    if check_compression:
      star_command.append('--readFilesCommand')
      star_command.append('gunzip')
      star_command.append('-c')

    subprocess.run(star_command)
    subprocess.run(['mv',os.path.join(star_dir,'Aligned.out.sam'),os.path.join(star_dir,(fastq_file_name + '.sam'))])
    subprocess.run(['mv',os.path.join(star_dir,'SJ.out.tab'),os.path.join(star_dir,(fastq_file_name + '.sj.tab'))])

  print ('Completed running STAR')

  print ('Sorting sam files into bams')

  # Should move the sorting below into a method that takes a dir as an argument
  sam_files = []
  for sam_file in glob.glob(star_dir + "/*.sam"):
    sam_files.append(sam_file)

  if not sam_files:
    raise IndexError('The list of sam files is empty, expected them in Star output dir. Star dir:\n%s' % star_dir)

  sorted_bam_files = []
  for sam_file in sam_files:
    sam_file_name = os.path.basename(sam_file)
    sam_temp_file_path = os.path.join(star_dir,(sam_file_name + ".tmp"))
    bam_sort_file_path = os.path.join(star_dir,re.sub('.sam','.bam',sam_file_name))

    if os.path.exists(bam_sort_file_path):
      print("Found an existing bam file, will not sort sam file again. Bam file:")
      print(bam_sort_file_path)

    else:
      print("Converting samfile into sorted bam file. Bam file:")
      print(bam_sort_file_path)
      subprocess.run(['samtools','sort','-@',str(num_threads),'-T',sam_temp_file_path,'-o',bam_sort_file_path,sam_file])



def run_subsample_script(fastq_file,fastq_file_pair,subsample_script_path):

  if fastq_file_pair:
    subprocess.run(['python3',subsample_script_path,'--fastq_file',fastq_file,'--fastq_file_pair',fastq_file_pair])
  else:
    subprocess.run(['python3',subsample_script_path,'--fastq_file',fastq_file])


def check_for_fastq_subsamples(fastq_file_list):
  # This should probably removed at some point as it is needlessly complicated
  # Would be better to just build into the previous step
  # Mainly just about making sure that if you have subsamples they're used and if you have pairs they're paired
  for idx,fastq_files in enumerate(fastq_file_list):
    fastq_file = fastq_files[0]
    subsample_file = fastq_file + ".sub"

    fastq_file_pair = ''
    subsample_file_pair = '' 
    if(len(fastq_files) == 2):
      fastq_file_pair = fastq_files[1]
      subsample_file_pair = fastq_file_pair + ".sub"
 
    # This bit will replace the list entry with a string, don't need a list after this function for each pair/file
    if os.path.exists(subsample_file):
      print("Found a subsampled file extension, will use that instead of the original file. Path:")
      print(subsample_file)
      fastq_file_list[idx] = subsample_file
    else:
      fastq_file_list[idx] = fastq_file

    # This bit just concats the paired file (or subsampled paired file) if it exists
    if os.path.exists(subsample_file_pair):
      print("Found a subsampled paired file extension, will use that instead of the original file. Path:")
      print(subsample_file_pair)
      fastq_file_list[idx] = subsample_file + ',' + subsample_file_pair
    elif fastq_file_pair:
      fastq_file_list[idx] = fastq_file + ',' + fastq_file_pair

    print("Entry at current index:")
    print(fastq_file_list[idx])

  return(fastq_file_list)


def run_minimap2_align(minimap2_path,paftools_path,main_output_dir,long_read_fastq_dir,genome_file,num_threads):

  if not minimap2_path:
    minimap2_path = 'minimap2'

  check_exe(minimap2_path)

  if not paftools_path:
    paftools_path = 'paftools.js'

  check_exe(paftools_path)

  minimap2_output_dir = create_dir(main_output_dir,'minimap2_output')

  genome_file_name = os.path.basename(genome_file)
  genome_file_index = (genome_file_name + '.mmi')
  minimap2_index_file = os.path.join(minimap2_output_dir,genome_file_index)
  minimap2_hints_file = os.path.join(minimap2_output_dir,'minimap2_hints.gff')

  fastq_file_list = []
  for fastq_file in glob.glob(long_read_fastq_dir + "/*.fastq"):
    fastq_file_list.append(fastq_file)

  for fastq_file in glob.glob(long_read_fastq_dir + "/*.fq"):
    fastq_file_list.append(fastq_file)

  if not fastq_file_list:
    raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % long_read_fastq_dir) 

  if not os.path.exists(minimap2_index_file):
    print ('Did not find an index file for minimap2. Will create now')
    subprocess.run([minimap2_path,'-t',str(num_threads),'-d',os.path.join(minimap2_index_file),genome_file])

  if not minimap2_index_file:
    raise IOError('The minimap2 index file does not exist. Expected path:\n%s' % minimap2_index_file)

  print ('Running minimap2 on the files in the long read fastq dir')
  for fastq_file_path in fastq_file_list:
    fastq_file_name = os.path.basename(fastq_file_path)
    sam_file = os.path.join(minimap2_output_dir,(fastq_file_name + '.sam'))
    bed_file = os.path.join(minimap2_output_dir,(fastq_file_name + '.bed'))
    bed_file_out = open(bed_file,'w+')
    print ("Processing %s" % fastq_file)
    subprocess.run([minimap2_path,'-t',str(num_threads),'--cs','--secondary=no','-ax','splice','-u','b',minimap2_index_file,fastq_file_path,'-o',sam_file])
    print("Creating bed file from SAM")
    subprocess.run([paftools_path,'splice2bed',sam_file],stdout=bed_file_out)
    bed_file_out.close()

  bed_to_gtf(minimap2_output_dir)

  print ('Completed running minimap2')


def bed_to_gtf(minimap2_output_dir):

  gtf_file_path = os.path.join(minimap2_output_dir,'annotation.gtf')
  gtf_out = open(gtf_file_path,"w+")
  exons_dict = {}
  gene_id = 1
  for bed_file in glob.glob(minimap2_output_dir + "/*.bed"):
    print("Converting bed to GTF:")
    print(bed_file)
    bed_in = open(bed_file)
    bed_lines = bed_in.readlines()
    for line in bed_lines:
      line = line.rstrip()
      elements = line.split('\t')
      seq_region_name = elements[0]
      offset = int(elements[1])
      hit_name = elements[3]
      strand = elements[5]
      block_sizes = elements[10].split(',')
      block_sizes = list(filter(None, block_sizes))
      block_starts = elements[11].split(',')
      block_starts = list(filter(None, block_starts))
      exons = bed_to_exons(block_sizes,block_starts,offset)
      transcript_line = [seq_region_name,'minimap','transcript',0,0,'.',strand,'.','gene_id "minimap_' + str(gene_id) +'"; ' + 'transcript_id "minimap_' + str(gene_id) + '"' ]
      transcript_start = None
      transcript_end = None
      exon_records = []       
      for i,exon_coords in enumerate(exons):
        if transcript_start is None or exon_coords[0] < transcript_start:
          transcript_start = exon_coords[0]

        if transcript_end is None or exon_coords[1] > transcript_end:
          transcript_end = exon_coords[1]
          
        exon_line = [seq_region_name,'minimap','exon',str(exon_coords[0]),str(exon_coords[1]),'.',strand,'.',
                     'gene_id "minimap_' + str(gene_id) +'"; ' + 'transcript_id "minimap_' + str(gene_id) + '"; exon_number "' + str(i+1) + '";']

        exon_records.append(exon_line)

      transcript_line[3] = str(transcript_start)
      transcript_line[4] = str(transcript_end)

      gtf_out.write("\t".join(transcript_line) + "\n")
      for exon_line in exon_records:
        gtf_out.write("\t".join(exon_line) + "\n")

      gene_id += 1

    gtf_out.close()


def bed_to_gff(input_dir,hints_file):

  gff_out = open(hints_file,"w+")
  exons_dict = {}
  for bed_file in glob.glob(input_dir + "/*.bed"):
    print("Processing file for hints:")
    print(bed_file)
    bed_in = open(bed_file)
    bed_lines = bed_in.readlines()
    for line in bed_lines:
      line = line.rstrip()
      elements = line.split('\t')
      seq_region_name = elements[0]
      offset = int(elements[1])
      hit_name = elements[3]
      strand = elements[5]
      block_sizes = elements[10].split(',')
      block_sizes = list(filter(None, block_sizes))
      block_starts = elements[11].split(',')
      block_starts = list(filter(None, block_starts))
      exons = bed_to_exons(block_sizes,block_starts,offset)
      for i,element in enumerate(exons):
        exon_coords = exons[i]
        exon_key = seq_region_name + ':' + exon_coords[0] + ':' + exon_coords[1] + ':' + strand
        if exon_key in exons_dict:
          exons_dict[exon_key][5] += 1
        else:
          gff_list = [seq_region_name,'CDNA','exon',exon_coords[0],exon_coords[1],1.0,strand,'.']
          exons_dict[exon_key] = gff_list

  for exon_key, gff_list in exons_dict.items():
    gff_list[5] = str(gff_list[5])
    gff_line = '\t'.join(gff_list) + '\tsrc=W;mul=' + gff_list[5] + ';\n'
    gff_out.write(gff_line)

  gff_out.close()

  sorted_hints_out = open((hints_file + '.srt'),'w+')
  subprocess.run(['sort','-k1,1','-k7,7','-k4,4','-k5,5',hints_file],stdout=sorted_hints_out)
  sorted_hints_out.close()

def bed_to_exons(block_sizes,block_starts,offset):
  exons = []
  for i,element in enumerate(block_sizes):
    block_start = offset + int(block_starts[i]) + 1
    block_end = block_start + int(block_sizes[i]) - 1

    if block_end < block_start:
      print('Warning: block end is less than block start, skipping exon')
      continue

    exon_coords = [str(block_start),str(block_end)]
    exons.append(exon_coords)

  return exons



# start gene g1
#1       AUGUSTUS        gene    1       33908   1       +       .       g1
#1       AUGUSTUS        transcript      1       33908   .       +       .       g1.t1
#1       AUGUSTUS        CDS     3291    3585    .       +       2       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        exon    3291    3585    .       +       .       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        CDS     11377   11510   .       +       1       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        exon    11377   11510   .       +       .       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        CDS     30726   30871   .       +       2       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        exon    30726   30871   .       +       .       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        CDS     32975   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        exon    32975   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        stop_codon      33500   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
#1       AUGUSTUS        tts     33908   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# protein sequence = [GGRGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEKKEKKEEEEEEEEEEEEEEEK
# EEGEGMEGRDIMGIRGKQKQPRKQPELSSSFPTFLLTQKTPYCPESIKLLLILRNNIQTIVFYKGPKGVINDWRKFKLESEDGDSIPPSKKEILRQMS
# SPQSRDDSKERMSRKMSIQEYELIHQDKEDESCLRKYRRQCMQDMHQKLSFGPRYGFVYELETGEQFLETIEKEQKVTTIVVNIYEDGVRGCDALNSS
# LACLAVEYPMVKFCKIKASNTGAGDRFSTDVLPTLLVYKGGELISNFISVAEQFAEEFFAVDVESFLNEYGLLPEREIHDLEQTNMEDEDIE]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 0
# CDS exons: 0/4
# CDS introns: 0/4
# 5'UTR exons and introns: 0/0
# 3'UTR exons and introns: 0/1
# hint groups fully obeyed: 0
# incompatible hint groups: 0
# end gene g1
###

def augustus_output_to_gtf(augustus_output_dir,augustus_genome_dir):

    gtf_file_path = os.path.join(augustus_output_dir,'annotation.gtf')
    gtf_out = open(gtf_file_path,'w+')
    record_count = 1
    for gff_file_path in glob.glob(augustus_genome_dir + "/*.aug"):
      gff_file_name = os.path.basename(gff_file_path)
      match = re.search(r'\.rs(\d+)\.re(\d+)\.',gff_file_name)
      start_offset = int(match.group(1))
      
      exon_number = 1
      current_exon_hints_total = 0
      current_exon_hints_match = 0
      current_intron_hints_total = 0
      current_intron_hints_match = 0
      current_record = []
      gff_in = open(gff_file_path,'r')
      line = gff_in.readline()
      while line:
        match = re.search(r'# CDS exons\: (\d+)\/(\d+)',line)
        if match:
          current_exon_hints_match = match.group(1)
          current_exon_hints_total = match.group(2)

        match = re.search(r'# CDS introns\: (\d+)\/(\d+)',line)
        if match:
          current_introns_hints_match = match.group(1)
          current_introns_hints_total = match.group(2)

        if re.search(r'# end gene',line):
          for output_line in current_record:
            gtf_out.write(output_line)
          current_record = []
          current_exon_hints_total = 0
          current_exon_hints_match = 0
          current_intron_hints_total = 0
          current_intron_hints_match = 0
          record_count += 1
          exon_number = 1

        values = line.split("\t")
        if len(values) == 9 and values[1] == "AUGUSTUS" and (values[2] == "transcript" or values[2] == "exon"):
          values[3] = str(int(values[3]) + (start_offset - 1))
          values[4] = str(int(values[4]) + (start_offset - 1))
          values[8] = 'gene_id "aug' + str(record_count) + '"; transcript_id "aug' + str(record_count) + '";'
          if values[2] == "exon":
            values[8] = values[8] + ' exon_number "' + str(exon_number) + '";'
            exon_number += 1
          values[8] = values[8] + "\n"
          current_record.append("\t".join(values))

        line = gff_in.readline()
      gff_in.close()
    gtf_out.close()



def run_augustus_predict(augustus_path,main_output_dir,masked_genome_file,num_threads):

  min_seq_length = 1000

  if not augustus_path:
    augustus_path = '/hps/nobackup2/production/ensembl/jma/src/Augustus/bin/augustus'

  bam2hints_path = '/homes/fergal/bin/bam2hints'
  bam2wig_path = '/homes/fergal/bin/bam2wig'
  wig2hints_path = '/homes/fergal/bin/wig2hints'
  check_exe(augustus_path)

  # Run bam2hints, bam2wig, wig2hints, then combine the hints into a single file
  # Multiprocess with all three steps in the MP as that would be fastest
  
  augustus_dir = create_dir(main_output_dir,'augustus_output')
  augustus_hints_dir = create_dir(augustus_dir,'hints')
  augustus_genome_dir = create_dir(augustus_dir,'genome_dir')
  augustus_evidence_dir = create_dir(augustus_dir,'evidence')
  augustus_hints_file = os.path.join(augustus_evidence_dir,'augustus_hints.gff')
  star_dir = os.path.join(main_output_dir,'star_output')
  minimap2_output_dir = os.path.join(main_output_dir,'minimap2_output')
  
  if(os.path.exists(star_dir)):
    print("Found a Star output dir, generating hints from any .sj.tab files")
    generate_hints(bam2hints_path,bam2wig_path,wig2hints_path,augustus_hints_dir,star_dir,num_threads)
    hints_out = open(augustus_hints_file,'w+')
    for gff_file in glob.glob(augustus_hints_dir + "/*.bam.hints.gff"):
      gff_in = open(gff_file,'r')
      line = gff_in.readline()
      while line:
        hints_out.write(line)
        line = gff_in.readline()
      gff_in.close()
    hints_out.close()
  
  seq_region_lengths = get_seq_region_lengths(genome_file,5000)
  slice_ids = create_slice_ids(seq_region_lengths,1000000,100000,5000)

  generic_augustus_cmd = [augustus_path,'--species=human','--UTR=on',('--extrinsicCfgFile=' + '/hps/nobackup2/production/ensembl/jma/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg')]

  pool = multiprocessing.Pool(int(num_threads))
  tasks = []

  for slice_id in slice_ids:
    pool.apply_async(multiprocess_augustus_id, args=(generic_augustus_cmd,slice_id,masked_genome_file,augustus_hints_file,augustus_genome_dir,))

  pool.close()
  pool.join()

  augustus_output_to_gtf(augustus_dir,augustus_genome_dir)
  

def create_slice_ids(seq_region_lengths,slice_size,overlap,min_length):
  if not slice_size:
    slice_size = 1000000

  if not overlap:
    overlap = 0

  if not min_length:
    min_length = 0

  slice_ids = []

  for region in seq_region_lengths:
    region_length = seq_region_lengths[region]
    if region_length < min_length:
      continue

    if region_length <= slice_size:
      slice_ids.append([region,1,region_length])
      continue

    start = 1
    end = start + slice_size - 1
    while end < region_length:
      start = start - overlap
      if start < 1:
        start = 1

      end = start + slice_size - 1
      if end > region_length:
        end = region_length

      if (end - start + 1) >= min_length:
        slice_ids.append([region,start,end])

      start = end + 1
        
  return slice_ids


def generate_hints(bam2hints_path,bam2wig_path,wig2hints_path,augustus_hints_dir,star_dir,num_threads):

  pool = multiprocessing.Pool(int(num_threads))
  for bam_file in glob.glob(star_dir + "/*.bam"):
    pool.apply_async(multiprocess_augustus_hints, args=(bam2hints_path,bam2wig_path,wig2hints_path,bam_file,augustus_hints_dir,))
  pool.close()
  pool.join()


def multiprocess_augustus_hints(bam2hints_path,bam2wig_path,wig2hints_path,bam_file,augustus_hints_dir):
  bam_file_name = os.path.basename(bam_file)
  print("Processing " + bam_file_name + " for Augustus hints")

  bam2hints_file_name = bam_file_name + '.hints.gff'
  bam2hints_file_path = os.path.join(augustus_hints_dir,bam2hints_file_name)
  bam2hints_cmd = [bam2hints_path,('--in=' + bam_file),('--out=' + bam2hints_file_path),'--maxintronlen=100000']
  print("bam2hints command:\n" + ' '.join(bam2hints_cmd))
  subprocess.run(bam2hints_cmd)

  
#  bam2wig_cmd = [bam2wig_path,'-D',augustus_hints_dir,bam_file]
#  print("bam2wig command:\n" + ' '.join(bam2wig_cmd))
#  subprocess.run(bam2wig_cmd)


  # wig2hints is odd in that it runs directly off STDIN and then just prints to STDOUT,
  # so the code below is implemented in steps as it's not good practice to use pipes and
  # redirects in a subprocess command
#  wig_file_name = re.sub('.bam','.wig',bam_file_name)
#  wig_file_path = os.path.join(augustus_hints_dir,wig_file_name)
#  wig_hints_file_name = (wig_file_name + '.hints.gff')
#  wig_hints_file_path =  os.path.join(augustus_hints_dir,wig_hints_file_name)
#  print("Writing wig file info to hints file:\n" + wig_hints_file_name)
#  wig2hints_out = open(wig_hints_file_path,'w+')
#  wigcat = subprocess.Popen(('cat',wig_file_path), stdout=subprocess.PIPE)
#  subprocess.run(wig2hints_path, stdin=wigcat.stdout, stdout=wig2hints_out)
#  wig2hints_out.close()


def multiprocess_augustus_id(cmd,slice_id,genome_file,hints_file,output_dir):

  region = slice_id[0]
  start = slice_id[1]
  end = slice_id[2]
  seq = get_sequence(region,start,end,1,genome_file,output_dir)

  region_fasta_file_name = region + ".rs" + str(start) + ".re" + str(end) + ".fa"
  region_fasta_file_path = os.path.join(output_dir,region_fasta_file_name)
  region_augustus_file_path = os.path.join(output_dir,(region_fasta_file_name + ".aug"))

  region_fasta_out = open(region_fasta_file_path,'w+')
  region_fasta_out.write(">" + region + "\n" + seq + "\n")  
  region_fasta_out.close()

  region_hints_file = create_slice_hints_file(region,start,end,hints_file,region_fasta_file_path)

  aug_out = open(region_augustus_file_path,'w+')

  augustus_forward = cmd.copy()
  augustus_forward.append(('--hintsfile=' + region_hints_file))
  augustus_forward.append('--strand=forward')
  augustus_forward.append(region_fasta_file_path)
  subprocess.run(augustus_forward,stdout=aug_out)

  augustus_backward = cmd.copy()
  augustus_backward.append(('--hintsfile=' + region_hints_file))
  augustus_backward.append('--strand=backward')
  augustus_backward.append(region_fasta_file_path)    
  subprocess.run(augustus_backward,stdout=aug_out)

  aug_out.close()


def create_slice_hints_file(region,start,end,hints_file,region_fasta_file_path):

  # Note this is trying to be memory and file efficient at the cost of speed
  # So files are only created as needed and the hints are being read line by line as written as needed
  # This comes with the downside of being slow, but it's only a very small amount of time relative
  # to how slow the step is in total. Given that this step in general eats up a low of memory, saving as much
  # as possible here is not a bad thing even if it's adding in an overhead by continuously reading the hints file

  region_hints_file_path = region_fasta_file_path + ".gff"
  hints_in = open(hints_file)
  hints_out = open(region_hints_file_path,"w+")
  hint_line = hints_in.readline()
  while hint_line:
    hint_line_values = hint_line.split("\t")
    if not len(hint_line_values) == 9:
      hint_line = hints_in.readline()
      continue

    hint_region = hint_line_values[0]
    hint_region_start = int(hint_line_values[3])
    hint_region_end = int(hint_line_values[4])    

    if hint_region == region and hint_region_start >= start and hint_region_end <= end:
      hint_line_values[3] = str(int(hint_line_values[3]) - (start - 1))
      hint_line_values[4] = str(int(hint_line_values[4]) - (start - 1))
      hints_out.write("\t".join(hint_line_values))
  
    hint_line = hints_in.readline()
  hints_in.close()
  hints_out.close()

  return region_hints_file_path


def run_cufflinks_assemble(cufflinks_path,cuffmerge_path,samtools_path,main_output_dir,genome_file,num_threads):

  max_cufflinks_threads = 6
  if num_threads > max_cufflinks_threads:
    print("Reducing threads to " + str(max_cufflinks_threads) + " for cufflinks to keep the memory footprint down")
    num_threads = max_cufflinks_threads

  if not cufflinks_path:
    cufflinks_path = shutil.which('cufflinks')
  check_exe(cufflinks_path)

  if not cuffmerge_path:
    cuffmerge_path = shutil.which('cuffmerge')
  check_exe(cuffmerge_path)

  if not samtools_path:
    samtools_path = shutil.which('samtools')
  check_exe(samtools_path)
  
  cufflinks_dir = create_dir(main_output_dir,'cufflinks_output')
  cuffmerge_dir = create_dir(cufflinks_dir,'merged_asm')
  cuffmerge_input_file = os.path.join(cufflinks_dir,'cufflinks_assemblies.txt')
  star_dir = os.path.join(main_output_dir,'star_output')

  if(os.path.exists(star_dir)):
    print("Found a Star output dir, will load sam file")

  sam_files = []
  for sam_file in glob.glob(star_dir + "/*.sam"):
    sam_files.append(sam_file)

  if not sam_files:
    raise IndexError('The list of sam files is empty, expected them in Star output dir. Star dir:\n%s' % star_dir)


  sorted_bam_files = []
  for sam_file in sam_files:
    sam_file_name = os.path.basename(sam_file)
    sam_temp_file_path = os.path.join(star_dir,(sam_file_name + ".tmp"))
    bam_sort_file_path = os.path.join(star_dir,re.sub('.sam','.bam',sam_file_name))

    if os.path.exists(bam_sort_file_path):
      print("Found an existing bam file, will not sort sam file again. Bam file:")
      print(bam_sort_file_path)

    else:
      print("Converting samfile into sorted bam file. Bam file:")
      print(bam_sort_file_path)
      subprocess.run(['samtools','sort','-@',str(num_threads),'-T',sam_temp_file_path,'-o',bam_sort_file_path,sam_file])

    sorted_bam_files.append(bam_sort_file_path)

  for sorted_bam_file in sorted_bam_files:
    sorted_bam_file_name = os.path.basename(sorted_bam_file)
    transcript_file_name = re.sub('.bam','.gtf',sorted_bam_file_name)
    skipped_file_name = re.sub('.bam','.skipped.gtf',sorted_bam_file_name)
    genes_fpkm_file_name = re.sub('.bam','.genes.fpkm',sorted_bam_file_name)
    isoforms_fpkm_file_name = re.sub('.bam','.isoforms.fpkm',sorted_bam_file_name)

    if os.path.exists(os.path.join(cufflinks_dir,transcript_file_name)):
      print("Found an existing cufflinks gtf file, will not overwrite. File found:")
      print(os.path.join(cufflinks_dir,transcript_file_name))
    else:
      print("Running cufflinks on: " + sorted_bam_file_name)
      print("Writing output to: " + os.path.join(cufflinks_dir,transcript_file_name))
      subprocess.run([cufflinks_path,'--output-dir',cufflinks_dir,'--num-threads',str(num_threads),sorted_bam_file])
      subprocess.run(['mv',os.path.join(cufflinks_dir,'transcripts.gtf'),os.path.join(cufflinks_dir,transcript_file_name)])
      subprocess.run(['mv',os.path.join(cufflinks_dir,'skipped.gtf'),os.path.join(cufflinks_dir,skipped_file_name)])
      subprocess.run(['mv',os.path.join(cufflinks_dir,'genes.fpkm_tracking'),os.path.join(cufflinks_dir,genes_fpkm_file_name)])
      subprocess.run(['mv',os.path.join(cufflinks_dir,'isoforms.fpkm_tracking'),os.path.join(cufflinks_dir,isoforms_fpkm_file_name)])

  # Now need to merge
  print("Creating cuffmerge input file: " + cuffmerge_input_file)

  # Note that I'm writing the subprocess this way because python seems to have issues with wildcards in subprocess.run and this
  # was the answer I found most often from googling
  gtf_list_cmd = 'ls ' + os.path.join(cufflinks_dir,'*.gtf') + ' | grep -v ".skipped." >' + cuffmerge_input_file
  gtf_list_cmd = subprocess.Popen(gtf_list_cmd,shell=True)
  gtf_list_cmd.wait()

  subprocess.run(['python2.7',cuffmerge_path,'-s',genome_file,'-p',str(num_threads),'-o',cuffmerge_dir,cuffmerge_input_file])


def run_stringtie_assemble(stringtie_path,samtools_path,main_output_dir,genome_file,num_threads):

  if not stringtie_path:
    stringtie_path = shutil.which('stringtie')
  check_exe(stringtie_path)

  if not samtools_path:
    samtools_path = shutil.which('samtools')
  check_exe(samtools_path)

  stringtie_dir = create_dir(main_output_dir,'stringtie_output')
  stringtie_merge_input_file = os.path.join(stringtie_dir,'stringtie_assemblies.txt')
  stringtie_merge_output_file = os.path.join(stringtie_dir,'annotation.gtf')
  star_dir = os.path.join(main_output_dir,'star_output')

  if(os.path.exists(star_dir)):
    print("Found a Star output dir, will load sam file")

  sorted_bam_files = []
  for bam_file in glob.glob(star_dir + "/*.bam"):
    sorted_bam_files.append(bam_file)

  if not sorted_bam_files:
    raise IndexError('The list of sorted bam files is empty, expected them in Star output dir. Star dir:\n%s' % star_dir)

  for sorted_bam_file in sorted_bam_files:
    sorted_bam_file_name = os.path.basename(sorted_bam_file)
    transcript_file_name = re.sub('.bam','.gtf',sorted_bam_file_name)
    transcript_file_path = os.path.join(stringtie_dir,transcript_file_name)

    if os.path.exists(transcript_file_path):
      print("Found an existing stringtie gtf file, will not overwrite. File found:")
      print(transcript_file_path)
    else:
      print("Running Stringtie on: " + sorted_bam_file_name)
      print("Writing output to: " + transcript_file_path)
      subprocess.run([stringtie_path,sorted_bam_file,'-o',transcript_file_path,'-p',str(num_threads),'-t','-a','15'])

  # Now need to merge
  print("Creating Stringtie merge input file: " + stringtie_merge_input_file)

  # Note that I'm writing the subprocess this way because python seems to have issues with wildcards in subprocess.run and this
  # was the answer I found most often from googling
  gtf_list_cmd = 'ls ' + os.path.join(stringtie_dir,'*.gtf') + ' | grep -v "annotation.gtf" >' + stringtie_merge_input_file
  gtf_list_cmd = subprocess.Popen(gtf_list_cmd,shell=True)
  gtf_list_cmd.wait()

  if os.path.exists(stringtie_merge_output_file):
    print("Found an existing stringtie merge file, will not overwrite. File found:")
    print(stringtie_merge_output_file)
  else:
    print("Merging Stringtie results. Writing to the following file:")
    print(stringtie_merge_output_file)
    # Note, I'm not sure stringtie merge actually uses threads, but it doesn't complain if -p is passed in
    subprocess.run([stringtie_path,'--merge','-p',str(num_threads),'-o',stringtie_merge_output_file,stringtie_merge_input_file,'-i'])


def run_scallop_assemble(scallop_path,stringtie_path,main_output_dir):

  if not scallop_path:
    scallop_path = shutil.which('scallop')
  check_exe(scallop_path)

  if not stringtie_path:
    stringtie_path = shutil.which('stringtie')
  check_exe(stringtie_path)

  scallop_dir = create_dir(main_output_dir,'scallop_output')
  stringtie_merge_input_file = os.path.join(scallop_dir,'scallop_assemblies.txt')
  stringtie_merge_output_file = os.path.join(scallop_dir,'annotation.gtf')
  star_dir = os.path.join(main_output_dir,'star_output')

  if(os.path.exists(star_dir)):
    print("Found a Star output dir, will load sam file")

  sorted_bam_files = []
  for bam_file in glob.glob(star_dir + "/*.bam"):
    sorted_bam_files.append(bam_file)

  if not sorted_bam_files:
    raise IndexError('The list of sorted bam files is empty, expected them in Star output dir. Star dir:\n%s' % star_dir)

  for sorted_bam_file in sorted_bam_files:
    sorted_bam_file_name = os.path.basename(sorted_bam_file)
    transcript_file_name = re.sub('.bam','.gtf',sorted_bam_file_name)
    transcript_file_path = os.path.join(scallop_dir,transcript_file_name)

    if os.path.exists(transcript_file_path):
      print("Found an existing scallop gtf file, will not overwrite. File found:")
      print(transcript_file_path)
    else:
      print("Running Scallop on: " + sorted_bam_file_name)
      print("Writing output to: " + transcript_file_path)
      subprocess.run([scallop_path,'-i',sorted_bam_file,'-o',transcript_file_path,'--min_flank_length','10'])


  # Now need to merge
  print("Creating Stringtie merge input file: " + stringtie_merge_input_file)

  # Note that I'm writing the subprocess this way because python seems to have issues with wildcards in subprocess.run and this
  # was the answer I found most often from googling
  gtf_list_cmd = 'ls ' + os.path.join(scallop_dir,'*.gtf') + ' | grep -v "annotation.gtf" >' + stringtie_merge_input_file
  gtf_list_cmd = subprocess.Popen(gtf_list_cmd,shell=True)
  gtf_list_cmd.wait()

  if os.path.exists(stringtie_merge_output_file):
    print("Found an existing stringtie merge file, will not overwrite. File found:")
    print(stringtie_merge_output_file)
  else:
    print("Merging Stringtie results. Writing to the following file:")
    print(stringtie_merge_output_file)
    # Note, I'm not sure stringtie merge actually uses threads and is very quick, so leaving out
    subprocess.run([stringtie_path,'--merge','-o',stringtie_merge_output_file,stringtie_merge_input_file])


def splice_junction_to_gff(input_dir,hints_file):
  
  sjf_out = open(hints_file,"w+")
  
  for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
    sjf_in = open(sj_tab_file)
    sjf_lines = sjf_in.readlines()
    for line in sjf_lines:
      elements = line.split('\t')
      strand = '+'
      # If the strand is undefined then skip, Augustus expects a strand
      if elements[3] == '0':
        continue
      elif elements[3] == '2':
        strand = '-'

      junction_length = int(elements[2]) - int(elements[1]) + 1
      if junction_length < 100:
        continue

      if not elements[4] and elements[7] < 10:
        continue
       
      # For the moment treat multimapping and single mapping things as a combined score
      score = float(elements[6]) + float(elements[7])
      score = str(score)
      output_line = [elements[0],'RNASEQ','intron',elements[1],elements[2],score,strand,'.',('src=W;mul=' + score + ';')]
      sjf_out.write('\t'.join(output_line) + '\n')

  sjf_out.close()


def model_builder(work_dir):

  star_output_dir = os.path.join(work_dir,'star_output')

  all_junctions_file = os.path.join(star_output_dir,'all_junctions.sj')
  sjf_out = open(all_junctions_file,"w+")

  for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
    sjf_in = open(sj_tab_file)
    sjf_lines = sjf_in.readlines()
    for line in sjf_lines:
      elements = line.split('\t')
      strand = '+'

#    my $slice_name = $eles[0];
#    my $start = $eles[1];
#    my $end = $eles[2];
#    my $strand = $eles[3];

      # If the strand is undefined then skip, Augustus expects a strand
      if elements[3] == '0':
        continue
      elif elements[3] == '2':
        strand = '-'

      junction_length = int(elements[2]) - int(elements[1]) + 1
      if junction_length < 100:
        continue

      if not elements[4] and elements[7] < 10:
        continue

      # For the moment treat multimapping and single mapping things as a combined score
      score = float(elements[6]) + float(elements[7])
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


def get_seq_region_lengths(genome_file,min_seq_length):
  current_header = ""
  current_seq = ""

  seq_regions = {}
  file_in = open(genome_file)
  line = file_in.readline()
  while line:
    match = re.search(r'>(.+)$',line)
    if match and current_header:
      if len(current_seq) > min_seq_length:
        seq_regions[current_header] = len(current_seq)
 
      current_seq = ""
      current_header = match.group(1)
    elif match:
      current_header = match.group(1)
    else:
      current_seq += line.rstrip()

    line = file_in.readline()

  if len(current_seq) > min_seq_length:
    seq_regions[current_header] = len(current_seq)

  return seq_regions


def run_finalise_geneset(main_output_dir,genome_file,seq_region_names,db_details,num_threads):
  final_annotation_dir = create_dir(main_output_dir,'annotation')
  db_connection = db_details.split(',')
  merged_gtf = os.path.join(main_output_dir,'annotation_to_finalise.gtf')
  file_out = open(merged_gtf, 'w')
  annotation_dirs = ['genblast_output','stringtie_output','scallop_output','minimap2_output']
  for annotation_dir in annotation_dirs:
    gtf_file = os.path.join(main_output_dir,annotation_dir,'annotation.gtf')
    if not os.path.exists(gtf_file):
      print("No annotation.gtf file found in " + annotation_dir + ", skipping")
      continue

    file_in = open(gtf_file)
    line = file_in.readline()
    while line:
      print(line.rstrip(),file=file_out)
      line = file_in.readline()
    file_in.close()
  file_out.close()

  pool = multiprocessing.Pool(int(num_threads))
  for seq_region_name in seq_region_names:
    pool.apply_async(multiprocess_finalise_geneset, args=(seq_region_name,merged_gtf,db_connection,final_annotation_dir,))

  pool.close()
  pool.join()


def multiprocess_finalise_geneset(seq_region_name,merged_gtf,db_connection,final_annotation_dir):

  output_dbname = db_connection[0]
  output_server = db_connection[1]
  output_port = db_connection[2]
  output_user = db_connection[3]
  output_pass = db_connection[4]

  dna_dbname = db_connection[0]
  dna_server = db_connection[1]
  dna_port = db_connection[2]
  dna_user = db_connection[5]

  finalise_cmd = ['perl','/homes/fergal/enscode/ensembl-common/scripts/process_transcriptomic_gtf.pl','-dbname',output_dbname,'-host',output_server,'-user',output_user,'-port',output_port,'-pass',output_pass,'-dna_dbname',dna_dbname,'-dna_host',dna_server,'-dna_user',dna_user,'-dna_port',dna_port,'-gtf_file',merged_gtf,'-specify_seq_region_name',seq_region_name,'-output_path',final_annotation_dir]

  print('Finalising ' + seq_region_name)
  print(' '.join(finalise_cmd))
  subprocess.run(finalise_cmd)


def get_sequence(seq_region,start,end,strand,fasta_file,output_dir):
  start -= 1
  bedtools_path = 'bedtools'

  # This creates a tempfile and writes the bed info to it based on whatever information
  # has been passed in about the sequence. Then runs bedtools getfasta. The fasta file
  # should have a faidx. This can be created with the create_faidx static method prior
  # to fetching sequence
  with tempfile.NamedTemporaryFile(mode='w+t', delete=False, dir=output_dir) as bed_temp_file:
    bed_temp_file.writelines(seq_region + "\t" + str(start) + "\t" + str(end))
    bed_temp_file.close()

  bedtools_command = [bedtools_path, 'getfasta','-fi', fasta_file,'-bed',bed_temp_file.name]
  bedtools_output = subprocess.Popen(bedtools_command, stdout=subprocess.PIPE)
  for idx,line in enumerate(io.TextIOWrapper(bedtools_output.stdout, encoding="utf-8")):
    if idx == 1:
      if strand == 1:
        sequence = line.rstrip()
      else:
        sequence = reverse_complement(line.rstrip())

  os.remove(bed_temp_file.name)
  return sequence


def reverse_complement(sequence):
  rev_matrix = str.maketrans("atgcATGC", "tacgTACG")
  return sequence.translate(rev_matrix)[::-1]



def seq_region_names(genome_file):
  region_list = []

  file_in = open(genome_file)
  line = file_in.readline()
  while line:
    match = re.search(r'>([^\s]+)',line)
    if match:
      region_name = match.group(1)
      if region_name == "MT":
        print ("Skipping region named MT")
        line = file_in.readline()
        continue
      else:
        region_list.append(match.group(1))
    line = file_in.readline()

  return region_list


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
  
  
def create_paired_paths(fastq_file_paths):
  path_dict = {}
  final_list = []

  for path in fastq_file_paths:
    match = re.search(r'(.+)_\d+\.(fastq|fq)',path)
    if not match:
      print("Could not find _1 or _2 at the end of the prefix for file. Assuming file is not paired:")
      print(path)
      final_list.append([path])
      continue

    prefix = match.group(1)
    if prefix in path_dict:
#      path_dict[prefix] = path_dict[prefix] + ',' + path
      path_dict[prefix].append(path)
    else:
      path_dict[prefix] = [path]

  for pair in path_dict:
    final_list.append(path_dict[pair])

  return(final_list)


def check_exe(exe_path):

  if not shutil.which(exe_path):
    raise OSError('Exe does not exist. Path checked: %s' % exe_path)

  

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--output_dir', help='Path where the output and temp files will write to. Uses current dir by default', required=False)
  parser.add_argument('--genome_file', help='Path to the fasta genome file', required=True)
  parser.add_argument('--num_threads', type=int, help='Number of threads to use', required=False)
  parser.add_argument('--run_masking', help='Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome', required=False)
  parser.add_argument('--red_path', help='Path to Red executable. See http://toolsmith.ens.utulsa.edu', required=False)
  parser.add_argument('--genblast_path', help='Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html', required=False)
  parser.add_argument('--convert2blastmask_path', help='Path to convert2blastmask executable', required=False)
  parser.add_argument('--makeblastdb_path', help='Path to makeblastdb executable', required=False)
  parser.add_argument('--run_genblast', help='Run GenBlast to align protein sequences', required=False)
  parser.add_argument('--protein_file', help='Path to a fasta file with protein sequences', required=False)
  parser.add_argument('--run_star', help='Run Star for short read alignment', required=False)
  parser.add_argument('--star_path', help='Path to Star for short read alignment', required=False)
  parser.add_argument('--max_reads_per_sample', nargs='?', const=0, type=int, help='The maximum number of reads to use per sample. Default=0 (unlimited)', required=False)
  parser.add_argument('--max_total_reads', nargs='?', const=0, type=int, help='The maximum total number of reads. Default=0 (unlimited)', required=False)
  parser.add_argument('--short_read_fastq_dir', help='Path to short read fastq dir for running with Star', required=False)
  parser.add_argument('--run_minimap2', help='Run minimap2 for long read alignment', required=False)
  parser.add_argument('--minimap2_path', help='Path to minimap2 for long read alignment', required=False)
  parser.add_argument('--paftools_path', help='Path to paftools for SAM to BED conversion', required=False)
  parser.add_argument('--long_read_fastq_dir', help='Path to long read fastq dir for running with minimap2', required=False)
  parser.add_argument('--run_augustus', help='Run Augustus with hints for gene/transcript prediction', required=False)
  parser.add_argument('--augustus_path', help='Path to Augustus', required=False)
  parser.add_argument('--run_cufflinks', help='Run Cufflinks on the results from the STAR alignments', required=False)
  parser.add_argument('--cufflinks_path', help='Path to Cufflinks', required=False)
  parser.add_argument('--cuffmerge_path', help='Path to Cuffmerge', required=False)
  parser.add_argument('--run_stringtie', help='Run Stringtie on the results from the STAR alignments', required=False)
  parser.add_argument('--run_scallop', help='Run Scallop on the results from the STAR alignments', required=False)
  parser.add_argument('--stringtie_path', help='Path to Stringtie', required=False)
  parser.add_argument('--scallop_path', help='Path to Scallop', required=False)
  parser.add_argument('--subsample_script_path', help='Path to gbiab subsampling script', required=False)
  parser.add_argument('--samtools_path', help='Path to subsampling script', required=False)
  parser.add_argument('--finalise_geneset', help='Used to finalise the gene set from the various GTF files generated', required=False)
  parser.add_argument('--db_details', help='A comma separated stinrg of dbname,host,port,user,pass', required=False)
  parser.add_argument('--run_sncrna', help='Search for sncRNA structures using Rfam and cmsearch', required=False)
  parser.add_argument('--annotate', help='Run a generalised annotation, will automatically check for input data and run tools based on that', required=False)
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
  run_genblast = args.run_genblast
  protein_file = args.protein_file
  run_star = args.run_star
  star_path = args.star_path
  short_read_fastq_dir = args.short_read_fastq_dir
  max_reads_per_sample = args.max_reads_per_sample
  max_total_reads = args.max_total_reads
  run_minimap2 = args.run_minimap2
  minimap2_path = args.minimap2_path
  paftools_path = args.paftools_path
  long_read_fastq_dir = args.long_read_fastq_dir
  run_augustus = args.run_augustus
  augustus_path = args.augustus_path
  run_cufflinks = args.run_cufflinks
  cufflinks_path = args.cufflinks_path
  cuffmerge_path = args.cuffmerge_path
  run_stringtie = args.run_stringtie
  run_scallop = args.run_scallop
  stringtie_path = args.stringtie_path
  scallop_path = args.scallop_path
  subsample_script_path = args.subsample_script_path
  samtools_path = args.samtools_path
  finalise_geneset = args.finalise_geneset
  db_details = args.db_details
  run_sncrna = args.run_sncrna
  annotate = args.annotate

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

  # If the annotate flag is set then we want to set a standardised set of analyses
  # The flag is passed in to each method to indicate that it's okay if data aren't present as it is a generalised run
  if annotate:
    run_masking = 1
    run_genblast = 1
    run_star = 1
    run_minimap2 = 1
    run_stringtie = 1
    run_scallop = 1
    finalise_geneset = 1

  # Collect a list of seq region names, most useful for multiprocessing regions
#  seq_region_names = seq_region_names(genome_file)  
#  for i in seq_region_names:
#    print(i)

  # Run masking
  if run_masking:
    print ("Running masking via Red")
    masked_genome_file = run_red(red_path,work_dir,genome_file)
    print ("Masked genome file: " + masked_genome_file)

  else:
    print ("Not running masking, presuming the genome file is softmasked")

  # Run GenBlast
  if run_genblast:
    print ("Running GenBlast")
    run_genblast_align(genblast_path,convert2blastmask_path,makeblastdb_path,work_dir,protein_file,masked_genome_file,num_threads)

  # Run STAR
  if run_star:
     print ("Running Star")
     run_star_align(star_path,subsample_script_path,work_dir,short_read_fastq_dir,genome_file,max_reads_per_sample,max_total_reads,num_threads)

  # Run minimap2
  if run_minimap2:
     print ("Running minimap2")
     run_minimap2_align(minimap2_path,paftools_path,work_dir,long_read_fastq_dir,genome_file,num_threads)

  # Run Augustus
  if run_augustus:
     print ("Running Augustus")
     run_augustus_predict(augustus_path,work_dir,masked_genome_file,num_threads)

  # Run Stringtie
  if run_stringtie:
     print ("Running Stringtie")
     run_stringtie_assemble(stringtie_path,samtools_path,work_dir,genome_file,num_threads)

  # Run Scallop
  if run_scallop:
     print ("Running Scallop")
     run_scallop_assemble(scallop_path,stringtie_path,work_dir)

  # Run Cufflinks
  if run_cufflinks:
     print ("Running Cufflinks")
     run_cufflinks_assemble(cufflinks_path,cuffmerge_path,samtools_path,work_dir,genome_file,num_threads)

  # Do some magic
  if finalise_geneset:
     print("Finalise geneset")
     run_finalise_geneset(work_dir,genome_file,seq_region_names,db_details,num_threads)


  if run_sncrna:
     print("Annotating sncRNAs")
     search_rfam(genome_file,None,None,None,work_dir,"insect",num_threads)

