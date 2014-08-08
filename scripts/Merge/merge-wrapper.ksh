#!/usr/bin/ksh
# $Id: merge-wrapper.ksh,v 1.15 2014-01-23 15:18:47 ak4 Exp $

########################################################################
# CONFIGURATION
#

# Read configuration from the file supplied.

if [[ ! -f "$1" ]]; then
  printf "Can not open file '%s' for reading (config)\n" "$1"
  exit 1
fi

. "$1"

#
# END CONFIGURATION
########################################################################

if [[ -n ${output_dir} && -d ${output_dir} ]]; then
  printf ">>> Output directory '%s' already exists!\n" ${output_dir}
  exit 1
fi

if [[ -z ${output_dir} ]]; then
  output_dir=$( mktemp -d merge_output.XXXX )
else
  if ! mkdir ${output_dir}; then
    printf ">>> Could not create directory '%s'\n" ${output_dir}
    exit 1
  fi
fi

printf ">>> Output directory is '%s'\n" ${output_dir}

bsub -q normal \
  -J "merge-run[1-${njobs}]%${concurrent}" \
  -oo "${output_dir}/merge-run-%I.out" \
  -eo "${output_dir}/merge-run-%I.err" \
  -M 1500 -R 'select[mem>1500]' -R 'rusage[mem=1500]' \
  -We 30 \
  perl ${ensembl_analysis_base}/scripts/Merge/merge.pl \
  --host_ensembl="${host_ensembl}" \
  --user_ensembl="${rouser}" \
  --password_ensembl="${ropassword}" \
  --database_ensembl="${database_ensembl}" \
  --host_havana="${host_havana}" \
  --user_havana="${rouser}" \
  --password_havana="${ropassword}" \
  --database_havana="${database_havana}" \
  --host_dna="${host_dna}" \
  --port_dna="${port_dna}" \
  --user_dna="${user_dna}" \
  --database_dna="${database_dna}" \
  --host_ccds="${host_ccds}" \
  --user_ccds="${rouser}" \
  --password_ccds="${ropassword}" \
  --database_ccds="${database_ccds}" \
  --host_output="${host_output}" \
  --user_output="${rwuser}" \
  --password_output="${rwpassword}" \
  --database_output="${database_output}" \
  --ensembl_include="${ensembl_include}" \
  --ensembl_exclude="${ensembl_include}" \
  --havana_include="${havana_include}" \
  --havana_exclude="${havana_exclude}" \
  --ensembl_tag="${ensembl_tag}" \
  --havana_tag="${havana_tag}" \
  --havana_gene_xref="${havana_gene_xref}" \
  --havana_transcript_xref="${havana_transcript_xref}" \
  --havana_translation_xref="${havana_translation_xref}" \
  --njobs="${njobs}" --job="\$LSB_JOBINDEX" $*

# Post-processing:  Copy all unprocessed Ensembl genes to the output
# database.

bsub -q normal \
  -J 'merge-copy-submit' \
  -w 'done(merge-run)' \
  -oo "${output_dir}/merge-copy-submit.out" \
  -eo "${output_dir}/merge-copy-submit.err" \
  <<MERGE_COPY_SUBMIT_END
cd "${output_dir}" || exit 1

awk '\$1 == "PROCESSED" { print \$2 }' merge-run-*.out |
sort -u -n -o genes-processed.txt

mysql -BN -h '${host_ensembl}' -u ensro -D '${database_ensembl}' \
  >genes-all.txt <<SQL_END
SELECT gene_id
FROM gene
JOIN seq_region USING (seq_region_id)
ORDER BY gene_id
SQL_END

diff genes-all.txt genes-processed.txt |
awk '\$1 == "<" { print \$2 }' |
sort -R -o genes-copy.txt

rm -f genes-copy-*
split -n l/${concurrent} genes-copy.txt genes-copy-

i=0
for list in genes-copy-*; do
  (( ++i ))
  sort -n -o \${list} \${list}
  bsub -q normal \
    -J "merge-copy-\${i}" \
    -oo "merge-copy-\${i}.out" \
    -eo "merge-copy-\${i}.err" \
    -M 500 -R 'select[mem>500]' -R 'rusage[mem=500]' \
    perl ${ensembl_analysis_base}/scripts/genebuild/copy_genes.pl \
    --file="\${list}" \
    --sourcehost='${host_ensembl}' \
    --sourceuser='${rouser}' \
    --sourcedbname='${database_ensembl}' \
    --targethost='${host_output}' \
    --targetuser='${rwuser}' \
    --targetpass='${rwpassword}' \
    --targetdbname='${database_output}' \
    --verbose
done
MERGE_COPY_SUBMIT_END

bsub -q small \
  -J 'merge-cleanup' \
  -w 'exit(merge-run) || ended(merge-copy-*)' \
  -oo "${output_dir}/merge-cleanup.out" \
  -eo "${output_dir}/merge-cleanup.err" \
  -M 100 -R 'select[mem>100]' -R 'rusage[mem=100]' \
  <<MERGE_CLEANUP_END
bkill -J 'merge-copy-submit' || exit 0
MERGE_CLEANUP_END
