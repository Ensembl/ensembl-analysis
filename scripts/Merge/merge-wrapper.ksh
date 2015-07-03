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
  -M 1500 -R"select[mem>1500] rusage[mem=1500]" \
  -We 30 \
  perl ${ensembl_analysis_base}/scripts/Merge/merge.pl \
  --host_secondary="${host_secondary}" \
  --user_secondary="${rouser}" \
  --password_secondary="${ropassword}" \
  --database_secondary="${database_secondary}" \
  --host_primary="${host_primary}" \
  --user_primary="${rouser}" \
  --password_primary="${ropassword}" \
  --database_primary="${database_primary}" \
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
  --secondary_include="${secondary_include}" \
  --secondary_exclude="${secondary_exclude}" \
  --primary_include="${primary_include}" \
  --primary_exclude="${primary_exclude}" \
  --secondary_tag="${secondary_tag}" \
  --primary_tag="${primary_tag}" \
  --primary_gene_xref="${primary_gene_xref}" \
  --primary_transcript_xref="${primary_transcript_xref}" \
  --primary_translation_xref="${primary_translation_xref}" \
  --njobs="${njobs}" --job="\$LSB_JOBINDEX" $*

# Post-processing:  Copy all unprocessed Ensembl genes to the output
# database.

bsub -q normal \
  -J 'merge-copy-submit' \
  -w 'done(merge-run)' \
  -oo "${output_dir}/merge-copy-submit.out" \
  -eo "${output_dir}/merge-copy-submit.err" \
  -M 100 -R"select[mem>100] rusage[mem=100]" \
  <<MERGE_COPY_SUBMIT_END
cd "${output_dir}" || exit 1

awk '\$1 == "PROCESSED" { print \$2 }' merge-run-*.out |
sort -u -n -o genes-processed.txt

if [ -n "${secondary_include}" ];then
    echo "INCLUDE"
    SQLQUERY=`echo ${secondary_include} | sed 's/,/","/g'`
    mysql -BN -h '${host_secondary}' -u ensro -D '${database_secondary}' \
      >genes-all.txt <<SQL_END
SELECT gene_id
FROM gene
JOIN seq_region USING (seq_region_id)
JOIN analysis USING (analysis_id)
WHERE logic_name IN ("\$SQLQUERY")
ORDER BY gene_id
SQL_END
else
    echo "NO INCLUDE"
    mysql -BN -h '${host_secondary}' -u ensro -D '${database_secondary}' \
      >genes-all.txt <<SQL_END
SELECT gene_id
FROM gene
JOIN seq_region USING (seq_region_id)
ORDER BY gene_id
SQL_END
fi

diff genes-all.txt genes-processed.txt |
awk '\$1 == "<" { print \$2 }' |
sort -R -o genes-copy.txt
if [ -n "${secondary_exclude}" ];then
    echo "EXCLUDE"
    SQLQUERY=`echo ${secondary_exclude} | sed 's/,/","/g'`
    mysql -BN -h '${host_secondary}' -u ensro -D '${database_secondary}' \
      >genes-excluded <<SQL_END
SELECT gene_id
FROM gene
JOIN seq_region USING (seq_region_id)
JOIN analysis USING (analysis_id)
WHERE logic_name IN ("\$SQLQUERY")
ORDER BY gene_id
SQL_END
    sort -n genes-copy.txt > genes-copy.tmp
    diff genes-copy.tmp genes-excluded |
    awk '\$1 == "<" { print \$2 }' |
    sort -R -o genes-copy.txt
fi

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
    -M 500 -R"select[mem>500] rusage[mem=500]" \
    perl ${ensembl_analysis_base}/scripts/genebuild/copy_genes.pl \
    --file="\${list}" \
    --sourcehost='${host_secondary}' \
    --sourceuser='${rouser}' \
    --sourcedbname='${database_secondary}' \
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
  -M 100 -R"select[mem>100] rusage[mem=100]" \
  <<MERGE_CLEANUP_END
bkill -J 'merge-copy-submit' || exit 0
MERGE_CLEANUP_END
