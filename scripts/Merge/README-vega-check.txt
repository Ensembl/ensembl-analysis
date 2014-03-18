
Example submission to LSF:

  bsub -M800000 -R"select[mem>800] rusage[mem=800]" \
  -oo vega_check_before.out "perl $GBP/sanity_scripts/vega_check.pl \
  -dbname mus_musculus_vega_fixed_nn -dbhost genebuildn -dnadbname mus_musculus_ensembl_nn \
  -dnadbhost genebuildn -coord_system toplevel -path GRCm38 -sql_output vega_check_before.sql"

Look at output:

  echo "count gene_biotype, transcript_biotype." ; grep "is not allowed" vega_check_before.out | awk '{print $9,$18}' | sort | uniq -c

Suggested SQL:

  cat $SCR9/vega_check_before.sql
update transcript set biotype='processed_transcript' where stable_id='OTTMUST000000nnnnn';
etc.


