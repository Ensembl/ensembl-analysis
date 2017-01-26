# Generic config for eg protein anno pipeline
WORK_DIR_FULL=/lustre/scratch109/ensembl/dm15/hive_humancdna_84_rerun_2/

export ENSGBSCRIPTS=/nfs/users/nfs_d/dm15/cvs_checkout_head/ensembl-personal/genebuilders/scripts
export ENSEMBL_REPO_ROOT=/nfs/users/nfs_d/dm15/enscode/
export EG_REPO_ROOT=/software/ensembl/genebuild/lib/eg_protein_annotation
PERL5LIB=""
PERL5LIB=${ENSEMBL_REPO_ROOT}/ensembl/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-hive/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-analysis/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-compara/scripts/pipeline
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-hive/scripts
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-production/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-killlist/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-killlist/scripts
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-pipeline/modules
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-pipeline/test_system/config
PERL5LIB=${PERL5LIB}:${EG_REPO_ROOT}/eg-proteinfeature/lib
PERL5LIB=${PERL5LIB}:${EG_REPO_ROOT}/eg-release/lib
PERL5LIB=${PERL5LIB}:${EG_REPO_ROOT}/eg-utils/lib
PERL5LIB=${PERL5LIB}:${EG_REPO_ROOT}/eg-pipelines/modules
PERL5LIB=${PERL5LIB}:${EG_REPO_ROOT}/String
PERL5LIB=${PERL5LIB}:${ENSEMBL_REPO_ROOT}/ensembl-io/modules
PERL5LIB=${PERL5LIB}:/software/ensembl/genebuild/lib/bioperl-1.2.3/
PATH=${PATH}:/software/ensembl/genebuild/usrlocalensemblbin
PATH=${PATH}:/software/ensembl/genebuild/usrlocalensembllib
PATH=${PATH}:${ENSEMBL_REPO_ROOT}/ensembl-compara/scripts/pipeline
PATH=${PATH}:${ENSEMBL_REPO_ROOT}/ensembl-hive/scripts


export PERL5LIB
export PATH

echo PERL5LIB is:
/software/ensembl/genebuild/bin/print_path $PERL5LIB
echo

echo PATH is:
/software/ensembl/genebuild/bin/print_path $PATH
echo

# Give the pipeline a name, e.g. human_prot_anno_e80
pipeline_name=human_cdna_84_rerun_2

# Full connection info for the hive pipeline db, needs write access
HIVE_HOST=genebuild11
HIVE_USER=ensadmin
HIVE_PASS=
HIVE_PORT=3306
HIVE_DBNAME=dm15_human_cdna_84_rerun_2

registry=${WORK_DIR_FULL}/registry.pm
echo Registry file for this run is: ${registry}

# Fill in core db details below
echo "{
package reg;" > $registry
echo "Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => 'ens-staging1',
    -port    => 3306,
    -user    => 'ensro',
    -dbname  => 'ensembl_production',
    -species => 'multi',
    -group   => 'production'
  ); " >> $registry
echo "1;
}" >> $registry

echo Using ${HIVE_DBNAME} on ${HIVE_HOST} as hive database

pipeline_dir=${WORK_DIR_FULL}/${pipeline_name}
echo Using ${pipeline_dir} as a temporary directory

hive_url=mysql://${HIVE_USER}:${HIVE_PASS}@${HIVE_HOST}:${HIVE_PORT}/${HIVE_DBNAME}
echo Hive url: ${hive_url}

# for data big files of the project
export ENSADMIN_PSW=''

