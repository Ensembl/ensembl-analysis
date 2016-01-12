##
# How to craete a UniProt database
##

# To create a unirpt database you need to have at least ensembl/modules, ensembl-hive/modules, ensembl-analysis/modules and BioPerl in your PERL5LIB
# we will use the Hive to create the databases so you will need to run the steps in screen
# You need to source this file
# Then run the init command printed
# Run sync command then the loop command given by the init script
#
# The script will check the version of UniProt and stops if you already have the database in your BASE_UNIPROT_PATH

HIVE_HOST=''
HIVE_USER=''
HIVE_PASS=''
HIVE_PORT=3306

BASE_UNIPROT_PATH="/data/blastdb/Ensembl"
export EMBL2FASTA_SCRIPT="$ENSCODE/ensembl-analysis/scripts/databases/embl2fasta.pl"
export PROCESS_ISOFORMS_SCRIPT="$ENSCODE/ensembl-analysis/scripts/databases/process_uniprot_isoforms.pl"

###
## You should only change values above this line.
## Parameters below this line should stay as they are
###
for S in "$EMBL2FASTA_SCRIPT" "$PROCESS_ISOFORMS_SCRIPT"; do
    if [ ! -e "$S" ]; then
        printf " \e[31m%s\e[0m does not exist\n" "$S"
    fi
done

for S in "ensembl-hive" "ensembl" "ensembl-analysis"; do
    if [ "`echo $PERL5LIB | sed \"s/$S\/modules//\"`" = "$PERL5LIB" ]; then
        printf " \e[31m%s\e[0m repository is not in your PERL5LIB\n" "$S"
    fi
done

export UNIPROT_VERSION=`wget -S --spider www.uniprot.org 2>&1 | grep 'X-UniProt-Release' | awk '{print $2}'`
export UNIPROT_DIR="$BASE_UNIPROT_PATH/uniprot_$UNIPROT_VERSION"
if [ -e "$UNIPROT_DIR" ]; then
    printf "The directory \e[31m%s\e[0m already exists:\n \e[32m-\e[0m this database has already been successfully created\n \e[31m-\e[0m delete the directory %s and start the pipeline\n\n" "$UNIPROT_DIR" "$UNIPROT_DIR"
fi
export UNIPROT_DATE=`wget -S --spider www.uniprot.org 2>&1 | grep 'Last-Modified' | sed 's/\s*Last-Modified:\s\+//'`
export UNIPROT_VERT_FILE="uniprot_vertebrata"
export VARSPLIC_FILE="uniprot_sprot_varsplic"
export UNIPROT_FILE='uniprot'
# Full connection info for the hive pipeline db, needs write access
pipeline_name="create_uniprot_$UNIPROT_VERSION"
HIVE_DBNAME="${USER}_$pipeline_name"
echo "Using ${HIVE_DBNAME} on ${HIVE_HOST} as hive database"

echo "Using $UNIPROT_DIR as working directory"

hive_url="mysql://${HIVE_USER}:${HIVE_PASS}@${HIVE_HOST}:${HIVE_PORT}/${HIVE_DBNAME}"
echo "Hive url: ${hive_url}"

init_pipe="ensembl-hive/scripts/init_pipeline.pl"
if [ -e "$ENSCODE/$init_pipe" ]; then
    init_pipe="$ENSCODE/$init_pipe"
fi
printf "\nInitiating the pipeline:\n \e[32m$\e[0m perl $init_pipe UniProtDB_conf -host ${HIVE_HOST} -port ${HIVE_PORT} -user ${HIVE_USER} -password ${HIVE_PASS} -dbname ${HIVE_DBNAME} -pipeline_name ${pipeline_name}\n"
unset init_pipe
