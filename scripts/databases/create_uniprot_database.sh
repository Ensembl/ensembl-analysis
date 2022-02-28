# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

##
# How to craete a UniProt database
##

# To create a UniProtKB database you need to have at least ensembl/modules, ensembl-hive/modules, ensembl-analysis/modules and BioPerl in your PERL5LIB
#
# The script will check the version of UniProt and stops if you already have the database in your BASE_UNIPROT_PATH
# Then run the pipeline

HIVE_HOST=''
HIVE_USER=''
HIVE_PASS=''
HIVE_PORT=3306
EHIVE_DRIVER="mysql" # This should not change unless you know what you are doing
KILL_HOST=''
KILL_USER=''
KILL_PASS=''
KILL_PORT=3306
KILL_DBNAME="gb_kill_list"

USAGE=0
INIT_PIPE=1
RUN_PIPE=1
BASE_UNIPROT_PATH="$BLASTDB_DIR/uniprot"
ENSEMBL_BASE=$ENSCODE
while getopts "s:d:h:u:p:P:IRj:k:l:m:n:" o; do
    case $o in
        d ) BASE_UNIPROT_PATH=$OPTARG;;
        s ) ENSEMBL_BASE=$OPTARG;;
        h ) HIVE_HOST=$OPTARG;;
        u ) HIVE_USER=$OPTARG;;
        p ) HIVE_PASS=$OPTARG;;
        P ) HIVE_PORT=$OPTARG;;
        j ) KILL_HOST=$OPTARG;;
        k ) KILL_USER=$OPTARG;;
        l ) KILL_PASS=$OPTARG;;
        m ) KILL_PORT=$OPTARG;;
        n ) KILL_DBNAME=$OPTARG;;
        I ) RUN_PIPE=0;;
        R ) INIT_PIPE=0;;
        * ) USAGE=1;;
    esac
done
export EMBL2FASTA_SCRIPT="$ENSEMBL_BASE/ensembl-analysis/scripts/databases/EMBL2fasta.pl"
export PROCESS_ISOFORMS_SCRIPT="$ENSEMBL_BASE/ensembl-analysis/scripts/databases/process_uniprot_isoforms.pl"

if [ $USAGE -eq 1 ]; then
  echo "You should not need to change anything, it should by default init the pipeline and start running it."
  echo "  d Destination directory, default is $BLASTDB_DIR/uniprot"
  echo "  s Base directory where your Perl API are, default is $ENSCODE"
  echo "  h Host server for the Hive pipeline database"
  echo "  u User name for the Hive pipeline database"
  echo "  p Password for the Hive pipeline database"
  echo "  P Port for the Hive pipeline database, default is 3306"
  echo "  j Host server for the kill list database, default to Hive host"
  echo "  k User name for the kill list database, default to Hive user"
  echo "  l Password for the kill list database, default to Hive password"
  echo "  m Port for the Hive pipeline database, default to Hive port"
  echo "  n Name of the kill list database, default is gb_kill_list"
  echo "  I To only run the init script"
  echo "  R To only run the pipeline"
  echo
  exit 1
fi

###
## You should only change values above this line.
## Parameters below this line should stay as they are
###
for S in "$EMBL2FASTA_SCRIPT" "$PROCESS_ISOFORMS_SCRIPT"; do
    if [ ! -e "$S" ]; then
        printf " \e[31m%s\e[0m does not exist\n" "$S"
    fi
done

for S in "ensembl-hive" "ensembl" "ensembl-analysis" "ensembl-killlist"; do
    if [ "`echo $PERL5LIB | sed \"s/$S\/modules//\"`" = "$PERL5LIB" ]; then
        printf " \e[31m%s\e[0m repository is not in your PERL5LIB\n" "$S"
    fi
done

if [ -z "$KILL_HOST" ]; then
  KILL_HOST="$HIVE_HOST";
fi

if [ -z "$KILL_PORT" ]; then
  KILL_PORT="$HIVE_PORT";
fi

if [ -z "$KILL_USER" ]; then
  KILL_USER="$HIVE_USER";
fi

if [ -z "$KILL_PASS" ]; then
  KILL_PASS="$HIVE_PASS";
fi

export UNIPROT_VERSION=`wget -S --spider www.uniprot.org 2>&1 | grep 'X-UniProt-Release' | awk '{print $2}'`
export UNIPROT_DIR="$BASE_UNIPROT_PATH/uniprot_$UNIPROT_VERSION"
if [ -e "$UNIPROT_DIR" ]; then
    printf "The directory \e[31m%s\e[0m already exists:\n \e[32m-\e[0m this database has already been successfully created\n \e[31m-\e[0m delete the directory %s and start the pipeline\n\n" "$UNIPROT_DIR" "$UNIPROT_DIR"
fi
export UNIPROT_DATE=`wget -S --spider www.uniprot.org 2>&1 | grep 'Last-Modified' | sed 's/\s*Last-Modified:\s\+//'`
# Full connection info for the hive pipeline db, needs write access
pipeline_name="create_uniprot_$UNIPROT_VERSION"
HIVE_DBNAME="${USER}_$pipeline_name"
echo "Using ${HIVE_DBNAME} on ${HIVE_HOST} as hive database"
echo "Using ${KILL_DBNAME} on ${KILL_HOST} as kill list database"

echo "Using $UNIPROT_DIR as working directory"

init_pipe="init_pipeline.pl"
if [ -e "$ENSEMBL_BASE/ensembl-hive/scripts/$init_pipe" ]; then
    init_pipe="$ENSEMBL_BASE/ensembl-hive/scripts/$init_pipe"
fi
if [ $INIT_PIPE -eq 1 ]; then
  perl $init_pipe Bio::EnsEMBL::Analysis::Hive::Config::UniProtDB_conf -host ${HIVE_HOST} -port ${HIVE_PORT} -user ${HIVE_USER} -password ${HIVE_PASS} -dbname ${HIVE_DBNAME} -pipeline_name ${pipeline_name} -kill -kill_list_host "$KILL_HOST" -kill_list_port "$KILL_PORT" -kill_list_user "$KILL_USER" -kill_list_pass "$KILL_PASS" -kill_list_dbname "$KILL_DBNAME"
  if [ $? -ne 0 ]; then
    echo "The init script failed"
    exit 1;
  fi
fi

export EHIVE_URL=$EHIVE_DRIVER://$HIVE_USER:$HIVE_PASS@$HIVE_HOST:$HIVE_PORT/$HIVE_DBNAME
echo "Hive URL: $EHIVE_URL"

beekeeper="beekeeper.pl"
if [ -e "$ENSEMBL_BASE/ensembl-hive/scripts/$beekeeper" ]; then
    beekeeper="$ENSEMBL_BASE/ensembl-hive/scripts/$beekeeper"
fi
if [ $RUN_PIPE -eq 1 ]; then
  perl $beekeeper -url $EHIVE_URL -loop -can_respecialize 1
  if [ $? -ne 0 ]; then
    echo "Beekeeper failed"
    exit 1;
  fi
fi
