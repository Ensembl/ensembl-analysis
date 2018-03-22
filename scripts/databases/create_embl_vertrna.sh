#!/bin/sh
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

BASE_ENA_FTP="ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/std"
VERTRNADIR="$BLASTDB_DIR/vertrna"
VERTRNA_FILE="embl_vertrna-1"
BLAST_TYPE="ncbi"
ENSEMBL_BASE=$ENSCODE

USAGE=0
INIT_PIPE=1
RUN_PIPE=1
EHIVE_DRIVER="mysql" # This should not change unless you know what you are doing
VERTRNA_VERSION=`wget -nv --spider ${BASE_ENA_FTP%std}/doc/Release_* |& grep File | sed 's/File .Release_\([0-9]\+\).*/\1/'`

while getopts "f:s:b:v:d:e:h:u:p:P:IR" o; do
    case $o in
        f ) BASE_ENA_FTP=$OPTARG;;
        s ) SCRIPT=$OPTARG;;
        b ) BLAST_TYPE=$OPTARG;;
        v ) VERTRNA_VERSION=$OPTARG;;
        d ) DEST_DIR=$OPTARG;;
        e ) ENSEMBL_BASE=$OPTARG;;
        h ) EHIVE_HOST=$OPTARG;;
        u ) HIVEDB_USER=$OPTARG;;
        p ) EHIVE_PASS=$OPTARG;;
        P ) EHIVE_PORT=$OPTARG;;
        I ) RUN_PIPE=0;;
        R ) INIT_PIPE=0;;
        * ) USAGE=1;;
    esac
done

if [ $USAGE -eq 1 ];then
  echo <<EOF
    By default the script will initialise the pipeline and run the pipeline.
    Available options are:
       -f FTP URL for for the ENA website, default is $BASE_ENA_FTP
       -s Script to parse the EMBL file into FASTA, ensembl-analysis/scripts/databases/embl2fasta.pl
       -b Blast type, ncbi or wu, default is ncbi
       -v Version of ENA
       -d Output directory, you will need to set the version, deafult is $BLASTDB_DIR/vertrna/[version]
       -e Directory where ensembl-analysis is
       -I If you only want to initialise the pipeline
       -R If you only want ot run the pipeline
        * ) USAGE=1;;
EOF
fi

if [ -z $DEST_DIR ];then
  DEST_DIR=$VERTRNADIR/$VERTRNA_VERSION
fi

if [ -n $SCRIPT ]; then
  if [ -e $SCRIPT ]; then
    export EMBL2FASTA_SCRIPT=$SCRIPT
  else
    echo "'$SCRIPT' does not exist"
    exit 1
  fi
fi

if [ -z $EHIVE_HOST ];then
  echo "Missing option -h"
  exit 1;
else
  export EHIVE_HOST=$EHIVE_HOST
fi

if [ -z $HIVEDB_USER ];then
  echo "Missing option -u"
  exit 1;
fi

if [ -z $EHIVE_PASS ];then
  echo "Missing option -p"
  exit 1;
else
  export EHIVE_PASS=$EHIVE_PASS
fi

if [ -z $EHIVE_PORT ];then
  echo "Missing option -P"
  exit 1;
else
  export EHIVE_PORT=$EHIVE_PORT
fi

PIPELINE_NAME="vertrna_$VERTRNA_VERSION"
if [ $INIT_PIPE -eq 1 ];then
  perl $ENSEMBL_BASE/ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Analysis::Hive::Config::VertRNA_conf -pipeline_name $PIPELINE_NAME -vertrna_dir $DEST_DIR -blast_type $BLAST_TYPE -vertrna_file $VERTRNA_FILE -user $HIVEDB_USER -vertrna_version $VERTRNA_VERSION -hive_force_init 1
  if [ $? -ne 0 ];then
    echo "Something wrong happened while initialising the pipeline"
    exit 2
  fi
fi

export EHIVE_URL=$EHIVE_DRIVER://$HIVEDB_USER:$EHIVE_PASS@$EHIVE_HOST:$EHIVE_PORT/${USER}_$PIPELINE_NAME

if [ $RUN_PIPE -eq 1 ];then
  perl $ENSEMBL_BASE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -sync
  if [ $? -ne 0 ];then
    echo "Something wrong happened while sync'ing the pipeline"
    exit 3
  fi
  perl $ENSEMBL_BASE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop -can_respecialize 1
  if [ $? -ne 0 ];then
    echo "Something wrong happened while running the pipeline"
    exit 4
  fi
fi
