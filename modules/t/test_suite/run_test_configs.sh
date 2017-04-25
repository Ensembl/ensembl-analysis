#!/usr/bin/env bash
## TODO: proper getopts
## Variables that can be modified
BIOPERL=$BIOPERL_LIB
NCBI_BLASTFMT="makeblastdb"
WU_BLASTFMT="xdformat"
BLAST_TYPE="ncbi"
INDICATE="indicate"
## Variables that should not be modified
ENA_BASE="http://www.ebi.ac.uk/ena/data/view/%s&display=fasta"
UNIPROT_BASE="http://www.uniprot.org/uniprot/%s.fasta"
MOUSEHOST='ensembldb.ensembl.org'
MOUSEDBNAME='mus_musculus_core_87_38'
MOUSEDBPORT=3306
MOUSEUSER="anonymous"

## Checking the variables to be able to connect
if [ -z "$EHIVE_HOST" ]; then
  echo "You are missing \$EHIVE_HOST"
  exit 1
fi
if [ -z "$EHIVE_PORT" ]; then
  echo "You are missing \$EHIVE_PORT"
  exit 1
fi
if [ -z "$EHIVE_USER" ]; then
  echo "You are missing \$EHIVE_USER"
  exit 1
fi
if [ -z "$EHIVE_ROUSER" ]; then
  echo "You are missing \$EHIVE_ROUSER"
  exit 1
fi
if [ -z "$EHIVE_PASS" ]; then
  echo "You are missing \$EHIVE_PASS"
  exit 1
fi

## Deleting created files
if [ -e "runtimestamps.txt" ]; then
  rm runtimestamps.txt
fi
if [ -d "rnaseq_test_out" ]; then
  rm -r rnaseq_test_out
fi
mkdir rnaseq_test_out

## Creating PERL5LIB and adding hive scripts to PATH
LOCALPERL5LIB="$PERL5LIB"
PERL5LIB=$PWD/../..
BACKDOTS="../../.."
if [ -d "$PWD/../../../../ensembl" ];then
  BACKDOTS="../../../.."
fi
BASEP5=$PWD/$BACKDOTS
export PERL5LIB=$PERL5LIB:$BASEP5/ensembl/modules:$BASEP5/ensembl-hive/modules:$BASEP5/ensembl-io/modules:$BASEP5/ensembl-pipeline/modules:$BASEP5/ensembl-compara/modules:$BASEP5/ensembl-killlist/modules:$BIOPERL
export PATH=$PATH:$BASEP5/ensembl-hive/scripts

## Functions
format_blast_db() {
# $1 should be prot or nucl
# $2 is the filename
  if [ "$BLAST_TYPE" == "ncbi" ];then
    $NCBI_BLASTFMT -dbtype $1 -parse_seqids -in $2 > /dev/null
    if [ "$?" -ne 0 ]; then
      echo "Something went wrong while formatting $2 with $NCBI_BLASTFMT"
      exit 3
    fi
  else
    if [ "$1" == "prot" ];then
      TYPE="T"
    else
      TYPE="F"
    fi
    $WU_BLASTFMT -p $TYPE $2 > /dev/null
    if [ "$?" -ne 0 ]; then
      echo "Something went wrong while formatting $2 with $WU_BLASTFMT"
      exit 3
    fi
  fi
}

database_checks() {
# $1 is the number of the exercise
  if [ $1 -eq 1 ];then
    DATA_DBNAME="test_workshop_rat_core"
    QUERIES=(
      "SELECT COUNT(*) FROM repeat_feature rf, analysis a WHERE rf.analysis_id = a.analysis_id AND a.logic_name = 'repeatmasker_repbase_rat'"
      "SELECT COUNT(*) FROM repeat_feature rf, analysis a WHERE rf.analysis_id = a.analysis_id AND a.logic_name = 'dust'"
      "SELECT COUNT(*) FROM repeat_feature rf, analysis a WHERE rf.analysis_id = a.analysis_id AND a.logic_name = 'trf'"
      "SELECT COUNT(*) FROM simple_feature sf, analysis a WHERE sf.analysis_id = a.analysis_id AND a.logic_name = 'eponine'"
      "SELECT COUNT(*) FROM simple_feature sf, analysis a WHERE sf.analysis_id = a.analysis_id AND a.logic_name = 'cpg'"
      "SELECT COUNT(*) FROM simple_feature sf, analysis a WHERE sf.analysis_id = a.analysis_id AND a.logic_name = 'trnascan'"
      "SELECT COUNT(*) FROM prediction_transcript"
      "SELECT COUNT(*) FROM prediction_exon"
    )
    RESULTS=(
      447
      851
      457
      35
      5
      53
      26
      240
    )
  elif [ $1 -eq 2 ];then
    DATA_DBNAME="test_workshop_rat_uniprot"
    QUERIES=(
      "SELECT COUNT(*) FROM gene"
      "SELECT COUNT(*) FROM transcript"
      "SELECT COUNT(*) FROM exon"
    )
    RESULTS=(
      7
      7
      48
    )
  elif [ $1 -eq 3 ];then
# Because I cannota get the correct data, it will wait for a test
#    DATA_DBNAME="test_workshop_rat_projection"
#    QUERIES=(
#      "SELECT COUNT(*) FROM gene"
#    )
#    RESULTS=(
#      1755
#    )
#  elif [ $1 -eq 4 ];then
    DATA_DBNAME="test_workshop_rat_layer"
    QUERIES=(
      "SELECT COUNT(*) FROM gene"
    )
    RESULTS=(
      1755
    )
  elif [ $1 -eq 5 ];then
    DATA_DBNAME="test_workshop_rat_genebuilder"
    QUERIES=(
      "SELECT COUNT(*) FROM gene"
    )
    RESULTS=(
      1755
    )
  fi

  EXIT_CODE=0
  for I in `seq 0 $((${#QUERIES[*]}-1))`; do
    RES=`mysql -u$EHIVE_ROUSER -h$EHIVE_HOST ${EHIVE_ROPASS:+-p$EHIVE_ROPASS} -P $EHIVE_PORT -D$DATA_DBNAME -NB -e "${QUERIES[$I]}"`
    if [ "$RES" != "${RESULTS[$I]}" ]; then
      printf "Test %d Query %d: '%s' \033[41mfailed\033[0m\n\t%s instead of %s\n" $1 $((I+1)) "${QUERIES[$I]}" "$RES" "${RESULTS[$I]}"
      EXIT_CODE=1
    fi
  done
  if [ $EXIT_CODE -eq 1 ]; then
    exit 1
  fi
}

## Preparing the data
## EX1 data
echo "Getting data for exercise 1"
CONTIG_FILE="data/assembly/contigs.fa"
if [ -e "$CONTIG_FILE" ]; then
  rm $CONTIG_FILE
fi

ACCESSIONS=(
  "AC133299.3"
  "AC106687.6"
  "AC126894.4"
  "AABR07073112.1"
);
for A in ${ACCESSIONS[*]}; do
  URL=`printf $ENA_BASE "$A"`
  wget -qq -O - "$URL" >> $CONTIG_FILE
  grep "$A" $CONTIG_FILE > /dev/null
  if [ "$?" -ne 0 ]; then
    echo "Something went wrong while downloading $A"
    exit 1
  fi
done

##EX2 data
echo "Getting data for exercise 2"
PROTEINS_FILE="data/uniprot_proteins/uniprot_proteins.fa"
if [ -e "$PROTEINS_FILE" ]; then
  rm $PROTEINS_FILE
fi

ACCESSIONS=(
"G3V7I5"
"B0BN16"
"Q4QR84"
"Q8WX77"
"H9GIP6"
"Q9P1Z9"
"H9GKT7"
"P28289"
"H9G4T0"
);
for A in ${ACCESSIONS[*]}; do
  URL=`printf $UNIPROT_BASE "$A"`
  wget -qq -O - "$URL" | sed '/^>/ s/^>[sptr]\{2\}|\([0-9A-Z]\+\)|[A-Z0-9]\+_\([A-Z]\+\).*SV=\([0-9]\+\)/>\1.\3 \2/' >> $PROTEINS_FILE
  grep "$A" $PROTEINS_FILE > /dev/null
  if [ "$?" -ne 0 ]; then
    echo "Something went wrong while downloading $A"
    exit 1
  fi
done
format_blast_db "prot" $PROTEINS_FILE

###EX3 data
#echo "Getting data for exercise 3"
#STABLEIDFILE="data/assembly/mouse_stable_ids.lst"
#TARGET_DBNAME="test_wrokshop_mouse_core"
#TEMPDATA=`mktemp`
#ACCESSIONS='"AC087559.69"'
#ksh scripts/clone_database.ksh -l -s "${MOUSEDBNAME}@${MOUSEHOST}:${MOUSEDBPORT}" -t "${TARGET_DBNAME}@${EHIVE_HOST}:${EHIVE_PORT}" -r ${MOUSEUSER} -w $EHIVE_USER -p${EHIVE_PASS}
#mysql -h $MOUSEHOST -u $MOUSEUSER -P $MOUSEDBPORT $MOUSEDBNAME -NB -e "SELECT * FROM dna d, seq_region sr WHERE d.seq_region_id = sr.seq_region_id AND sr.name IN (${ACCESSIONS})" > ${TEMPDATA}
#mysql -h $EHIVE_HOST -u $EHIVE_USER -p$EHIVE_PASS -P $EHIVE_PORT -e "LOAD DATA INFILE '$TEMPDATA' INTO TABLE dna"
#perl scripts/genebuild/copy_genes.pl --sourcehost $MOUSEHOST -sourceport $MOUSEDBPORT --sourceuser $MOUSEUSER --sourcedbname $MOUSEDBNAME --dnahost $MOUSEHOST --dnauser $MOUSEUSER --dnaport $MOUSEDBPORT --dnadbname $MOUSEDBNAME --targethost $EHIVE_HOST --targetuser $EHIVE_USER --targetpass $EHIVE_PASS --targetdbname ${TARGET_DBNAME} --file $STABLEIDFILE --stable_id

##EX4 data
echo "Getting data for exercise 4"
PROTEINS_FILE="data/rna_seq/blast_database/uniprot_blastdb.fa"
if [ -e "$PROTEINS_FILE" ]; then
  rm $PROTEINS_FILE
fi

ACCESSIONS=(
"D4ADA5"
"G3V7I5"
"B0BN16"
"Q4QR84"
"F1LXD9"
"Q9R1R4"
"Q6IMZ5"
"Q56A27"
);
for A in ${ACCESSIONS[*]}; do
  URL=`printf $UNIPROT_BASE "$A"`
  wget -qq -O - "$URL" | sed '/^>/ s/^>[sptr]\{2\}|\([0-9A-Z]\+\).*SV=\([0-9]\+\)/>\1.\2/' >> $PROTEINS_FILE
  grep "$A" $PROTEINS_FILE > /dev/null
  if [ "$?" -ne 0 ]; then
    echo "Something went wrong while downloading $A"
    exit 1
  fi
done
format_blast_db "prot" $PROTEINS_FILE
$INDICATE -d "`dirname $PROTEINS_FILE`" -f "`basename $PROTEINS_FILE`" -i "`dirname $PROTEINS_FILE`/index" -p singleWordParser > /dev/null
if [ "$?" -ne 0 ]; then
  echo "Something went wrong while indexing $PROTEINS_FILE with $INDICATE"
  exit 2
fi

echo "Start:" > runtimestamps.txt
date >> runtimestamps.txt
CONFIGS=(
  "GenomePreparation_conf.pm"
  "UniProtAlignment_conf.pm"
  "RNASeq_conf.pm"
  )
#  "ProjectTranscripts_conf.pm" Ex3
#  "FinalGeneSet_conf.pm"

for I in `seq 1 ${#CONFIGS[*]}`; do
  init_pipeline.pl "${CONFIGS[$((I-1))]}" -pipeline_name test_workshop_ex${I} -pipe_dbname test_workshop_ex${I}_pipe -dna_dbname test_workshop_rat_core -hive_force_init 1 -drop_databases 1
  if [ "$?" -gt 0 ]; then
    echo "Problem while initialising exercise $I"
    exit 2
  fi
  beekeeper.pl -url mysql://$EHIVE_USER:$EHIVE_PASS@$EHIVE_HOST:$EHIVE_PORT/test_workshop_ex${I}_pipe -sync
  if [ "$?" -gt 0 ]; then
    echo "Problem while synchronising exercise $I"
    exit 3
  fi
  beekeeper.pl -url mysql://$EHIVE_USER:$EHIVE_PASS@$EHIVE_HOST:$EHIVE_PORT/test_workshop_ex${I}_pipe -loop -sleep 0.2 -local -can_respecialize 1
  if [ "$?" -gt 0 ]; then
    echo "Problem while running exercise $I"
    exit 4
  fi
  echo "Pipeline $I end:" >> runtimestamps.txt
  date >> runtimestamps.txt

  database_checks $I
done

# Settings variable back to normal
export PERL5LIB=$LOCALPERL5LIB
export PATH=${PATH%$BASEP5/ensembl-hive/scripts}
