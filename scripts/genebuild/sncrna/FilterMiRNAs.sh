#!/bin/sh

usage(){
	echo -e "\n[ GB MiRNA pipeline Filter ]\n\n

	DESCRIPTION:
		Given a list of coordinates for DNA aligned features mapping to putative stem-loop structures,
		this script derives additional metrics and applies additional filters using user-supplied 
		random forest classifier. 

	INPUT FILES:
		List of DAFs in BED format
		List of repeat features in BED format
		Genomic reference is FASTA format
		Pre-built RFC model and scaler

	PARAMETERS:

	Mandatory:
		-d DAFS in BED format
		-r repeats in BED format
		-g genome reference in FASTA format
		-m rfc model pickle
		-s pre-processing scaler pickle
		-w working_dir
		-c mysql connection string (eg. mysql-ens-genebuild-prod-1:\$port:\$dbname:\$user:\$password)

	Optional:
		-D turn on delete identified stem-loops from DB flag

	sh FilterMiRNAs.sh <options> ";
	exit 1
}

CODEBASE=$(dirname "$0")
DELETE_FLAG=false
DCK=false
RCK=false
GCK=false
MCK=false
SCK=false
WCK=false

while getopts ":w:d:r:g:m:s:c:G" opt; do
        case $opt in
                d)      if [ ! -e "$OPTARG" ]; then
            								echo "(E) You have not provided a valid file containing coordinates of DNA align features, exiting"
						            		exit 1
						            fi
						            DCK=true; DAFS="$OPTARG";;						
                        
                r)      if [ ! -e "$OPTARG" ]; then
								            echo "(E) You have not provided a valid file containing coordinates of repeats, exiting"
								            exit 1
						            fi
                        
						            RCK=true; REPEATS="$OPTARG";;						
                g)      if [ ! -e "$OPTARG" ]; then
								            echo "(E) You have not provided a valid file containing genome sequence in FASTA format, exiting"
								            exit 1
						            fi
                        
						            GCK=true; GENOME_FASTA="$OPTARG";;						
                m)      if [ ! -e "$OPTARG" ]; then
								            echo "(E) You have not provided a valid path for the RFC model, exiting"
								            exit 1
						            fi
                        
						            MCK=true; MODEL="$OPTARG";;
                s)      if [ ! -e "$OPTARG" ]; then
								            echo "(E) You have not provided a valid path for the pre-processing scaler, exiting"
								            exit 1
						            fi
                        
						            SCK=true; SCALER="$OPTARG";;
				        c)      DELETE_FLAG=true; DBSTRING="$OPTARG";;
                w)      if [ ! -d "$OPTARG" ]; then
                            echo "(I) The path provided is not a directory, using current directory"
                            WORKING_DIR=`pwd`
                        fi
                        
                        WCK=true; WORKING_DIR="$OPTARG";;                        
                G)      DELETE_FLAG=true;; 
                h)      usage;;
                \?)
                   echo -e "Invalid option: -$OPTARG\n"
                         usage
                         ;;
                :)
                   echo "Option -$OPTARG requires an argument!"
                        usage
                        exit 1
                         ;;
        esac
done

if ! $DCK || ! $RCK || ! $WCK || ! $SCK || ! $MCK || ! $GCK
then
	echo "(E) Some mandatory input parameters are missing; please check parameters and try again..exiting!"
	usage
	exit 1
fi


# derive additional metrics required by the model
annotateBed -i $DAFS -files $REPEATS | \
	bedtools nuc -fi $GENOME_FASTA -bed stdin -s | \
		cut -f1-11,13 | \
		grep -v '#' > ${WORKING_DIR}/annotated_dafs.tsv 

# check for RNAFold results
if [ ! -e ${WORKING_DIR}/rna_fold_results.txt ]
then
  echo "(E) Missing ${WORKING_DIR}/rna_fold_results.txt, please generate to continue, exiting! "
fi 

if [ ! -e ${WORKING_DIR}/identified_mirnas.bed ]
then
    echo "(E) Missing ${WORKING_DIR}/identified_mirnas.bed, please dump identified stem-loops from DB to continue, exiting! "
fi

# run rfc model
python ${CODEBASE}/FilterDafs.py $MODEL $SCALER $WORKING_DIR

if $DELETE_FLAG
then
	IFS=: read HOST PORT DBNAME DBUSER DBPASSWORD <<< "$DBSTRING"

	perl ${CODEBASE}/delete_genes.pl \
		-dbhost $HOST \
		-dbport $PORT \
		-dbuser $DBUSER \
		-dbpass $DBPASSWORD \
		-dbname $DBNAME \
		-idfile ${WORKING_DIR}/mirnas_to_delete.txt
fi

