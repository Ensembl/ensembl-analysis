#!/bin/sh

REPEATS=`[ -e $1 ]; echo $?`
GENOME_FASTA=`[ -e $2 ]; echo $?`
DBSTRING=$3 # expects string like: host;dbname;port;user eg. ensembldb.ensembl.org;homo_sapiens_core_90_38;3306;anonymous
WORKING_DIR=$4

mkdir -p ${WORKING_DIR}/data_dir ${WORKING_DIR}/output_dir ${WORKING_DIR}/models

OLDIFS=$IFS
IFS=: read DBHOST DBNAME DBPORT DBUSER <<< "$DBSTRING"
IFS=$'\n'

TRANSCRIPT_BED=${WORKING_DIR}/data_dir/repeats.bed

# prepare annotations
if (( $(echo $REPEATS '>' 0 ) ))
then
        mysql -u $DBUSER -h $DBHOST -P $DBPORT -D $DBNAME -NB -e \
                "select concat('chr', s.name) as chrom, r.seq_region_start,
                        r.seq_region_end, r.seq_region_strand from repeat_feature r,
                        seq_region s where r.seq_region_id = s.seq_region_id ;" > $TRANSCRIPT_BED
else
        cat $1 > $TRANSCRIPT_BED
fi

if (( $(echo $GENOME_FASTA '>' 0 ) ))
then
    perl ${ENSCODE}/ensembl-analysis/scripts/sequence_dump.pl \
          -dbhost $DBHOST \
          -dbuser $DBUSER \
          -dbname $DBNAME \
          -dbport $DBPORT \
          -coord_system_name toplevel \
          -output_dir $WORKING_DIR/data_dir \
          -mask -mask_repeat RepeatMask \
          -filename  $WORKING_DIR/data_dir/genome.fasta \
          -header rnaseq \
          -softmask -onefile -padded_human_Chr_Y -species homo_sapiens

        #cat ${WORKING_DIR}/*.fa > $WORKING_DIR/genome.fasta
        sed -i 's/>/>chr/g' $WORKING_DIR/data_dir/genome.fasta
else
        cp $2 $WORKING_DIR/data_dir/genome.fasta
fi
