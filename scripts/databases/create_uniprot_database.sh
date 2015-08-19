#!/bin/sh

# Script to download the Uniprot db, classification info and create BLAST dbs
# off PE level and classification
#
# The first section of the code will set all the files downloading
# The second section processes taxonomic divisions
# The third section sets up BLAST dbs based off PE levels and taxonomic divisions
# The final section processes the main Uniprot db files, these usually take longest
# to download.
#
# The code tries to progress in the most logical order of file processing in terms
# of the amount of time taken to download each file
#
# IF SOMETHING BREAKS:
# In this case the first thing to do is to check what section it broke in
# The most important thing is to check that all the files had downloaded successfully
# If they didn't you can put in an if statement to skip over the ones that did complete
# in the downloads section.
# If they did all download successfully you can just comment out/delete the code up
# to the relevant section.
# There are some conditionals in the code that check if various job names have completed
# successfully. If you've ran the jobs independently of the script without using the
# correct job name, then you can just comment out the conditional if you know the files
# are ready
#
# NOTE: all PE3 classification downloads have been disabled as the query has stopped
# completing successfully. All the code is still present, you need to use -t


CURRENT=$PWD
DESTDIR=$PWD
EMB2FASTA_SCRIPT="`dirname $0`/embl2fasta.pl"

BASEFTP="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"
TAXFTP="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/"

DOWNLOAD=0
UNZIP=0
MAIN=0
INDICATE=0
PELEVEL=0
PELEVEL3=0
FULL_ENTRY_LOC=0

while getopts ":o:s:dumiptlha" o; do
    case $o in
        o ) DESTDIR=$OPTARG;;
        a ) AUTO=1;;
        d ) DOWNLOAD=1;;
        u ) UNZIP=1;;
        m ) MAIN=1;;
        i ) INDICATE=1;;
        p ) PELEVEL=1;;
        t ) PELEVEL3=1;;
        l ) FULL_ENTRY_LOC=1;;
        s ) SCRIPTS=$OPTARG;;
        * ) USAGE=1;;
    esac
done

if [ -n "$USAGE" ];then
    echo "You need to specify at least an output directory using -o"
    echo "Optional options are:"
    echo "  -d Download the main files"
    echo "  -a Fetch the UniProt version, no need to add uniprot_YYYY_XX to your path, just the root"
    echo "  -u Unzip the downloaded files"
    echo "  -m Process the main database for blast"
    echo "  -i Index the vertebrate proteins with indicate"
    echo "  -p Process the PE12 level files"
    echo "  -s Specify the path to the script embl2fasta.pl"
    echo "  -t Process PE3 files, DO NOT USE at the moment"
    exit 0
fi

if [ "$DOWNLOAD" -eq 0 ] && [ "$UNZIP" -eq 0 ] && [ "$MAIN" -eq 0 ] && [ "$INDICATE" -eq 0 ] && [ "$PELEVEL" -eq 0 ]; then
    DOWNLOAD=1
    UNZIP=1
    MAIN=1
    INDICATE=1
    PELEVEL=1
fi

DIE=0
echo "Checking your binaries..."
for B in "/software/ensembl/bin/indicate" "xdformat"
    do
        which $B > /dev/null
        if [ "$?" -eq 1 ];then
            DIE=1
            echo "  $B is not in your PATH"
        fi
done
if [ ! -e "$EMB2FASTA_SCRIPT" ];then
    DIE=1
    echo "Could not find $EMB2FASTA_SCRIPT"
fi
if [ "$DIE" -eq 1 ];then
    echo "Failed"
    exit 11
else
    echo "Done"
fi

UNIPROT_VERSION=`wget -S --spider www.uniprot.org 2>&1 | grep 'X-UniProt-Release' | awk '{print $2}'`
if [ "$AUTO" -eq 1 ]; then
  DESTDIR="$DESTDIR/uniprot_$UNIPROT_VERSION"
fi

if [ ! -e "$DESTDIR" ]; then
    mkdir -m 775 $DESTDIR
    if [ "$?" -ne 0  ]; then
        echo "Could not create $DESTDIR"
        exit 10
    fi
elif [ "$AUTO" -eq 1 ];then
    printf "The directory already exists and you are using -a:\n - this database has already been successfully created\n - delete the directory %s and retry\n - rerun the script without -a\n" "$DESTDIR"
    exit 0
fi

UNIPROT_DATE=`wget -S --spider www.uniprot.org 2>&1 | grep 'Last-Modified' | sed 's/\s*Last-Modified:\s\+//'`
COUNT=0


FASTA=".fasta"
DAT=".dat"
GZ=".gz"

tax_files=(
            'uniprot_sprot_human'
            'uniprot_sprot_mammals'
            'uniprot_sprot_rodents'
            'uniprot_sprot_vertebrates'
            'uniprot_trembl_human'
            'uniprot_trembl_mammals'
            'uniprot_trembl_rodents'
            'uniprot_trembl_vertebrates'
            'uniprot_trembl_unclassified'
          )

main_files=(
             'uniprot_sprot.fasta'
             'uniprot_trembl.fasta'
           )

# If you had a new file to the pe levels, you may need to add it here to
pelevel_files=(
                'uniprot_PE12_vertebrata_initial'
                'uniprot_PE12_vertebrata_frag_initial'
                'uniprot_PE12_nonvert_initial'
                'uniprot_PE12_nonvert_frag_initial'
              )

pelevel3_files=(
        'uniprot_PE3_nonvert_initial'
        'uniprot_PE3_nonvert_frag_initial'
        'uniprot_PE3_vertebrata_initial'
        'uniprot_PE3_vertebrata_frag_initial'
               )

# Extensions to make job names more distinct
exten1="_updl"
exten2="_gzip"
exten3="_conv"

# Download log directory:
mkdir -p $DESTDIR/download_logs
mkdir -p $DESTDIR/gzip_logs
mkdir -p $DESTDIR/convert_logs

# Now let's go in the DESTDIR to make things simple
exit_clean () {
    cd $CURRENT
    exit $1
}

cd $DESTDIR
#######################################################################
#
# Begin file download code
#
#######################################################################

# This section of the code will set the following downloading:
#
# From: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete
# uniprot_sprot.fasta.gz
# uniprot_trembl.fasta.gz
#
# From: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
# uniprot_sprot_human.dat.gz
# uniprot_sprot_mammals.dat
# uniprot_sprot_rodents.dat
# uniprot_sprot_vertebrates.dat
# uniprot_trembl_human.dat
# uniprot_trembl_mammals.dat
# uniprot_trembl_rodents.dat
# uniprot_trembl_unclassified.dat
# uniprot_trembl_vertebrates.dat
#
# From direct web query of uniprot.org
# uniprot_PE12_nonvert
# uniprot_PE12_nonvert_frag
# uniprot_PE12_vertebrata
# uniprot_PE12_vertebrata_frag
# uniprot_PE3_nonvert
# uniprot_PE3_nonvert_frag
# uniprot_PE3_vertebrata
# uniprot_PE3_vertebrata_frag
#
# Note: There is obviously redundancy here. The main example is between the
# taxonomic division files and the main db files. The sensible thing to do
# would to have simply download all taxonomic division files, not just vert
# and then there would be no need for the main db files. However the problem
# is that the taxonomic division files are only in embl format at the moment
# and the nonvert ones are far too big to download. So presently it is fastest
# to download the vert tax files in embl format and convert them, so that we
# can run indicate on them while downloading the main files in fasta format.
# If we fix indicate then the taxonomic division files and step will be
# removed altogether.
#
# Requirements for executing this code:
# None
#
# Requirements for exiting this code:
# None
#
# Outcomes of executing this code:
# All listed files are bsubbed for download with jobnames used in later checkpoints


if [ "$DOWNLOAD" -eq 1 ]; then
# Main uniprot files
    for F in ${main_files[@]};
        do
            jobname=$F$exten1
            echo "Downloading $F$GZ..."
            bsub -q long -o $DESTDIR/download_logs/$F.dl.out -e $DESTDIR/download_logs/$F.dl.err -J "$jobname" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' wget -O "$DESTDIR/$F$GZ" -q "$BASEFTP/$F$GZ"
        done

# Files for taxonomic divisions excluding nonvert
    for F in ${tax_files[@]};
        do
            jobname=$F$DAT$exten1
            echo "Downloading $F$DAT..."
            bsub -o $DESTDIR/download_logs/$F.dl.out -e $DESTDIR/download_logs/$F.dl.err -J "$jobname" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' wget -O "$DESTDIR/$F$DAT$GZ" -q "$TAXFTP/$F$DAT$GZ"
        done


# Download PE12 vert/nonvert/frag/nonfrag if you add or change a file here be sure that it is properly naemed in pelevel_files
    echo "Downloading PE12 vert and nonvert, with frag and nonfrag"
    bsub -o $DESTDIR/download_logs/PE12_vert.out -e $DESTDIR/download_logs/PE12_vert.err -J "PE12_vert_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+AND+taxonomy%3aCraniata+AND+fragment:no&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE12_vertebrata_initial$GZ"
    bsub -o $DESTDIR/download_logs/PE12_vert_frag.out -e $DESTDIR/download_logs/PE12_vert_frag.err -J "PE12_vert_frag_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+AND+taxonomy%3aCraniata+AND+fragment:yes&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE12_vertebrata_frag_initial$GZ"
    bsub -o $DESTDIR/download_logs/PE12_nonvert.out -e $DESTDIR/download_logs/PE12_nonvert.err -J "PE12_nonvert_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+NOT+taxonomy%3aCraniata+AND+fragment%3ano&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE12_nonvert_initial$GZ"
    bsub -o $DESTDIR/download_logs/PE12_nonvert_frag.out -e $DESTDIR/download_logs/PE12_nonvert_frag.err -J "PE12_nonvert_frag_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+NOT+taxonomy%3aCraniata+AND+fragment%3ayes&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE12_nonvert_frag_initial$GZ"

# NOTE: I'm disaabling this section for the moment as the query for PE3 nonvert
# usually fails to download correctly. It is not essential.
# Download PE3 vert/nonvert/frag/nonfrag
    if [ "$PELEVEL3" -eq 1 ];then
        echo "Downloading PE3 vert and nonvert, with frag and nonfrag"
        bsub -o $DESTDIR/download_logs/PE3_vert.out -e $DESTDIR/download_logs/PE3_vert.err -J "PE3_vert_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22inferred+from+homology%22)+AND+taxonomy%3aCraniata+AND+fragment:no&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE3_vertebrata_initial$GZ"
        bsub -o $DESTDIR/download_logs/PE3_vert_frag.out -e $DESTDIR/download_logs/PE3_vert_frag.err -J "PE3_vert_frag_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22inferred+from+homology%22)+AND+taxonomy%3aCraniata+AND+fragment:yes&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE3_vertebrata_frag_initial$GZ"
        bsub -o $DESTDIR/download_logs/PE3_nonvert.out -e $DESTDIR/download_logs/PE3_nonvert.err -J "PE3_nonvert_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22inferred+from+homology%22)+NOT+taxonomy%3aCraniata+AND+fragment%3ano&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE3_nonvert_initial$GZ"
        bsub -o $DESTDIR/download_logs/PE3_nonvert_frag.out -e $DESTDIR/download_logs/PE3_nonvert_frag.err -J "PE3_nonvert_frag_updl" -G ensembl-genebuild -M 500 -R 'select[mem>500] rusage[mem=500]' "wget -q -O - \"http://www.uniprot.org/uniprot/?query=(existence%3a%22inferred+from+homology%22)+NOT+taxonomy%3aCraniata+AND+fragment%3ayes&compress=yes&format=fasta\" > $DESTDIR/uniprot_PE3_nonvert_frag_initial$GZ"
    fi

fi

#######################################################################
#
# End file download code
#
#######################################################################


#######################################################################
#
# Begin taxonomic division processing and indicate code
#
#######################################################################


# This section of the code bsub jobs involved in the processing of the taxonomic
# division files. For each taxonomic division file in the previous section two
# jobs are bsubbed. The first unzips the dat file, with the condition that the
# dat file has finished downloading. The second converts the dat file from embl
# to fasta format, with the condition that the dat file has been unzipped.
# The code will wait after bsubbing these jobs to ensure all the taxonomic division
# files have finished downloading. Once they have it will concat the final set of
# fasta files into a single file called uniprot_all_vert_indicate.
# It then runs indicate on this file, with the results placed in uniprot_index
#
# Requirements for executing this code:
# None
#
# Requirements for exiting this code:
# All taxonomic division files must have been downloaded, processed and indicate
# must have been run successfully
#
# Outcomes of executing this code:
# When finished indicate will have been run so that we can classify all vertebrate
# proteins during the genebuild. All processing of taxonomic files will be complete.

if [ "$UNZIP" -eq 1 ];then

# Get taxonomic division files in EMBL format:
    for F in ${tax_files[@]};
        do
            echo "Processing taxomonic division $F$DAT..."

            if [ -e "$DESTDIR/$F$DAT" ];then
                echo "$F$DAT is already unzipped"
                OPTIONS_LSF_C2=""
            else
                if [ "$DOWNLOAD" -eq 1 ];then
                    OPTIONS_LSF_C1="-w done($F$DAT$exten1)"
                else
                    OPTIONS_LSF_C1=""
                fi
                jobname=$F$DAT$exten2
                OPTIONS_LSF_C2="-w done($F$DAT$exten2)"
                bsub -o $DESTDIR/gzip_logs/$F.unzip.out -e $DESTDIR/gzip_logs/$F.unzip.err -J "$jobname" $OPTIONS_LSF_C1 -M 500 -R 'select[mem>500] rusage[mem=500]' "gunzip $DESTDIR/$F$DAT$GZ"
            fi

            if [ -e "$DESTDIR/$F$FASTA" ];then
                echo "$DESTDIR/$F$DAT already converted"
            else
                jobname=$F$DAT$exten3
                bsub -o $DESTDIR/convert_logs/$F.conv.out -e $DESTDIR/convert_logs/$F.conv.err -J "$jobname" $OPTIONS_LSF_C2 -M 500 -R 'select[mem>500] rusage[mem=500]' "perl $EMB2FASTA_SCRIPT $DESTDIR/$F$DAT"
            fi
        done


# Concat the fasta files into one for indicate
# At this point all convert jobs must be done

# I've decided on an interactive bsub at this point as a kind of checkpoint.
# I could just keep bsubbing with dependencies, but that will get very hard to debug
# So nothing moves forward here till the converts are done on the tax divisions

# WARNING: if your /bin/sh is not bash, you may have

# Maybe a bit ugly but it removes the need of dependencies in the jobs
    IDX=$((${#tax_files[*]}-1))
    ITERATION=0

#    for I in ${tax_files[@]};do
    while [ "${#DONE[*]}" -ne "${#tax_files[@]}" ]; do
        if [ "$ITERATION" -eq 100 ];then
            break
        fi
        for I in `seq 0 $IDX`;do
            if [ -e "${tax_files[$I]}$FASTA" ];then
                if [ -e "${tax_files[$I]}$DAT" ];then
                    if [ -z "${HASH[$I]}" ];then
                        HASH[$I]=`grep -c '^SQ   SEQUENCE' ${tax_files[$I]}$DAT`
                    fi
                    FASTACOUNT=`grep -c '>' ${tax_files[$I]}$FASTA`;
                    if [ "$FASTACOUNT" -eq "${HASH[$I]}" ];then
                        DONE[$I]=1
                        echo "${tax_files[$I]} is good"
                    else
                        echo "Fasta: $FASTACOUNT <=> Dat: ${HASH[$I]}"
                    fi
                fi
            fi
        done
        ITERATION=$(($ITERATION+1))
        sleep 2m
    done

    bsub -I -J "all_tax_div_concat" -M 500 -R 'select[mem>500] rusage[mem=500]' "cat ${tax_files[@]/%/$FASTA} > $DESTDIR/uniprot_all_vert_indicate"

fi

if [ "$INDICATE" -eq 1 ];then
# Formatting for OBDA Index: Similarity,...
bsub -M 2000 -R "select[mem>2000] rusage[mem=2000]" -I "/software/ensembl/bin/indicate -d $DESTDIR -f uniprot_all_vert_indicate --index $DESTDIR/uniprot_index  --parser singleWordParser 2>&1 | grep Error; if [ \"\$?\" = 2 ]; then exit \$?; elif [ \"\$?\" = 1 ]; then exit 0; else exit 1; fi"

if [ "$?" -ne 0 ];then
    echo "Failed to format uniprot with indicate"
    exit_clean 7
fi

# Changing permissions again
chmod -R g+w $DESTDIR/uniprot_index
fi

COUNTV=`grep -c \> $DESTDIR/uniprot_all_vert_indicate`
echo "xdformatting uniprot_vertebrate..."
bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/uniprot_vertebrata -t uniprot_all_vert_indicate -v 'uniprot_$UNIPROT_VERSION' -d '$UNIPROT_DATE' $DESTDIR/uniprot_all_vert_indicate > $DESTDIR/xdformat.log"
if [ "$?" -ne 0 ];then
    echo "Failed to format uniprot with xdformat"
    exit_clean 8
fi

XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
if [ "$XCOUNT" -eq "$COUNTV" ];then
    echo "xdformat successful"
    rm $DESTDIR/xdformat.log
else
    echo "xdformat failed"
    exit_clean 9
fi
rm ${tax_files[@]/%/$FASTA}
if [ "$?" -ne 0 ];then
    echo "Failed to delete files for indicate: ${tax_files[@]/%/$FASTA}"
fi

rm ${tax_files[@]/%/$DAT}
if [ "$?" -ne 0 ];then
    echo "Failed to delete files for indicate: ${tax_files[@]/%/$DAT}"
fi

#######################################################################
#
# End taxonomic division processing and indicate code
#
#######################################################################


#######################################################################
#
# Begin classification processing code
#
#######################################################################

# At this point the code will wait until all the classifaction files that were set
# downloading through queries of uniprot.org have completed downloading. Once they
# have the code begins executing. The first thing is a count of the sequences each
# of the PE12 files. Then xdformat is run on each PE12 file to create a BLAST db
# and a check is done against the log file to ensure the sequence counts match.
# After this the process is repeated on PE3
#
# Requirements for executing this code:
# All classification files must have finished downloading
#
# Requirements for exiting this code:
# All xdformats of the files must complete successfully
#
# Outcomes of executing this code:
# Once this code exits all classifcation files will have been downloaded from the
# uniprot website through a web query. All files will have been turned into a BLAST
# by by running xdformat on them


# NOTE: This is the original bsub wait from when PE3 was also downloaded. This is
# disabled as PE3 is not currently downloaded. If the PE3 download is reinstated
# then the below bsub can replace the current one
#echo "Waiting for PE level divisions to download..."
#bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' -w "done("PE12_vert_updl") && done("PE12_vert_frag_updl") && done("PE12_nonvert_updl") && done("PE12_nonvert_frag_#updl") && done("PE3_vert_updl") && done("PE3_vert_frag_updl") && done("PE3_nonvert_updl") && done("PE3_nonvert_frag_updl")" echo "...PE level divisions downloaded"

# I don't understand why, but this failed on memory once, so I've put memory in.
# Seems ridiculous though
if [ "$PELEVEL" -eq 1 ];then
    if [ "$DOWNLOAD" -eq 1 ];then
        echo "Waiting for PE level divisions to download..."
        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' -w "done("PE12_vert_updl") && done("PE12_vert_frag_updl") && done("PE12_nonvert_updl") && done("PE12_nonvert_frag_updl")" echo "...PE level divisions downloaded"
    else
        PELEVEL_ERROR=0
        for F in ${pelevel_files[*]};do
            if [ -e "$F$GZ" ];then
                echo "$F$GZ has been downloaded"
            else
                echo "You're missing $F$GZ"
                PELEVEL_ERROR=1
            fi
        done
        if [ "$PELEVEL_ERROR" -eq 1 ];then
            exit_clean 1
        fi
    fi

    for F in ${pelevel_files[*]};do
        echo "Processing $F"
        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' gunzip "$F$GZ"
        if [ "$?" -ne 0 ];then
            echo "Failed to unzip $F$GZ. File may be corrupted due to early query termination. Exiting"
            exit_clean 1
        fi

        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "sed -r -e 's/^>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]).*/>\1.\4/' -e '/^[^>]/ s/O/K/g' $F > ${F%_initial}"
        if [ "$?" -ne 0 ]; then
            echo "Failed to parse $F. Uniprot accession format may have changed. Exiting"
            exit_clean 1
        fi

        echo "Successfully unzipped and parsed $F$GZ"
        rm $DESTDIR/$F
    done



# Run on PE12
    COUNTV=`grep -c \> $DESTDIR/uniprot_PE12_vertebrata`
    COUNTVF=`grep -c \> $DESTDIR/uniprot_PE12_vertebrata_frag`
    COUNTNV=`grep -c \> $DESTDIR/uniprot_PE12_nonvert`
    COUNTNVF=`grep -c \> $DESTDIR/uniprot_PE12_nonvert_frag`


    echo "xdformatting PE12..."
    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE12_vertebrata -t uniprot_PE12_vertebrata -v 'uniprot_$UNIPROT_VERSION' -d '$UNIPROT_DATE' $DESTDIR/uniprot_PE12_vertebrata > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$COUNTV" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE12_vertebrata_wfrag -t 'uniprot_PE12_vertebrata,uniprot_PE12_vertebrata_frag' -v 'uniprot_$UNIPROT_VERSION' -d '$UNIPROT_DATE' $DESTDIR/uniprot_PE12_vertebrata $DESTDIR/uniprot_PE12_vertebrata_frag > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$(($COUNTV+$COUNTVF))" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE12 -t 'uniprot_PE12_vertebrata,uniprot_PE12_nonvert' -v 'uniprot_$UNIPROT_VERSION' -d '$UNIPROT_DATE' $DESTDIR/uniprot_PE12_vertebrata $DESTDIR/uniprot_PE12_nonvert > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$(($COUNTV+$COUNTNV))" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE12_wfrag -t 'uniprot_PE12_vertebrata,uniprot_PE12_vertebrata_frag,uniprot_PE12_nonvert,uniprot_PE12_nonvert_frag' -v 'uniprot_$UNIPROT_VERSION' -d '$UNIPROT_DATE' $DESTDIR/uniprot_PE12_vertebrata $DESTDIR/uniprot_PE12_vertebrata_frag $DESTDIR/uniprot_PE12_nonvert $DESTDIR/uniprot_PE12_nonvert_frag > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$(($COUNTV+$COUNTVF+$COUNTNV+$COUNTNVF))" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi
fi

if [ "$PELEVEL3" -eq 1 ];then
# NOTE: all the unzips and conversions for PE3 are disables below
# Code for unzipping and then running the perl one liner to parse uniprot_PE3_nonvert_initial.gz
    if [ "$DOWNLOAD" -eq 1 ];then
        echo "Waiting for PE level divisions to download..."
        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' -w "done("PE3_vert_updl") && done("PE3_vert_frag_updl") && done("PE3_nonvert_updl") && done("PE3_nonvert_frag_updl")" echo "...PE level divisions downloaded"
    else
        PELEVEL3_ERROR=0
        for F in ${pelevel3_files[*]};do
            if [ -e "$F$GZ" ];then
                echo "$F$GZ has been downloaded"
            else
                echo "You're missing $F$GZ"
                PELEVEL3_ERROR=1
            fi
        done
        if [ "$PELEVEL3_ERROR" -eq 1 ];then
            exit_clean 1
        fi
    fi

    for F in ${pelevel3_files[*]};do
        echo "Processing $F"
        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' gunzip "$F$GZ"
        if [ "$?" -ne 0 ];then
            echo "Failed to unzip $F$GZ. File may be corrupted due to early query termination. Exiting"
            exit_clean 1
        fi

        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "sed -r -e 's/^>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]).*/>\1.\4/' -e '/^[^>]/ s/O/K/g' $F > ${F%_initial}"
        if [ "$?" -ne 0 ]; then
            echo "Failed to parse $F. Uniprot accession format may have changed. Exiting"
            exit_clean 1
        fi

        echo "Successfully unzipped and parsed $F$GZ"
        rm $DESTDIR/$F
    done


# NOTE: PE3 xdformats are disabled below
# Run on PE3
    COUNTV=`grep -c \> $DESTDIR/uniprot_PE3_vertebrata`
    COUNTVF=`grep -c \> $DESTDIR/uniprot_PE3_vertebrata_frag`
    COUNTNV=`grep -c \> $DESTDIR/uniprot_PE3_nonvert`
    COUNTNVF=`grep -c \> $DESTDIR/uniprot_PE3_nonvert_frag`

    echo "xdformatting PE3..."
    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE3_vertebrata $DESTDIR/uniprot_PE3_vertebrata > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$COUNTV" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE3_vertebrata_wfrag $DESTDIR/uniprot_PE3_vertebrata $DESTDIR/uniprot_PE3_vertebrata_frag > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$(($COUNTV+$COUNTVF))" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE3 $DESTDIR/uniprot_PE3_vertebrata $DESTDIR/uniprot_PE3_nonvert > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$(($COUNTV+$COUNTNV))" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE3_wfrag $DESTDIR/uniprot_PE3_vertebrata $DESTDIR/uniprot_PE3_vertebrata_frag $DESTDIR/uniprot_PE3_nonvert $DESTDIR/uniprot_PE3_nonvert_frag > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$(($COUNTV+$COUNTVF+$COUNTNV+$COUNTNVF))" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 9
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE123 $DESTDIR/uniprot_PE12_vertebrata $DESTDIR/uniprot_PE12_nonvert $DESTDIR/uniprot_PE3_vertebrata $DESTDIR/uniprot_PE3_nonvert > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    else
        rm $DESTDIR/xdformat.log
    fi

    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -o $DESTDIR/PE123_wfrag $DESTDIR/uniprot_PE12_vertebrata $DESTDIR/uniprot_PE12_vertebrata_frag $DESTDIR/uniprot_PE12_nonvert $DESTDIR/uniprot_PE12_nonvert_frag $DESTDIR/uniprot_PE3_vertebrata $DESTDIR/uniprot_PE3_vertebrata_frag $DESTDIR/uniprot_PE3_nonvert $DESTDIR/uniprot_PE3_nonvert_frag > $DESTDIR/xdformat2.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 8
    else
        rm $DESTDIR/xdformat2.log
    fi

fi
#######################################################################
#
# End classification processing code
#
#######################################################################


#######################################################################
#
# Begin main db files processing code
#
#######################################################################

# This is the final section of code. It first checks the release stats
# for the expected sequence counts. After that it waits to ensure the
# two main db files have downloaded. Once they have it processes them,
# checking the counts are correct and building the entry_loc file. If
# these are successful it will then concatenate both files together and
# run xdformat on the concatenated file. Once this is complete everything
# is done.
#
# Requirements for executing this code:
# uniprot_trembl.fasta.gz and uniprot_trembl.fasta.gz must have downloaded.
#
# Requirements for exiting this code:
# The sequence count must match the stats page on the website, the entry.loc
# file has to be created and xdformat must run successfully
#
# Outcomes of running this code:
# The uniprot BLAST db will be complete and everything should be finished.

if [ "$MAIN" -eq 1 ];then
    RCOUNT=0
# Getting the total number of sequences
    for C in "http://web.expasy.org/docs/relnotes/relstat.html" "http://www.ebi.ac.uk/uniprot/TrEMBLstats/"
        do
            RCOUNT=$(( RCOUNT += `wget -q -O - "$C" | perl -ne 'if(/(\d+)\s+sequence\s+entries/) {print $1}'` ))
        done

# Now wait for the download of the main db files to finish, this should take longest
    if [ "$DOWNLOAD" -eq 1 ];then
        echo "Waiting for main db files to complete downloading..."
        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' -w "done("uniprot_sprot.fasta_updl") && done("uniprot_trembl.fasta_updl")" echo "...main db files downloaded"
    else
        MAINF_ERROR=0
        for F in ${main_files[@]};do
            if [ -e "$DESTDIR/$F$GZ" ];then
                echo "$F$GZ has been downloaded"
            else
                echo "$F$GZ is missing"
                MAINF_ERROR=1
            fi
        done
        if [ "$MAINF_ERROR" -eq 1 ];then
            exit_clean 3
        fi
    fi


# Formatting the header
    for F in ${main_files[@]};
        do
            echo "Processing $F..."

            bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "gunzip -c $DESTDIR/$F$GZ | perl -ne 'if (/>/) { s/>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]).*/>\$1.\$4 \$3/} else {s/O/K/g} print;' > $DESTDIR/$F"
            if [ "$?" -ne 0 ];then
                echo "Failed to unzip the file $F$GZ and clean the header"
                exit_clean 1
            fi


            FCOUNT=`grep -c \> "$DESTDIR/$F"`
            COUNT=$(( COUNT += $FCOUNT ))
            grep \> $DESTDIR/$F | awk "{print \$1, \"${F}\"}" | sed -e 's/>//;s/uniprot_sprot.fasta/STD/;s/uniprot_trembl.fasta/PRE/' > $DESTDIR/$F.entry_loc
            ECOUNT=`wc -l "$DESTDIR/$F.entry_loc" | awk '{print $1}'`
            if [ "$FCOUNT" -ne "$ECOUNT" ];then
                echo "$DESTDIR/$F.entry_loc creation has failed, rerun:"
                echo "grep \> $DESTDIR/$F | awk \"{print \\\$1, \\\"${F}\\\"}\" | sed -e 's/>//;s/uniprot_sprot.fasta/STD/;s/uniprot_trembl.fasta/PRE/' > $DESTDIR/$F.entry_loc"
            fi
            echo "done"
        done

    if [ "$RCOUNT" -eq "$COUNT" ];then
        echo "Files downloaded successfully"
    else
        echo "From the stats: $RCOUNT"
        echo "From your files: $COUNT"
        exit_clean 3
    fi

# Concatenating the files
    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "cat $DESTDIR/uniprot_sprot.fasta $DESTDIR/uniprot_trembl.fasta > $DESTDIR/uniprot"
    if [ "$?" -ne 0 ];then
        echo "Failed to concatenate uniprot files"
        exit_clean 11
    fi

    COUNT=`grep -c \> $DESTDIR/uniprot`
# Concatenating the entry_loc files
    if [ "$FULL_ENTRY_LOC" -eq 1 ]; then
        bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "cat $DESTDIR/uniprot_sprot.fasta.entry_loc $DESTDIR/uniprot_trembl.fasta.entry_loc > $DESTDIR/entry_loc"
        if [ "$?" -ne 0 ];then
            echo "Failed to concatenate entry_loc files"
            exit_clean 11
        fi
        ECOUNT=`wc -l $DESTDIR/entry_loc | awk '{print $1}'`
        if [ "$ECOUNT" -ne "$COUNT" ];then
            echo "From the entry_loc file: $ECOUNT"
            echo "From your files: $COUNT"
            echo "You need to re-create the entry_loc file"
        fi
    else
        mv $DESTDIR/uniprot_sprot.fasta.entry_loc $DESTDIR/entry_loc
        echo "You only have SwissProt accession in your entry_loc file"
    fi

    if [ "$RCOUNT" -eq "$COUNT" ];then
        echo "Concatenation successful"
        rm $DESTDIR/uniprot_sprot.* $DESTDIR/uniprot_trembl.*
    else
        echo "Concatenation failed"
        exit_clean 4
    fi


# Formatting for blast
    echo "Formating with xdformat..."
    bsub -I -M 500 -R 'select[mem>500] rusage[mem=500]' "xdformat -p -t uniprot -v 'uniprot_$UNIPROT_VERSION' -d '$UNIPROT_DATE' $DESTDIR/uniprot > $DESTDIR/xdformat.log"
    if [ "$?" -ne 0 ];then
        echo "Failed to format uniprot with xdformat"
        exit_clean 5
    fi

    XCOUNT=`cat $DESTDIR/xdformat.log | grep 'written:' | awk '{print $6}' | sed 's/,//g'`
    if [ "$XCOUNT" -eq "$COUNT" ];then
        echo "xdformat successful"
        rm $DESTDIR/xdformat.log
    else
        echo "xdformat failed"
        exit_clean 6
    fi

    echo "Your uniprot databases are done!"
fi


#######################################################################
#
# End main db files processing code
#
#######################################################################


#######################################################################
#
# Begin changing permissions for the group
#
#######################################################################

chmod -R g+w $DESTDIR

#######################################################################
#
# End changing permissions for the group
#
#######################################################################
