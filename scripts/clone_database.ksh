#!/bin/ksh
# $Id: clone_database.ksh,v 1.2 2013-11-06 16:55:45 ak4 Exp $

self="$0"
self_base="${self##*/}"

function parse_dbarg
{
  # Arguments:
  #     $1 a name of a parameter to be used for output.
  #     $2 is a database specification on the form "database@host:port".
  #
  # This function parses the database specification (the ":port" part is
  # optional) and returns a string which can be eval'ed and which then
  # sets the '.dbname', '.dbhost' and '.dbport' fields of the compound
  # variable whose name is given as the first argument.

  print $2 | tr '@:' '  ' |
  awk -vparam="$1" '
    NF == 3 {
      printf("%s=( dbname='%s' dbhost='%s' typeset -i dbport=%d )",
        param, $1, $2, $3);
    }
    NF == 2 {
      printf("%s=( dbname='%s' dbhost='%s' typeset -i dbport=3306 )",
        param, $1, $2);
    }'
}

function usage
{
  cat <<USAGE_END
USAGE:
  ${self_base} [-f|-u] -s database@host[:port] -t database@host[:port]

  ${self_base} -s database@host[:port] -o filename

  ${self_base} -h

Creates a new Core database and copies important tables from a reference
(source) database to it.  Alternatively, updates important tables in an
existing database with the ones from a reference (source) database.

  -s    ("source") Source database on the form "database@host:port",
        where the ":port" part is optional (default port is 3306).

  -t    ("target") Target database on the form "database@host:port",
        where the ":port" part is optional (default port is 3306).

  -f    ("force") Force creation of the target database by dropping it
        if it already exists.  The -f switch is ignored if -u is also
        used.

  -l    ("lock_table") Do NOT lock the table when dumping. This can be
        dangerous but useful if your read-only user has no LOCK TABLES
        priviledge. Default is to lock the tables.

  -u    ("update") Do not try to create or recreate the target database
        but update it with the important tables from the source
        database.  In this mode, the -f is ignored.
        NOTE:  The important tables are dropped and recreated.

  -o    ("output") Instead of sending the output to a new database, dump
        it to a file.  The filename must be given as an argument to this
        flag.  In this mode, only the -s flag will be used.

  -r    ("read-only") user name of the source database should be read-only
        and without password (or change the script)

  -w    ("read-write") user name of the target database

  -P    ("password") password of the user of the target database

  -h    ("help") Help, displays this text and exists.

Important tables:
  analysis
  analysis_description
  assembly
  assembly_exception
  attrib_type
  coord_system
  external_db
  meta
  misc_set
  seq_region
  seq_region_attrib
  seq_region_mapping

Example:

  ${self_base} -u -s ak4_chrysemys_picta_ref@genebuild1 \\
    -t ak4_chrysemys_picta_exonerate@genebuild2 \\
    -r ro_user -w rw_user -P pass

USAGE_END
}

opt_update=0
opt_force=0
opt_output=''
opt_locktable='true'

while getopts 's:t:uflo:r:w:P:h' opt; do
  case ${opt} in
    s)  eval $( parse_dbarg 'opt_source' ${OPTARG} )   ;;
    t)  eval $( parse_dbarg 'opt_target' ${OPTARG} )   ;;
    u)  opt_update=1    ;;
    f)  opt_force=1     ;;
    l)  opt_locktable='false'     ;;
    o)  opt_output=${OPTARG}    ;;
    r)  opt_ro_user=${OPTARG}    ;;
    w)  opt_rw_user=${OPTARG}    ;;
    P)  opt_rw_pass=${OPTARG}    ;;
    h)  usage; exit 0   ;;
    *)  usage; exit 1   ;;
  esac
done

if [[ -z ${opt_source.dbhost} || -z ${opt_source.dbname} ]]; then
  print -u2 "ERROR: Incomplete source specification."
  usage; exit 1
elif [[ -z ${opt_output} && ( -z ${opt_target.dbhost} || -z ${opt_target.dbname} ) ]]; then
  print -u2 "ERROR: Incomplete target specification."
  usage; exit 1
fi
if [[ -z ${opt_ro_user} ]]; then
  print -u2 "ERROR: Missing read-only user."
  usage; exit 1
fi
if [[ -z ${opt_rw_user} ]]; then
  print -u2 "ERROR: Missing read-write user."
  usage; exit 1
fi
if [[ -z ${opt_rw_pass} ]]; then
  print -u2 "ERROR: Missing read-write password."
  usage; exit 1
fi

set -o pipefail
if [[ -z ${opt_output} ]]; then
  if (( !opt_update )); then
    # Create the new database:
    if (( opt_force )); then
      mysql --host=${opt_target.dbhost} --port=${opt_target.dbport} \
        --user=${opt_rw_user} --password=${opt_rw_pass} \
        --verbose --execute="DROP DATABASE IF EXISTS ${opt_target.dbname}"
    fi

    mysql --host=${opt_target.dbhost} --port=${opt_target.dbport} \
      --user=${opt_rw_user} --password=${opt_rw_pass} \
      --verbose --execute="CREATE DATABASE ${opt_target.dbname}" || exit 1

    # Dump table definitions from the source database and apply them to
    # the new target database.
    mysqldump --host=${opt_source.dbhost} --port=${opt_source.dbport} \
      --user=${opt_ro_user} \
      --no-data --lock-tables=${opt_locktable} \
      ${opt_source.dbname} |
    mysql --host=${opt_target.dbhost} --port=${opt_target.dbport} \
      --user=${opt_rw_user} --password=${opt_rw_pass} \
      --database=${opt_target.dbname}
  fi

  output_cmd="mysql --host=${opt_target.dbhost}
  --port=${opt_target.dbport}
  --user=${opt_rw_user} --password=${opt_rw_pass}
  --database=${opt_target.dbname}"

  print "WRITING OUTPUT TO DATABASE ${opt_target.dbname}"
else
  output_cmd="cat >${opt_output}"

  print "WRITING OUTPUT TO FILE ${opt_output}"
fi

# Dump important tables from the source database and load them into the
# target database:
mysqldump --host=${opt_source.dbhost} --port=${opt_source.dbport} \
  --user=${opt_ro_user} \
  --verbose --lock-tables=${opt_locktable} \
  ${opt_source.dbname} \
  analysis \
  analysis_description \
  assembly \
  assembly_exception \
  attrib_type \
  coord_system \
  external_db \
  meta \
  misc_set \
  seq_region \
  seq_region_attrib \
  seq_region_synonym \
  karyotype \
  mapping_set \
  data_file \
  seq_region_mapping | $( eval ${output_cmd} )
