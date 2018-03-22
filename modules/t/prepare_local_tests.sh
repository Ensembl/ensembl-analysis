#!/usr/bin/env bash
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

CURRENT_DIR=`dirname $0`
if [ ! -e "$CURRENT_DIR/MultiTestDB.conf" ];then
  echo "You should create a file $CURRENT_DIR/MultiTestDB.conf containing connection details to a database:"
  cat <<EOF
{
  port => 3306,
  user => EHIVE_USER,
  pass => EHIVE_PASS,
  host => HOST,
  driver => 'mysql',
}
EOF
  echo
fi

BASEDIR="modules/t/test-genome-DBs/homo_sapiens/core"
ENSEMBLDIR="../ensembl/$BASEDIR"
if [ ! -e "$ENSEMBLDIR" ];then
  if [ -n "$PERL5LIB" ];then
    for D in `echo $PERL5LIB | sed 's/:/\n/g'`; do
      if [ "$D" != "${D/ensembl\/modules}" ]; then
        ENSEMBLDIR=$D
      fi
    done
  else
    printf "\033[31mPERL5LIB is not set\033[0m\n"
    exit 1
  fi
fi

if [ ! -e "$BASEDIR" ]; then
  mkdir -p "$BASEDIR"
fi

for F in ${ENSEMBLDIR}/*; do
# We also want the SQLite table in case we start testing it too
  cp -r "$F" "$BASEDIR"
done
