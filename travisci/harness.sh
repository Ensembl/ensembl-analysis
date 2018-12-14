#!/bin/bash
export PERL5LIB=$PWD/bioperl-live:$PWD/ensembl/modules:$PWD/ensembl-variation/modules:$PWD/ensembl/modules:$PWD/ensembl-external/modules:$PWD/modules:$PWD/scripts:$PWD/scripts/buildchecks:$PWD/ensembl-compara/modules:$PWD/ensembl-funcgen/modules:$PWD/ensembl-killlist/modules:$PWD/ensembl-pipeline/scripts:$PWD/ensembl-pipeline/modules:$PWD/ensembl-hive/modules:$PWD/ensembl-io/modules:$PWD/bioperl-run/lib:$PWD/ensembl-56/modules:$PWD/GIFTS/modules:$PWD/ensembl-test/modules

export WORK_DIR=$PWD

if [ "$DB" = 'mysql' ]; then
    (cd modules/t && ln -sf MultiTestDB.conf.mysql MultiTestDB.conf)
elif [ "$DB" = 'sqlite' ]; then
    (cd modules/t && ln -sf MultiTestDB.conf.SQLite MultiTestDB.conf)
else
    echo "Don't know about DB '$DB'"
    exit 1;
fi

echo "Running test suite"
echo "Using $PERL5LIB"
rt=0
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
else
  # just test the basic syntax for all the scripts and modules - start by renaming example modules
  # We need to fake some configuration files
  # Some of these config files should probably be remove from the module
  PERL_IMPORT='sub import { my ($callpack) = caller(0); my $pack = shift; my @vars = @_ ? @_ : keys(%Config); return unless @vars; eval "package $callpack; use vars qw(" . join(" ", map { "\$".$_ } @vars) . ")"; die $@ if $@; foreach (@vars) { if (defined $Config{ $_ }) { no strict "refs"; *{"${callpack}::$_"} = \$Config{ $_ }; } else { die "Error: Config: $_ not known\n"; } } }'
  echo '1;' > modules/Bio/EnsEMBL/Analysis/Config/HaplotypeProjection.pm
  # Funcgen config
  printf "@RUNNABLE_CONFIG = ();\n\$ANALYSIS_WORK_DIR = '%s';\n\$ANALYSIS_INPUT_DIR = '%s';\n\$ANALYSIS_TARGET_DIR = '%s';%s\n1;\n" "$PWD" "$PWD" "$PWD" "$PERL_IMPORT" > $PWD/runnable_config.pm
  perl -c $PWD/runnable_config.pm
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  cp $PWD/scripts/RNASeq/setup_rnaseq_pipeline_config.pm_example $PWD/modules/setup_rnaseq_pipeline_config.pm
  # We need to fake this module but it may be cleaned better later
  sed 's/Solexa2Genes/Solexa2GenesLiteNew/' $PWD/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/Solexa2Genes.pm.example > $PWD/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/Solexa2GenesLiteNew.pm
  find $PWD/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  find $PWD/scripts -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  find $PWD/ensembl-pipeline/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  find $PWD/scripts -type f -name "*.pl" | xargs -i perl -c {}
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
# We avoid the Finished directory at the moment
#  Exonerate2Array.pm as it is a FuncGen module
#  ExonerateRefinedCloneEnds.pm as we have a newer module for the clone ends
#  SWEmbl.pm, InputSet.pm has been removed from ensembl-funcgen and it doesn't affect us

  M=( "Bio/EnsEMBL/Analysis/RunnableDB/ExonerateRefinedCloneEnds.pm" \
  "Bio/EnsEMBL/Analysis/RunnableDB/ExonerateClones.pm" \
  "Bio/EnsEMBL/Analysis/RunnableDB/Exonerate2Array.pm" \
  "Bio/EnsEMBL/Analysis/RunnableDB/Funcgen/SWEmbl.pm" \
  "Bio/EnsEMBL/Analysis/RunnableDB/FilterGenes.pm" )
  printf "\e[31mWe will not test:\e[0m\n - \e[33m%s\e[0m\n" "Annacode modules"
  for S in `seq 0 $((${#M[@]}-1))`; do
      printf " - \e[33m%s\n\e[0m" "${M[$S]}"
      RES=${RES}" ! -name `basename ${M[$S]}`"
  done
  if [ "$OLDPERL" = 'true' ];then
    M=( "Bio/EnsEMBL/Analysis/Hive/Config/TSLsAppris_conf.pm" \
    "Bio/EnsEMBL/Analysis/Hive/Config/UniProtDB_conf.pm" \
    "Bio/EnsEMBL/Analysis/Hive/RunnableDB/JiraTicket.pm" \
    "Bio/EnsEMBL/Analysis/Hive/RunnableDB/HiveLoadPDBProteinFeatures.pm" \
    "Bio/EnsEMBL/Analysis/Runnable/MakePDBProteinFeatures.pm" \
    "Bio/EnsEMBL/Analysis/Hive/Config/HiveRelCoDuties_conf.pm" )
    for S in `seq 0 $((${#M[@]}-1))`; do
        printf " - \e[33m%s\n\e[0m" "${M[$S]}"
        RES=${RES}" ! -name `basename ${M[$S]}`"
    done
  fi
  find $PWD/modules -type f -name "*.pm" ! -path "*Finished*" `echo "$RES"` | xargs -i perl -c {}
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  export EHIVE_ROOT_DIR=$PWD/ensembl-hive
  perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi

fi
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
