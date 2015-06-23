#!/bin/bash
export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/ensembl-external/modules:$PWD/modules:$PWD/scripts:$PWD/scripts/buildchecks:$PWD/ensembl-compara/modules:$PWD/ensembl-funcgen/modules:$PWD/ensembl-killlist/modules:$PWD/ensembl-pipeline/scripts:$PWD/ensembl-pipeline/modules:$PWD/ensembl-hive/modules:$PWD/bioperl-live:$PWD/bioperl-run/lib:$PWD/ensembl-56/modules

export WORK_DIR=$PWD

echo "Running test suite"
echo "Using $PERL5LIB"
if [ "$COVERALLS" = 'true' ]; then
  export PERL5LIB=$PERL5LIB:$PWD/ensembl-test/modules
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
else
  # just test the basic syntax for all the scripts and modules - start by renaming example modules
  # We need to fake some configuration files
  # Some of these config files should probably be remove from the module
  PERL_IMPORT='sub import { my ($callpack) = caller(0); my $pack = shift; my @vars = @_ ? @_ : keys(%Config); return unless @vars; eval "package $callpack; use vars qw(" . join(" ", map { "\$".$_ } @vars) . ")"; die $@ if $@; foreach (@vars) { if (defined $Config{ $_ }) { no strict "refs"; *{"${callpack}::$_"} = \$Config{ $_ }; } else { die "Error: Config: $_ not known\n"; } } }'
  echo '1;' > modules/Bio/EnsEMBL/Analysis/Config/HaplotypeProjection.pm
  echo '1;' > modules/Bio/EnsEMBL/Analysis/Config/ExonerateClones.pm
#  echo '1;' > modules/Bio/EnsEMBL/Analysis/Config/ExonerateRefinedCloneEnds.pm
  echo "%Config = ( FILTER_CONFIG_BY_LOGIC => { DEFAULT => {}});$PERL_IMPORT;1;" > modules/Bio/EnsEMBL/Analysis/Config/FilterGenes.pm
  # Funcgen config
  printf "@RUNNABLE_CONFIG = ();\n\$ANALYSIS_WORK_DIR = '%s';\n\$ANALYSIS_INPUT_DIR = '%s';\n\$ANALYSIS_TARGET_DIR = '%s';%s\n1;\n" "$PWD" "$PWD" "$PWD" "$PERL_IMPORT" > $PWD/runnable_config.pm
  perl -c $PWD/runnable_config.pm
  cp $PWD/scripts/RNASeq/setup_rnaseq_pipeline_config.pm_example $PWD/modules/setup_rnaseq_pipeline_config.pm
  # We need to fake this module but it may be cleaned better later
  sed 's/Solexa2Genes/Solexa2GenesLiteNew/' $PWD/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/Solexa2Genes.pm.example > $PWD/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/Solexa2GenesLiteNew.pm
  find $PWD/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  find $PWD/scripts -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  find $PWD/ensembl-pipeline/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  find $PWD/scripts -type f -name "*.pl" | xargs -i perl -c {}
# We avoid the Finished directory at the moment
#  Exonerate2Array.pm as it is a FuncGen module
#  ExonerateRefinedCloneEnds.pm as we have a newer module for the clone ends
#  M=( "Bio/EnsEMBL/Analysis/RunnableDB/Exonerate2Array.pm" "Bio/EnsEMBL/Analysis/RunnableDB/ExonerateRefinedCloneEnds.pm" "Bio/EnsEMBL/Analysis/RunnableDB/FilterGenes.pm" )
  M=( "Bio/EnsEMBL/Analysis/RunnableDB/ExonerateRefinedCloneEnds.pm" )
  ARRAY=`seq 0 $((${#M[@]}-1))`
  printf "\e[31mWe will not test:\n - %s\n" "Annacode modules" 
  for S in $ARRAY; do
      printf " - %s\n" "${M[$S]}"
      N[$S]=`basename ${M[$S]}`
  done
  printf "\e[0m\n"
  find $PWD/modules -type f -name "*.pm" ! -path "*Finished*" `for I in $ARRAY; do RES=${RES}" ! -name ${N[$I]}"; done; echo "$RES"` | xargs -i perl -c {}
fi
rt=$?
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
