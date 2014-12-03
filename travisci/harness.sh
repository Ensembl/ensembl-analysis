#!/bin/bash
export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/ensembl-external/modules:$PWD/modules:$PWD/scripts:$PWD/scripts/buildchecks:$PWD/ensembl-compara/modules:$PWD/ensembl-funcgen/modules:$PWD/ensembl-killlist/modules:$PWD/ensembl-pipeline/scripts:$PWD/ensembl-pipeline/modules

echo "Running test suite"
echo "Using $PERL5LIB"
if [ "$COVERALLS" = 'true' ]; then
  export PERL5LIB=$PERL5LIB:$PWD/ensembl-test/modules
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
else
  # just test the basic syntax for all the scripts and modules - start by renaming example modules
# We need to fake some configuration files
  echo '1;' > modules/Bio/EnsEMBL/Analysis/Config/HaplotypeProjection.pm
  cp $PWD/scripts/RNASeq/setup_rnaseq_pipeline_config.pm_example $PWD/modules/setup_rnaseq_pipeline_config.pm
  find $PWD/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  find $PWD/scripts -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  find $PWD/ensembl-pipeline/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  find $PWD/scripts -type f -name "*.pl" | xargs -i perl -c {}
# We avoid the Finished directory at the moment
  find $PWD/modules -type f -name "*.pm" ! -path "*Finished*" | xargs -i perl -c {}
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
