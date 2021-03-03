#!/usr/bin/env bash


# exit on any error
set -e


# update the following values
################################################################################
COMMON_NAME="African pygmy mouse"
SCIENTIFIC_NAME="Mus minutoides"
ASSEMBLY_ACCESSION="GCA_902729485.1"

SCHEMA_VERSION="103"
################################################################################


# copy rnaseq merge files to production
################################################################################
COMMON_NAME_UNDERSCORES="${COMMON_NAME// /_}"
SCIENTIFIC_NAME_UNDERSCORES="${SCIENTIFIC_NAME// /_}"
SCIENTIFIC_NAME_UNDERSCORES_LOWER_CASE="${SCIENTIFIC_NAME_UNDERSCORES,,}"

ANNOTATION_NAME="${COMMON_NAME_UNDERSCORES}-${SCIENTIFIC_NAME_UNDERSCORES}-${ASSEMBLY_ACCESSION}"

DATA_DIRECTORY_NAME="${ANNOTATION_NAME}"

ANNOTATIONS_DATA_ROOT="/hps/nobackup2/production/ensembl/genebuild/production"
ANNOTATION_DATA_DIRECTORY="${ANNOTATIONS_DATA_ROOT}/${DATA_DIRECTORY_NAME}/${SCIENTIFIC_NAME_UNDERSCORES_LOWER_CASE}/${ASSEMBLY_ACCESSION}"

RAPID_RELEASE_FTP_ROOT="/nfs/production/panda/ensembl/production/ensemblftp/rapid-release/species"

RNASEQ_DIRECTORY="${RAPID_RELEASE_FTP_ROOT}/${SCIENTIFIC_NAME_UNDERSCORES}/${ASSEMBLY_ACCESSION}/rnaseq"

mkdir --parents --verbose "$RNASEQ_DIRECTORY"
chmod g+w "$RNASEQ_DIRECTORY"

rsync --recursive --copy-links --times --devices --specials --whole-file --verbose --human-readable --progress --stats "${ANNOTATION_DATA_DIRECTORY}"/rnaseq/merge/* "${RNASEQ_DIRECTORY}/"
################################################################################


# create symlink for web
################################################################################
ASSEMBLY_ACCESSION_LOWER_CASE="${ASSEMBLY_ACCESSION,,}"
ASSEMBLY_ACCESSION_REMOVE_UNDERSCORE="${ASSEMBLY_ACCESSION_LOWER_CASE/_/}"
ASSEMBLY_ACCESSION_REPLACE_DOT="${ASSEMBLY_ACCESSION_REMOVE_UNDERSCORE/\./v}"

CORE_DATABASE_NAME="${SCIENTIFIC_NAME_UNDERSCORES_LOWER_CASE}_${ASSEMBLY_ACCESSION_REPLACE_DOT}_core_${SCHEMA_VERSION}_1"

SPECIES_PRODUCTION_NAME=$(gb1 $CORE_DATABASE_NAME --skip-column-names -e "SELECT meta_value FROM meta WHERE meta_key='species.production_name';")
ASSEMBLY_NAME=$(gb1 $CORE_DATABASE_NAME --skip-column-names -e "SELECT meta_value FROM meta WHERE meta_key='assembly.name';")

RAPID_RELEASE_DATA_FILES_ROOT="/nfs/production/panda/ensembl/production/ensemblftp/rapid-release/web/data_files/"

DATA_FILES_BASE_DIRECTORY="${RAPID_RELEASE_DATA_FILES_ROOT}/${SPECIES_PRODUCTION_NAME}"
DATA_FILES_DIRECTORY="${DATA_FILES_BASE_DIRECTORY}/${ASSEMBLY_NAME}"

mkdir --parents --verbose "$DATA_FILES_DIRECTORY"
chmod --recursive g+w "$DATA_FILES_BASE_DIRECTORY"

ln --symbolic --verbose "$RNASEQ_DIRECTORY" "$DATA_FILES_DIRECTORY"
################################################################################
