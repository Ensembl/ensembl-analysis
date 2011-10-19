# Conventions:
#  - use lower case and underscores
#  - internal ids are integers named tablename_id
#  - same name is given in foreign key relations

################################################################################
#
# First clean the tables... :D
#
DROP TABLE IF EXISTS track_evidence, analysis_run, evidence, evidence_coord, \
    input_seq, reasons, dbs;


################################################################################
#
# Table structure for table 'track_evidence'
#

CREATE TABLE track_evidence (

  evidence_id                INT(10) UNSIGNED NOT NULL,
  analysis_run_id            SMALLINT(3) UNSIGNED NOT NULL,
  reason_id                  SMALLINT(3) UNSIGNED NOT NULL DEFAULT 0,
  input_id                   VARCHAR(100) NOT NULL,
#  is_current                 BOOLEAN NOT NULL DEFAULT 1,

  KEY evidence_idx (evidence_id),
  KEY analysis_run_idx (analysis_run_id),
  KEY reason_idx (reason_id),
  KEY input_idx (input_id)

) COLLATE=latin1_swedish_ci ENGINE=InnoDB;


################################################################################
#
# Table structure for table 'analysis_run'
#

CREATE TABLE analysis_run (

  analysis_run_id            SMALLINT(3) UNSIGNED NOT NULL AUTO_INCREMENT,
  analysis_id                SMALLINT(3) UNSIGNED NOT NULL,
  run_date                   DATETIME DEFAULT '0000-00-00 00:00:00' NOT NULL,
  output_db_id               VARCHAR(5) NOT NULL,
  input_db_id                VARCHAR(5) NOT NULL,

  PRIMARY KEY (analysis_run_id),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;



################################################################################
# Uniquely define each of the input protein/mRNA alignments by its
# hit_name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand
#
# Table structure for table 'evidence'
#
#
#CREATE TABLE evidence (
#  evidence_id                 INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
#  input_seq_id                INT(10) UNSIGNED NOT NULL,
#  seq_region_id               INT(10) UNSIGNED NOT NULL,
#  seq_region_start            INT(10) UNSIGNED NOT NULL,
#  seq_region_end              INT(10) UNSIGNED NOT NULL,
#  seq_region_strand           TINYINT(1) NOT NULL,
#
#  PRIMARY KEY (evidence_id),
#  KEY input_seq_idx (input_seq_id),
#  KEY seq_region_idx (seq_region_id)
#
#) COLLATE=latin1_swedish_ci ENGINE=InnoDB;




################################################################################
# Uniquely define each of the input protein/mRNA alignments by its
# hit_name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand
#
# Table structure for table 'evidence'
#

CREATE TABLE evidence (
  evidence_id                 INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  input_seq_id                INT(10) UNSIGNED NOT NULL,
  is_aligned                  ENUM('u', 'y', 'n') DEFAULT 'u',

  PRIMARY KEY (evidence_id),
  KEY input_seq_idx (input_seq_id)

) COLLATE=latin1_swedish_ci ENGINE=InnoDB;




################################################################################
# Uniquely define each of the input protein/mRNA alignments by its
# hit_name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand
#
# Table structure for table 'evidence_coord'
#

CREATE TABLE evidence_coord (
  evidence_id                 INT(10) UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) NOT NULL,

  PRIMARY KEY (evidence_id),
  KEY seq_region_idx (seq_region_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;




################################################################################
# A list of protein or mRNA accessions that are present in the input files
#
# Table structure for table 'inputseq'
#

CREATE TABLE input_seq (

  input_seq_id                INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  hit_name                    VARCHAR(40) NOT NULL,
  external_db_id              SMALLINT UNSIGNED,
  molecule_type               ENUM('PROTEIN','MRNA','EST'),
  submission_date             DATE DEFAULT '0000-00-00' NOT NULL, 

  PRIMARY KEY (input_seq_id),
  KEY hit_namex (hit_name),
  KEY molecule_typex (molecule_type),
  KEY external_db_idx (external_db_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;




################################################################################
# A list of reason why it was accepted or rejected
#
# Table structure for table 'reasons'
#

CREATE TABLE reasons (

  reason_id                 SMALLINT(3) UNSIGNED NOT NULL,
  code						CHAR(20) NOT NULL,
  reason                    VARCHAR(60) NOT NULL,

  PRIMARY KEY (reason_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;




################################################################################
# A list of databases used during the process
#
# Table structure for table 'dbs'
#

CREATE TABLE dbs (

  db_id                 TINYINT(3) UNSIGNED NOT NULL AUTO_INCREMENT,
  db_name               VARCHAR(40) NOT NULL,
  instance              VARCHAR(20) NOT NULL,

  PRIMARY KEY (db_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
