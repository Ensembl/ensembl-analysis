
CREATE TABLE assembly (

  assembly_id                 INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gca                         VARCHAR(14) NOT NULL,
  species_id                  INT(10) NOT NULL,

  PRIMARY KEY (assembly_id),

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE species (

  species_id                  INT(10) NOT NULL AUTO_INCREMENT,
  taxon_id                    INT(10) UNSIGNED NOT NULL,
  common_name                 VARCHAR(40) NOT NULL,
  group_name                  VARCHAR(40) NOT NULL,

  PRIMARY KEY (species_id),

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE repeat_sequence (

  repeat_sequence_id          INT(10) NOT NULL AUTO_INCREMENT,
  repeat_class_id             INT(10) NOT NULL,
  species_id                  INT(10) UNSIGNED NOT NULL,
  assembly_id                 INT(10) UNSIGNED NOT NULL,

  PRIMARY KEY (repeat_sequence_id),

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE repeat_class (

  repeat_class_id             INT(10) NOT NULL AUTO_INCREMENT,
  repeat_name                 VARCHAR(255) NOT NULL,
  repeat_class                VARCHAR(100) NOT NULL,
  repeat_type                 VARCHAR(40) NOT NULL,
  repeat_sequence             LONGTEXT NOT NULL,

  PRIMARY KEY (repeat_class_id),
  KEY name (repeat_name),
  KEY class (repeat_class),
  KEY type (repeat_type),

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
