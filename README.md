# Ensembl-analysis

The EnsEMBL Genome Annotation System


## Requirements

### Softwares

We use Linuxbrew to install all our software you will need to tap cask from https://github.com/Ensembl/homebrew-cask before being able to install `ensembl/cask/genebuild-annotation`.

You will need to ask for a licence for some of the software we use:
- crossmatch, part of phred/phrap: used by RepeatMasker, you could use NCBI blast instead
- trf
- repbase: you could use RepeatModeler to generate a repeat library. However the generated library could mask some genes coding for some protein families


### Perl EnsEMBL repositories you need to have

We recommend that you clone all the repositories into one directory
| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-hive | default | https://github.com/Ensembl/ensembl-hive.git |
| ensembl-compara | release/98 | https://github.com/Ensembl/ensembl-compara.git |
| ensembl-production | default | https://github.com/Ensembl/ensembl-production.git |
| ensembl-taxonomy | default | https://github.com/Ensembl/ensembl-taxonomy.git |
| ensembl-orm | default | https://github.com/Ensembl/ensembl-orm.git |
| ensembl-killlist | default | https://github.com/Ensembl/ensembl-killlist.git |
| ensembl-datacheck | default | https://github.com/Ensembl/ensembl-datacheck.git |
| ensembl-metadata | default | https://github.com/Ensembl/ensembl-metadata.git |
| ensembl-io | default | https://github.com/Ensembl/ensembl-io.git |
| ensembl-variation | default | https://github.com/Ensembl/ensembl-variation.git |
| core_meta_updates | default | https://github.com/Ensembl/core_meta_updates.git |

For each of these repository, you will need to install their dependencies using the cpanfile provided in their Git repositories

You can use the [Ensembl git commands](https://github.com/Ensembl/ensembl-git-tools) and run the following command to clone the repositories
```
git ensembl --clone genebuild
```

### Python EnsEMBL repositories you need to have

| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl-genes | default | https://github.com/Ensembl/ensembl-genes.git |

### Python virtual environment

You will need to create two virtual environment:
- `genebuild` using the requirements.txt file; it needs to be activated for the pipeline to run
- `genebuild-mirna` using the requirements\_p36\_ncrna.txt file; it will be used directly from the analysis, no need to activate it

### Shell environment

If you are not part of the Ensembl Genebuild team, you will need to set some shell environment variables to avoid having to provide the information to the configuration files. We will assume you are using your home directory
| Variable | Value | Hive configuration parameter | Description |
|----------|-------|------------------------------|-------------|
| ENSCODE | $HOME | -enscode\_root\_dir | Directory path where you cloned all the Perl repositories |
| ENSEMBL\_SOFTWARE\_HOME | $HOME | -software\_base\_path | Directory where pyenv, plenv and linuxbrew are installed |
| LINUXBREW\_HOME | $HOME/.linuxbrew | -linuxbrew\_home\_path | Base directory for your Linuxbrew installation |
| PYTHONPATH | $HOME/ensembl-genes/ensembl\_genes:$HOME/ensembl-hive/wrappers/python3/ | | It needs to be set until the package can be installed properly |
| BLASTDB\_DIR | $HOME | | It will be used to find the path to the entry\_loc file which is the list of accession from SwissProt which would be located at $HOME/uniprot/2021\_03/entry\_loc |

### MySQL

We currently use MySQL databases to store our data. To avoid having to do many changes to the configuration files we recommend having one read-only user and one read-write user. It is also better to use different servers for keeping the eHive pipeline database, the DNA database and the "data" databases.


## Running the EnsEMBL Genome Annotation System

There is a main configuration file, `Bio::EnsEMBL::Analysis::Hive::Config::Genome_annotation_conf`, which will generate a set of pipelines to:
- load the DNA into the `dna_db` database
- mask repeats
- align multiple sets of data to the genome to produce gene models depending on the data available
- select the best transcripts for each loci
- produce the finalised databases for an Ensembl release
The whole system is explained in more details below

### Initialising the pipeline

You will need to activate the genebuild virtual environment
```
pyenv activate genebuild
```

#### Filling the main configuration automatically

If you are operating within an environment prepared for Ensembl with the assembly registry you can use the `$ENSCODE/ensembl-analysis/scripts/genebuild/create_annotation_configs.pl`.

You would need to edit `$ENSCODE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/genome_annotation.ini`

Then you can run
```
perl $ENSCODE/ensembl-analysis/scripts/genebuild/create_annotation_configs.pl --config_file $ENSCODE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/genome_annotation.ini --assembly_registry_host <host_name> --assembly_registry_port <port>
```

#### Filling the main configuration manually

If you're setup do not work with the create\_annotation\_configs.pl script, you would need to edit `$ENSCODE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/Genome_annotation_conf.pm` and fill in any information according to you environment

Then you can run
```
perl $ENSCODE/ensembl-hive/scripts/init_Pipeline.pl Bio::EnsEMBL::Analysis::Hive::Config::Genome_annotation_conf [extra parameters]
```

### Running the pipeline

To start the pipeline you need the URL to your pipeline database which will be provided when running the init\_Pipeline.pl script. If you initialised the pipeline automatically, you need to look at the command file created in your `working_dir` directory to retrieve the information.
```
export EHIVE_URL=mysql://readwrite_user:password@host:port/dbname
```

You can now start the pipeline with
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop
```

If you only want to run some analyses, you can run
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop -analyses_pattern 1..5
```

### Monitoring the pipeline

#### GuiHive

To follow the pipeline steps, it is better to use GuiHive, a graphical interface to ensembl-hive, which allows you to change parameters, debug your problems and much more https://github.com/Ensembl/guiHive

#### What to do when the main pipeline fails

You should first look at the job tab to know the reason of the failure
* Insufficient memory: you can either use a different resource or add a new one more suited to your needs
* Error in the code: I'm afraid you will need to do proper debugging

Once you are happy with your fix, you would need to reset the jobs with
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -reset_failed_jobs
```
and restart the pipeline

#### How can I debug a job

By default ensembl-hive redirect all output to `/dev/null` unless you used some logging parameters.

You will need to run the problematic job with runWorker. First you will need to retrieve the job id using GuiHive or the pipeline database. Then you can run
```
perl $ENSCODE/ensembl-hive/scripts/runWorker.pl -url $EHIVE_URL -debug 1 -job_id XX
```

Using a higher value for `-debug` is usually not useful as it is mostly seen as a boolean flag.

#### What do to when a subpipeline fails

First you would need to go on the job tab and do a middle click/right click on the `guihive` link. You would need to insert the password twice, this is normal and not a scam.

Then from the "GuiHive" of the sub pipeline, you need to make your changes/debugging the same way as a normal pipeline. Once this is done, you can simply reset the main pipeline and restart the main pipeline


## What is the difference between a "main" pipeline and a "sub" pipeline

We use the main/sub terminology to define the different part of the system. We could have used parent/child
* The "main" pipeline will run multiple analyses and will initialise at least one sub pipeline and run this sub pipeline
* The "sub" pipeline will run multiple analyses and can potentially be the "main" pipeline of a different pipeline. A sub pipeline can be run on its own

A pipeline needs three analyses to be called a "main" pipeline:
* `create_<sub pipeline name>_jobs`
* `initialise_<sub pipeline name>`
* `run_<sub pipeline name>`

### create\_\<sub pipeline name\>\_jobs

The main reason for this analysis is to provide an easy access to the sub pipeline guiHive and the pipeline database details

It is usually a `Bio::EnsEMBL::Hive::RunnableDB::JobFactory` where the input list will have at least one element with four key value pair. The reason for this analysis it
  - ehive\_url: the sub pipeline database connection URI, `mysql://rw\_user:password@host:port/database`, useful for debugging as it can be used
  - external\_link: the sub pipeline guiHive URL
  - meta\_pipeline\_db: the sub pipeline database connection hash, `{'-dbname' => 'database','-driver' => 'mysql','-host' => 'host','-pass' => 'password','-port' => 3306,'-user' => 'rw_user'}`
  - pipeline\_name: the name of the sub pipeline

**Example analysis config**
```
{
  -logic_name => 'create_rnaseq_db_pipeline_job',
  -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
  -parameters => {
    column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link'],
    inputlist => [[$rnaseq_db_pipe_db, $rnaseq_db_pipe_url, 'rnaseq_db_'.$self->o('production_name'), $rnaseq_db_guihive]],
    guihive_host => $self->o('guihive_host'),
    guihive_port => $self->o('guihive_port'),
  },
  -rc_name => 'default',
  -max_retry_count => 0,
  -flow_into => {
    2 => ['initialise_rnaseq_db'],
  }
},
```

### initialise\_\<sub pipeline name\>

This analysis will initialise the sub pipeline.

It must be a `Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit`, it will simply generate the commandline which will be run with `$ENSCODE/ensembl-hive/scripts/init_pipeline.pl`

There is a set of parameters which are needed:
* hive\_config: the name of the config to use
* databases: list of keys which contains the connection details of any database to be used such as 'dna\_db'. The key value pairs should not be in `extra_parameters`.
* extra\_parameters: any parameter to be provided to the initialisation script. `species_name => 'homo_sapiens'` would be used as `-species_name homo_sapiens`
* metadata\_pipeline\_db: the connection details of the pipeline database, this would usually be provided by the `create_<sub pipeline name>_jobs` analysis

There is one special parameter which is not required but can be helpful, `commandline_params`, which will not process the value associated. It can be used for re-initialising a sub pipeline with something like `-hive_force_init 1` or `-hive_force_init 1 -drop_databases 1` if the pipeline has already been started/run and you want to drop any databases created by the pipeline.

Arrays and hashes cannot be easily passed to the init script on the command line. Also Hive does not overwrite arrays, it always add elements to the parameter if the array already exist.

**Example analysis config**
```
{
  -logic_name => 'initialise_rnaseq_db',
  -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit',
  -parameters => {
    hive_config => $self->o('hive_rnaseq_db_config'),
    databases => ['rnaseq_db', 'rnaseq_refine_db', 'rnaseq_blast_db', 'dna_db'],
    rnaseq_db => $self->o('rnaseq_db'),
    rnaseq_refine_db => $self->o('rnaseq_refine_db'),
    rnaseq_blast_db => $self->o('rnaseq_blast_db'),
    dna_db => $self->o('dna_db'),
    enscode_root_dir => $self->o('enscode_root_dir'),
    extra_parameters => {
      output_path => $self->o('output_path'),
      user_r => $self->o('user_r'),
      dna_db_host => $self->o('dna_db_host'),
      dna_db_port => $self->o('dna_db_port'),
      databases_host => $self->o('databases_host'),
      databases_port => $self->o('databases_port'),
      release_number => $self->o('release_number'),
      production_name => $self->o('production_name'),
      species_name => $self->o('species_name'),
      registry_host => $self->o('registry_host'),
      registry_port => $self->o('registry_port'),
      registry_db => $self->o('registry_db'),
      assembly_name => $self->o('assembly_name'),
      assembly_accession => $self->o('assembly_accession'),
    },
  },
  -rc_name      => 'default',
  -max_retry_count => 0,
  -flow_into => {
    1 => ['run_rnaseq_db'],
  },
},
```

### run\_\<sub pipeline name\>

This analysis will run the sub pipeline

It must be a `Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun`, it will run beekeeper in a loop. If the job is a retry, it will first reset all failed jobs in the sub pipeline and start beekeeper.

There is only one parameter to be passed, `beekeeper_script`, it is not required but recommended.

**Example analysis config**
```
{
  -logic_name => 'run_rnaseq_db',
  -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun',
  -parameters => {
    beekeeper_script => $self->o('hive_beekeeper_script'),
  },
  -rc_name      => 'default',
  -max_retry_count => 1,
},
```


## The different parts of the EnsEMBL Genome Annotation System

### The main pipeline

The main pipeline will query ENA to retrieve the possible short read and long read accessions which would be used in the corresponding sub pipelines. Then start each sub pipeline below in the order they appear in this document.

By default, ENA is queried using the NCBI taxonomy id of the species, but you can provide either a single project id or a list of project ids to `study_accession`

It will create a registry file which will be used for the whole genome alignment and the DataChecks

It will provide stats on the assembly and on the repeat masking which will be sent to the email provided.

Before starting the alignments sub pipelines, it will reset a set of arrays in the "Transcript selection" sub pipeline to allow any related sub pipeline to provide the database connection details some analyses will need.

### Loading the assembly

#### What it does

It will create the core database which is referenced as `dna_db`, `reference_db` and sometimes `core_db` in configuration files. It will load a set of static tables. It will process the assembly report file to load the sequence names and synonyms and the dna. If UCSC sequence accessions exist in the report, they will be loaded as sequence synonyms

When there is a known mitochondrial sequence, it will load the mitochondrial dna and genes. Finally, it will dump the genome in a fasta file and index it for faidx access.

#### Notifications

None

#### Caveats

If RefSeq has not created a sister assembly (GCF\_\*), we will not load any mitochondrial data even if a RefSeq mitochondrial annotation exists. Because the mitochodrial sequence is not used in the annotation process it can be checked before the gene set is released with a query on the NCBI website, searching for a mitocondrial annotation from RefSeq for the species of interest.

### RefSeq annotation import

#### What it does

It checks on the assembly report if a corresponding RefSeq assembly exists. It will then load the RefSeq sequence accessions as sequence synonyms and the RefSeq gene set using their GFF3 annotation file

The gene set is not used for annotation purpose. It can be used to compare the gene set generated and Ensmbl users appreciate the possibilty of looking at both gene sets in one location.

#### Notifications

None

#### Caveats

None

### Repeat masking the genome

#### What it does

It will run RepeatMasker using RepBase with the closest clade library and it will run dustmasker and TRF.

It will verify the presence of a RepeatModeler library file and run RepeatMasker with this library if needed.

It will run Red, a different repeat masking program if `replace_repbase_with_red_to_mask` is set to 1.

#### Notifications

None

#### Caveats

Red and RepeatMasker with a RepeatModeler library may mask the genome where protein gene families could be. However in the absence of a RepBase entry for the species of interest, they will be helpful.

### Whole genome alignment against a high quality assembly with LastZ

#### What it does

It uses an Ensembl Compara pipeline to do the whole genome alignment of your species genome against a high quality assembly. We usually use human GRCh38 for all species except the rodents where we would use mouse GRCm39. You can choose a different species for the source. You will need to check the quality of the assembly and the quality of the annotation.

You will need to use a Compara master database.

It is possible to avoid running this part by setting `skip_projection` to 1.

The results of this sub pipeline will be used for the projection sub pipeline.

#### Notifications

None

#### Caveats

The pipeline provides the best result when the two species are not too divergent. For example, it does not provide good information when used with any fish if human is the source.

### Projecting a high quality annotation on the species of interest

#### What it does

It will use the whole genome alignment of the LastZ sub pipeline to project the high quality source annotation onto the genome of interest. It will also use CESAR2.0 to project the high quality source annotation onto the genome of interest. Then we select the best model for each projected gene based on coverage and identity. The selected model can be used later in the transcript selection sub pipeline.

#### Notifications

* proj\_sel database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

#### Caveats

None

### Aligning proteins from related species

#### What it does

We align a controlled set of proteins grouped by species/family/clade depending on the clade of the species of interest using GenBlast. The models are classified based on their coverage and identity. The models with the highest coverage and the highest identity will be preferably used when selecting the transcripts.

The controlled sets are found in modules/Bio/EnsEMBL/Analysis/Hive/Config/UniProtCladeDownloadStatic.pm

#### Notifications

* genblast\_nr database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

#### Caveats

None

### Aligning proteins and cDNAs from the species of interest

#### What it does

We align the cDNAs specific to the species onto the genome using Exonerate. First we align all possible sequences, they will be used later to add UTRs. Then we only align the full length cDNAs using the first alignments to reduce the search space and using the reported CDS start and CDS end. When internal stops are found, the model is removed unless we found only one. In this case we replace the stop codon with an intron to avoid any assembly error to interfere.

We also align the species specific proteins retrieved from UniProt with protein existence (PE) level 1 and 2 (protein evidence and transcript evidence). We use PMatch to reduce the search space then Exonerate and GeneWise in two different analyses. We ran a final Exonerate on the whole genome as back up.

We align the species specific seleno proteins and mark the internal stop with an attribute

The final step will select the best model for each protein of a locus with the following ranking:
* cdna
* seleno protein
* edited cdna
* protein
* backup protein

#### Notifications

* bt database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

* cdna database to transcript selection sub pipeline
  - run\_utr\_addition
  - run\_utr\_addition\_10GB
  - run\_utr\_addition\_30GB
  - utr\_memory\_failover

#### Caveats

None

### Creating an IG/TR annotation

#### What it does

It uses a known set of human IG/TR genes and align them using GenBlast with customised parameters to restrict the size of the introns. It collapses the models created at each locus

#### Notifications

* igtr database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

#### Caveats

None

### Generating a short non coding gene set

#### What it does

It aligns sequence from the RFam database to the genome using Blast and uses the Infernal software suite to assess the short non coding genes.

Sequences from miRBase are aligned to the genome and are filtered using a 

#### Notifications

None

#### Caveats

The `genebuild-mirna` virtual environment is needed but do not need to be activacted as we use the full path to Python.

### Short read alignments

#### What it does

We use STAR to align the set of Illumina short reads downloaded. Then we generate models with Scallop based on the STAR alignments. We determine the coding potential of the models using blast on a sub set of UniProt proteins, PE 1 and 2. For non vertebrate we also use CPC2 and Samba as UniProt does not have enough data. The models are collapsed by loci for the transcript selection pipeline.

#### Notifications

* rnalayer\_nr database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

* rnalayer\_nr database to transcript selection sub pipeline
  - run\_utr\_addition
  - run\_utr\_addition\_10GB
  - run\_utr\_addition\_30GB
  - utr\_memory\_failover

* rnalayer\_nr database to homology rnaseq sub pipeline
  - genblast\_rnaseq\_support
  - genblast\_rnaseq\_support\_himem

* scallop\_blast database to main pipeline
  - initialise\_rnaseq

#### Caveats

The fastq files downloaded will be deleted to free disk space as it can easily be above the terabyte

### Homology intron check

#### What it does

It will use all the splice sites found by the short read alignment step and compare them with the introns found in the protein models from GenBlast. The models will be ranked based on the completeness of the intron support.

#### Notifications

* gb\_rnaseq\_nr database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

#### Caveats

It is not run when there is no short read data

### Long read alignments

#### What it does

We use Minimap2 to align the long reads to the genome and the resulting models are collapsed to help reduce the number of isoforms when there is a slight difference in final exons for example. It also try to find the correct canonical splice sites using the redundancy. The default data type is PacBio but ONT can be used too.

#### Notifications

* lrfinal database to transcript selection sub pipeline
  - create\_toplevel\_slices
  - split\_slices\_on\_intergenic
  - layer\_annotation

* lrfinal database to transcript selection sub pipeline
  - run\_utr\_addition
  - run\_utr\_addition\_10GB
  - run\_utr\_addition\_30GB
  - utr\_memory\_failover

#### Caveats

The fastq files downloaded will be deleted to free disk space as it can easily be above the terabyte

### Transcript selection

#### What it does

It will look at all the protein coding data generated by sub pipelines and using a layer mechanism it will decide which model to use for each loci. It will also add the UTR regions to the models based on the transcriptomic data sets; cDNA, RNA-seq and long reads. It will trim the UTR regions when necessary.

It will look for pseudogenes. The criterion are single exon models which have a good blast hit with a multi exon model or if the model is within an LTR region. Multi exon models with very short introns are also seen as pseudogenes.

Long non coding RNA will be labelled but we need transcriptomic data.

Short non coding RNA will be copied over unless they overlap protein coding genes.

Some cleaning is done to remove low quality models which haven't been labelled as pseudogenes or lncRNA or when more information was needed about the gene set.

#### Notifications

None

#### Caveats

None

### Gene set finalisation

#### What it does

We copy the final gene set into the core database while making a last check for readthroughs and bad UTRs.

The stable id needs to be generated, either with the simple script when it's a new species or with the stable id mapping pipeline when it is an update to the assembly. If it is an update to the assembly, the `mapping_required` parameter of run\_stable\_ids should be set to 1. Once the mapping has been done you can add a `skip_analysis` parameter to the same analysis and set it to 1 to restart the pipeline.

Then the database is prepared to be ready for handover.

#### Notifications

None

#### Caveats

The external db ids of the supporting evidences are set using a script which when it breaks can brak many things. So if the analysis load\_external\_db\_ids\_and\_optimise\_af\_tables fails, you have to be very careful rerunning the job.

The keys of the meta table are set for RapidRelease. If you want to handover for the Main website, you will need to manually do some tweaks

### \_otherfeatures\_ database creation

#### What it does

It will create an empty database from the core database and copy the genes from the following databases:
* cdna
* refseq
* lrinitial

The database is then prepared for handover.

#### Notifications

None

#### Caveats

When this pipeline is started after a successful core handover and the hopefully unlikely modification have been done on the database from the pipeline, the meta table should not need any updates

### \_rnaseq\_ database creation

#### What it does

It will copy the database provided by the short read pipeline and sort the dna\_align\_feature table which contains the introns. This is necessary to speed up the MySQL queries as the table contains millions of rows.
It will the assign the external db id for protein evidences and make sure the data\_file table which contains the filename of the BAM/BigWig files. It will upload the web\_data and analyses descriptions to the Production database for the configuration matrix to work as expected on the website.
It will generate BigWig files from the tissue BAM files and generate a README and a md5sum file for the FTP. And finally prepare the database for handover.

#### Notifications

None

#### Caveats

When this pipeline is started after a successful core handover and the hopefully unlikely modification have been done on the database from the pipeline, the meta table should not need any updates
