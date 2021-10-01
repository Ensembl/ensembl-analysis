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
| ensembl-variation | default | https://github.com/Ensembl/ensembl-variation.git |
| ensembl-production | default | https://github.com/Ensembl/ensembl-production.git |
| ensembl-taxonomy | default | https://github.com/Ensembl/ensembl-taxonomy.git |
| ensembl-orm | default | https://github.com/Ensembl/ensembl-orm.git |
| ensembl-killlist | default | https://github.com/Ensembl/ensembl-killlist.git |
| ensembl-datacheck | default | https://github.com/Ensembl/ensembl-datacheck.git |
| ensembl-metadata | default | https://github.com/Ensembl/ensembl-metadata.git |
| ensembl-io | default | https://github.com/Ensembl/ensembl-io.git |

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
- `genebuild` using the requirements.txt file
- `genebuild-mirna` using the requirements\_p36\_ncrna.txt file

### Shell environment
If you are not part of the Ensembl Genebuild team, you will need to set some shell environment variables to avoid having to provide the information to the configuration files. We will assume you are using your home directory
| Variable | Value | Hive configuration parameter | Description |
|----------|-------|------------------------------|-------------|
| ENSCODE | $HOME | -enscode\_root\_dir | Directory path where you cloned all the Perl repositories |
| ENSEMBL\_SOFTWARE\_HOME | $HOME | -software\_base\_path | Directory where pyenv, plenv and linuxbrew are installed |
| LINUXBREW\_HOME | $HOME/.linuxbrew | -linuxbrew\_home\_path | Base directory for your Linuxbrew installation |
| PYTHONPATH | $HOME/ensembl-genes/ensembl\_genes | It needs to be set until the package can be install properly |
| BLASTDB\_DIR | $HOME | It will be used to find the path to the entry\_loc file which is the list of accession from SwissProt which would be located at $HOME/uniprot/2021\_03/entry\_loc |

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

#### Filling the main configuration manually
You would need to edit `$ENSCODE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/Genome_annotation_conf.pm` and fill in any information according to you environment

Then you can run
```
perl $ENSCODE/ensembl-hive/scripts/init_Pipeline.pl Bio::EnsEMBL::Analysis::Hive::Config::Genome_annotation_conf [extra parameters]
```

#### Filling the main configuration automatically
If you are operating within an environment prepared for Ensembl with the assembly registry you can use the `$ENSCODE/ensembl-analysis/scripts/genebuild/create_annotation_configs.pl`.

You would need to edit `$ENSCODE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/genome_annotation.ini`

Then you can run
```
perl $ENSCODE/ensembl-analysis/scripts/genebuild/create_annotation_configs.pl --config_file $ENSCODE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/genome_annotation.ini --assembly_registry_host <host_name> --assembly_registry_port <port>
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

## The different parts of the EnsEMBL Genome Annotation System
### The main pipeline
The main pipeline will query ENA to retrieve the possible short read and long read accessions which would be used in the corresponding sub pipelines. Then start each sub pipeline below in the order they appear in this document.

By default, ENA is queried using the NCBI taxonomy id of the species, but you can provide either a single project id or a list of project ids to `study_accession`

It will create a registry file which will be used for the whole genome alignment and the DataChecks

It will provide stats on the assembly and on the repeat masking which will be sent to the email provided.

### Loading the assembly

### RefSeq annotation import

### Repeat masking the genome

### Whole genome alignment against a high quality assembly with LastZ

### Projecting a high quality annotation on the species of interest

### Aligning proteins from related species

### Aligning proteins and cDNAs from the species of interest

### Creating an IG/TR annotation

### Generating a short non coding gene set

### Short read alignments

### Long read alignments

### Transcript selection

### Gene set finalisation

### \_otherfeatures\_ database creation

### \_rnaseq\_ database creation
