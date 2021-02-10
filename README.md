# Ensembl-analysis

[![Build Status](https://travis-ci.com/Ensembl/ensembl-analysis.svg?branch=dev/hive_master)][travis]

[travis]: https://travis-ci.com/Ensembl/ensembl-analysis

The Ensembl Annotation System provide all the Perl modules needed to annotate vertebrate genomes.


## Download

To clone the Ensembl Annotation System, use the following command:

```
git clone https://github.com/Ensembl/ensembl-analysis.git
```


## API requirements

In order to use the Ensembl Annotation System, an installation of [BioPerl 1.6.924 core modules](https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip) (bioperl-live) is required.

You will also need the following repositories:
* https://github.com/Ensembl/ensembl
* https://github.com/Ensembl/ensembl-ehive
* https://github.com/Ensembl/ensembl-compara
* https://github.com/Ensembl/ensembl-variation
* https://github.com/Ensembl/ensembl-io
* https://github.com/Ensembl/ensembl-taxonomy
* https://github.com/Ensembl/ensembl-killist
* https://github.com/Ensembl/ensembl-production
* https://github.com/Ensembl/ensembl-datacheck
* https://github.com/Ensembl/ensembl-orm

A guide for installing all Ensembl APIs and their respective prerequisites is available here:
http://www.ensembl.org/info/docs/api/api_installation.html

Once all dependencies are installed, you should run the following command to install Perl modules:

```
cpanm --installdeps --with-recommends --notest
```


## Softwares

In order to run our pipeline you will need to install the following softwares:
* exonerate 0.9.0
* geneWise
* GenBlast
* NCBI Blast
* BWA
* samtools
* minimap2
* RepeatMasker
* RepeatModeler
* RepeatDetector
* dust
* TRF
* lastZ
* Cesar2.0
* ensc-annotation-tools
* infernal
* Genscan
* tRNAscan-SE


## Contributions

If you wish to contribute to this repository or any Ensembl repository, please refer to [our contribution guide](https://github.com/Ensembl/ensembl/blob/master/CONTRIBUTING.md).


## Contact us
Please email comments or questions to the public Ensembl developers list at

<http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at

<http://www.ensembl.org/Help/Contact>.
