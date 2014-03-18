$Id: README-practical.txt,v 1.3 2014-01-28 15:22:55 ak4 Exp $

This is a document outlining the practicalities around running the new
merge code.

The new merge code (I'll drop the word "new" from this from now on)
consists of three pieces:

  1)  merge.pl, the main script.
      This is the actual Perl script that does performs the merge.

  2)  merge-wrapper.ksh, the convenience wrapper.
      This is a KornShell script (ksh, a variant of bash) that reads the
      configuration (see below) and submits the merge.pl script to the
      farm as a job array.  It also submits a post-processing job and a
      cleanup job (more about this later).

  3)  merge.conf, the configuration.
      This is a file that has bash syntax (no spaces around '='
      etc.) and that defines a number of variables that makes up the
      configuration of the merge (database names, user names, passwords
      etc.).  The wrapper script will source it.

I will now write about these file individually in more detail, and then
I will show how to run the code.

------------------------------------------------------------------------
THE MAIN SCRIPT
------------------------------------------------------------------------

You should hopefully not have to run the main script directly, ever.
It has far too many possible command line options.  This is what the
wrapper script is for.  It may however be useful to have a look at the
various options, so run the script with the "--help" option at least
once.

The main script implements the algorithm described in the file
"README-algorithm.txt", more or less (you know how these things are,
documentation never stays up to date with the code, but I've done my
best to insure that it's valid).

It performs the merge between the Havana gene set and the Ensembl gene
set, or whatever takes the place of these two.  You might, for example,
at some point want to run the merge between RefSeq and Ensembl in which
case all textual references to "Havana" obviously should be changed to
"RefSeq" below.

The main thing to remember is that the Havana gene set will always,
and without exception, be copied from the Havana input database to the
output database in full.  There will never be any Havana genes not
copied.

The Ensembl gene set, on the other hand, will be incorporated (merged
or copied) into the Havana set depending on overlap with Havana genes.
This means that some Ensembl genes probably won't get copied into the
output database.  This can happen if the Ensembl gene isn't overlapping
with any Havana gene, or if the algorithm decides that the Ensembl model
should be deleted for one reason or other.  Ensembl genes that doesn't
overlap any Havana genes needs to be copied separately, outside of the
main script.  This is what the post-processing job that the wrapper
script is submitting is doing (see below).

Apart from supplementing Havana genes and transcripts with Ensembl
annotation, the script also adds xrefs to all written genes, transcripts
and translations.  These are simply the original Havana stable IDs and
makes it easy for Havana to trace back to their identifiers once we've
run stable ID mapping (which should happen at some point after running
the merge).

The other modification that may go unnoticed unless I mention it
here is the tagging of analysis logic names.  The logic names of all
supporting features (exon, transcript and intron support) as well
as the logic names associated with the genes, transcripts, exons
and translations will be tagged with either "_ensembl" or "_havana"
depending on origin.  These tags may be changed using the --ensembl_tag
and --havana_tag options.  These options also affect the "source" of
Ensembl and Havana genes and transcripts and the source set for merged
genes and transcripts (these will get a source value of the two tags
concatenated with "_", e.g. "ensembl_havana").

------------------------------------------------------------------------
THE CONVENIENT WRAPPER SCRIPT
------------------------------------------------------------------------

The wrapper script will read the configuration file and invoke the main
script as an LSF job array of configurable length with the options
supplied therein.

Since there always will be the possibility that there are Ensembl genes
that are not overlapping with any Havana gene, these will have to be
copied separately.  The wrapper script will do this by submitting a
post-processing script that will parse the log files of the main script,
looking for the gene dbIDs of the processed Ensembl genes.  From that it
figures out what genes have not been processed and goes off to extract
the gene dbIDs of those.  It uses the standard copy_genes.pl script to
copy these unprocessed genes in parallel to the output database.

Note that these gene are copied verbatim (as-is) and thus will not get
the tagging (of logic names etc.) that other genes will get.

There is a third job that gets submitted by the wrapper script.  It is
a small cleanup job that I added because I noticed that the job that
copies the gene still ran even if the main job array failed or was
killed.  It simply kills the copy job if the job array fails.

There is one problem with the wrapper script and its job logic and that
is that if the main array fails (due to it exceeding memory resources
etc.) the two other jobs, the copy job and the cleanup job, will still
hang around and will have to be killed manually.  If these are the only
two jobs pending, then you may expect to find that one or several of the
job array jobs were terminated by LSF.

------------------------------------------------------------------------
THE CONFIGURATION FILE
------------------------------------------------------------------------

The configuration file is a shell script that simply assigns values to
a set of variables.  These variables are later passed on by the wrapper
script to the main Perl script by means of command line options.

The configuration file is a bash-like shell script, so bash syntax rules
needs to be followed.

Have a look at the merge.conf file and the wrapper script,
merge-wrapper.ksh, to see how the values of the variables are
transferred to the main script, merge.pl.

------------------------------------------------------------------------
HOW TO RUN THE MERGE
------------------------------------------------------------------------

First make a copy of the configuration file (or modify it directly, it's
up to you) and fill it out.  The options are well documented in the
configuration file, and even more so in the documentation provided via
the --help option to the main Perl script.

Make sure that the output database exists and is in the correct state.
It should be an empty clone of the Ensembl database (or Havana database,
at this stage in the process, the "essential tables" should have been
synchronised already).  You may use the clone_database.ksh script to
create a suitable empty database.  This script lives in

  ensembl-personal/genebuilders/scripts/

Run it with the -h option to get usage information.

The configuration file contains an option to name an output directory.
This is the directory will be used for logging output.  The merge
process does not dump any tables, so the size should hopefully stay
limited.  For the recent merge of zebrafish (for release 75), the
directory ended up using 20M, for a human test merge, it ended up at
between 70M and 95M.

The output directory should not exist already, so don't create it.
Delete it or rename it if it exists from an old run.

Then we set the code running:

  ./merge-wrapper.ksh ./merge.conf

Notice that you will have to give a full path (not necessarily an
absolute path though) to the configuration file.

You may then follow the progress of the merge by means of bjobs and by
looking at the log files written into the output directory.

If you need to kill the job, use bkill as usual.  Make sure to remove or
rename the output directory and to re-clone the output database before
running the wrapper script again.  Output is written to the output
database pretty much immediately, so you probably shouldn't assume that
you were quick enough to kill the job before any output had been written
if you accidentally set it off running.  The clone_database.ksh has a
useful -f flag that will drop an already existing database if you need
to re-clone the output database.
