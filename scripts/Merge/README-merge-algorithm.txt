$Id: README-algorithm.txt,v 1.1 2014-01-22 14:27:24 ak4 Exp $

A cluster of overlapping Ensembl and Havana genes is considered.  Each
Havana gene is tested against all Ensembl genes on transcript level.

The process of deciding what models to merge and what models to copy
does not change any data.  It only creates two lists of candidates for
merging and copying.

For each Ensembl gene:
  For each Ensembl transcript:

    For each Havana gene:

      The Havana gene might be on the opposite strand from the Ensembl
      gene.  In that case skip the comparison between these two genes
      completely.

      For each Havana transcript:

        If the Havana transcript is a single exon transcript:
          Compare the single Havana exon with the set of exons from the
          Ensembl transcript.

          The comparison yields true if both transcripts are single exon
          transcripts with identical start and stop coordinates, but if
          that is not the case, and both transcripts are coding, the
          coding regions are compared in the same way.

          Special case:  If either single exon transcript is a stop
          codon longer at the end, then the comparison yields true.

        Else (it is multi-exon):
          Compare the set of Havana introns with the set of introns from
          the Ensembl transcript.

          The comparison yields true if the two sets are identical in
          number of introns and if all introns have pairwise identical
          start and stop coordinates.

        If the comparison in either branch above yields true:
          Save the Ensembl transcript as a candidate for merging into
          this particular Havana transcript.

        Else (intron structure not same):
          Examine the two sets of exons for any overlap.

          If there is overlap:
            Save the Ensembl transcript as a candidate for copying into
            this particular Havana gene.

            Special case: If the Ensembl gene is a single exon RNA gene
            and the Havana gene is a multi exon coding gene, the Ensembl
            model will not be copied.

We end up with a list of Ensembl transcripts to merge into Havana
transcripts ("merge candidates") and another list of Ensembl transcripts
to copy into Havana genes ("copy candidates").  At this point, an
Ensembl transcript may theoretically occur multiple times in each list.

Go through the list of "merge candidates" and merge these Ensembl
transcripts with the corresponding Havana transcript.  We allow for
merging one Ensembl transcript into multiple Havana transcripts.

The merging process has two special cases:

* If the Havana gene has a reference error, the Ensembl transcript will
  be copied instead.  In addition, the Havana gene biotype will be
  promoted to whatever the Ensembl transcript biotype is.

* If the Ensembl transcript model is a CCDS model and if the Havana
  transcript biotype is different from the Ensembl transcript biotype,
  the model will be copied instead.

Go through the list of "merge candidates" (again).  This time save
all un-merged "sibling" transcripts (transcripts of the same gene) as
candidates for copying into the Havana gene of the transcript that the
merged "sibling" was merged into.  This expands the list of Ensembl
transcripts to copy into Havana genes from above.

Since the list of "copy candidates" may contain the same Ensembl
transcript multiple times (and some of these may be conflicting in what
Havana gene the transcript should be copied into), we need to create a
new list of distinct and definite "copy candidates".

We go through the list of "copy candidates" and for each candidate (that
has not been already merged) we calculate the absolute exon overlap
(in bases) between the Ensembl transcript and the particular Havana
transcript that yielded the candidate.  If we find a candidate whose
Ensembl transcript has been seen before, we keep the candidate with the
largest overlap.  This also takes care of conflicts where an Ensembl
transcript has been found to overlap Havana transcripts in multiple
Havana genes.

Go through the list of definite "copy candidates" (considering any CCDS
models first since the copy procedure might modify the biotype of Havana
genes when encountering these, see below) and if the Ensembl transcript
has complete start and stop codons, we copy it into the Havana gene
found to have the greatest overlap.  If the Ensembl transcript has
incomplete start or stop codon, we save this "copy candidate" into a new
list of candidates that may possibly be ignored.

Go through the list of "copy candidates" to possibly ignore and see if
any other Ensembl transcript from the same Ensembl gene has previously
been either merged into a Havana transcript or copied into a Havana
gene.  If so, the "copy candidate" is finally ignored.  If not, it is
copied.  This may copy more than one Ensembl transcript with incomplete
start or stop codon from an Ensembl gene into a Havana gene if no other
Ensembl transcript has previously been merged or copied from that
particular Ensembl gene.

When copying, there are a few special cases for how biotypes and
translations are handled:

* If the Havana gene has a reference error and the Ensembl transcript is
  coding, the Havana gene biotype is promoted to the Ensembl transcript
  biotype.

* If the Havana gene is not coding but the Ensembl transcript is:

  - If the Ensembl transcript is a CCDS transcript, the Havana gene
    biotype is promoted to the Ensembl transcript biotype.

  - If the Ensembl transcript is not a CCDS transcript, the Ensembl
    translation is removed and the Ensembl biotype is demoted to the
    Havana gene biotype.

The cluster of Havana genes is finally stored into the output database
and the next cluster of Havana and Ensembl genes is considered.

As a post-processing step, the Ensembl genes that were not merged with
any Havana transcript or copied into any Havana gene are identified and
copied verbatim.
