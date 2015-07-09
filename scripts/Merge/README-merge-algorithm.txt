$Id: README-algorithm.txt,v 1.1 2014-01-22 14:27:24 ak4 Exp $

A cluster of overlapping Secondary (usually Ensembl) and Primary (usually Havana) genes is considered.  Each
Primary gene is tested against all Secondary genes on transcript level.

The process of deciding what models to merge and what models to copy
does not change any data.  It only creates two lists of candidates for
merging and copying.

For each Secondary gene:
  For each Secondary transcript:

    For each Primary gene:

      The Primary gene might be on the opposite strand from the Secondary
      gene.  In that case skip the comparison between these two genes
      completely.

      For each Primary transcript:

        If the Primary transcript has a genome patch truncated attribute:
          # "REVERSE AND PARTIAL" MERGE
          Swap Primary and Secondary transcripts.
          
          If the exons of the Secondary transcript are a subset of the
          exons of the Primary transcript (start and end are allowed to differ):
            Save the Secondary transcript as a candidate for merging into
            this particular Primary transcript.
          
          Else (exon structure not same):
            Examine the two sets of exons for any overlap.

            If there is overlap:
              Save the Secondary transcript as a candidate for copying into
              this particular Primary gene.

              Special case: If the Secondary gene is a single exon RNA gene
              and the Primary gene is a multi exon coding gene, the Secondary
              model will not be copied.

        Elsif the Primary transcript is a single exon transcript:
          # "NORMAL" MERGE
          Compare the single Primary exon with the set of exons from the
          Secondary transcript.

          The comparison yields true if both transcripts are single exon
          transcripts with identical start and stop coordinates, but if
          that is not the case, and both transcripts are coding, the
          coding regions are compared in the same way.

          Special case:  If either single exon transcript is a stop
          codon longer at the end, then the comparison yields true.

        Else (it is multi-exon):
          Compare the set of Primary introns with the set of introns from
          the Secondary transcript.

          The comparison yields true if the two sets are identical in
          number of introns and if all introns have pairwise identical
          start and stop coordinates.

        If the comparison in either branch above yields true:
          Save the Secondary transcript as a candidate for merging into
          this particular Primary transcript.

        Else (intron structure not same):
          Examine the two sets of exons for any overlap.

          If there is overlap:
            Save the Secondary transcript as a candidate for copying into
            this particular Primary gene.

            Special case: If the Secondary gene is a single exon RNA gene
            and the Primary gene is a multi exon coding gene, the Secondary
            model will not be copied.

We end up with a list of Secondary transcripts to merge into Primary
transcripts ("merge candidates") and another list of Secondary transcripts
to copy into Primary genes ("copy candidates").  At this point, an
Secondary transcript may theoretically occur multiple times in each list.

Go through the list of "merge candidates" and merge these Secondary
transcripts with the corresponding Primary transcript.  We allow for
merging one Secondary transcript into multiple Primary transcripts.

The merging process has two special cases:

* If the Primary gene has a reference error, the Secondary transcript will
  be copied instead.  In addition, the Primary gene biotype will be
  promoted to whatever the Secondary transcript biotype is.

* If the Secondary transcript model is a CCDS model and if the Primary
  transcript biotype is different from the Secondary transcript biotype,
  the model will be copied instead.

Go through the list of "merge candidates" (again).  This time save
all un-merged "sibling" transcripts (transcripts of the same gene) as
candidates for copying into the Primary gene of the transcript that the
merged "sibling" was merged into.  This expands the list of Secondary
transcripts to copy into Primary genes from above.

Since the list of "copy candidates" may contain the same Secondary
transcript multiple times (and some of these may be conflicting in what
Primary gene the transcript should be copied into), we need to create a
new list of distinct and definite "copy candidates".

We go through the list of "copy candidates" and for each candidate (that
has not been already merged) we calculate the absolute exon overlap
(in bases) between the Secondary transcript and the particular Primary
transcript that yielded the candidate.  If we find a candidate whose
Secondary transcript has been seen before, we keep the candidate with the
largest overlap.  This also takes care of conflicts where an Secondary
transcript has been found to overlap Primary transcripts in multiple
Primary genes.

Go through the list of definite "copy candidates" (considering any CCDS
models first since the copy procedure might modify the biotype of Primary
genes when encountering these, see below) and if the Secondary transcript
has complete start and stop codons, we copy it into the Primary gene
found to have the greatest overlap.  If the Secondary transcript has
incomplete start or stop codon, we save this "copy candidate" into a new
list of candidates that may possibly be ignored.

Go through the list of "copy candidates" to possibly ignore and see if
any other Secondary transcript from the same Secondary gene has previously
been either merged into a Primary transcript or copied into a Primary
gene.  If so, the "copy candidate" is finally ignored.  If not, it is
copied.  This may copy more than one Secondary transcript with incomplete
start or stop codon from an Secondary gene into a Primary gene if no other
Secondary transcript has previously been merged or copied from that
particular Secondary gene.

When copying, there are a few special cases for how biotypes and
translations are handled:

* If the Primary gene has a reference error and the Secondary transcript is
  coding, the Primary gene biotype is promoted to the Secondary transcript
  biotype.

* If the Primary gene is not coding but the Secondary transcript is:

  - If the Secondary transcript is a CCDS transcript, the Primary gene
    biotype is promoted to the Secondary transcript biotype.

  - If the Secondary transcript is not a CCDS transcript, the Secondary
    translation is removed and the Secondary biotype is demoted to the
    Primary gene biotype.

The cluster of Primary genes is finally stored into the output database
and the next cluster of Primary and Secondary genes is considered.

As a post-processing step, the Secondary genes that were not merged with
any Primary transcript or copied into any Primary gene are identified and
copied verbatim.
