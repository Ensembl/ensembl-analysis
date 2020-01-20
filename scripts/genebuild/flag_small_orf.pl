use warnings;
use strict;
use feature 'say';


use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $dbname = '';
my $user   = 'ensadmin';
my $host   = $ENV{GBS6};
my $port   = $ENV{GBP6};
my $pass   = 'ensembl';


my $dna_dbname = '';
my $dna_user   = 'ensadmin';
my $dna_host   = $ENV{GBS6};
my $dna_port   = $ENV{GBP6};
my $dna_pass   = 'ensembl';

my $options = GetOptions ("user|dbuser|u=s"      => \$user,
                          "host|dbhost|h=s"      => \$host,
                          "port|dbport|P=i"      => \$port,
                          "dbname|db|D=s"    => \$dbname,
                          "dbpass|pass|p=s" => \$pass,
                          "dnadbname=s" => \$dna_dbname,);

my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $dna_user,
  -host    => $dna_host,
  -dbname  => $dna_dbname,
  -pass    => $dna_pass);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

$db->dnadb($dna_db);

my $gene_adaptor = $db->get_GeneAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

my $genes = $gene_adaptor->fetch_all();

my $assembly_name = $db->get_MetaContainer->single_value_by_key('assembly.default');

my $small_orf_cutoff = 100;
my $small_intron = 75;
#say "Looping through genes...";
foreach my $gene (@$genes) {
#  say $gene->slice->name;
#  die;
  unless($gene->biotype eq 'protein_coding') {
    next;
  }

  my $single_exon_only = 1;
  my $transcripts = $gene->get_all_Transcripts;

  foreach my $transcript (@{$transcripts}) {
    my $introns = $transcript->get_all_Introns;
    if(scalar(@$introns) > 1) {
      $single_exon_only = 0;
      last;
    }

    my $intron = shift(@$introns);
    if($intron && $intron->length >= $small_intron) {
      $single_exon_only = 0;
      last;
    } #elsif($intron && $intron->length < $small_intron) {
      #say "Found a potential pseudogene: ".$gene->start.":".$gene->end.":".$gene->strand;
      #my $cds_start = $transcript->coding_region_start;
      #my $cds_end = $transcript->coding_region_end;
      #my $cds_slice_name = "primary_assembly:".$assembly_name.":".$transcript->seq_region_name.":".$cds_start.":".$cds_end.":1";
      #say $cds_slice_name;
      #my $slice = $slice_adaptor->fetch_by_name($cds_slice_name);
      #my $seq =  $slice->seq();
      #if($gene->strand == -1) {
      #  $seq = reverse($seq);
      #  $seq =~ tr/atgcATGC/tacgTACG/;
      #}
      #say translate($seq);
      #say "Original translation: ";
      #say $transcript->translate->seq;
    #}

    if($transcript->translation->length >= $small_orf_cutoff) {
      $single_exon_only = 0;
      last;
    }
  }

  if($single_exon_only) {
#    $gene->biotype('small_orf');
    $gene_adaptor->remove($gene);
#    my $transcript =  shift(@$transcripts);
#    say ">".$transcript->dbID;
#    say $transcript->translate->seq;
  }
}

exit;


sub translate {
  my ($seq) = @_;

  my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=> 'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC' =>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

  my $protein='';
  my $codon;
  for(my $i=0; $i<(length($seq)-2); $i+=3) {
    $codon = substr($seq,$i,3);
    $protein .= $g{$codon};
  }

  return($protein);
}
