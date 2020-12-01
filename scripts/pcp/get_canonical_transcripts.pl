use warnings;
use strict;

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $coord_system = 'toplevel';
my $dna_dbname = 'homo_sapiens_core_101_38';
my $dna_user   = 'ensro';
my $dna_host   = 'mysql-ens-mirror-1.ebi.ac.uk';
my $dna_port   = '4240';
my $dna_pass;

my $options = GetOptions ("user|dbuser|u=s"	 => \$dna_user,
                          "host|dbhost|h=s"	 => \$dna_host,
                          "port|dbport|P=i"	 => \$dna_port,
                          "dbname|db|D=s"    => \$dna_dbname,
                          "dbpass|pass|p=s" => \$dna_pass);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $dna_user,
  -host    => $dna_host,
  -dbname  => $dna_dbname,
  -pass    => $dna_pass);

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all($coord_system);
my $gene_adaptor = $db->get_GeneAdaptor();

for my $slice (@$slices) {
	my $genes = $slice->get_all_Genes();
	foreach my $gene (@$genes) {
		if ($gene->biotype eq "protein_coding") {
			my $transcript = $gene->canonical_transcript;
			print ">", $transcript->seq_region_name, ":", $transcript->seq_region_start, "-", $transcript->seq_region_end, "\n", $transcript->seq->seq, "\n";
		}
	}
}