#!usr/bin/perl

use strict;
use warnings;
use Utilities;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;

print ("Uses cleared\n");

#open (CARS_T,'<DataFiles/t_ids_w_intropolis.txt') || die "Could not open file";#CARS_Refseq_TIDs.txt
open (CARS_SYMS,'<DataFiles/CARS_Gene_Syms_TIDs.txt') || die "Could not open file";
open (RS_SYMS,'DataFiles/RefSeq_Gene_Syms.txt') || die "Could not open file";
#open (RS_T,'DataFiles/t_ids_w_intropolis');
#open (INT_RS,'>DataFiles/RF_introns.txt');

my $dbname = 'homo_sapiens_core_94_38';
my $dbuser = 'ensro';
my $dbpass = '';
my $dbhost = 'mysql-ensembl-mirror.ebi.ac.uk';
my $dbport = '4240';

my $rfsdbname = 'refseq_human';
my $rfsdbuser = 'ensro';
my $rfsdbhost = 'mysql-ens-havana-prod-1.ebi.ac.uk';
my $rfsdbport = '4581';

my $coredb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $dbhost,
  -port => $dbport,
  -user => $dbuser,
  -dbname => $dbname,
  -pass => $dbpass,
);

my $rfsdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $rfsdbhost,
  -port => $rfsdbport,
  -user => $rfsdbuser,
  -dbname => $rfsdbname,
  -pass => $dbpass,
);

print ("Connections to DBs established ...\n");

my $rsa = $rfsdb->get_SliceAdaptor();
my $rsta = $rfsdb->get_TranscriptAdaptor();
my $ta = $coredb->get_TranscriptAdaptor();

print ("Adaptors retreived ...\n");

my %coding_exons;
my $count = 0;
my @coding_transcripts;
my @allowed_biotypes = ('protein_coding','TR_V_gene','TR_J_gene','TR_D_gene','TR_C_gene','IG_D_gene','IG_J_gene','IG_V_gene','IG_C_gene');

# my $slices = $rsa->fetch_all('toplevel');
# print "\n***slice", @$slices[1],"\n";
# my $slice_count = scalar(@$slices);
# my $processing_count = 0;
#
# foreach my $slice (@$slices) {
#   if (index($slice->name,"chromosome:GRCh38:22")!= -1){
#     $processing_count++;
#     print "Processing slice: ".$slice->name." (".$processing_count."/".$slice_count.")";
#     my $genes = $slice->get_all_Genes();
#     foreach my $generef (@$genes) {
#       my $gene = $generef->stable_id();
#       #print $generef->external_name(),"\n";
#       $gene_count++;
#       # print "Gene Stable ID :", $generef->stable_id(), "\n";
#       my @transcripts = @{$generef->get_all_Transcripts()};
#
#       foreach my $element (@transcripts) {
#         if (grep $_ eq ($element->biotype), @allowed_biotypes) { # check the transcript has a biotype in the allowed list of biotypes
#           push(@coding_transcripts, $element);
#         }
#       }
#
#       unless ($coding_exons{$gene}) {
#         #print "Gene Stable ID :", $gene, "\n";
#
#         foreach my $transcript (@coding_transcripts) {
#           #my $trans = $transcript->stable_id;
#           my @introns = @{$transcript->get_all_Introns()};
#
#           print ("Introns retreived\n");
#
#           my $totalIntrons = scalar (@introns);
#           print ("Number of introns retreived : ", $totalIntrons ,"\n");
#
#           for ( @introns){
#               my $intronStart = $_->seq_region_start()-1;
#               #print ("Intron start : ", $intronStart, "\n");
#               print INT_RS $intronStart, "\t";
#
#               my $intronEnd = $_->seq_region_end();
#               #print ("Intron end : ", $intronEnd, "\n");
#               print INT_RS $intronEnd, "\n";
#           }
#         }
#       }
#     }
#   }
# }
my $transcript;
#****************************!# 
# while (<CARS_T>){
#      my $t_stableID = $_;
#      chomp $t_stableID;
#      $transcript = $ta->fetch_by_stable_id($t_stableID);
#      print $transcript->get_Gene()->external_name(), "\t";
#      print "$t_stableID\n";
    
     # my @introns = @{$transcript->get_all_Introns()};
#
#      print ("Introns retreived ... \n");
#
#      my $totalIntrons = scalar (@introns);
#      print ("Number of introns retreived : ", $totalIntrons ,"\n");
#
#      for ( @introns){
#        my $intronStart = $_->seq_region_start()-1;
#        #print ("Intron start : ", $intronStart, "\n");
#        print INT_RS $intronStart, "\t";
#
#         my $intronEnd = $_->seq_region_end();
#         #print ("Intron end : ", $intronEnd, "\n");
#         print INT_RS $intr onEnd, "\n";
#     }
# }
#    print ("Transcript belongs to chromosome : ", $transcript->seq_region_name(),"\n");
#
#    print ("transcript created\nTranscript start : ", $transcript->seq_region_start());

#****************************!# 
my %carsHash;
my $h_key;
my $h_value;
my @line;

while (<CARS_SYMS>){
    @line = split(/\s+/,$_);
    $h_key = $line[0];
    $h_value = $line[1];

    $carsHash{$h_key} = $h_value;
}

close (CARS_SYMS);

my $c = 0;
my $rs_trans;
#my %rsHash;
my @cars_m_trans;
my @rs_m_trans;

while (<RS_SYMS>){
    chomp $_;
    @line = split(/\s+/,$_);
    #$rsHash{$line[0]} = $line[1];
    if (exists $carsHash{$line[0]}){
        # print ++$c, ". ", $line[0], "\n";
        
        $transcript = $ta->fetch_by_stable_id($carsHash{$line[0]});
        my @introns = @{$transcript->get_all_Introns()};
        
        $rs_trans = $rsta->fetch_by_stable_id($line[1]);
        my @intronsRef = @{$rs_trans->get_all_Introns()};
        
        # print "size of ref introns: ", scalar(@intronsRef), "\n";
        #foreach my $i1, $i2 (@introns, @intronsRef)
        my $s1 = scalar(@introns);
        my $s2 = scalar(@intronsRef);
        my $max = ($s1 < $s2) ? $s1 : $s2;
        #print $s1,"\t",$s2,"\n"; 
        if($s1==$s2){
            my $all_match = 0;
            for my $it (0.. $max-1){
                #if (($i1->seq_region_start()-1 == $i2->seq_region_start()-1) && ($i1->seq_region_end()-1 == $i2->seq_region_end()-1))
                #print "$it . ref seq introns : $intronsRef[$it]\n";
                if (($introns[$it]->seq_region_start()-1 == $intronsRef[$it]->seq_region_start()-1) && ($introns[$it]->seq_region_end()== $intronsRef[$it]->seq_region_end())){
                    # print "Cars trans : $carsHash{$line[0]}\nRefSeq Trans : $line[1]\n";
                    # print $introns[$it]->seq_region_start(), "\n";
                    # print $intronsRef[$it]->seq_region_start(), "\n\n";
                    # print $introns[$it]->seq_region_end(), "\n";
                    # print $intronsRef[$it]->seq_region_end(), "\n\n";
                    $all_match++;
                }
            }
            if ($all_match == $s1 && $s1 >= 1){
                push (@cars_m_trans, $carsHash{$line[0]});
                push (@rs_m_trans, $line[1]);
            }
        }
    }
}

my @final_cars_trans = uniq(@cars_m_trans);
my @final_rs_trans = uniq(@rs_m_trans);

print scalar(@final_cars_trans),"\n", scalar(@final_rs_trans), "\n";

for (my $it = 0; $it<scalar(@final_cars_trans); $it++ ){
   print $final_cars_trans[$it],"\t", $final_rs_trans[$it],"\n";
}

close (RS_SYMS);
#close (INT_RS);

exit;

