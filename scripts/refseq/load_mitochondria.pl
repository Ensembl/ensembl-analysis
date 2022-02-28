#!/usr/bin/env perl
#
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 Synopsis

load_mitochondria.pl

=head1 Description

Parses mitochondrial genes out of genbank format NC_ file  and writes them to the database specified.
Can also load the sequence of the chromosome, the corresponding assembly entries and the appropriate attribute types if needed.

=head1 Example command line

NOTE : Do NOT use the --contig flag to supply a name like **AY172335** or **AY172335.1.16299**.
The script will fail.

       perl ensembl-pipeline/scripts/DataConversion/mitochondria/load_mitochondria.pl \
        -dbhost HOST -dbuser USER -dbport 9999 -dbpass PASS\
        -dbname DATABASE_NAME  \
        -chromosome  MT          \
        -name        MT \
        -scaffold    NC_005089   \
        -clone       AY172335  \
        -toplevel    chromosome  \
        -gene_type   test\
        -trna_type   test-Mt-tRNA \
        -rrna_type   test-Mt-tRNA \
        -logic_name  gmap_ncbi  \
        -genbank_file      $BASE/downloads/sequence.gb


=head1 Config

All configuration is done through MitConf.pm

=cut

use warnings ;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::SeqRegionSynonym;
use Bio::EnsEMBL::IO::Parser::Genbank;

use Time::localtime;

use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace throw warning);

use constant {
    INFO_TYPE           => 'DIRECT',
    INFO_TEXT           => 'Imported from GeneBank',
    REFSEQ_PEP          => 'RefSeq_peptide',
    REFSEQ_XPEP         => 'RefSeq_peptide_predicted',
    NCBI_EXTERNAL_DB_ID => 700,
    CHECK_CODON_TABLE   => 2
};

my %EXTERNAL_DB = ( 'GeneID'               => 'EntrezGene',
                    'UniProtKB/Swiss-Prot' => 'Uniprot/SWISSPROT',
                  );
my $help;
my @genes;;
my $count;
#use scaffold or supercontig:
my $scaffold = 'scaffold';
my $contig = 'contig';
my $clone = 'clone';
my %opt;
my $path = undef;
my $non_interactive = 0;
my $MIT_DBHOST;
my $MIT_DBUSER;
my $MIT_DBPASS;
my $MIT_DBPORT;
my $MIT_DBNAME;
my $MIT_GENBANK_FILE='';
my $MIT_LOGIC_NAME;
my $MIT_SOURCE_NAME;
my $MIT_NAME;
my $MIT_TOPLEVEL;
my $MIT_GENE_TYPE;
my $MIT_TRNA_TYPE;
my $MIT_RRNA_TYPE;
my $MIT_DB_VERSION;
my $MIT_DEBUG;
my $MIT_DB_FILE;
my $MIT_SCAFFOLD_SEQNAME;
my $MIT_CLONE_SEQNAME;
my $MIT_CONTIG_SEQNAME;


 # options submitted with commandline override MitConf.pm

GetOptions( \%opt,
            '-h|help',
            'dbhost=s',
            'dbuser=s',
            'dbpass=s',
            'dbport=i',
            'dbname=s',
            'logic_name=s',     # analysis.logic_name for new created genes
            'contig=s',         # MIT_CONTIG_SEQNAME
            'chromosome=s',     # MIT_CHROMOSOME_SEQNAME
            'scaffold=s',       # MIT_SCAFFOLD_SEQNAME
            'clone=s',          # MIT_CLONE_SEQNAME
            'toplevel=s',       # MIT_TOPLEVEL
            'gene_type=s',
            'trna_type=s',
            'rrna_type=s',
            'source=s',         # Allow to set the gene.source by commandline option | MIT_SOURCE_NAME
            'name=s',
            'genbank_file=s',
            'download!',
            'noxref!',
            'path=s',
            'accession=s',
            'non_interactive!',
             );
            # or &usage();

if ( $opt{path}){
  print "you specify path: ", $opt{path},"\n";
  $path =  $opt{path};
}

if (    $opt{dbhost} && $opt{dbuser}
     && $opt{dbname} && $opt{dbpass}
     && $opt{dbport} ) {
  $MIT_DBHOST = $opt{dbhost};
  $MIT_DBUSER = $opt{dbuser};
  $MIT_DBPASS = $opt{dbpass};
  $MIT_DBPORT = $opt{dbport};
  $MIT_DBNAME = $opt{dbname};
}

$MIT_GENBANK_FILE = $opt{genbank_file} if $opt{genbank_file};
$MIT_LOGIC_NAME   = $opt{logic_name}   if $opt{logic_name};
$MIT_SOURCE_NAME  = $opt{source}       if $opt{source};
$MIT_NAME         = $opt{name}         if $opt{name};
$MIT_TOPLEVEL     = $opt{toplevel}     if $opt{toplevel};
$MIT_GENE_TYPE    = $opt{gene_type}    if $opt{gene_type};
$MIT_TRNA_TYPE    = $opt{trna_type}    if $opt{trna_type};
$MIT_RRNA_TYPE    = $opt{rrna_type}    if $opt{rrna_type};
$MIT_DB_VERSION   = $opt{accession}    if $opt{accession};
$MIT_SCAFFOLD_SEQNAME  = $opt{scaffold}     if $opt{scaffold};
$MIT_CONTIG_SEQNAME    = $opt{contig}       if $opt{contig};
$MIT_CLONE_SEQNAME     = $opt{clone}        if $opt{clone};
$non_interactive       = 1                  if $opt{non_interactive};

unless (    $MIT_DBHOST && $MIT_DBUSER
         && $MIT_DBNAME && $MIT_GENBANK_FILE
         && !$help ) {
  warning( "Can't run without MitConf.pm values: "
      . "MIT_DBHOST $MIT_DBHOST "
      . "MIT_DBUSER $MIT_DBUSER "
      . "MIT_DBNAME $MIT_DBNAME "
      . "MIT_DBPASS $MIT_DBPASS "
      . "MIT_GENBANK_FILE $MIT_GENBANK_FILE "
      . "MIT_LOGIC_NAME $MIT_LOGIC_NAME "
      . "MIT_NAME $MIT_NAME "
      . "MIT_TOPLEVEL $MIT_TOPLEVEL "
      . "MIT_GENE_TYPE $MIT_GENE_TYPE "
      . "MIT_TRNA_TYPE $MIT_TRNA_TYPE "
      . "MIT_RRNA_TYPE $MIT_RRNA_TYPE " );
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}

if (exists $opt{download}) {
  my $base_cmd = 'wget -q "http://eutils.ncbi.nlm.nih.gov/entrez/eutils';
  my $ftp_cmd = $base_cmd.'/efetch.fcgi?db=nuccore&retmod=text&rettype=gb&id='.$MIT_DB_VERSION.'"';
  if (system($ftp_cmd.' -O '.$MIT_GENBANK_FILE)) {
    die("Could not retrieve the genbank file\n");
  }
}

############################
#PARSE GENBANK FILE FIRST
my $genbank_parser = Bio::EnsEMBL::IO::Parser::Genbank->open($MIT_GENBANK_FILE);
$genbank_parser->next;

####################
#Open db for writing

my $output_db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( '-host'     => $MIT_DBHOST,
                                      '-user'     => $MIT_DBUSER,
                                      '-pass'     => $MIT_DBPASS,
                                      '-dbname'   => $MIT_DBNAME,
                                      '-port'     => $MIT_DBPORT,
                                      '-no_cache' => 1, );

########################
# Check taxon is correct

my $meta_container = $output_db->get_MetaContainer();
if ($genbank_parser->get_taxon_id != $meta_container->get_taxonomy_id) {
  # comment out if you are sure you have the correct MT. it can happen
  throw('Your taxon ids differ: '.$meta_container->get_taxonomy_id.' and '.$genbank_parser->get_taxon_id);
}

###########################################
# write chromosome if it is not already there

my $dbe_adaptor = $output_db->get_DBEntryAdaptor;
my $slice_adaptor = $output_db->get_SliceAdaptor;
my $attrib_adaptor = $output_db->get_AttributeAdaptor;
my $slice = $slice_adaptor->fetch_by_region('toplevel', $MIT_NAME);
my $slices =  &get_chromosomes($genbank_parser, $output_db);

if ($slice){
  print "Found chromosome ".$slice->name."\n" ;
}  else {
  print "There is no chromosome in $MIT_DBNAME called $MIT_NAME\nHave found the following slices to load:\n";
  foreach my $cs (keys %$slices){
    print "$cs    \t".$slices->{$cs}->name."\n";
  }
  print "Do you want to load the chromosome and "
      . "the associated assembly entries from the genbank file?(Y/N) ";
  my $answer;
  if($non_interactive) {
    $answer = "y";
  } else {
    $answer = <>;
    chomp $answer;
  }

  if ($answer eq "y" or $answer eq "Y"){
    &load_chromosomes($slices, $output_db, $genbank_parser->get_sequence );
    $slice = $slice_adaptor->fetch_by_region('toplevel',$MIT_NAME);
  }
  else {
    print "Ok gene load aborted.\n";
    exit 0;
  }
}
if ($slice->is_chromosome) {
    my $max_attrib = 0;
    my $answer;
    my ($attrib) = @{$slice->get_all_Attributes('karyotype_rank')};
    if (defined $attrib) {
        print "Your mitochondrion already have a karyotype rank: ", $attrib->value, "\n"
            , "Shall we remove it and insert a new rank? (Y/N) ";
        if($non_interactive) {
          $answer = "y";
         } else {
          $answer = <>;
          chomp $answer;
         }

        if ($answer eq 'y' or $answer eq 'Y') {
            $attrib_adaptor->remove_on_Slice($slice, $attrib);
        }
    }
    foreach my $chromosome (@{$slice->adaptor->fetch_all('chromosome')}) {
        foreach my $attrib (@{$chromosome->get_all_Attributes('karyotype_rank')}) {
            $max_attrib = $attrib->value if ($max_attrib < $attrib->value);
        }
    }
    if ($max_attrib) {
      push my @rank,
        Bio::EnsEMBL::Attribute->new(
                      -CODE        => 'karyotype_rank',
                      -VALUE       => ++$max_attrib );
      print "Karyotype rank for the mitochondrion: $max_attrib\n";
      $attrib_adaptor->store_on_Slice( $slice, \@rank ) unless (defined $answer and ($answer eq 'n' or $answer eq 'N'));
    }
    else {
      print "NOT adding karyotype rank as it would be 1\n";
    }
}

#########################################
#Check that there are not genes in the MT chromosome before uploading the new ones.

my @mt_genes = @{$slice->get_all_Genes};

if (@mt_genes && scalar(@mt_genes) > 0){
  print "There are already ",scalar(@mt_genes)," genes in $MIT_NAME chromosome\n";

  my $g_answer;
  if($non_interactive) {
    $g_answer = "y";
  } else {
    print "Do you want to remove them?(Y/N)\n";
    $g_answer = <>;
    chomp $g_answer;
  }
  if ($g_answer eq "y" or $g_answer eq "Y"){
    my $gene_adaptor = $output_db->get_GeneAdaptor;
    foreach my $mt_gene(@mt_genes){
      print "Removing gene: ",$mt_gene->dbID,"\n";
      $gene_adaptor->remove($mt_gene);
    }
    print "Genes removed succesfully, moving to new genes load\n";
  }else{
    my $load_answer;
    if($non_interactive) {
      $load_answer = "y";
    } else {
      print "You choose not to remove the genes\n"
         . "Do you want to keep loading the MT genes? "
         . "(This may create duplicated entries)(Y/N)?\n";
      $load_answer = <>;
      chomp $load_answer;
    }
    if ($load_answer eq "y" or $load_answer eq "Y"){
      print "Loading genes without removing existing ones, Thanks\n";
    }else{
      print "You choose to abort gene loading. Program exiting\n";
      exit 0;
    }
  }
}

#########################################
#Fetch analysis object

my $logic_name;

if ($MIT_LOGIC_NAME){
  $logic_name = $MIT_LOGIC_NAME;
} else {
  $logic_name = 'mt_genbank_import';
  print "Cannot find MIT_LOGIC_NAME - using standard logic name ensembl\n" if $MIT_DEBUG;
}
if (!$MIT_DB_FILE) {
  warning("MIT_DB_FILE not defined");
}
if (!$MIT_DB_VERSION) {
  warning("DB_VERSION not defined");
}

my $ensembl_analysis = $output_db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if(!defined $ensembl_analysis){
  #croak "analysis $logic_name not found\n";
  $ensembl_analysis =
    Bio::EnsEMBL::Analysis->new( -logic_name => 'mt_genbank_import', -db_version => $MIT_DB_VERSION );
  if (defined $MIT_DB_FILE) {
    $ensembl_analysis->db_file($MIT_DB_FILE);
  }
  print "You have no ".$logic_name." defined creating new object\n";
}


 #########################################################
 # Create gene objects from array of hashes
 # contining parsed coords
 # index 0 contains the accession
 # index 1 contains source data (taxon and sequence file)
 # subsequent indecies contain all the other annotations
 #

my $codon_table_stored = 0;
my $features = $genbank_parser->get_features;
my $length = scalar(@$features);
my $entry_count = 1;
for (my $index = 1; $index < $length; $index++) {
    my $feature = $features->[$index];
    next unless ($feature->{header} eq 'CDS' or $feature->{header} =~ /[tr]RNA/);
    printf "ENTRY %d\n========\n", $entry_count;
    if ($feature->{note}->[0] =~ /frameshift/i) {
        warning('There is a frameshift in '.$feature->{product}->[0].", the gene has more than one exon!\n");
    }
    if (exists $feature->{note}->[0] and $feature->{note}->[0] =~ /tRNAscan-SE/) {
        warning('Skipping '.$feature->{product}->[0].' : '.$feature->{note}->[0]."\n");
        next;
    }
    if (exists $feature->{trans_splicing}) {
        # some trans_splicing features could be represented by exons
        # but other would require exons within the same transcript on different strands
        # which is not supported by the current schema so all trans_splicing are skipped
        warning('Skipping '.$feature->{product}->[0].' : '.$feature->{trans_splicing}->[0]."\n");
        next;
    }

    #############
    # TRANSCRIPTS

    my $transcript = new Bio::EnsEMBL::Transcript;
    my $exon_number = 0;
    my $start_exon;
    my $end_exon;
    my $strand = 1;
    $strand = -1 if ($feature->{position} =~ /complement/);
    my @exons;
    while ($feature->{position} =~ /(\d+)\.\.\>?(\d+)/gc) {
        my $exon = new Bio::EnsEMBL::Exon;
        $exon->start($1);
        $exon->end($2);
        $exon->strand($strand);
        $exon->slice($slice);

        if ($feature->{header} !~ /CDS/) {
          $exon->phase(-1);
          $exon->end_phase(-1);
        } else {
          $exon->phase(0);
          $exon->end_phase(($exon->end - $exon->start + 1)%3);
        }

        if ($strand == 1) {
          push(@exons,$exon);
        } else {
          unshift(@exons,$exon);
        }
        $start_exon = $exon if ($exon_number == 0);
        $end_exon = $exon;
        ++$exon_number;
    }

    foreach my $exon (@exons) {
        $transcript->add_Exon($exon);
    }

    if ($strand == -1) {
        $transcript->start_Exon($end_exon);
        $transcript->end_Exon($start_exon);
    } else {
        $transcript->start_Exon($start_exon);
        $transcript->end_Exon($end_exon);
    }

    my $type;
    my $gene = new Bio::EnsEMBL::Gene;
    if ($feature->{header} =~ /(\w)RNA/) {
        if ($1 eq 't') {
            $type = $MIT_TRNA_TYPE;
        }
        elsif ($1 eq 'r') {
            $type = $MIT_RRNA_TYPE;
        }
        else {
            $type = 'UNKNOWN';
            warning('Unknow type for '.$feature->{header}."\n");
        }
    }
    elsif ($feature->{header} eq 'CDS') {
        #############
        # TRANSLATION

        # Make and store the mitochondrial codon usage attribute
        my $aa = $output_db->get_AttributeAdaptor();
        if (not $codon_table_stored) {
            if ($feature->{transl_table}) {
                $codon_table_stored = $feature->{transl_table}->[0];
            } else {
                # genbank's transl_table default is 1 if it's not present
                $codon_table_stored = 1;
            }
            push my @codonusage,
              Bio::EnsEMBL::Attribute->new(-CODE        => 'codon_table',
                                           -NAME        => 'Codon Table',
                                           -DESCRIPTION => 'Alternate codon table',
                                           -VALUE       => $codon_table_stored);
            $aa->store_on_Slice($slice,\@codonusage);
        } elsif ($feature->{transl_table}) {
            if ($feature->{transl_table}->[0] ne $codon_table_stored) {
                throw('Translation table is '.$feature->{transl_table}->[0].' instead of previously stored '.$codon_table_stored."\n");
            }
        }

        my $translation = new  Bio::EnsEMBL::Translation(
                             -START_EXON => $transcript->start_Exon(),
                             -END_EXON   => $transcript->end_Exon(),
                             -SEQ_START  => 1,
                             -SEQ_END    => $transcript->end_Exon()->length(),
                             );

        my %h_dbentry;
        foreach my $dbentry (@{$dbe_adaptor->fetch_all_by_name($feature->{gene}->[0])}) {
            next unless $dbentry->info_type eq 'UNMAPPED';
            push(@{$h_dbentry{$dbentry->primary_id}}, $dbentry);
        }
        if (!exists $opt{noxref}) {
            foreach my $db_xref (@{$feature->{db_xref}}) {
              my ($key, $value) = split(':', $db_xref);
                next if ($key eq 'GI');
                if (exists $h_dbentry{$value}) {
                    # add dbentry to transcript only if they are of unmapped type
                    warning('Xref already exists for '.$value."\t".$feature->{gene}->[0]."\n");
                    foreach my $dbentry (@{$h_dbentry{$value}}) {
                        $dbentry->info_type(INFO_TYPE);
                        $dbentry->info_text(INFO_TEXT);
                        $gene->add_DBEntry($dbentry);
                    }
                }
                else {
                    if (exists $EXTERNAL_DB{$key}) {
                        warning('Adding xref '.$value. ' from '.$EXTERNAL_DB{$key});
                        my $db_entry = Bio::EnsEMBL::DBEntry->new(
                          -adaptor     => $dbe_adaptor,
                          -primary_id  => $value,
                          -display_id  => $feature->{gene}->[0],
                          -description => $feature->{product}->[0],
                          -dbname      => $EXTERNAL_DB{$key},
                          -info_type   => INFO_TYPE,
                          -info_text   => INFO_TEXT,
                        );
                        $gene->add_DBEntry($db_entry);
                    }
                    else {
                        warning('Not using xref '.$value. ' from '.$key);
                    }
                }
            }
            if (exists $feature->{protein_id}) {
                my $dbname = $feature->{protein_id}->[0] =~ /^NP/ ? REFSEQ_PEP : REFSEQ_XPEP;
                my $dbxentry = Bio::EnsEMBL::DBEntry->new(
                  -adaptor     => $dbe_adaptor,
                  -primary_id  => $feature->{protein_id}->[0],
                  -display_id  => $feature->{gene}->[0],
                  -description => $feature->{product}->[0],
                  -dbname      => $dbname,
                  -info_type   => INFO_TYPE,
                  -info_text   => INFO_TEXT,
                );
                $translation->add_DBEntry($dbxentry);
            }
        }
        if (exists $feature->{note}->[0] and $feature->{note}->[0] =~ /frameshift/) {
            warning('There is a frameshift in: '.$feature->{gene}->[0]."\nNote: ".$feature->{note}->[0]);
        }
        if (exists $feature->{transl_except}) {
            if ($feature->{transl_except}->[0] =~ /pos:(\d+),aa:(\w+)/) {
                my $alt_seq;
                my $pos = $1;
                if ($2 eq 'TERM') {
                    $alt_seq = 'AA' if ($feature->{note}->[0] =~ /TAA/);
                    warning("Adding SeqEdit for TAA stop codon completion by 3' AA residues addition.");
                    my $seq_edit = Bio::EnsEMBL::SeqEdit->new(
                          -CODE    => '_rna_edit',
                          -NAME    => 'rna_edit',
                          -DESC    => 'RNA edit',
                          -START   => $pos,
                          -END     => $pos+1,
                          -ALT_SEQ => $alt_seq
                          );
                    $translation->add_Attributes($seq_edit->get_Attribute());
                }
            }
        }
        if (exists $feature->{fragment}) {
            warning('The gene '.$feature->{gene}->[0]." is fragmented, a methionine will be added!\n");
        }
        $transcript->translation($translation);
        if ($transcript->translate()->seq() =~ /\*/) {
          throw("Stop codon found in translation ".$transcript->translate()->seq()." Transcript: ".$transcript->stable_id()." ".$transcript->seq_region_start()." ".$transcript->seq_region_end());
        }
        if ($transcript->translate()->seq() =~ /^[^M]/) {
          warning("Adding SeqEdit for non-methionine start codon in translation ".$transcript->translate()->seq());
          my $seqedit = Bio::EnsEMBL::SeqEdit->new(
                -CODE    => 'amino_acid_sub',
                -NAME    => 'Amino acid substitution',
                -DESC    => 'Some translations have been manually curated for amino acid substitiutions. For example a stop codon may be changed to an amino acid in order to prevent premature truncation, or one amino acid can be substituted for another.',
                -START   => 1,
                -END     => 1,
                -ALT_SEQ => 'M'
                );
          $transcript->translation()->add_Attributes($seqedit->get_Attribute());
        }
        $type = $MIT_GENE_TYPE;
    }
    eval {
      $gene->biotype($type);
      $gene->analysis($ensembl_analysis);
      $gene->description($feature->{product}->[0]);
      $transcript->biotype($type);
      $transcript->source($MIT_SOURCE_NAME);
      $transcript->analysis($ensembl_analysis);
      $gene->add_Transcript($transcript);
      $gene->canonical_transcript($transcript);
      $gene->source($MIT_SOURCE_NAME);
      $count++;
    };
    if ($@){
      print "Error: $@\n";
      exit;
    }
    printf "\t%-12s %s\n\t%-12s %s\n\t%-12s %5d\n\t%-12s %5d\n\t%-12s %5d\n***********************************************\n\n", 'Description:', $gene->description, 'Biotype:', $gene->biotype, 'Start:',  $gene->start, 'End:', $gene->end, 'Strand:', $gene->strand;
    ++$entry_count;
    push @genes,$gene;
  }


print "Have ".scalar(@genes)." gene objects\n";

print " LOADING : Do you want to load them into the db ? (Y/N) ";
  my $answer;
  if($non_interactive) {
    $answer = "y";
  } else {
    $answer = <>;
    chomp $answer;
  }

  if ($answer eq "y" or $answer eq "Y"){
    &load_db($output_db,\@genes);
  }
else{
  exit 0;
}


######
#TEST

exit 0 unless $MIT_DEBUG;

print "\n\n################################################################
# Testing gene load\n\n" if $MIT_DEBUG;
my $new_genes_adaptor = $output_db->get_GeneAdaptor;
my @new_genes = @{$new_genes_adaptor->fetch_all_by_Slice($slice)};
if (!scalar(@new_genes)) {
  throw("No genes loaded");
}

foreach my $new_gene (@new_genes){
  print ref($new_gene)."\t with dbID ".$new_gene->dbID."\n" ;
  foreach my $new_transcript (@{$new_gene->get_all_Transcripts()}) {
    print "Transcript:\n".$new_transcript->seq->seq."\n";
    if ($new_transcript->translation) {
      print "Translation:\n".$new_transcript->translate()->seq()."\n" if $MIT_DEBUG;
    }
  }
}

exit 0;

sub load_db(){
  my ($output_db,$genes)=@_;
  print "storing genes\n" ;
  #############
  # STORE GENES

  my $gene_adaptor = $output_db->get_GeneAdaptor;
  foreach my $gene(@{$genes}){
    print "Loading gene ",$gene,"\t" if $MIT_DEBUG;
    my $dbid = $gene_adaptor->store($gene);
    print "dbID = $dbid\n" if $MIT_DEBUG;
    my $stored_gene = $output_db->get_GeneAdaptor->fetch_by_dbID($dbid);
  }
  my $geneset_date = sprintf("%d-%02d", (localtime->year()+1900), localtime->mon);
  my $last_geneset_update = $meta_container->single_value_by_key('genebuild.last_geneset_update');
  if ($last_geneset_update) {
      $meta_container->update_key_value('genebuild.last_geneset_update', $geneset_date);
      print STDOUT "Meta key last_geneset_update was $last_geneset_update changed to $geneset_date\n";
  }
  return 0;
}

################################
# Get the sequence if requested

sub get_chromosomes {
  my ($genbank_parser,$output_db,) = @_;
  my %slices;

  my %assembly;
  $assembly{$MIT_TOPLEVEL} = $MIT_NAME;
  if ($MIT_SCAFFOLD_SEQNAME) {
      $assembly{$scaffold} = $MIT_SCAFFOLD_SEQNAME;
  }
  else {
      $assembly{$scaffold} = $genbank_parser->get_sequence_name;
  }
  if ($MIT_CLONE_SEQNAME) {
      $assembly{$clone} = $MIT_CLONE_SEQNAME;
  }
  else {
    my $comments = join(' ', @{$genbank_parser->get_raw_comment});
    if ($comments =~ /reference\s+sequence\s+[a-z ]+([A-Z]{1,2}\d+)./) {
      $assembly{$clone} = $1;
    }
  }
  if ($MIT_CONTIG_SEQNAME) {
      $assembly{$contig} = $MIT_CONTIG_SEQNAME;
  }
  else {
    my $comments = join(' ', @{$genbank_parser->get_raw_comment});
    if ($comments =~ /reference\s+sequence\s+[a-z ]+([A-Z]{1,2}\d+)./) {
      $assembly{$contig} = $1;
    }
  }

  my $csa = $output_db->get_CoordSystemAdaptor();
  my $sa  = $output_db->get_SliceAdaptor();
  # Get all coord systems in the database:
  # Make a slice for each coord system

  foreach my $cs (@{$csa->fetch_all()}) {
    my $name = $cs->name;
    print STDERR $name, ':', $cs->is_sequence_level, "\n";
    $name =  'top_level' if ($cs->name eq $MIT_TOPLEVEL);
    $name =  'seq_level' if ($cs->is_sequence_level);
    if ($assembly{$cs->name}){
      $slices{$name}  = Bio::EnsEMBL::Slice->new
      (
        -coord_system      => $cs,
        -start             => 1,
        -end               => $genbank_parser->get_length,
        -strand            => 1,
        -seq_region_name   => $assembly{$cs->name},
        -seq_region_length => $genbank_parser->get_length,
      );
    }
  }

  # Die before storing anything unless you have sequences that are top level and seq level
  # Unless you only have one coord system in which case you set it to both seq and top level
  die "Cannot find seq_level coord system" unless $slices{'seq_level'};
  die "Cannot find top_level coord system $MIT_TOPLEVEL"
    unless (    scalar( keys %slices ) > 1 && $slices{'top_level'}
             or scalar( keys %slices ) == 1 );

return \%slices;

}

sub load_chromosomes {
  my ($slices, $output_db, $seq_ref) = @_;
  my $sa = $output_db->get_SliceAdaptor();
  my $aa = $output_db->get_AttributeAdaptor();
  # Store slices
  # add the sequence if the coord system is contig
  # add the top level attribute if the coord system is chromosome

  # Make top level seq attribute
  push my @toplevel,
    Bio::EnsEMBL::Attribute->new(
                    -CODE        => 'toplevel',
                    -NAME        => 'Top Level',
                    -DESCRIPTION => 'Top Level Non-Redundant Sequence Region',
                    -VALUE       => 1 );

  foreach my $cs (  sort keys %$slices ) {
    my $slice = $slices->{$cs};
    print "Slice " . $slice->name . "\n" if $MIT_DEBUG;
    if ( $cs eq 'seq_level' ) {
      $sa->store( $slice, \$seq_ref );
      print "Storing seqlevel \n" if $MIT_DEBUG;
      # If only have 1 coord systen it needs to be both seq_level
      # and top level
      if ( scalar( keys %$slices ) == 1 ) {
        $aa->store_on_Slice( $slice, \@toplevel );
      }
      next;
    }
    print "Storing slice \n" if $MIT_DEBUG;
    $sa->store( $slice );
    if ( $cs eq 'top_level' ) {
      $aa->store_on_Slice( $slice, \@toplevel );
    }
  }

  # if you only have 1 coordsystem dont need an assembly
  return 0 if ( scalar( keys %$slices ) == 1 );

  # load the assembly
  # Load a chromosome - scaffold entry in the asembly, if these
  # coord stestems exist

  if ( $slices->{'top_level'} && $slices->{$scaffold} ) {
    print "Making assembly for chromosome vs scaffold\n" if $MIT_DEBUG;
    &load_assembly( $slices->{'top_level'}->get_seq_region_id,
                    $slices->{'top_level'}->start,
                    $slices->{'top_level'}->end,
                    $slices->{$scaffold}->get_seq_region_id,
                    $slices->{$scaffold}->start,
                    $slices->{$scaffold}->end,
                    1,
                    $output_db );
  }

  if ( scalar( keys %$slices ) == 2 ) {
      my $seqregionsynonym = Bio::EnsEMBL::SeqRegionSynonym->new(
        -seq_region_id  => $slices->{'top_level'}->get_seq_region_id,
        -external_db_id => NCBI_EXTERNAL_DB_ID,
        -synonym        => 'MT');
      my $srs= $output_db->get_SeqRegionSynonymAdaptor();
      $srs->store($seqregionsynonym);
  }
  # Load assemby tables for each other coord system vs seq_level

  foreach my $cs ( keys %$slices ) {
    print "Slice " . $slices->{$cs}->name . "\n" if $MIT_DEBUG;
    next if ( $cs eq 'seq_level' );
    print "Making assembly for $cs vs seq level\n" if $MIT_DEBUG;
    &load_assembly( $slices->{$cs}->get_seq_region_id,
                    $slices->{$cs}->start,
                    $slices->{$cs}->end,
                    $slices->{'seq_level'}->get_seq_region_id,
                    $slices->{'seq_level'}->start,
                    $slices->{'seq_level'}->end,
                    1,
                    $output_db );
  }
  return 0;
} ## end sub load_chromosomes

##################################################################
# Do the sql statement to load the values into the assembly table

sub load_assembly {
  my ( $chr_id, $chr_start, $chr_end,
       $contig, $contig_start, $contig_end,
       $contig_ori, $db ) = @_;

  if ( !$contig ) {
    die "contig id must be defined for this to work\n";
  }
  my $sql = "insert into assembly(asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) values(?, ?, ?, ?, ?, ?, ?)";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute( $chr_id, $chr_start, $chr_end, $contig, $contig_start,
                 $contig_end, $contig_ori );
  return 0;
}

##################################################
# Example of genbank entry


#     CDS             join(13818..13986,16435..16470,18954..18991,20508..20984,
#                     21995..22246,23612..23746,25318..25342,26229..26701)
#                     /gene="COX1"
#                     /locus_tag="Q0045"
#                     /EC_number="1.9.3.1"
#                    /note="Subunit I of cytochrome c oxidase, which is the
#                     terminal member of the mitochondrial inner membrane
#                     electron transport chain; one of three
#                     mitochondrially-encoded subunits;
#                     go_component: respiratory chain complex IV (sensu
#                     Eukaryota) [goid GO:0005751] [evidence IPI] [pmid
#                     1331058];
#                     go_function: cytochrome-c oxidase activity [goid
#                     GO:0004129] [evidence IDA] [pmid 1331058];
#                     go_process: aerobic respiration [goid GO:0009060]
#                     [evidence IMP] [pmid 9724417]"
#                     /codon_start=1
#                     /evidence=experimental
#                     /transl_table=3
#                     /product="Cox1p"
#                     /protein_id="NP_009305.1"
#                     /db_xref="GI:6226519"
#                     /db_xref="SGD:S000007260"
#                     /db_xref="GeneID:854598"
#                     /translation="MVQRWLYSTNAKDIAVLYFMLAIFSGMAGTAMSLIIRLELAAPG
#                     SQYLHGNSQLFNVLVVGHAVLMIFFLVMPALIGGFGNYLLPLMIGATDTAFPRINNIA
#                     FWVLPMGLVCLVTSTLVESGAGTGWTVYPPLSSIQAHSGPSVDLAIFALHLTSISSLL
#                     GAINFIVTTLNMRTNGMTMHKLPLFVWSIFITAFLLLLSLPVLSAGITMLLLDRNFNT
#                     SFFEVSGGGDPILYEHLFWFFGHPEVYILIIPGFGIISHVVSTYSKKPVFGEISMVYA
#                     MASIGLLGFLVWSHHMYIVGLDADTRAYFTSATMIIAIPTGIKIFSWLATIHGGSIRL
#                     ATPMLYAIAFLFLFTMGGLTGVALANASLDVAFHDTYYVVGHFHYVLSMGAIFSLFAG
#                     YYYWSPQILGLNYNEKLAQIQFWLIFIGANVIFFPMHFLGINGMPRRIPDYPDAFAGW
#                     NYVASIGSFIATLSLFLFIYILYDQLVNGLNNKVNNKSVIYNKAPDFVESNTIFNLNT
#                     VKSSSIEFLLTSPPAVHSFNTPAVQS"
