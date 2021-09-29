#!/usr/bin/env perl
#
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Time::localtime;

use Getopt::Long;

use constant {
    INFO_TYPE           => 'DIRECT',
    INFO_TEXT           => 'Imported from GeneBank',
    REFSEQ_PEP          => 'RefSeq_peptide',
    REFSEQ_XPEP         => 'RefSeq_peptide_predicted',
    NCBI_EXTERNAL_DB_ID => 700,
    REFSEQ_EXTERNAL_DB_ID => 1830,
    INSDC_EXTERNAL_DB_ID => 50710,
};

my %EXTERNAL_DB = ( 'GeneID'               => 'EntrezGene',
                    'UniProtKB/Swiss-Prot' => 'Uniprot/SWISSPROT',
                  );
my $help;
my @genes;
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
            'xref!',
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

my $ftp_cmd = 'wget -q "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmod=text&rettype=%s&id=%s"';
if (exists $opt{download}) {
  if (system(sprintf($ftp_cmd, 'gb', $MIT_DB_VERSION).' -O '.$MIT_GENBANK_FILE)) {
    throw("Could not retrieve the genbank file");
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
                                      );

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

if (!$slice) {
  my $slices = &get_chromosomes($genbank_parser, $output_db);
  my $sequence;
  if ($genbank_parser->{record}->{_tax} eq 'CON') {
    my $data = $genbank_parser->{record}->{_raw_contig};
    $data =~ s/^join\(//;
    # I'm using length-1 to remove the last bracket
    my @fragments = split(',', substr($data, 0, length($data)-1));
    foreach my $fragment (@fragments) {
      if (substr($fragment, 0, 3) eq 'gap') {
        my ($length) = $fragment =~ /(\d+)/;
        $sequence .= 'N'x$length;
      }
      else {
        my ($rev, $accession, $start, $end) = $fragment =~ /(complement\()?([^:]+):(\d+)\.\.(\d+)/;
        my $length = $end-$start+1;
        open(my $fasta_file, sprintf($ftp_cmd, 'fasta', $accession).' -O - |') or throw("Could not fetch fasta for $accession");
        my $fasta_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($fasta_file);
        $fasta_parser->next;
        my $tmp_sequence = substr($fasta_parser->getSequence, $start-1, $length);
        if ($rev) {
          $tmp_sequence = reverse($tmp_sequence);
          $tmp_sequence =~ tr/ATGC/TACG/;
        }
        $sequence .= $tmp_sequence;
      }
    }
  }
  else {
    $sequence = $genbank_parser->get_sequence;
  }
  &load_chromosomes($slices, $output_db, $sequence);
  $slice = $slice_adaptor->fetch_by_region('toplevel',$MIT_NAME);
}

# Trying to add karyotype_rank attribute if we have chromosomes
my $max_attrib = 0;
foreach my $chromosome (@{$slice->adaptor->fetch_all_karyotype}) {
  foreach my $attrib (@{$chromosome->get_all_Attributes('karyotype_rank')}) {
    $max_attrib = $attrib->value if ($max_attrib < $attrib->value);
  }
}
if ($max_attrib) {
  my ($attrib) = @{$slice->get_all_Attributes('karyotype_rank')};
  if ($attrib and $attrib->value != $max_attrib) {
    $attrib_adaptor->remove_on_Slice($slice, $attrib);
  }
  unless ($attrib and $attrib->value != $max_attrib) {
    $attrib = Bio::EnsEMBL::Attribute->new(
        -CODE        => 'karyotype_rank',
        -VALUE       => ++$max_attrib
        );
    $attrib_adaptor->store_on_Slice($slice, [$attrib]);
  }
}
else {
  warning("NOT adding karyotype rank as it would be 1");
}

#########################################
#Check that there are not genes in the MT chromosome before uploading the new ones.

my $mt_genes = $slice->get_all_Genes;

if (@$mt_genes && scalar(@$mt_genes) > 0){
  warning("There are already ",scalar(@$mt_genes)," genes in $MIT_NAME chromosome");
  my $gene_adaptor = $output_db->get_GeneAdaptor;
  foreach my $mt_gene(@$mt_genes){
    warning('Removing gene: '.$mt_gene->dbID);
    $gene_adaptor->remove($mt_gene);
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
    if (exists $feature->{note}->[0] and $feature->{note}->[0] =~ /frameshift/i) {
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
        if (exists $opt{xref}) {
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

&load_db($output_db,\@genes);

exit 0;

=head2 load_db

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor database to store the data into
 Arg [2]    : Arrayref of Bio::EnsEMBL::Gene
 Description: Load the genes into the database and update the genebuild.last_geneset_update
              meta key if needed
 Returntype : Boolean 0 if successful
 Exceptions : None

=cut

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

=head2 get_chromosomes

 Arg [1]    : Bio::EnsEMBL::IO::Parser::Genbank Genbank file parser
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor database to query
 Description: Create the MT slice(s) based on the default coordinate system
              and prepare the synonyms
 Returntype : Hashref of Bio::EnsEMBL::Slice
 Exceptions : Throws if it cannot find at least one coordinate system
              Throws if it cannot find a sequence level coordinate system

=cut

sub get_chromosomes {
  my ($genbank_parser,$output_db,) = @_;
  my %slices;

  my %assembly;
  my $original_entry;
  foreach my $line (@{$genbank_parser->get_raw_comment}) {
    if ($line =~ /reference\s+sequence\s+[a-z ]+([A-Z]{1,2}\d+)./) {
      $original_entry = $1;
      last;
    }
  }

  my $csa = $output_db->get_CoordSystemAdaptor();
  # Get all coord systems in the database:
  # Make a slice for each coord system

  my $region_name;
  my $has_toplevel = 0;
  my $has_seqlevel = 0;
  foreach my $cs (@{$csa->fetch_all()}) {
    my @synonyms;
    if ($cs->rank == 1) {
      $region_name = 'MT';
      $has_toplevel = 1;
      push(@synonyms, $genbank_parser->get_sequence_name);
      if ($cs->is_sequence_level) {
        push(@synonyms, $original_entry);
      }
    }
    else {
      $region_name = $original_entry;
    }
    if ($cs->is_sequence_level) {
      $has_seqlevel = 1;
    }
    $slices{$cs->name} = Bio::EnsEMBL::Slice->new(
      -coord_system      => $cs,
      -start             => 1,
      -end               => $genbank_parser->get_length,
      -strand            => 1,
      -seq_region_name   => $region_name,
      -seq_region_length => $genbank_parser->get_length,
    );
    if (@synonyms) {
      $slices{$cs->name}->{_gb_synonyms} = \@synonyms;
    }
  }
  throw('Could not find toplevel and/or seqlevel coord_system') unless ($has_toplevel and $has_seqlevel);

  return \%slices;
}

=head2 load_chromosomes

 Arg [1]    : Hashref of Bio::EnsEMBL::Slice
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor the database to store the data into
 Arg [3]    : String the DNA sequence of the mitochondrion
 Description: Load the regions, the dna and the attributes for the mitochondrion into
              a core database
              If it is a chromosome assembly, it will add the karyotype_rank attribute
              for the mitochondrion to be the last sequence on the diagram
 Returntype : Boolean 0 if successful
 Exceptions : Throws if you have more than one slice to load but no assembly information

=cut

sub load_chromosomes {
  my ($slices, $output_db, $seq_ref) = @_;
  my $sa = $output_db->get_SliceAdaptor();
  my $aa = $output_db->get_AttributeAdaptor();
  my $srs= $output_db->get_SeqRegionSynonymAdaptor();
  # Store slices
  # add the sequence if the coord system is contig
  # add the top level attribute if the coord system is chromosome

  my $cs_version;
  foreach my $slice (values %$slices) {
    print "Slice " . $slice->name . "\n" if $MIT_DEBUG;
    if ( $slice->coord_system->is_sequence_level) {
      $sa->store( $slice, \$seq_ref );
    }
    else {
      $sa->store( $slice );
    }
    if ($slice->coord_system->rank == 1) {
      $cs_version = $slice->coord_system->version;
      my $toplevel_attribute = Bio::EnsEMBL::Attribute->new(
        -CODE        => 'toplevel',
        -NAME        => 'Top Level',
        -DESCRIPTION => 'Top Level Non-Redundant Sequence Region',
        -VALUE       => 1,
      );
      $aa->store_on_Slice( $slice, [$toplevel_attribute]);
      if ($slice->seq_region_name eq 'MT') {
        if (exists $slice->{_gb_synonyms}) {
          foreach my $synonym (@{$slice->{_gb_synonyms}}) {
            my $external_db_id = INSDC_EXTERNAL_DB_ID;
            if ($synonym =~ /N\w_\d+/) {
              $external_db_id = REFSEQ_EXTERNAL_DB_ID;
            }
            my $seqregionsynonym = Bio::EnsEMBL::SeqRegionSynonym->new(
              -seq_region_id  => $slice->get_seq_region_id,
              -external_db_id => $external_db_id,
              -synonym        => $synonym);
            $srs->store($seqregionsynonym);
          }
        }
      }
      else {
        my $seqregionsynonym = Bio::EnsEMBL::SeqRegionSynonym->new(
          -seq_region_id  => $slice->get_seq_region_id,
          -external_db_id => undef,
          -synonym        => 'MT');
        $srs->store($seqregionsynonym);
      }
    }
  }

  # if you have more than 1 slice you need to load an assembly
  # It will use the assembly_mapping keys as they should be present
  if (scalar(keys %$slices) > 1) {
    my $assembly_mappings = $output_db->get_MetaContainer->list_value_by_keys('assembly.mapping');
    if (@$assembly_mappings) {
      foreach my $assembly_mapping (@$assembly_mappings) {
        my @coord_systems = split(':', $assembly_mapping);
        if (@coord_systems == 2) {
          my ($asm_name) = $coord_systems[0] =~ /^(\w+)[|#]$cs_version$/;
          my ($cmp_name) = $coord_systems[1] =~ /^(\w+)([|#]$cs_version)?$/;
          if ($asm_name and $cmp_name) {
            &load_assembly(
              $slices->{$asm_name}->get_seq_region_id,
              $slices->{$asm_name}->start,
              $slices->{$asm_name}->end,
              $slices->{$cmp_name}->get_seq_region_id,
              $slices->{$cmp_name}->start,
              $slices->{$cmp_name}->end,
              1,
              $output_db
            );
          }
          else {
            warning("$assembly_mapping does not seem to be for $cs_version");
          }
        }
      }
    }
    else {
      throw('Something went wrong, you have '.scalar(keys %$slices).' slices but no "assembly.mapping" in the meta table');
    }
  }
  return 0;
} ## end sub load_chromosomes

##################################################################
# Do the sql statement to load the values into the assembly table

=head2 load_assembly

 Arg [1]    : Int object dbID
 Arg [2]    : Int object start
 Arg [3]    : Int object end
 Arg [4]    : Int component dbID
 Arg [5]    : Int component start
 Arg [6]    : Int component end
 Arg [7]    : Int component direction compared to object
 Description: Load the assembly table for mapping purposes
 Returntype : Boolean 0 if successful
 Exceptions : None

=cut

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
