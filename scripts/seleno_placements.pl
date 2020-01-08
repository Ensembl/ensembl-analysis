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

use warnings ;
use strict;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Runnable::Genewise;
use Bio::EnsEMBL::SeqEdit;
use Text::Wrap;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Getopt::Long;

# Need to have a fasta file of selenocystein-containing proteins as input
#my $file = '/ecs2/scratch2/fsk/mouse2/data/selenocysteins.sp';

# my $temp_seq_file = '/ecs2/scratch2/fsk/mouse2/data/tmp/seleno_';

# The exonerate rule has to align the proteins in the above file to the genome.
# Output format is via the ryo option: %S %pi %ql %tl %g
# my $exonerate_output_file = '/ecs2/scratch2/fsk/mouse2/data/exonerate.out';

# We need to fetch all the sequence in the db around the exonerate alignments,
# so we need a reference db:

#optional: a seperate DNA db
#dna db

my $protein_residue_window = 10;

my $input_id;
my $verbose = 0; 

my %opt ;   

$opt{dbport} = 3306;
$opt{dnaport} = 3306;
$opt{dnauser} = 'ensro' ;  

&GetOptions( 
            \%opt, 
            'write!'     , 
            'input_id=s' , 
	    'verbose!' , 
            'dbhost=s' , 
            'dbuser=s' , 
            'dbpass=s' , 
            'dbport=s' , 
            'dbname=s' , 
            'dnahost=s' , 
            'dnauser=s' , 
            'dnapass=s' , 
            'dnaport=s' , 
            'dnadbname=s' , 
            'skip_ids=s@' ,    # ids to skip for mouse it is P63300 
            'exonerate_file=s' , # ouput of exonerate  
            'sp_file=s',   # swissprott-file with selenocystein-protein-sequences 
            'temp_seq_file=s',          
           ); 


# for mouse, skip id P63300 'cause it causes trouble 
my @ids_to_skip = map {split/,/} @{$opt{skip_ids}} ; 

$verbose = $opt{verbose} if $opt{verbose}   ; 
my $blessed_seq_io =
  Bio::SeqIO->new(
		  '-file'   => "<$opt{sp_file}",
		  '-format' => 'swiss',
		 );

my %sec_location_hash;
my %sequence_hash;
my @sequence_array;
my %alignment_hash;

print "Reading SwissProt file: $opt{sp_file} \n" if($verbose);
# soak up the input protein ids to run with from the command line;
# gather the sequences into a temporary fasta file,
# and store the locations of the selenocysteins in an array into %sec_location_hash 
# -- and store the array keyed by sequence name into %sequence_hash.

read_swissprot_file_for_sec_positions_and_sequence(
         $blessed_seq_io, \@sequence_array, \%sequence_hash, \%sec_location_hash);


print "Reading exonerate output file: \n" if($verbose);
#Read the exonerate alignment output. This has been produced using the 
#ryo string:%S %pi %ql %tl %g, which means the fields will be:
#queryid, qstart, qend, qorient, targetid, tstart, tend, tstrand, score, rank, score, pid, gene orient
# Place the best alignment for each peptide into %alignment_hash

find_best_exonerate_output($opt{exonerate_file} , \%alignment_hash);

# now we have the best alignment position AND the location(s) of the SEL
# for each peptide. Need a db adaptor to do seq fetching

print "creating db adaptor\n" if($verbose);
my $db; 

  $db    = Bio::EnsEMBL::DBSQL::DBAdaptor->new (
						-user   => $opt{dbuser},
						-dbname => $opt{dbname},
						-pass   => $opt{dbpass},
						-host   => $opt{dbhost},
						-port   => $opt{dbport},
					       );



if($opt{dnadbname} ) { 
  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new (
						   -user   => $opt{dnauser},
						   -dbname => $opt{dnadbname},
						   -host   => $opt{dnahost},
						   -port   => $opt{dnaport},
						  ); 
  $db->dnadb($dnadb) ; 
 }  


print "created db adaptor: ".ref($db)."\n" if($verbose);

#realign the 20 residues around the selenocysteine back to
#the little piece of the genomic sequence from the match.
#Check that we have a stop codon where we expect it.
# Fix the stop codon -> cysteine.
#Then run genewise against the the fixed genomic sequence.

align_protein_subseqs_against_genomic_sequence(
  $db, 
  \%alignment_hash, 
  \%sequence_hash, 
  \%sec_location_hash, 
  $protein_residue_window
);


sub read_swissprot_file_for_sec_positions_and_sequence{
  my($blessed_seq_io, $sequence_array, $sequence_hash, $sec_location_hash) = @_;
  while(my $sequence = $blessed_seq_io->next_seq()) {
    if($opt{input_id} && !($opt{input_id} eq $sequence->accession)){
      next;
    }

    push @{$sequence_array}, $sequence;
    $sequence_hash->{$sequence->accession} = $sequence;
    if (!$sequence->accession || $sequence->accession eq "unknown"){
      throw "Your protein sequence file is not in swissprot format"; 
    } 
    print ">".$sequence->accession."\n";
    print $sequence->seq."\n";
    #
    # this only works if you use a newer bioperl version, ie. bioperl-live
    #
    my @features = $sequence->get_SeqFeatures;
    my @sec_array = ();
    my $counter = 0;
    foreach my $feature(@features){
      if(!($feature->primary_tag eq 'SE_CYS')){
        next;
      }
      $sec_array[$counter] = $feature->start;
      $counter++;
    }
    $sec_location_hash->{$sequence->accession} = \@sec_array;
  }

  print "done gathering protein sequences\n";
}



sub find_best_exonerate_output {
  my ($exonerate_output_file, $alignment_hash) = @_;
  
  open(EXONERATE, "<$exonerate_output_file") or die 
    "can find exonerate output file: $exonerate_output_file\n";

  print "processing exonerate lines\n";
  #store the best alignment we can see for each input protein id
  LINE:while(<EXONERATE>){
    chomp;
    #ryo string:%S %pi %ql %tl %g, which means the fields will be:
    my (
        $queryid, $qstart, $qend, $qorient, 
        $targetid, $tstart, $tend, $tstrand, 
        $score, $pid, $querylength, $tlength, $gene_orient
       ) = split;

    ##for mouse, this protein/gene caused trouble... 
    for my $skip_id  ( @ids_to_skip ) { 
     if($queryid eq $skip_id ){
       print STDERR "\nSkipping ID $skip_id !!\n";
       next LINE;
     }
    }

#    if($queryid eq 'P63300'){
#      print STDERR "\nSkipping P63300!\n";
#      next LINE;
#    }

    my @align_vector = 
      (
       $queryid, $qstart, $qend, $qorient,
       $targetid, $tstart, $tend, $tstrand,
       $score, $pid, $querylength, $tlength, $gene_orient
      ); 
     
     my $count = ($targetid =~ tr/://) ; 

     next LINE unless $count==5 ;  

    #Specific hack to move this protein from chr 1 to chr 19
    
    if(exists $alignment_hash->{$queryid}){
      #specific hack for a wrong placement onto chr 1 when it should go onto chr 19
      my @old_align_vector = @{$alignment_hash->{$queryid}};
      my $old_score = $old_align_vector[8];
      if ($old_score < $score){
        print "already seen this queryid - replacing score for $queryid with $score on $targetid\n";
        $alignment_hash->{$queryid} = \@align_vector; 
      }
    }else{
      #print $targetid ."\n";
      print "setting score for $queryid\n";
      $alignment_hash->{$queryid} = \@align_vector; 
    }
    #print "value of hashref: ".ref($alignment_hash{$queryid})."\n";
  }
}

sub align_protein_subseqs_against_genomic_sequence{
  my ($db, $alignment_hash, $sequence_hash, $sec_location_hash, $protein_residue_window) = @_;
  
  my $slice_adaptor = $db->get_SliceAdaptor;
  
  foreach my $key(keys %alignment_hash){ 
    if($opt{input_id} && !($key eq $opt{input_id})){
      next;
    }

    my @align_vector = @{$alignment_hash->{$key}};
    print "align vector: @align_vector\n" if($verbose);
    my (
        $queryid, $qstart, $qend, $qorient, 
        $targetid, $tstart, $tend, $tstrand, 
        $score, $pid, $querylength, $tlength, $gene_orient
    ) = @align_vector;
  
    my @split_out_name = split (/:/, $targetid);
    my $coord_system = $split_out_name[0];
    my $version = $split_out_name[1];
    my $chromosome_name = $split_out_name[2];
    if($tstart > $tend){
      my $temp = $tstart;
      $tstart = $tend;
      $tend = $temp;
    }

    #Exonerate alignments are 'inbetween' - have to increment by one to get actual genomic coords.
    $qstart+=1;
    $qend+=1;
    $tstart+=1;
    $tend+=1;
    
    print "Fetching seq from : $chromosome_name  $tstart -> $tend\n";
    my $slice = $slice_adaptor->fetch_by_region($coord_system, $chromosome_name, $tstart-1000, $tend+1000);
    my $sequence_string = $slice->get_repeatmasked_seq(undef, 1)->seq;
    
    #print "First 100 chars of sequence: ".substr($sequence_string, 0, 100)."\n";

    my $edited_sequence_string = uc($sequence_string); # (the string we will alter, changing stop->cys)
    my $filename_stem = $opt{temp_seq_file}.$queryid;
    
    #write the genomic sequence
    my $genome_seq = new Bio::Seq(-seq=>$sequence_string, -accession=>$slice->name, -id=>$slice->name);

    my $genome_seq_filename = $filename_stem."_genomic.fa";
    my $bio_seqio = new Bio::SeqIO(-format=>'Fasta', -file=>">".$genome_seq_filename);
    $bio_seqio->write_seq($genome_seq);

    #write the pieces of the protein which are around the sec locations in the protein
    # Then exonerate those pieces against the genomic sequence
    # Then interpret the exonerate output to find the sec location in the genomic (should be a stop codon)
    my $original_protein = $sequence_hash->{$queryid};
    my $long_protein_seq = $original_protein->seq;
    #print "full protein seq: $long_protein_seq\n";
    
    my @sec_pos_array = @{$sec_location_hash->{$queryid}};
    my %sec_location_to_stop_codon_location_map = {};

    #
    #There may be a number of locations of selenocysteine in the protein: the
    # array has all of them. For each one we have to re-align against the genome
    # fragment, and do an edit on the genomic sequence
    foreach my $sec_location(@sec_pos_array){
      my $stop_codon_location = 
        get_stop_codon_location($filename_stem, $tstrand, $long_protein_seq, $queryid, $sec_location, $genome_seq_filename, $sequence_string, $genome_seq);

      $sec_location_to_stop_codon_location_map{$sec_location} = $stop_codon_location;

      if($stop_codon_location >=1){
        #print "Editing\n";
        if($tstrand eq '+'){
          #print "On the outside we get: ".substr($edited_sequence_string, $stop_codon_location - 1, 3)."\n";
          substr($edited_sequence_string, $stop_codon_location - 1, 3) = 'TGT';
          #print "After editing, it's  ".substr($edited_sequence_string, $stop_codon_location - 1, 3)."\n";
        }elsif($tstrand eq '-'){
          #print "On the outside we get: ".substr($edited_sequence_string, $stop_codon_location - 1, 3)."\n";
          substr($edited_sequence_string, $stop_codon_location - 1, 3) = 'ACA';
          #print "After editing, it's  ".substr($edited_sequence_string, $stop_codon_location - 1, 3)."\n";
        }else{
          #print "Not editing\n";
        }
      }
    }

    my $edited_genomic_bioseq = new Bio::Seq(-seq=>$edited_sequence_string, -accession=>$slice->name, -id=>$slice->name);
    my $protein_bioseq = new Bio::Seq(-seq=>$long_protein_seq, -accession=>$queryid, -id=>$queryid);

    my $reverse = undef;
    #print "Target strand is $tstrand \n";
    if($tstrand eq '-'){
      #print "reversing genewise\n";
      $reverse = 1; 
    }else{
      $reverse = undef; 
    }

    my $genewise_runnable = 
      new Bio::EnsEMBL::Pipeline::Runnable::Genewise(
        $edited_genomic_bioseq,
        $protein_bioseq,
        $slice,
        undef,
        $reverse 
      );
      
    $genewise_runnable->run;
    my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('selenocysteine');
    my $full_slice = $db->get_SliceAdaptor->fetch_by_region($coord_system, $chromosome_name);
    
    foreach my $old_gene($genewise_runnable->output){

      my $gene = new Bio::EnsEMBL::Gene;
      
      $gene->stable_id($queryid."--test");
      $gene->version(1);
      $gene->slice($slice);
      
      

      foreach my $old_transcript(@{$old_gene->get_all_Transcripts}){

        my $total_length = 0;
        my $counter = 0;
	my $have_evidence = 0;

	#convert transcript supporting feature if necessary from
	#Bio::EnsEMBL::FeaturePair to Bio::EnsEMBL::DnaPepAlignFeature
	my @all_supp_features;
        foreach my $e (@{$old_transcript->get_all_Exons}) {
	  push @all_supp_features, @{$e->get_all_supporting_features};
	}
        foreach my $f (@{$old_transcript->get_all_supporting_features}) {
	  if($f->isa("Bio::EnsEMBL::BaseAlignFeature")){
	    $have_evidence = 1;
	    $f->analysis($analysis);
	    $f->slice($slice);
	  }
        }
	if(!$have_evidence){
	  print "\nworking on transcript supporting evidence.\n";
	  #remove any non-conform evidence
	  $old_transcript->{_supporting_evidence} = [];
	  if (scalar(@all_supp_features)) {
	    my $daf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@all_supp_features);
	    $daf->analysis($analysis);
	    $daf->slice($slice);
	    $old_transcript->add_supporting_features($daf);
	  }
	  else{
	    throw(" Can't add transcript_supporting_evidence: no features! ");
	  }
	}

        foreach my $e (@{$old_transcript->get_all_Exons}) {
          foreach my $f (@{$e->get_all_supporting_features}) {
            $f->analysis($analysis);
            $f->slice($slice);
          }
          $counter++;
          $total_length += abs(($e->start - $e->end));
        }

        my $transcript = set_stop_codon($old_transcript);
        $gene->add_Transcript($transcript);

        $transcript->stable_id($queryid."--test");
        $transcript->version(1);
        $transcript->slice($slice);

        my $translation = $transcript->translation;
        $translation->stable_id($queryid."--test");
        $translation->version(1);

        $gene->analysis($analysis);
        $gene->transfer($full_slice);

        if($opt{write}){
          $db->get_GeneAdaptor->store($gene);
          print "stored gene: ".$gene->dbID."\n";
        }
        foreach my $transcript(@{$gene->get_all_Transcripts}){
          my $translation = $transcript->translation;
          foreach my $sec_location(@sec_pos_array){
            my $total_nucleotides_before_stop_codon=0;
            my $total_residues_before_stop_codon=0;

            my $stop_codon_location = $sec_location_to_stop_codon_location_map{$sec_location};
            #print "finding new sec location for genomic stop codon: $stop_codon_location\n";

            EXON:
            foreach my $exon(@{$transcript->get_all_Exons}){
              #print "looking at exon: ".$exon->start."-".$exon->end.":".$exon->strand."\n";
              if($stop_codon_location > $exon->start && $stop_codon_location < $exon->end){
                if($exon->strand > 0){
                  #print "stop in + exon: adding ".($stop_codon_location - $exon->start)."\n";
                  $total_nucleotides_before_stop_codon+= $stop_codon_location - $exon->start;
                }else{
                  #print "stop in - exon: adding ".(abs($stop_codon_location + 2 - $exon->end))."\n";
                  $total_nucleotides_before_stop_codon+= abs($stop_codon_location + 2 - $exon->end);
                }
                last EXON;
              }else{
                 #print "stop outside this exon: adding ".(abs($exon->end - $exon->start) + 1)."\n";
                 $total_nucleotides_before_stop_codon+=abs($exon->start - $exon->end) + 1;
              }
            }
           
            $total_residues_before_stop_codon = $total_nucleotides_before_stop_codon / 3;
            print "counted: $total_nucleotides_before_stop_codon, giving $total_residues_before_stop_codon residues\n";
            
            my $seq_edit = 
              Bio::EnsEMBL::SeqEdit->new(
                -CODE    => '_selenocysteine',
                -NAME    => 'Selenocysteine',
                -DESC    => 'Selenocysteine',
                -START   => $total_residues_before_stop_codon+1,
                -END     => $total_residues_before_stop_codon+1,
                -ALT_SEQ => 'U'
              );
  
             my $attribute = $seq_edit->get_Attribute();
             my $translation = $transcript->translation();
             my $attribute_adaptor = $db->get_AttributeAdaptor();
            $attribute_adaptor->store_on_Translation($translation, [$attribute]) if($opt{write});
	    print "Added seqedit at position: ($total_residues_before_stop_codon+1)\n" if($verbose);
          }
        }

      }
    }
  }
}

sub run_exonerate{
  my ($query_seq, $target_seq) = @_;
  my %output_hash;
  
  my $command = 
    "exonerate-0.9.0 --querytype protein --targettype dna ".
    "--bestn 1 ".
    "--ryo \"%S %V\\n\" ".
    "--query $query_seq --target $target_seq ".
    "--showalignment false --showvulgar false --showsugar false";
    
  print "\n>>Exonerate command: $command\n";

  open( EXONERATE_COMMAND, "$command |" ) || die("Error running exonerate $!");
  while(<EXONERATE_COMMAND>){
    if(/Command/){
      next;
    }
    
    if(/Hostname/){
      next;
    }
    
    if(/completed/){
      next;
    }
    
    my (
      $queryid, $qstart, $qend, $qorient, 
      $targetid, $tstart, $tend, $tstrand, 
      $score, @align_components
    ) = split;
    
    #print "REALIGN: $queryid, $qstart, $qend, $qorient,$targetid, $tstart, $tend, $tstrand,$score -- @align_components\n";
    
    $output_hash{queryid}=$queryid;
    $output_hash{qstart}=$qstart;
    $output_hash{qend}=$qend;
    $output_hash{qorient}=$qorient;
    $output_hash{targetid}=$targetid;
    $output_hash{tstart}=$tstart;
    $output_hash{tend}=$tend;
    $output_hash{tstrand}=$tstrand;
    $output_hash{score}=$score;
    $output_hash{align_components}=\@align_components;
    
    last;
  }
  close (EXONERATE_COMMAND);

  return \%output_hash;
}

sub get_stop_codon_location{
  my($filename_stem, $tstrand, $long_protein_seq, $queryid, $sec_location, $genome_seq_filename, $sequence_string, $genome_seq) = @_; 
  #print "Considering ".$queryid."\n";
  # I subseq out 10 residues before the sec and 10 residues after
  my $protein_seq = new Bio::Seq(-seq=>$long_protein_seq, -accession=>$queryid, -id=>$queryid);
  my $prot_end_coord = $sec_location+$protein_residue_window;
  my $prot_start_coord = $sec_location-$protein_residue_window;
  
  if($sec_location+$protein_residue_window > $protein_seq->length){
    $prot_end_coord = $protein_seq->length;
  }
  
  if($sec_location-$protein_residue_window < 1){
    $prot_start_coord = 1;
  }
    
  #print "prot subseq start: $prot_start_coord - $prot_end_coord\n";
    
  my $protein_subseq =
    new Bio::Seq(
      -seq=>
        $protein_seq->subseq(
          $prot_start_coord, 
          $prot_end_coord
        ),
      -accession=>$queryid, 
      -id=>$queryid
    );
    
  my $protein_seq_filename = $filename_stem."_protein_".$sec_location.".fa";
  #print "location: $sec_location\n";
  #print "subseq: ".$protein_subseq->seq."\n";
  my $bio_seqio = new Bio::SeqIO(-format=>'Fasta', -file=>">".$protein_seq_filename);
  $bio_seqio->write_seq($protein_subseq);
  
  #
  #Execute exonerate on the subsequences we've just written
  my $output_hashref = run_exonerate($protein_seq_filename, $genome_seq_filename);
  
  my $target_start = $output_hashref->{tstart};
  my $target_end = $output_hashref->{tend};
  
  #A reverse-strand hit will have end less than start. We only want
  # the start.
  #print  "STRAND: $tstrand\n";
  
  if($target_start > $target_end){
    my $tmp = $target_start;
    #print "switching target limits\n"; 
    $target_start = $target_end;
    $target_end = $tmp;
  }
  
  my $query_start = $output_hashref->{qstart};
  my $query_end = $output_hashref->{qend};
  
  if($query_start > $query_end){
    my $tmp = $query_start;
    #print "switching query limits\n"; 
    $query_start = $query_end;
    $query_end = $tmp;
  }

  #print "EXON offsets: query: $query_start - $query_end target:$target_start - $target_end\n";

  $query_start+=1;
  $query_end+=1;
  $target_start+=1;
  $target_end+=1;

  #print "GENOMIC limits: query: $query_start - $query_end target:$target_start - $target_end\n";
  #print "Genomic subsequence: ".substr($sequence_string, $target_start - 1, $target_end - $target_start + 1)."\n";

  my $genomic_offset;

  if($tstrand eq '+'){
    #You know the sec is at position 11: so there are 30 bp
    $genomic_offset = $target_start+(($protein_residue_window-($query_start-1)) * 3);
  }elsif($tstrand eq '-'){
    #print "neg strand prot offset: target start: $target_start: prot subseq length: ".$protein_subseq->length." prt res wind: ".$protein_residue_window."\n";
    #$genomic_offset = $target_start + 3*($protein_subseq->length - ($protein_residue_window +1));
    $genomic_offset = $target_start + 3*($query_end - ($protein_residue_window+1+1));
  }
  
  #print "sec offset: $genomic_offset\n";
    
  #Now verify that the triplet at the genomic offset is indeed a stop:
  my $stop_seq = $genome_seq->subseq($genomic_offset, $genomic_offset+2);
  
  if(($tstrand eq '+') && uc($stop_seq) eq 'TGA'){
    #print "STOP:(+strand)\n";
    return $genomic_offset; 
  }elsif(($tstrand eq '-') && uc($stop_seq) eq 'TCA'){
    #print "STOP:(-strand)\n";
    return $genomic_offset;
  }else{
    #print "NOT STOP codon? :".$stop_seq."\n";
    return -1;
  }
}
