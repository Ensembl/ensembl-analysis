
=head1 LICENSE

 Copyright [2019-2020] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Star

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Star->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence. Star is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionGroup;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(set_alignment_supporting_features);
# use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs);
use File::Spec;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


# create slices only for projection db. 

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  my $projection_dna_dba = $self->hrdb_get_dba($self->param('projection_dna_db'));
  
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');
  $self->hrdb_set_con($projection_dna_dba,'projection_dna_db');

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  my $projection_gene_dba = $self->hrdb_get_dba($self->param('projection_gene_db'));
  
  $source_gene_dba->dnadb($source_dna_dba);
  $target_gene_dba->dnadb($target_dna_dba);
  $projection_gene_dba->dnadb($projection_dna_dba);
  
  $self->hrdb_set_con($source_gene_dba,'source_gene_db');
  $self->hrdb_set_con($target_gene_dba,'target_gene_db');
  $self->hrdb_set_con($projection_gene_dba,'projection_gene_db');

# $name  = 'chromosome:NCBI34:X:1000000:2000000:1';
  # my $slice_array = ['primary_assembly:CHM13_T2T_v1.0:1:1:248387497:1']; # $self->param('iid');
  my $slice_array = ['1']; 
  
  my $slice_adaptor_proj = $projection_gene_dba->get_SliceAdaptor();
  my $sequence_adaptor_proj = $projection_gene_dba->get_SequenceAdaptor();

  my $slice_adaptor = $source_gene_dba->get_SliceAdaptor();
  my $sequence_adaptor = $source_gene_dba->get_SequenceAdaptor();

  # print "DEBUG::yohooo\n" . $projection_gene_dba; 

  my $source_transcript_id_hash = {};
  my $input_genes = [];
  my $projection_genes = [];
  foreach my $slice_name (@$slice_array) {
  	print "DEBUG:: $slice_name \n"; 
    # my $slice = $slice_adaptor->fetch_by_region('chromosome', $slice_name,);
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $slice_name, 1, 2000000 );
    
    my $genes = $slice->get_all_Genes();
    
    # my $slice_projection = $slice_adaptor_proj->fetch_by_name($slice_name);
    my $slice_projection = $slice_adaptor_proj->fetch_by_region( 'chromosome', $slice_name );
    # print "DEBUG::slice_projection::" . $slice_projection . "\n"; 
    
    $projection_genes = $slice_projection->get_all_Genes(); 
    say "Processing from projection db: ".scalar(@$projection_genes)." genes";
    if ( scalar(@$projection_genes) < 1) {
    	die; 
    }
    my $filtered_genes = [];
    # Put the filtering in a subroutine at some point if it gets more complex
    # At the moment the filtering is just on readthroughs
    say "Got ".scalar(@$genes)." unfiltered genes";
    foreach my $gene (@$genes) {
      my $is_readthrough = 0;
      my $transcripts = $gene->get_all_Transcripts();
      foreach my $transcript (@$transcripts) {
        my $source_transcript_id = $transcript->dbID();
        $source_transcript_id_hash->{$source_transcript_id} = $transcript;      	
        
        if($is_readthrough) {
          last;
        }

        my $attributes = $transcript->get_all_Attributes();
        foreach my $attribute (@{$attributes}) {
          if($attribute->value eq 'readthrough') {
            $is_readthrough = 1;
            last;
          }
        } # foreach my $attribute
      } # End foreach my $transcript
      
      
      ### BK: I have readthroughs in MY DB 
      # unless($is_readthrough) {
        push(@$filtered_genes,$gene);
      #}
    } # End foreach my $gene
    push(@$input_genes,@$filtered_genes);
  }

  say "Processing ".scalar(@$input_genes)." genes";
  
  my $sorted_input_genes = [sort { $a->slice->name() cmp $b->slice->name() or
                                   $a->start() <=> $b->start() or
                                   $a->end() <=> $b->end() }  @{$input_genes}];

  my $parent_gene_id_hash = {};
  for(my $i=0; $i<scalar(@$sorted_input_genes); $i++) {
    my $gene = ${$sorted_input_genes}[$i];
    my $gene_id = $gene->dbID();
    my $gene_stable_id = $gene->stable_id();
    my $gene_version = $gene->version();
    my $gene_biotype = $gene->biotype();


    my $transcripts = $gene->get_all_Transcripts();

    foreach my $transcript (@$transcripts) {
      my $transcript_id = $transcript->dbID();
      my $transcript_stable_id = $transcript->stable_id();
      my $transcript_version = $transcript->version();
      my $biotype = $transcript->get_Biotype();
      my $biotype_group = $biotype->biotype_group();
      my $is_canonical = $transcript->is_canonical();
      my $source = $transcript->source();

      $parent_gene_id_hash->{$transcript_id}->{'gene_id'} = $gene_id;
      $parent_gene_id_hash->{$transcript_id}->{'gene_stable_id'} = $gene_stable_id;
      $parent_gene_id_hash->{$transcript_id}->{'gene_version'} = $gene_version;
      $parent_gene_id_hash->{$transcript_id}->{'gene_biotype'} = $gene_biotype;
      $parent_gene_id_hash->{$transcript_id}->{'transcript_stable_id'} = $transcript_stable_id;
      $parent_gene_id_hash->{$transcript_id}->{'transcript_version'} = $transcript_version;
      $parent_gene_id_hash->{$transcript_id}->{'biotype_group'} = $biotype_group;
      $parent_gene_id_hash->{$transcript_id}->{'is_canonical'} = $is_canonical;
      $parent_gene_id_hash->{$transcript_id}->{'source'} = $source;
      
      # print "DEBUG:: $transcript_id \n"; 
    }

    my $slice = $gene->slice();
    my $stable_id = $gene->stable_id.".".$gene->version;


    my $strand = $gene->strand;
  } #close foreach sorted_input_gene
  say "Processing from projection db 2 : ".scalar(@$projection_genes)." genes";

  $self->process_results($parent_gene_id_hash, $source_transcript_id_hash, $projection_genes);

} 



sub process_results {
  my ($self, $parent_gene_ids, $source_transcript_id_hash, $input_genes) = @_;


  my $output_genes = $input_genes ;  # all genes for this region...  $runnable->output();
  # my $parent_gene_ids = $self->parent_gene_ids();
  my $final_gene_hash = {};
  
  # foreach my $key (keys %$parent_gene_ids) {
  #  print "DEBUG:: $key \n"; 
  # }
  
    
  foreach my $gene (@$output_genes) {
    my $transcripts = $gene->get_all_Transcripts();
    print "DEBUG:: gene_info: " . $gene->dbID() . "\n";
    foreach my $transcript (@$transcripts) {
    print "DEBUG:: transcript_info: " . $transcript->dbID() . "_" . $transcript->seq_region_start() . "_" . $transcript->stable_id() .  "\n";

      unless($parent_gene_ids->{$transcript->stable_id()}) {

        print "DEBUG:NOTDEF:WARNING\n"; 
      	next; # THIS IS VERY DANGEROUS !!!!!!
        $self->throw("The following mapped transcript stable id was not found in the initial list of dbIDs: ".$transcript->stable_id());
      }
      my $parent_gene_id = $parent_gene_ids->{$transcript->stable_id()}->{'gene_id'};
      my $biotype_group = $parent_gene_ids->{$transcript->stable_id()}->{'biotype_group'};
      my $is_canonical = $parent_gene_ids->{$transcript->stable_id()}->{'is_canonical'};
      my $source = $parent_gene_ids->{$transcript->stable_id()}->{'source'};

#       print "DEBUG:: " . $parent_gene_ids->{$transcript->stable_id()} . "\n"; 

      unless($final_gene_hash->{$parent_gene_id}) {
        $final_gene_hash->{$parent_gene_id} = [];
      }

      my $transcript_id = $transcript->stable_id();
      my $source_transcript = $source_transcript_id_hash->{$transcript_id};
      if (!defined($source_transcript)) {
      	print "DEBUG:NOTDEF\n"; 
      	next; 
      }
# print "DEBUG:: $source_transcript exist \n"; 
      my ($coverage,$percent_id) = (0,0);
      # ($coverage,$percent_id) = align_nucleotide_seqs($source_transcript->seq->seq(),$transcript->seq->seq());
#      if($coverage >= 95 && $percent_id > 90) {
      $transcript->biotype($source_transcript->biotype());
      $transcript->created_date($coverage);
      $transcript->modified_date($percent_id);
      my $cov_string = "cov: ".$coverage." perc_id: ".$percent_id;
      $transcript->description($cov_string);

#      } else {
#        $transcript->biotype($source_transcript->biotype()."_weak");
#      }

      $transcript->source($source);

      push(@{$final_gene_hash->{$parent_gene_id}},$transcript);
    }
  } # End foreach my $gene

  my $final_genes = [];
  print "DEBUG::start the output loop\n"; 
  foreach my $gene_id (keys(%{$final_gene_hash})) {
    
    my $transcripts = $final_gene_hash->{$gene_id};
    
    print "DEBUG::gene_id output: $gene_id \n"; 
    # use Data::Dumper; 
    # print "DEBUG::transcripts: Dumper: " .  Dumper($transcripts) . "\n";
    # my $gene = Bio::EnsEMBL::Gene->new(-analysis => 'final_group');
    my $gene = Bio::EnsEMBL::Gene->new();
    my $analysis_new = Bio::EnsEMBL::Analysis->new( -logic_name => 'final_group');
    
    $gene->analysis($analysis_new);

    # Should change this at some point to have a separate has for gene meta data to be tidy
    # if (!defined($parent_gene_ids->{${$transcripts}[0]}) ) { 
    #	print "DEBUG::NEXT::WARNING::gene_id output\n"; 
    #	next ;
    # }
    
    if (!defined(${$transcripts}[0]->stable_id())) { 
      print "DEBUG::NEXT!!! NOT EXIST STABLE_ID \n "; 
      next; 	
    } 
    
    my $parent_gene_id = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_id'}; 
    my $parent_gene_stable_id = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_stable_id'};

    my $parent_gene_version = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_version'};
    my $parent_gene_biotype = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_biotype'};

    foreach my $transcript (@$transcripts) {
      my $parent_transcript_stable_id = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_stable_id'};
      my $parent_transcript_version = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_version'};
      $transcript->stable_id($parent_transcript_stable_id);
      $transcript->version($parent_transcript_version);
      $gene->add_Transcript($transcript);
      my $exons = $transcript->get_all_Exons;
      
    }

    $gene->{'parent_gene_id'} = $parent_gene_id;
    $gene->stable_id($parent_gene_stable_id);
    $gene->version($parent_gene_version);
    $gene->biotype($parent_gene_biotype);
    print_gene($gene);
    print "DEBUG:: Another gene in final output\n"; 
    push(@$final_genes, $gene); 
  }
  
  
  print "DEBUG:: final output is ready to be stored \n"; 

  return($self->output($final_genes) ) ;
}





sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_gene_db');
  my $output_gene_adaptor = $output_dba->get_GeneAdaptor;
  my $output_genes = $self->output();
  foreach my $output_gene (@$output_genes) {
    $output_gene_adaptor->store($output_gene);
  }

  return 1;
}


sub print_gene {
  my ($gene) = @_;

  print STDERR "GENE: ", join(' ', $gene->display_id, $gene->seq_region_name, $gene->seq_region_start, $gene->seq_region_end, $gene->strand, $gene->biotype), "\n";
  foreach my $dbe (@{$gene->get_all_DBEntries}) {
    print STDERR "  DBE: ", $dbe->dbname, ' ', $dbe->description, "\n";
  }
  foreach my $attribute (@{$gene->get_all_Attributes}) {
    print STDERR '  ATTRIBUTE: ', $attribute->code, ' ', $attribute->value, "\n";
  }
  foreach my $transcript (@{$gene->get_all_Transcripts}) {
    print STDERR " TRANSCRIPT: ", join(' ', $transcript->display_id, $transcript->seq_region_name, $transcript->seq_region_start, $transcript->seq_region_end, $transcript->strand, $transcript->biotype), "\n";
    foreach my $dbe (@{$transcript->get_all_DBEntries}) {
      print STDERR "   DBE: ", $dbe->dbname, ' ', $dbe->description, "\n";
    }
    foreach my $attribute (@{$transcript->get_all_Attributes}) {
      print STDERR '   ATTRIBUTE: ', $attribute->code, ' ', $attribute->value || 'NULL', "\n";
    }
    if ($transcript->translation) {
      foreach my $attribute (@{$transcript->translation->get_all_Attributes}) {
        print STDERR '    ATTRIBUTE: ', $attribute->code, ' ', $attribute->value || 'NULL', "\n";
      }
      print STDERR '   TRANSLATION: ', $transcript->translation->display_id, ' ', $transcript->translation->start, ' ', $transcript->translation->end, ' ', $transcript->translation->{phase} || 'NULL', "\n", $transcript->translation->seq, "\n";
      print STDERR "   CDNA: \n", $transcript->translateable_seq, "\n";
    }
    
    if (@{$transcript->get_all_Exons} < 1) {
    	die('NO EXONS???');
    }
    foreach my $exon (@{$transcript->get_all_Exons}) {
      print STDERR "  EXON: ", join(' ', $exon->display_id, $exon->seq_region_name, $exon->seq_region_start, $exon->seq_region_end, $exon->strand, $exon->phase, $exon->end_phase), "\n";
    }
  }
}



1;