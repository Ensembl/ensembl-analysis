package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudopipeParseOutput;

use warnings;
use strict;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

sub fetch_input {
  my ($self) = @_;
  $self->create_analysis;

  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));

  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    my $dna_dba = $self->hrdb_get_dba($self->param_required('dna_db'));
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $pgene_type = $self->param_required('pseudogene_type');
  open(my $fh, $self->param_required('pseudogene_file')) or die "Could not open file '$self->param_required('pseudogene_file'))' $!";
  while(my $pgene = <$fh>){
    print "LINE: ".$pgene."\n";
    my @pgene_info = split /\t/, $pgene;

    my $seq_region_name = $pgene_info[0];

    print "SRN: ".$seq_region_name."\n";
    my $slice_adaptor = $output_dba->get_SliceAdaptor;
#    my $slice = $slice_adaptor->fetch_by_region('primary_assembly', $seq_region_name);
    my $slice = $slice_adaptor->fetch_by_region('toplevel', $seq_region_name);
    my $strand = 1;
    if ($pgene_info[3] eq '-'){
      $strand = -1;
    }

    my $gene_start = $pgene_info[1];
    my $gene_end = $pgene_info[2];
    my $source = "ensembl";

#    my $slice_start = $gene_start - 500;
#    my $slice_end = $gene_end + 500;

#    my $coordsys = $dna_dba->get_CoordSystemAdaptor->fetch_by_name("primary_assembly");
#    my $coordsys_version = $coordsys->version;
#    my $slice_name = "primary_assembly:".$coordsys_version.":".$seq_region_name.":".$slice_start.":".$slice_end.":".$strand;
#    my $slice = $dna_dba->get_SliceAdaptor->fetch_by_name($slice_name);

    my $biotype = "pseudogene";
    if ($pgene_type eq "processed"){
      $biotype = "processed_pseudogene";
    }
    my $analysis = Bio::EnsEMBL::Analysis->new(
                                                -logic_name => 'pseudopipe',
                                                -module => 'HivePseudopipeParseOutput',
                                              );

    my $exon_locs = $pgene_info[13];
    if ($exon_locs =~ m/com/) {
      $exon_locs =~ s/com//;
      #flip the exons?
    }
    $exon_locs =~ s/\(//g;
    $exon_locs =~ s/\)//g;
    my @exon_locs = split / /, $exon_locs;

    my @exons;
    for my $coord_pair (@exon_locs) {
      my @coords = split /\.{2}/, $coord_pair;
      my $exon_start = int($coords[0]);
      my $exon_end = int($coords[1]);

      print "COORDS: ".$exon_start." ".$exon_end." ".$slice->seq_region_name." ".$slice->seq_region_start." ".$slice->seq_region_end."\n";

      my $exon = new Bio::EnsEMBL::Exon(
					-START     => $exon_start,
					-END       => $exon_end,
					-STRAND    => $strand,
					-SLICE     => $slice,
					-ANALYSIS  => $analysis,
				       );
      print "EXON: ".$exon->start()." ".$exon->end()." ".$exon->slice->seq_region_name." ".$exon->slice->seq_region_start." ".$exon->slice->seq_region_end."\n";
      push(@exons, $exon);
    }

    my $transcript = Bio::EnsEMBL::Transcript->new(
					 -EXONS    => \@exons,
					 -ANALYSIS => $analysis,
					 -BIOTYPE  => $biotype,
                                                  );

    print "TRANSCRIPT: ".$transcript->start()." ".$transcript->end()." ".$transcript->slice->seq_region_name." ".$transcript->slice->seq_region_start." ".$transcript->slice->seq_region_end."\n";

    calculate_exon_phases($transcript, 0);
    my @transcripts = ($transcript);

    my $new_gene = Bio::EnsEMBL::Gene->new(
	      -TRANSCRIPTS => \@transcripts,
	      -CANONICAL_TRANSCRIPT => $transcript,
	      -START     => $transcript->start(),
              -END       => $transcript->end(),
              -STRAND    => $transcript->strand(),
              -SLICE     => $transcript->slice(),
              -ANALYSIS  => $analysis,
	      -BIOTYPE   => $biotype,
      );

    print "GENE: ".$new_gene->start()." ".$new_gene->end()." ".$new_gene->slice->seq_region_name." ".$new_gene->slice->seq_region_start." ".$new_gene->slice->seq_region_end."\n";

    $self->output([$new_gene]);
  }
}

sub write_output {
  my ($self) = @_;
  my $genes = $self->output;

  print "IN WRITE OUTPUT!\n";

  # write genes out to a different database from the one we read genes from.
  my $out_dba = $self->hrdb_get_con('output_db');

  my $gene_adaptor = $out_dba->get_GeneAdaptor;
  foreach my $gene ( @{$genes} ) {
    print "FINAL GENE: ".$gene->start()." ".$gene->end()." ".$gene->biotype()."  SLICE: ".$gene->slice->seq_region_name." ".$gene->slice->seq_region_start." ".$gene->slice->seq_region_end."\n";

    empty_Gene($gene);
    $gene_adaptor->store($gene);
  }

  return 1;
}

1;
