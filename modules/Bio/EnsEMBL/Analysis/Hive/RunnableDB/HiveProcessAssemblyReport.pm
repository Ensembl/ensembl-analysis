=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport;

use strict;
use warnings;

use Net::FTP;
use Time::Piece;
use File::Fetch;
use File::Temp;
use File::Spec::Functions qw(catfile splitpath catdir);
use File::Path qw(make_path);
use Digest::MD5;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;

use Bio::EnsEMBL::Hive::Utils qw(destringify);

use Bio::EnsEMBL::IO::Parser::Fasta;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Attribute;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               toplevel_as_sequence_levels => 1,
               _report_name => '#assembly_accession#_#assembly_name#_assembly_report.txt',
               _md5checksum_name => 'md5checksums.txt',
               _genome_file_name => '#assembly_accession#_#assembly_name#_genomic.fna',
               _genome_zip_ext => '.gz',
               _molecule_relations => {
                 'assembled-molecule' => 'chromosome',
                 'unplaced-scaffold' => 'scaffold',
                 'unlocalized-scaffold' => 'scaffold',
                 'alt-scaffold' => 'scaffold',
               },
               ftp_user => 'anonymous',
               ftp_pass => undef,
               load_non_nuclear => 0,
               _agp_branch => 3,
               _coord_systems => {},
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    
    # if any DNA slice sequence length is greater than _MAX_SLICE_LENGTH base pairs,
    # all DNA slice sequences will be cut into _SUBSLICE_LENGTH-base-pair long sequences
    _MAX_SLICE_LENGTH => 950000000,  # 950 million base pairs
    _SUBSLICE_LENGTH => 10000000,    # 10 million base pairs
    _exceeded_max_slice_length => 0, # it will be set to 1 when any DNA slice sequence length is greater than _MAX_SLICE_LENGTH
    
    toplevel_as_sequence_levels => 1,
    _report_name => '#assembly_accession#_#assembly_name#_assembly_report.txt',
    _md5checksum_name => 'md5checksums.txt',
    _genome_file_name => '#assembly_accession#_#assembly_name#_genomic.fna',
    _genome_zip_ext => '.gz',
    _molecule_relations => {
      'assembled-molecule' => 'chromosome',
      'unplaced-scaffold' => 'scaffold',
      'unlocalized-scaffold' => 'scaffold',
      'alt-scaffold' => 'scaffold',
    },
    ftp_user => 'anonymous',
    ftp_pass => undef,
    load_non_nuclear => 0,
    _agp_branch => 3,
    _coord_systems => {},
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Retrieve external db ids for INSDC, RefSeq and UCSC. Retrieve the
              attribute information for 'toplevel' and 'karyotype_rank'.
              Download the report file using 'full_ftp_path' and '_report_name'
 Returntype : None
 Exceptions : Throws if 'target_db' is not set
              Throws if 'assembly_accession' is not set
              Throws if 'assembly_name' is not set
              Throws if 'full_ftp_path' is not set

=cut

sub fetch_input {
  my ($self) = @_;

  my $db = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($db, 'target_db');
  my $assembly_accession = $self->param_required('assembly_accession');
  my $assembly_name = $self->param_required('assembly_name');
  my %external_db_ids;
  my $external_db_adaptor = $db->get_DBEntryAdaptor;
  foreach my $external_db ('INSDC', 'RefSeq_genomic', 'UCSC') {
    $external_db_ids{$external_db} = $external_db_adaptor->get_external_db_id($external_db, undef, 1);
  }
  $self->param('_external_db_ids', \%external_db_ids);
  my ($id, $code, $name, $desc) = @{$db->get_AttributeAdaptor->fetch_by_code('toplevel')};
  $self->param('toplevel_attribute', Bio::EnsEMBL::Attribute->new(
    -code => $code,
    -name => $name,
    -description => $desc,
    -value => 1,
  ));
  ($id, $code, $name, $desc) = @{$db->get_AttributeAdaptor->fetch_by_code('karyotype_rank')};
  $self->param('karyotype_rank', {
    -code => $code,
    -name => $name,
    -description => $desc,
  });
  my $report_dir;
  if ($self->param_is_defined('output_path')) {
    $report_dir = $self->param('output_path');
  }
  else {
    $report_dir = File::Temp->newdir;
    $self->param('output_path', $report_dir);
  }
  if (!-d $report_dir) {
    make_path($report_dir);
  }
  my $fetcher = File::Fetch->new(uri => $self->param_required('full_ftp_path').'/'.$self->param('_report_name'));
  $fetcher->fetch(to => $report_dir);
  $fetcher = File::Fetch->new(uri => $self->param_required('full_ftp_path').'/'.$self->param('_md5checksum_name'));
  $fetcher->fetch(to => $report_dir);
  if ($self->param('toplevel_as_sequence_levels')) {
    $fetcher = File::Fetch->new(uri => $self->param_required('full_ftp_path').'/'.$self->param('_genome_file_name').$self->param('_genome_zip_ext'));
    $fetcher->fetch(to => $report_dir);
  }
}


=head2 run

 Arg [1]    : None
 Description: Process the Genbank report file of the assembly to know which sequences to
              load in the database and retrieve some meta data attached
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $file = catfile($self->param_required('output_path'), $self->param('_report_name'));
  my $one_system = $self->param('toplevel_as_sequence_levels');
  my $synonym_id = $self->param('_external_db_ids');
  my $molecule_matcher = $self->param('_molecule_relations');
  my $load_non_nuclear = $self->param('load_non_nuclear');
  my $karyotype_rank_data = $self->param('karyotype_rank');
  my $no_chromosome = 0;
  my $common_name;
  my $strain;
  my @slices;
  my @dirs;
  my @chromosomes;
  my @agp_files;
  my @fastas;
  my $karyotype_rank = 0;
  open(RH, $file) || $self->throw("Could not open $file");
  while (my $line = <RH>) {
    $line =~ s/\R$//;
    if (!$one_system and $line =~ /^##/) {
      my @data = split("\t", substr($line, 3));
      if ($data[2]) {
        if ($data[0] =~ /^GC/) {
          if ($load_non_nuclear and $data[2] eq 'non-nuclear') {
            push(@dirs, 'non-nuclear');
          }
          elsif ($data[2] ne 'non-nuclear') {
            $data[2] =~ s/ /_/;
            push(@dirs, $data[2]);
          }
        }
      }
    }
    elsif ($line =~ /^#\s+([^:]+):\s+(.+)/) {
      if ($1 eq 'Assembly level') {
        if ($2 eq 'Chromosome') {
          $self->param('chromosomes_present', 1);
        }
        else {
          $no_chromosome = 1;
        }
        if ($2 eq 'Contig') {
          $self->param('toplevel_as_sequence_levels', 1);
          $one_system = 1;
        }
      }
      elsif ($1 eq 'WGS project') {
        $self->param('wgs_project', $2);
      }
      elsif ($1 eq 'Date') {
        my $assembly_date = $2;
        $self->param('assembly_date', sprintf("%d-%02d", $assembly_date =~ /(\d+)-(\d+)/));
      }
      elsif ($1 eq 'Taxid') {
        $self->param('taxon_id', $2);
      }
      elsif ($1 eq 'Organism name') {
        $common_name = $2;
        $self->param('common_name', $1) if ($common_name =~ /\(([^)]+)\)/);
        $common_name =~ s/\s+\(.+//;
        $self->param('scientific_name', $common_name);
      }
      elsif ($1 eq 'Infraspecific name') {
        $strain = $2;
        $self->param('strain', $1) if ($strain =~ /[a-zA-Z]+\=(.+)/);
        $strain =~ s/\=.+//;
        $self->param('strain_type', $strain);
      }
      elsif ($1 eq 'RefSeq assembly accession') {
        if ($self->param_is_defined('assembly_refseq_accession')) {
          $self->warning('RefSeq assembly accession are different, using the one from the file')
            if ($2 ne $self->param('assembly_refseq_accession'));
        }
        $self->param('assembly_refseq_accession', $2);
      }
    }
    elsif ($line =~ /^\w/) {
      my @data = split("\t", $line);
      if (!$load_non_nuclear and $data[7] eq 'non-nuclear') {
        if ($data[1] eq 'assembled-molecule' and $data[6] =~ /^NC/) {
          if ($self->param_is_defined('mt_accession')) {
            $self->warning('MT accession are different, using the one from the file')
              if ($data[6] ne $self->param('mt_accession'));
          }
          $self->param('mt_accession', $data[6]);
        }
        else {
          $self->warning($data[0].' is a non nuclear sequence');
        }
      }
      else {
        my $seq_region_name = $data[4];
        my $karyotype_attribute;
        my $coord_system;
        if ($data[1] eq 'assembled-molecule') {
          $seq_region_name = $data[2];
          push(@chromosomes, $seq_region_name);
          $karyotype_attribute = Bio::EnsEMBL::Attribute->new(
            %$karyotype_rank_data,
            -value => ++$karyotype_rank,
          );
        }
        elsif ($data[7] ne 'Primary Assembly') {
          $seq_region_name = $data[0];
        }
        if ($one_system) {
          $coord_system = $self->get_coord_system('primary_assembly');
        }
        else {
          $coord_system = $self->get_coord_system($molecule_matcher->{$data[1]}, $no_chromosome);
        }
        my $slice = Bio::EnsEMBL::Slice->new(
          -seq_region_name => $seq_region_name,
          -start => 1,
          -end => $data[8] eq 'na' ? 1 : $data[8],
          -coord_system => $coord_system,
        );
        # This is not great but the easiest
        $slice->{karyotype_rank} = $karyotype_attribute if ($karyotype_attribute);
        $slice->{_gb_insdc_name} = $data[4];
        if ($data[4] eq $seq_region_name) {
          push(@{$slice->{add_synonym}}, [$data[0], undef]);
        }
        else {
          push(@{$slice->{add_synonym}}, [$data[4], $synonym_id->{INSDC}]);
        }
        if ($data[6] ne 'na') {
          push(@{$slice->{add_synonym}}, [$data[6], $synonym_id->{RefSeq_genomic}]);
        }
        if ($data[9] ne 'na') {
          push(@{$slice->{add_synonym}}, [$data[9], $synonym_id->{UCSC}]);
        }

        # if the maximum slice length is exceeded for any slice, all the slices will be cut later on
        # and an internal coord_system will be used
        if ($slice->length() > $self->param('_MAX_SLICE_LENGTH')) {
          $self->param('_exceeded_max_slice_length',1);
        }
        push(@slices, $slice);
      }
    }
  }
  close(RH) || $self->throw("Could not close $file");
  open(MH, catfile($self->param('output_path'), $self->param('_md5checksum_name'))) || $self->throw('Cannot open md5checksum file');
  my %checksum;
  while(my $line = <MH>) {
    $line =~ s/\R$//;
    $line =~ /^(\S+)\s+\S+\/([^\/]+)$/;
    $checksum{$2} = $1;
  }
  close(MH) || $self->throw('Cannot close md5checksum file');
  if ($one_system) {
    $self->say_with_header('Checking genome file');
    my $file = catfile($self->param('output_path'), $self->param('_genome_file_name').$self->param('_genome_zip_ext'));
    my $digest = Digest::MD5->new();
    open(my $fh, $file) || $self->throw('Could not open '.$file);
    binmode $fh;
    my $md5sum = $digest->addfile($fh)->hexdigest;
    close($fh) || $self->throw('Could not close '.$file);
    $self->throw('MD5 not as expected '.$checksum{$self->param('_genome_file_name').$self->param('_genome_zip_ext')}.' but '.$md5sum)
      unless ($md5sum eq $checksum{$self->param('_genome_file_name').$self->param('_genome_zip_ext')});
    my ($output) = $file =~ /^(\S+)\.\w+$/;
    anyuncompress $file => $output
      or $self->throw("anyuncompress failed: $AnyUncompressError");
    unlink $file;
  }
  else {
    if (@dirs) {
      my $wgs_project = $self->param_is_defined('wgs_project') ? $self->param('wgs_project') : undef;
      my ($type, $hostname, $base_ftp_dir) = $self->param('full_ftp_path') =~ /^(\w+:\/\/)([^\/]+)\/(\S+)/;
      my $client = Net::FTP->new($hostname);
      $client->login($self->param('ftp_user'), $self->param('ftp_pass'));

      my $structure_dir = $base_ftp_dir.'/'.$self->param('assembly_accession').'_'.$self->param('assembly_name').'_assembly_structure';
      my $base_url = "$type$hostname/$structure_dir";
      if ($client->cwd($structure_dir)) {
        foreach my $dir (@dirs) {
          $client->cwd($dir) || $self->throw("Could not go into $dir");
          foreach my $inode ($client->ls) {
            if ($client->cwd($inode)) {
              foreach my $agp ($client->ls('AGP')) {
                $agp =~ s/AGP\///;
                push(@agp_files, {
                  url => "$base_url/$dir/$inode/AGP/$agp",
                  md5sum => $checksum{$agp},
                  output_dir => catdir('#output_path#', $dir),
                  wgs_project => $wgs_project,
                });
              }
              if (@chromosomes) {
                if ($client->ls("FASTA/chr".$chromosomes[0].'.fna.gz')) {
                  foreach my $fasta ($client->ls('FASTA')) {
                    $fasta =~ s/FASTA\///;
                    push(@fastas, {
                      url => "$base_url/$dir/$inode/FASTA/$fasta",
                      md5sum => $checksum{$fasta},
                      output_dir => catdir('#output_path#', 'FASTA'),
                    });
                  }
                  $self->throw('Could not gat all the fasta files') unless (@chromosomes == @fastas or @chromosomes+1 == @fastas);
                }
              }
              if ($client->ls('FASTA/unplaced.scaffold.fna.gz')) {
                push(@fastas, {
                  url => "$base_url/$dir/$inode/FASTA/unplaced.scaffold.fna.gz",
                  md5sum => $checksum{'unplaced.scaffold.fna.gz'},
                  output_dir => catdir('#output_path#', 'FASTA'),
                });
              }
            }
            if ($inode eq 'component_localID2acc' or $inode eq 'scaffold_localID2acc') {
              push(@agp_files, {
                url => "$base_url/$dir/$inode",
                md5sum => $checksum{$inode},
                output_dir => catdir('#output_path#', $dir),
                uncompress => 0,
              });
            }
          }
        }
      }
      else {
        $self->warning('No assembly structure directory found in ftp');
      }
    }
    if (@agp_files) {
      $self->param('agp_files', \@agp_files);
    }
    else {
      $self->throw('Could not find any AGP files');
    }
    if (@fastas) {
      $self->param('fasta_files', \@fastas);
    }
    else {
      $self->throw('Could not find any FASTA files for checking');
    }
  }
  $self->output(\@slices);
}


=head2 write_output

 Arg [1]    : None
 Description: Write the Bio::EnsEMBL::Slice object corresponding to the toplevel
              sequences to the database.It will also add some meta data from the
              report file when they are present
              On branch '_branch_to_flow_to' (2), send FASTA files to be downloaded
              if needed
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $db = $self->hrdb_get_con('target_db');

  if ($self->param('_exceeded_max_slice_length')) {
    # the assembly.mapping meta_key is required to be able to fetch any sequence
    my $mc = $db->get_MetaContainer();
    $mc->store_key_value('assembly.mapping','primary_assembly:'.$self->param('assembly_name').'|ensembl_internal:'.$self->param('assembly_name'));

    $self->get_coord_system('ensembl_internal');
    $self->get_coord_system('primary_assembly')->{sequence_level} = 0;
  }

  my $coord_system_adaptor = $db->get_CoordSystemAdaptor;
  foreach my $coord_system (values %{$self->param('_coord_systems')}) {
    $coord_system_adaptor->store($coord_system);
  }

  my $slice_adaptor = $db->get_SliceAdaptor;
  my $attribute_adaptor = $db->get_AttributeAdaptor;
  my $toplevel_attribute = $self->param('toplevel_attribute');
  my $toplevel_as_sequence_levels = $self->param('toplevel_as_sequence_levels');
  my $genome_file;
  if ($toplevel_as_sequence_levels) {
    $genome_file = Bio::EnsEMBL::IO::Parser::Fasta->open(catfile($self->param('output_path'), $self->param('_genome_file_name')));
  }
  my %sequences;
  foreach my $slice (@{$self->output}) {
    if ($toplevel_as_sequence_levels) {
      my $insdc_name = $slice->{_gb_insdc_name};
      my $seq;
      if (exists $sequences{$insdc_name}) {
        $seq = $sequences{$insdc_name};
        delete $sequences{$insdc_name};
      } else {
        while ($genome_file->next) {
          $genome_file->getHeader =~ /^(\S+)/;
          if ($insdc_name eq $1) {
            $seq = uc($genome_file->getSequence);
            last;
          } else {
            $sequences{$1} = uc($genome_file->getSequence);
          }
        }
      }
      if ($slice->length <= 1) {
        my %new_slice = %$slice;
        $new_slice{end} = length($seq);
        $new_slice{seq_region_length} = $new_slice{end};
        $slice = Bio::EnsEMBL::Slice->new_fast(\%new_slice);
      }

      if ($self->param('_exceeded_max_slice_length')) {
        $slice_adaptor->store($slice);
        $self->cut_and_store_subslices($slice_adaptor,$slice,\$seq);
      } else {
        $slice_adaptor->store($slice,\$seq);
      }
    } else {
      $slice_adaptor->store($slice);
    }

    if (exists $slice->{add_synonym}) {
      foreach my $synonym_data (@{$slice->{add_synonym}}) {
        $slice->add_synonym(@$synonym_data);
      }
      $slice_adaptor->update($slice);
    }

    if (exists $slice->{karyotype_rank}) {
      $attribute_adaptor->store_on_Slice($slice,[$slice->{karyotype_rank}]);
    }
    $attribute_adaptor->store_on_Slice($slice,[$toplevel_attribute]);
  }
  my $display_name;
  my $common_name = $self->param('common_name');
  my $meta_adaptor = $db->get_MetaContainerAdaptor;
  my $date = localtime->strftime('%Y-%m-Ensembl');
  $meta_adaptor->store_key_value('genebuild.start_date', $date);
  $meta_adaptor->store_key_value('assembly.date', $self->param('assembly_date'));
  $meta_adaptor->store_key_value('species.scientific_name', $self->param('scientific_name'));
  $display_name = $self->param('scientific_name');
  $display_name =~ s/^(\w)/\U$1/;
  $meta_adaptor->store_key_value('species.common_name', $self->param('common_name'));
  $meta_adaptor->store_key_value('species.taxonomy_id', $self->param('taxon_id'));
  $meta_adaptor->store_key_value('assembly.accession', $self->param('assembly_accession'));
  $meta_adaptor->store_key_value('assembly.default', $self->param('assembly_name'));
  $meta_adaptor->store_key_value('assembly.name', $self->param('assembly_name'));
  $meta_adaptor->store_key_value('assembly.web_accession_source', 'NCBI');
  $meta_adaptor->store_key_value('assembly.web_accession_type', 'INSDC Assembly ID');
  $common_name =~ s/^(\w)/\U$1/;
  if ($self->param_is_defined('strain')) {
    $meta_adaptor->store_key_value('species.strain', $self->param('strain'));
    $meta_adaptor->store_key_value('strain.type', $self->param('strain_type'));
    $display_name .= ' ('.$self->param('strain').')';
  }
  else {
    $meta_adaptor->store_key_value('species.strain', 'reference');
    $meta_adaptor->store_key_value('strain.type', 'strain');
    $display_name .= ' ('.$common_name.')';
  }
  $display_name .= ' - '.$self->param('assembly_accession');
  $meta_adaptor->store_key_value('species.display_name', $display_name);

# Not sure it will properly add the new values, hopefully it will and not cause problems
  my $job_params = destringify($self->input_job->input_id);
  if ($self->param_is_defined('mt_accession')) {
    $job_params->{mt_accession} = $self->param('mt_accession');
  }
  if ($self->param_is_defined('assembly_refseq_accession')) {
    $job_params->{assembly_refseq_accession} = $self->param('assembly_refseq_accession');
  }
  if (!$toplevel_as_sequence_levels) {
    $self->dataflow_output_id($self->param('agp_files'), $self->param('_agp_branch'));
    $self->dataflow_output_id($self->param('fasta_files'), $self->param('_branch_to_flow_to'));
  }
  $self->dataflow_output_id($job_params, 'MAIN');
}


=head2 get_coord_system

 Arg [1]    : String type, should be either 'primary_assembly', 'chromosome', 'scaffold' or 'ensembl_internal'.
              'ensembl_internal' represents the sequence level which is created to deal with the assemblies having
              any sequence exceeding the maximum slice length.
 Arg [2]    : Boolean, false by default, set to true if scaffold is the highest coordinate system
 Description: Create the Bio::EnsEMBL::CoordSystem object based on Arg[1], cache the object to avoid
              duplication
 Returntype : Bio::EnsEMBL::CoordSystem
 Exceptions : None

=cut

sub get_coord_system {
  my ($self, $type, $no_chromosome) = @_;

  if (!exists $self->param('_coord_systems')->{$type}) {
    my $rank = 1;
    $rank = 2 if (($type eq 'scaffold' or $type eq 'ensembl_internal') and !$no_chromosome); 

    my $seq_level = $self->param('toplevel_as_sequence_levels');
    if ($type eq 'ensembl_internal') {
      $seq_level = 1;
    }

    my $coord_system;
    if ($type eq 'ensembl_internal') {
      $coord_system = Bio::EnsEMBL::CoordSystem->new(
        -name => $type,
        -rank => $rank,
        -default => 1,
        # note version will be NULL as ensembl_internal is not part of the official assembly
        -sequence_level => $seq_level
      );
    } else {
      $coord_system = Bio::EnsEMBL::CoordSystem->new(
        -name => $type,
        -rank => $rank,
        -default => 1,
        -version => $self->param('assembly_name'),
        -sequence_level => $seq_level
      );
    }
    $self->param('_coord_systems')->{$type} = $coord_system;
  }
  return $self->param('_coord_systems')->{$type};
}

=head2 cut_and_store_subslices

 Arg [1]    : SliceAdaptor type.
 Arg [2]    : Slice type.
 Arg [3]    : string ref.
 Description: Cut the slice Slice into _SUBSLICE_LENGTH-base-pair long slices and store them all into the seq_region, assembly and dna tables of the SliceAdaptor adaptor.
              It returns the seq_region_id assigned to the slice.
 Returntype : int
 Exceptions : None

=cut

sub cut_and_store_subslices {
  my ($self,$slice_adaptor,$slice,$seq) = @_;

  my $slice_seq_region_id = $slice_adaptor->get_seq_region_id($slice);

  # store all subslices with dna
  my $i = $slice->sub_Slice_Iterator($self->param('_SUBSLICE_LENGTH'));
  my $subslice_index = 0;
  while ($i->has_next()) {
    $subslice_index++;
    my $subslice = $i->next();
    my $subseq = substr($$seq,$subslice->start()-1,$self->param('_SUBSLICE_LENGTH'));
    my $subseq_length = length($subseq);
    my $new_subslice = Bio::EnsEMBL::Slice->new(-coord_system => $self->get_coord_system('ensembl_internal'),
                                                -start => 1,
                                                -end => $subseq_length,
                                                -strand => 1,
                                                -seq_region_name => $subslice->seq_region_name()."_".$subslice_index,
                                                -seq_region_length => $subseq_length,
                                                -adaptor => $slice_adaptor);

    my $new_subslice_seq_region_id = $slice_adaptor->store($new_subslice,\$subseq);

    # store the subslice as a component of the whole slice
    my $sql = "INSERT INTO assembly(asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) VALUES(?, ?, ?, ?, ?, ?, ?)";
    my $sth = $slice_adaptor->dbc()->prepare($sql);
    my $new_subslice_start_on_slice = ($subslice_index-1)*$self->param('_SUBSLICE_LENGTH')+1;
    $sth->execute($slice_seq_region_id,
                  $new_subslice_start_on_slice,
                  $new_subslice_start_on_slice+$subseq_length-1,
                  $new_subslice_seq_region_id,
                  1,
                  $subseq_length,
                  1);
  }
}

1;
