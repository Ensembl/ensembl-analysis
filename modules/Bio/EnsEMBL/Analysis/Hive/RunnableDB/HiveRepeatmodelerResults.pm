=head1 LICENSE

Copyright [2018-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatmodelerResults

=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatmodelerResults;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;

use Bio::EnsEMBL::Hive::Utils qw(destringify);
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  my $path_to_genomic = $self->param_required('path_to_genomic_fasta');
  my $path_to_store_assembly_libs = $self->param_required('path_to_assembly_libs');
  my $path_to_store_species_libs  = $self->param_required('path_to_species_libs');
  my $consensi_files_to_process = $self->param_required('min_consensi_files');
  my $assembly_registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(%{$self->param_required('assembly_registry_db')});
  $self->hrdb_set_con($assembly_registry_dba,'assembly_registry_db');
}


sub run {
  my ($self) = @_;

  my $full_path = catfile($self->param_required('path_to_genomic_fasta'),'output');
  my $gca = $self->param_required('iid');
  my ($chain,$version) = $self->split_gca($gca);

  my $assembly_output_path = catfile($self->param_required('path_to_assembly_libs'),$chain,$version);
  my $output_file = catfile($assembly_output_path,$gca.".repeatmodeler.fa");
  if(-e $output_file) {
    $self->throw("Found an existing repeatmodeler file on the assembly path. Will not overwrite. Path:\n".$output_file);
  }

  $self->process_assembly_files($full_path,$output_file,$assembly_output_path,$gca);
  $self->update_species_file($gca,$output_file);
}


sub write_output {
  my ($self) = @_;

  my $job_params = destringify($self->input_job->input_id());
  my $total_consensi_files_processed = $self->total_consensi_files_processed();
  $job_params->{'run_count'} = $total_consensi_files_processed;
  $self->dataflow_output_id($job_params,1);
}


sub process_assembly_files {
  my ($self,$full_path,$output_file,$assembly_output_path,$gca) = @_;
  my $consensi_files_to_process = $self->param_required('min_consensi_files');
  unless(-d $assembly_output_path) {
    if(system('mkdir -p '.$assembly_output_path)) {
      $self->throw("Could not make assembly output dir. Path used:\n".$assembly_output_path);
    }
  }

  my $assembly_run_hash;
  my $total_consensi_files_processed = 0;
  my $family_file = "";
  for(my $i=0; $i<$consensi_files_to_process; $i++) {
    my $run_dir = catfile($full_path,$i);
    unless(-d $run_dir) {
      $self->warning("Warning, did not find the following output dir, will skip:\n".$run_dir);
      next;
    }
    my $cmd = 'find '.$run_dir.' -maxdepth 1 -name RM*';
    my @repeatmodeler_subdirs = `$cmd`;
    unless(scalar(@repeatmodeler_subdirs)) {
      $self->warning("Warning, did not find a repeatmodeler output dir in the run dir. Run dir used:\n$run_dir");
    }

     foreach my $subdir (@repeatmodeler_subdirs) {
      chomp($subdir);
      my $run_lib_file = catfile($subdir,'consensi.fa');
      my $family_lib_file = catfile($subdir,'families.stk');
      unless(-e $run_lib_file) {
        $self->warning("Did not find a consensi.fa lib on path:\n".$run_lib_file);
        next;
      }
      unless(-e $family_lib_file) {
        $self->warning("Did not find a families.stk lib on path:\n".$family_lib_file);
        next;
      }
      if (-e $family_lib_file){
        $family_file = $assembly_output_path . '/families.stk_'.$i;
        say "output is $family_file";
      }
      my $assembly_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($run_lib_file);
      while($assembly_parser->next()) {
        my $sequence = $assembly_parser->getSequence;
        my $header = $assembly_parser->getHeader;
        unless($assembly_run_hash->{$sequence}) {
          $assembly_run_hash->{$sequence} = $header;
        }
      }

      my $copy_cmd = 'cp '.$run_lib_file.' '.catfile($assembly_output_path,'consensi.fa_'.$i);
      if(system($copy_cmd)) {
        $self->throw("Failed to copy the run consensi.fa.classified to the assembly lib dir. Commandline useed:\n".$copy_cmd);
      }
      $copy_cmd = 'cp '.$family_lib_file.' '.catfile($assembly_output_path,'families.stk_'.$i);
      if(system($copy_cmd)) {
        $self->throw("Failed to copy the run families.stk to the assembly lib dir. Commandline useed:\n".$copy_cmd);
      }
      $total_consensi_files_processed++;
    } # end foreach my $subdir
  } # end for(my $i=0; $i<$consensi; $i++)
  unless(open(OUT_ASSEMBLY,">".$output_file)) {
    $self->throw("Could not open an ouput file for writing in assembly storage dir. Path used:\n".$output_file);
  }

  foreach my $seq (keys(%{$assembly_run_hash})) {
    my $header = $assembly_run_hash->{$seq};
    say OUT_ASSEMBLY ">".$header;
    say OUT_ASSEMBLY $seq;
  }

  close OUT_ASSEMBLY;
  
  #looping through to compress families-classified.stk file
  my $cmd = 'find '.$assembly_output_path.' -maxdepth 1 -name families.stk*';
  my @families_subdirs = `$cmd`;
  unless(scalar(@families_subdirs)) {
    $self->warning("Warning, did not find a repeatmodeler output dir in the run dir. Run dir used:\n$assembly_output_path");
  }

  foreach my $subdir (@families_subdirs) {
    chomp($subdir);
    say "Compressing family file now...";
    my $compressed_family_lib = catfile($assembly_output_path,$gca.'.families.stk.gz');
    my $cmd = "gzip -cvf $subdir > $compressed_family_lib";
    say "line is $cmd";
    `$cmd`;
  }

  unless($total_consensi_files_processed >= $self->param_required('min_consensi_files')) {
    $self->throw("Only found ".$total_consensi_files_processed." files, required ".$self->param_required('min_consensi_files'));
  }

  $self->total_consensi_files_processed($total_consensi_files_processed);
}


sub update_species_file {
  my ($self,$gca,$output_file) = @_;

  my $assembly_registry_dba = $self->hrdb_get_con('assembly_registry_db');
  my $species_name = $assembly_registry_dba->fetch_species_name_by_gca($gca);
  $species_name = lc($species_name);
  $species_name =~ s/ +$//; # This is an issue in the assembly registry db that needs fixing
  $species_name =~ s/ +/\_/g;

  unless($species_name) {
    $self->throw("Could not find species name for the GCA in the assembly registry. GCA used: ".$gca);
  }

  my $path_to_store_species_libs  = $self->param('path_to_species_libs');
  my $species_dir_path = catfile($path_to_store_species_libs,$species_name);
  unless(-d $species_dir_path) {
    if(system('mkdir -p '.$species_dir_path)) {
      $self->throw("Could not make species output dir. Path used:\n".$species_dir_path);
    }
  }

  my $species_file = $species_dir_path."/".$species_name.".repeatmodeler.fa";
  unless(-e $species_file) {
    if(system('cp '.$output_file.' '.$species_file)) {
      $self->throw("Could not copy output files to species output dir. Path copy was attemted to:\n".$output_file);
    }
  } else {
    my $updated_species_hash = $self->update_lib($species_file,$output_file);
    open(OUT_SPECIES,">".$species_file);
    foreach my $seq (keys(%{$updated_species_hash})) {
      my $header = $updated_species_hash->{$seq};
      say OUT_SPECIES ">".$header;
      say OUT_SPECIES $seq;
    }
    close OUT_SPECIES;
  }
}


sub update_lib {
  my ($self,$species_file,$output_file) = @_;

  my $species_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($species_file);
  my $assembly_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($output_file);

  my $updated_species_hash;
  while($species_parser->next()) {
    my $sequence = $species_parser->getSequence;
    my $header = $species_parser->getHeader;
    unless($updated_species_hash->{$sequence}) {
      $updated_species_hash->{$sequence} = $header;
    }
  }

  while($assembly_parser->next()) {
    my $sequence = $assembly_parser->getSequence;
    my $header = $assembly_parser->getHeader;
    unless($updated_species_hash->{$sequence}) {
      $updated_species_hash->{$sequence} = $header;
    }
  }

  return($updated_species_hash);
}


sub split_gca {
  my ($self,$chain_version) = @_;
  unless($chain_version =~ /^(GCA\_\d{9})\.(\d+)$/) {
    $self->throw("Could not parse versioned GCA. GCA used: ".$chain_version);
  }

  my $chain = $1;
  my $version = $2;

  return($chain,$version);
}

sub total_consensi_files_processed {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_total_consensi_files_processed',$val);
  }

  return($self->param('_total_consensi_files_processed'));
}
1;
