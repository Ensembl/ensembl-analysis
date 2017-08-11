use warnings;
use strict;
use feature 'say';
use Net::FTP::File;
use Data::Dumper;
my $config_file = $ARGV[0];

my $general_hash = {};
my @species = ();

open(IN,$config_file);
while(<IN>) {
    my $line = $_;
    if($line =~ /^\[.+\]$/) {
		push(@species,$line);
    }elsif($line =~ /(.+)\=(.+)\:(.+)/){
		my $key = $1;
		my $server =$2;
		my $port = $3;
		my $tmp_hash = {};
		$tmp_hash->{'server'} = $server;
		$tmp_hash->{'port'} = $port;
		$general_hash->{$key} = $tmp_hash;
    }elsif($line =~ /(.+)\=([^:]+)\n/) {
		my $key = $1;
		my $value = $2;
		$general_hash->{$key} = $value;
    }elsif($line eq "\n") {
	# Skip
    }else {
		say "Line format not recognised. Skipping line:\n".$line;
    }
}
print Dumper ($general_hash);
close IN;

unless(-e $general_hash->{'work_dir'}) {
    system("mkdir -p ".$general_hash->{'work_dir'});
}

open(LOOP_CMD,">".$general_hash->{'work_dir'}."/beekeeper_cmds.txt");
open(CLEAN_CMD,">".$general_hash->{'work_dir'}."/clean_dir_dbs.sh");
system("chmod +x ".$general_hash->{'work_dir'}."/clean_dir_dbs.sh");
foreach my $row (@species) {
  unless($row =~ /^\[([^,]+),([^,]+),([^,]+),([^,]+)\]$/) {
      die "Issue parsing the following row:\n".$row."\nExpected format:\n[repeatmasker_library_name,uniprot_set_name,assembly_ftp_link,refseq_ftp_link]";
  }

  my $repeatmasker_library = $1;
  my $uniprot_set = $2;
  my $ftp_path = $3;
  my $refseq_ftp_path = $4;

  unless($ftp_path =~ /\/genomes\/genbank\/[^\/]+\/([^\/]+)\/all_assembly_versions\/([^\/]+)\//) {
    die "Failed to parse the following line:\n".$ftp_path."\n\nExpected a link to the main dir for the species. For example:\n".
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Propithecus_coquereli/all_assembly_versions/GCA_000956105.1_Pcoq_1.0/";
  }

  my $species_name = $1;
  $species_name = lc($species_name);

  # Cleaning commands:
 # my $clean_commands = clean_commands($general_hash,$species_name);
 # say CLEAN_CMD $clean_commands;

  my $gca_and_assembly_version = $2;
  unless($gca_and_assembly_version =~ /^(GCA\_\d+\.\d+)\_(.+)/) {
    die "Couldn't parse the assembly version out of the ftp link. Tried to parse off end of:\n".$gca_and_assembly_version;
  }
  my $gca = $1;
  my $assembly_version = $2;
  
  my $output_path = $general_hash->{'work_dir'}."/".$species_name."/".$assembly_version."/";
  print "\n\nDEBUG:: ".$output_path."\n\n".$general_hash->{'work_dir'}."\nEND\n";
  my $output_hash = {};
  $output_hash->{'gca'} = $assembly_version;
  $output_hash->{'assembly_version'} = $assembly_version;

  my $cmd = "wget -P ".$output_path." ".$ftp_path."/".$gca_and_assembly_version."_assembly_report.txt";
  my $return = system($cmd);
  if($return) {
      die "wget for assembly report failed. Commandline used:\n".$cmd;
  }

  open(ASSEMBLY_REPORT,$output_path."/".$gca_and_assembly_version."_assembly_report.txt");
  my @assembly_report = <ASSEMBLY_REPORT>;
  close ASSEMBLY_REPORT;

  my $taxon_id;
  my $assembly_level;
  my %wgs_code = ();
  my $assembly_refseq_accession;
  foreach my $assembly_line (@assembly_report) {
      chomp $assembly_line;

      if($assembly_line =~ /^#/) {
	  if($assembly_line =~ /Taxid\:[^\d]+(\d+)/) {
	      $taxon_id = $1;
	  } elsif($assembly_line =~ /Assembly level: ([a-zA-Z]+)/) {
	      $assembly_level = $1;
	  } elsif($assembly_line =~ /RefSeq assembly accession: (.*)/){
	  	  $assembly_refseq_accession = $1;
	  	  $assembly_refseq_accession =~ s/\s+$//;
	  }
	  next;
      }
      
      my @assembly_columns = split("\t",$assembly_line);
      my $accession = $assembly_columns[4];
      unless($accession && ($accession =~ /(^[A-Z]{4})/ || $accession =~ /^gb\|([A-Z]{4})/)) {
	  next;
      }

      my $code = $1;
      if(exists($wgs_code{$code})) {
	  $wgs_code{$code}++;
      } else {
	  $wgs_code{$code} = 1;
      }
  }

  unless($taxon_id) {
      die "Failed to find and parse 'Taxid' line from report file. File used:\n".$gca_and_assembly_version."_assembly_report.txt";
  }

  unless($assembly_level) {
      die "Failed to find and parse 'Assembly level' line from report file. File used:\n".$gca_and_assembly_version."_assembly_report.txt";
  }

  $assembly_level = lc($assembly_level);
  unless($assembly_level eq 'scaffold' || $assembly_level eq 'chromosome') {
      die "Parsed assembly level from report file but it was not 'scaffold' or 'chromosome'. Level found:\n".$assembly_level;
  }

  my $chromosomes_present = 0;
  if($assembly_level eq 'chromosome') {
      $chromosomes_present = 1;
  }

  # This section deals with either getting 0, 1 or many potential wgs codes. Getting 1 is ideal
  my @code_types = keys(%wgs_code);
  if(scalar(@code_types == 0)) {
      die "Failed to parse any potential 4 letter wgs code from the accessions in ".$gca_and_assembly_version."_assembly_report.txt";
  }

  if(scalar(@code_types > 1)) {
    die "Failed because of multiple potential 4 letter wgs code from the accessions in ".$gca_and_assembly_version."_assembly_report.txt, ".
        " codes parsed:\n".@code_types;
  }

  my $wgs_id = $code_types[0];

  my $full_ftp_path = $ftp_path."/".$gca_and_assembly_version."_assembly_structure";

# Is contigs_source 'NCBI' or 'ENA'? 
# Check for existence of files in the ftp using lftp command
my $half_wgs = substr(lc($wgs_id), 0, 2);
my $ncbi_suffix = 'wgs.'.$wgs_id.'.*.fsa_nt.gz';
my $ena_suffix = $half_wgs.'/'.$wgs_id.'*';
my $ncbi_lftp = `lftp -u anonymous,password -e "ls /genbank/wgs/$ncbi_suffix ;quit" ftp.ncbi.nih.gov`;
my $ena_lftp = `lftp -u anonymous,password -e "ls /pub/databases/ena/wgs_fasta/$ena_suffix ;quit" ftp.ebi.ac.uk`;

my $contigs_source = '';

if(!$ncbi_lftp){
    if(!$ena_lftp){
        die "Failed to find contigs at NCBI or ENA\n";
    }
    else{
        $contigs_source = 'ENA';
    }
}
else{
    $contigs_source = 'NCBI';
}
# End set contigs_source

# variables from input_config
  $output_hash->{'user_r'} = $general_hash->{'user_r'};
  $output_hash->{'user_w'} = $general_hash->{'user_w'};
  $output_hash->{'password'} = $general_hash->{'password'};
  $output_hash->{'port'} = $general_hash->{'port'};
  $output_hash->{'dbowner'} = $general_hash->{'farm_user_name'};
  $output_hash->{'email_address'} = $general_hash->{'farm_user_name'}.'@ebi.ac.uk';
  $output_hash->{'release_number'} = $general_hash->{'release_number'};
  $output_hash->{'genebuilder_id'} = $general_hash->{'genebuilder_id'};
  $output_hash->{'farm_user_name'} = $general_hash->{'farm_user_name'};
  
# servers/ports
  $output_hash->{'pipe_db_port'} = $general_hash->{'pipe_server'}->{'port'};
  $output_hash->{'pipe_db_server'} = $general_hash->{'pipe_server'}->{'server'};
  $output_hash->{'killlist_db_port'} = $general_hash->{'pipe_server'}->{'port'};
  $output_hash->{'killlist_db_server'} = $general_hash->{'pipe_server'}->{'server'};
  $output_hash->{'dna_db_port'} = $general_hash->{'core_server'}->{'port'};
  $output_hash->{'dna_db_server'} = $general_hash->{'core_server'}->{'server'};
  $output_hash->{'output_db_port'} = $general_hash->{'output_server'}->{'port'};
  $output_hash->{'output_db_server'} = $general_hash->{'output_server'}->{'server'};
  $output_hash->{'projection_source_db_name'} = '';#$general_hash->{''};
  
 # paths 
  $output_hash->{'enscode_root_dir'} = $general_hash->{'enscode_dir'};
  $output_hash->{'uniprot_blast_db_path'} = $general_hash->{'uniprot_blast_db_path'};
  $output_hash->{'unigene_blast_db_path'} = $general_hash->{'unigene_blast_db_path'};
  $output_hash->{'vertrna_blast_db_path'} = $general_hash->{'vertrna_blast_db_path'};

  $output_hash->{'output_path'} = $output_path;
  $output_hash->{'genome_file'} = $output_path.'/genome_dumps/'.$species_name.'_softmasked_toplevel.fa';
  $output_hash->{'contigs_source'} = $contigs_source;
  $output_hash->{'assembly_name'} = $assembly_version;
  $output_hash->{'assembly_accession'} = $gca;
  $output_hash->{'assembly_refseq_accession'} = $assembly_refseq_accession;
  $output_hash->{'gca'} = $gca;
  $output_hash->{'gca_and_assembly_version'} = $gca_and_assembly_version;
  $output_hash->{'taxon_id'} = $taxon_id;
  $output_hash->{'assembly_level'} = $assembly_level;
  $output_hash->{'chromosomes_present'} = $chromosomes_present;
  $output_hash->{'species_name'} = $species_name;
  $output_hash->{'production_name'} = $species_name;
  $output_hash->{'full_ftp_path'} = $full_ftp_path;
  $output_hash->{'wgs_id'} = $wgs_id;
  $output_hash->{'coord_system_version'} = $assembly_version;
  $output_hash->{'repeatmasker_library'} = $repeatmasker_library;
  $output_hash->{'repeatmasker_species'} = $repeatmasker_library;
  $output_hash->{'uniprot_set'} = $uniprot_set;

	print Dumper $output_hash;

  create_config($output_hash);

  chdir($output_path);
  $cmd = "init_pipeline.pl Genome_annotation_conf.pm";
  my $result = `$cmd`;
  unless($result =~ /beekeeper.+\-sync/) {
      die "Failed to run init_pipeline for ".$species_name."\nCommandline used:\n".$cmd;
  }

  my $sync_command = $&;
  $return = system($sync_command);
  if($return) {
      die "Failed to sync the pipeline for ".$species_name."\nCommandline used:\n".$cmd;
  }

  my $loop_command = $sync_command;
  $loop_command =~ s/sync/loop \-sleep 0.5/;
  say LOOP_CMD $loop_command;
}
close LOOP_CMD;
close CLEAN_CMD;

exit;

sub create_config {
    my ($output_hash) = @_;

    my $output_path = $output_hash->{'output_path'};

    my $cmd = "cp /homes/leanne/development/annotation_configs_automation/config_template.pm ".$output_path.'/Genome_annotation_conf.pm.tmp';
    print $cmd."\n\n";
    my $return = system($cmd);
    if($return) {
	die "Failed to copy parent config. Commandline used:\n".$cmd;
    }

    my $conf_file = "";
    open(IN,$output_path.'/Genome_annotation_conf.pm.tmp');
    while(<IN>) {
	$conf_file .= $_;
    }

  # write to output
  
  while((my $key, my $value) = each %$output_hash){
  	$conf_file =~ s/('$key' +\=\> +)''/$1'$value'/;

  }

  # If the path has been set this will delete a commented out input id in the config file and make the
  # branch of analyses for the refseq import run. Note that in the future there will probably be a better
  # way of doing this
    if($output_hash->{'refseq_ftp_path'} =~ /^ftp/) {
    $conf_file =~ s/##download_refseq_gff##//;
    }

    open(FINAL_CONF,">".$output_path.'/Genome_annotation_conf.pm');
    print FINAL_CONF $conf_file;
    close FINAL_CONF;

 # $cmd = "rm ".$output_path.'/Genome_annotation_conf.pm.tmp';
    system($cmd);
}

sub clean_commands {
    my ($general_hash,$species_name) = @_;
    my $cmds = "";
  # Clean species dir:
    $cmds .= "rm -rf ".$general_hash->{'work_dir'}."/".$species_name."\n";

    my $user = $general_hash->{'user_w'};
    my $pass = $general_hash->{'password'};
    my $port = $general_hash->{'port'};
    my $farm_user = $general_hash->{'farm_user_name'};
    my @db_names = ('core','pipe','genblast','refseq','cdna','ncrna');
    my$db_type = '';
    foreach my $db_name (@db_names) {
	if ($db_name eq 'core' || $db_name eq 'pipe'){
	    $db_type = $db_name;
	} 
	else{
	    $db_type = 'output'
	}
	
	$cmds .= "mysql -u".$user." -p".$pass." -P".$port." -h".$general_hash->{$db_type.'_server'}." -e 'drop database ".$farm_user."_".$species_name."_".$db_name."'\n";
    }

    return($cmds);
}
