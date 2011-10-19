=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils - Utilities
for the tracking system

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils qw(is_evidence_stored get_date setup_pipeline cleanup_meta_tracking);

=head1 DESCRIPTION

  Utilities for the tracking system, checking if the input sequence exists
  in the tracking system, format the date of an entry, setting up the pipeline
  and cleaning up the pipeline.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils;

use strict;
use vars qw(@ISA  @EXPORT_OK);

use constant TREMBL_DB  => 2000;
use constant SWISS_DB   => 2200;
use constant REFSEQ_DNA => 1800;
use constant REFSEQ_RNA => 1801;
use constant REFSEQ_AA  => 1810;
use constant EMBL  => 700;

use LWP::Simple;
use XML::Simple; 
use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;
use Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Database;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;
use Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Reason;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_database_connection_parameters_by_string);
use Data::Dumper;

@ISA = qw(Exporter);

@EXPORT_OK = qw( is_evidence_stored
              get_date
              setup_pipeline
              unlock_tracking
              cleanup_meta_tracking
              get_id_trembl
              get_id_swiss
              get_id_embl
              get_id_dna
              get_id_rna
              get_id_aa );


=head2 is_evidence_stored

 Arg [1]    : string, name of the hit
 Arg [2]    : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor object
 Arg [3]    : string, the molecule type: PROTEIN, MRNA, EST
 Example    : use Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils qw(is_evidence_stored);
              is_evidence_stored($name, $input_seq_adaptor, $molecule_type);
 Description: Check if the protein has an entry in the input_seq table, create one if not
 Returntype : int, 1 if the evidence is stored, dbID of the stored object otherwise
 Exceptions : 


=cut

sub is_evidence_stored {
    my ($name, $input_seq_adaptor, $molecule_type ) = @_;

    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    if ($name !~ /\./ or ! $input_seq_adaptor->fetch_by_hit_name($name)) {
        my $date;
        my $db;
        my $version;
        if ($molecule_type eq 'PROTEIN') {
            if ($name =~ s/^([A-Z0-9]{6})(\.\d)?$/$1/o) {
                my $cmd = 'mfetch -f lau -i acc:'.$name;
                open INF, "$cmd |" || die ('Could not execute mfetch');
                while (<INF>) {
                    if (/(UniProt\S+)/) {
                        $date = &get_date($_);
                        $db = ($_ =~ /Swiss-Prot/o) ? &SWISS_DB : &TREMBL_DB;
                    }
                    elsif (/sequence version (\d+)/) {
                        $version = $1;
                    }
                    else {
                        last;
                    }
                }
                close(INF);
                $name .= '.'.$version;
            }
            else {
                my $efetch_url = $base .'efetch.fcgi?db=protein&id='.$name.'&rettype=gp&retmode=xml';
                my $efetch_out = get($efetch_url);
                my $xmlfile = new XML::Simple;
                my $parsed_xml = $xmlfile->XMLin($efetch_out);
                my $seq = $parsed_xml->{'GBSeq'};
                $name   = $seq->{'GBSeq_accession-version'};
                $date = &get_date($seq->{'GBSeq_create-date'});
                $db = &REFSEQ_AA;
            }
        }
        else {
            if ($name =~ s/^([A-Z]{2}_\S+)/$1/o) {
                my $efetch_url = $base .'efetch.fcgi?db=nucleotide&id='.$name.'&rettype=gb&retmode=xml';
                my $efetch_out = get($efetch_url);
                my $xmlfile = new XML::Simple;
                my $parsed_xml = $xmlfile->XMLin($efetch_out);
                my $seq = $parsed_xml->{'GBSeq'};
                $name   = $seq->{'GBSeq_accession-version'};
                $date = &get_date($seq->{'GBSeq_create-date'});
                $db = &REFSEQ_RNA;
            }
        }
        my $input_seq = Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq->new(
            -hit_name        => $name,
            -molecule_type   => $molecule_type,
            -submission_date => $date,
            -external_db_id  => $db
            );
        return $input_seq_adaptor->store($input_seq);
    }
    return 1;
}


=head2 get_date

 Arg [1]    : string, a date in the format DD-MMM-YYYY (03-JAN-2009)
 Example    : &get_date($line);
 Description: Change the date format from DD-MMM-YYYY to YYYY-MM-DD
 Returntype : a string
 Exceptions :  None


=cut

sub get_date {
    my $line = shift;

    my ($day, $month, $year) = $line =~ /(\d{2})-(\w{3})-(\d{4})/o;
    my %h_months = (
        'JAN' => '01',
        'FEB' => '02',
        'MAR' => '03',
        'APR' => '04',
        'MAY' => '05',
        'JUN' => '06',
        'JUL' => '07',
        'AUG' => '08',
        'SEP' => '09',
        'OCT' => '10',
        'NOV' => '11',
        'DEC' => '12',
        );
    return $year.'-'.$h_months{uc($month)}.'-'.$day;
}


=head2 setup_pipeline

 Arg [1]    : $db, Bio::EnsEMBL::DBSQL::DBAdaptor
 Arg [2]    : $default_runnabledb_path, string
 Arg [3]    : $queue_config, listref from BatchQueue.pm
 Arg [4]    : $ra_analyses_to_run, listref of Bio::EnsEMBL::Analysis::Analysis
 Arg [5]    : $rh_to_keep, hashref of string, logic_name of analysis that we want
               to keep the same analysis run
 Example    : Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils::setup_pipeline($db, $default_runnabledb_path, $queue_config, $ra_analyses_to_run, $rh_to_keep);
 Description: Add the meta key 'tracking.analysis' for each analysis we want to track
 Returntype : 
 Exceptions : throw if not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor object


=cut

sub setup_pipeline {
    my ($db, $default_runnabledb_path, $queue_config, $ra_analyses_to_run, $rh_to_keep) = @_;

    throw('Should be a Bio::EnsEMBL::DBSQL::DBAdaptor object') unless ($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    my $track_db = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
            -dbconn   => $db->dbc
            );

    my $analysisrun_adaptor = $track_db->get_AnalysisRunAdaptor;
    foreach my $analysis (@{$ra_analyses_to_run}) {
        my $module_name = $analysis->module;
        my $logic_name = $analysis->logic_name;
        if (! $module_name) {
            warning('It is not a perl module that run this analysis '.$logic_name);
            next;
        }
        if (exists $rh_to_keep->{$analysis}) {
            next if ($track_db->has_meta_value_by_key('tracking.analysis', $analysis->dbID));
        }
        if ($track_db->has_meta_value_by_key('tracking.analysis', $analysis->dbID)) {
            throw('The analysis '.$analysis->logic_name.' is probably already running but you did not ask to keep the run!');
        }
        $track_db->store_meta_key_value('tracking.analysis', $analysis->dbID);
        my $runnable_path = $default_runnabledb_path;
        foreach my $anal (@{$queue_config}) {
            if ($anal->{'logic_name'} eq $logic_name) {
                $runnable_path = $anal->{'runnabledb_path'}
                if (exists $anal->{'runnabledb_path'});
                last;
            }
        }
        require "$runnable_path/$module_name.pm";
        $runnable_path =~ s'\/'::'g;
        $module_name = $runnable_path.'::'.$module_name;
        my $module = $module_name->new(
                -db => $db,
                -analysis => $analysis,
                -input_id => 'input_id'
                );
#This is hacking!!! depending on the module the name of the input/output DB change...
        my @a_input_dbs = qw(PAF_SOURCE_DB GENE_SOURCE_DB QUERYSEQS);
        my @a_output_dbs = qw(OUTPUT_DB TARGET_DB OUTDB);
        my $current_analysis = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
                -analysis_id => $analysis->dbID,
                -input_dbs   => _get_databases($module, \@a_input_dbs, $track_db),
                -output_dbs  => _get_databases($module, \@a_output_dbs, $track_db)
                );
        $analysisrun_adaptor->store($current_analysis);
    }
}

=head2 check_trackevidences

 Arg [1]    : 
 Example    : $;
 Description: 
 Returntype : 
 Exceptions : 


=cut

sub check_trackevidences {
    my ($rulemanager, $default_runnabledb_path, $queue_config, $ra_analyses_to_run) = @_;

    my $db = $rulemanager->db;
    foreach my $analysis (@{$ra_analyses_to_run}) {
        my $logic_name = $analysis->logic_name;
        foreach my $anal (@{$queue_config}) {
            if ($anal->{'logic_name'} eq $logic_name) {
                if (exists $anal->{'check_tracks'}) {
                    my $runnable_path = exists $anal->{'runnabledb_path'} ? $anal->{'runnabledb_path'} : $default_runnabledb_path;
                    my $module_name = $anal->{'check_tracks'};
                    my $module_path = $module_name;
                    $module_path =~ s'::'/'g;
                    require $module_path.'.pm';
                    $analysis->module($module_name);
                    $logic_name .= '_check';
                    $analysis->logic_name($logic_name);
                    my $module = $module_name->new(
                            -db => $db,
                            -analysis => $analysis,
                            -input_id => 'input_id'
                            );
                    my $job = $rulemanager->create_and_store_job($logic_name, $analysis);
                    $job->batch_runRemote();
                }
            }
        }
    }
}

=head2 _get_databases

 Arg [1]    : $module, a Bio::EnsEMBL::Analysis::RunnableDB object
 Arg [2]    : $ra_dbs, a listref of string representing the databases to look for
 Arg [3]    : $track_db, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor object
 Example    : $ra_input_db_id = _get_databases($module, $ra_dbs, $track_db);
 Description: Get the database used for the analysis according to the list of possible databases 
 Returntype : a listref of Bio::EnsEMBL::Analysis::EvidenceTracking::Database
 Exceptions : 


=cut

sub _get_databases {
    my ($module, $ra_dbs, $track_db) = @_;

    my @a_ids;
    my $database_adaptor = $track_db->get_DatabaseAdaptor;
    foreach my $db (@{$ra_dbs}) {
        if ($db eq 'QUERYSEQS') {
            my $edb = Bio::EnsEMBL::Analysis::EvidenceTracking::Database->new(
                    -db_name => $module->GENOMICSEQS,
                    -instance => 'file'
                    );
            push(@a_ids, $edb);
        }
        else {
            my %h_database;
            if($module->can($db)) {
                my $input_db = $module->$db;
                if (ref($input_db) ne 'HASH') {
                    %h_database = %{get_database_connection_parameters_by_string($input_db)};
                }
                else {
                    %h_database = %{$input_db};
                }
                my $edb = Bio::EnsEMBL::Analysis::EvidenceTracking::Database->new(
                        -db_name => $h_database{-dbname},
                        -instance => $h_database{-host}
                        );
                $database_adaptor->store($edb) unless $edb->is_stored($track_db);
                push(@a_ids, $edb);
            }
        }
    }
    if (! @a_ids) {
        my %h_database = %{get_database_connection_parameters_by_string('REFERENCE_DB')};
        my $edb = Bio::EnsEMBL::Analysis::EvidenceTracking::Database->new(
                -db_name => $h_database{-dbname},
                -instance => $h_database{-host}
                );
        $database_adaptor->store($edb) unless $edb->is_stored($track_db);
        push(@a_ids, $edb);
    }
    return \@a_ids;
}


=head2 cleanup_meta_tracking

 Arg [1]    : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor object
 Example    : Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils::cleanup_meta_tracking($db);
 Description: Update or remove the tracking.analysis meta_key when the analysis is done
 Returntype : 
 Exceptions : throw if not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor object


=cut

sub cleanup_meta_tracking {
    my $db = shift;

    throw('Should be a Bio::EnsEMBL::DBSQL::DBAdaptor object') unless ($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    my $dba = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
        -dbconn => $db->dbc
        );
    my @meta_keys = $dba->get_meta_values_by_key('tracking.analysis');
    while (my $meta_value = shift @meta_keys) {
        if($dba->is_analysis_done($meta_value)) {
            if ($dba->has_meta_value_by_key('tracking.done', $meta_value)) {
                $dba->remove_meta_key('tracking.analysis', $meta_value);
            }
            else {
                $dba->update_meta_key_by_value('tracking.analysis', 'tracking.done', $meta_value) ;
            }
        }
    }
}


=head2 unlock_tracking

 Arg [1]    : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor object
 Example    : Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils::unlock_tracking($db);
 Description: Remove the tracking.analysis meta_key for all analysis
 Returntype : 
 Exceptions : throw if not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor object


=cut

sub unlock_tracking {
    my $db = shift;

    throw('Should be a Bio::EnsEMBL::DBSQL::DBAdaptor object') unless ($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    my $dba = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
        -dbconn => $db->dbc
        );
    my @meta_keys = $dba->get_meta_values_by_key('tracking.analysis');
    while (my $meta_value = shift @meta_keys) {
        $dba->remove_meta_key('tracking.analysis', $meta_value);
    }
}

=head2 update_from_list

 Arg [1]    : $ra_to_update, listref of String, the the ids of the sequence to update
 Arg [2]    : $db, Bio::EnsEMBL::DBSQL::DBAdaptor object, the pipeline database
 Arg [3]    : $code, String, the reason code 
 Example    : update_from_list($ra_to_update, $db, $code);
 Description: Update the reason of the sequences given in parameters. It write directly
              in the database!
 Returntype : void
 Exceptions : 


=cut

sub update_from_list {
    my ($ra_to_update, $db, $reason_code) = @_;

    my $dba = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
        -dbconn => $db->dbc);
    my $dummy_analysisrun = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
        -dbID => 0);
    my $reason = $dba->get_ReasonAdaptor->fetch_by_code($reason_code);
    my $evidencetrack_adaptor = $dba->get_EvidenceTrackAdaptor;
    my $inputseq_adaptor = $dba->get_InputSeqAdaptor;
    foreach my $acc (@{$ra_to_update}) {
        my $evidence = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
            -input_seq => $inputseq_adaptor->fetch_by_hit_name($acc),
            -is_aligned => 'n');
        my $evidence_track = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
            -input_id => 'inputseq',
            -evidence => $evidence,
            -analysis_run => $dummy_analysisrun,
            -reason => $reason);
        $evidencetrack_adaptor->store($evidence_track);
    }
}

=head2 update_from_files

 Arg [1]    : $db, object, the pipeline database
 Arg [2]    : $code, string, the reason code
 Arg [3]    : $fileA, string, path to a file
 Arg [4]    : $fileB, string, path to a file
 Example    : update_from_files($db, $code, $fileA, $fileB);
 Description: Update the reason for a list of sequence that
              are missing in $fileA compare to $fileB
              $fileA must have the lowest number of lines
 Returntype : void
 Exceptions : 


=cut

sub update_from_files {
    my ($db, $reason_code, $file1, $file2) = @_;

    open(FR, $file1) || die('Could not open '.$file1);
    my @a_lines = <FR>;
    close(FR);
    my %h_ids;
    while(my $line = shift @a_lines) {
        if ($line =~ s/^>(\S+).*/$1/) {
            chomp $line;
            $h_ids{$line} = $line;
        }
    }
    my @a_ids;
    open(FR, $file2) || die('Could not open '.$file2);
    while(my $line = <FR>) {
        if ($line =~ s/^>(\S+).*/$1/) {
            chomp $line;
            if (!exists $h_ids{$line}) {
                push(@a_ids, $line);
            }
        }
    }
    close(FR);
    &update_from_list(\@a_ids, $db, $reason_code);
}

=head2 get_id_trembl

 Example    : my $ext_id = get_id_trembl;
 Description: Getter for the external_db id of TrEMBL
 Returntype : Int
 Exceptions : 


=cut

sub get_id_trembl {
    return &TREMBL_DB;
}

=head2 get_id_swiss

 Example    : my $ext_id = get_id_swiss;
 Description: Getter for the external_db id of Swiss-Prot
 Returntype : Int
 Exceptions : 


=cut

sub get_id_swiss {
    return &SWISS_DB;
}

=head2 get_id_dna

 Example    : my $ext_id = get_id_dna;
 Description: Getter for the external_db id of TrEMBL
 Returntype : Int
 Exceptions : 


=cut

sub get_id_dna {
    return &REFSEQ_DNA;
}

=head2 get_id_embl

 Example    : my $ext_id = get_id_embl;
 Description: Getter for the external_db id of ENA
 Returntype : Int
 Exceptions : 


=cut

sub get_id_embl {
    return &EMBL;
}

=head2 get_id_rna

 Example    : my $ext_id = get_id_rna;
 Description: Getter for the external_db id of RefSeq DNA
 Returntype : Int
 Exceptions : 


=cut

sub get_id_rna {
    return &REFSEQ_RNA;
}

=head2 get_id_aa

 Example    : my $ext_id = get_id_aa;
 Description: Getter for the external_db id of RefSeq protein
 Returntype : Int
 Exceptions : 


=cut

sub get_id_aa {
    return &REFSEQ_AA;
}


1;
