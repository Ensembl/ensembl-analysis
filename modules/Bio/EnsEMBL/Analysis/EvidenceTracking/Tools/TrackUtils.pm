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

Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils;

use strict;
use warnings;
use constant TREMBL_DB  => 2000;
use constant SWISS_DB   => 2200;
use constant REFSEQ_RNA => 1810;
use constant REFSEQ_AA  => 1810;

use LWP::Simple;
use XML::Simple; 
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;
use Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Database;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_database_connection_parameters_by_string);

=head2 is_evidence_stored

 Arg [1]    : string, name of the hit
 Arg [2]    : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor object
 Arg [3]    : string, the molecule type: PROTEIN, MRNA, EST
 Example    : Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils::is_evidence_stored($name, $input_seq_adaptor, $molecule_type);
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
                        $db = ($_ =~ /Swiss-Prot/o) ? &SWISS_DB : &TREMBL_DB;
                    }
                    elsif (/sequence version (\d+)/) {
                        $version = $1;
                    }
                    else {
                        $date = &get_date($_);
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
        my $input_seq = new Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq (
            -hit_name        => $name,
            -molecule_type   => $molecule_type,
            -submission_date => $date,
            -external_db_id  => $db
            );
        return $input_seq_adaptor->store($input_seq);
    }
#    return $_->dbID;
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
    return $year.'-'.$h_months{$month}.'-'.$day;
}

=head2 

 Arg [1]    : 
 Example    : $;
 Description: 
 Returntype : 
 Exceptions : 


=cut

sub setup_pipeline {
    my ($db, $default_runnabledb_path, $queue_config, $ra_analyses_to_run, $rh_to_keep) = @_;

    print STDERR "In setup\n";
    my $track_db = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
            -dbconn   => $db->dbc
            );

    print STDERR "Adaptor created\n";
    my $analysisrun_adaptor = $track_db->get_AnalysisRunAdaptor;
    my $runnable_path = $default_runnabledb_path;
    $runnable_path =~ s'/'::'g;
    foreach my $analysis (@{$ra_analyses_to_run}) {
    print STDERR 'meta fill', "\n";
        print STDERR $analysis->logic_name, "\n";
        if (exists $rh_to_keep->{$analysis}) {
            next if ($track_db->has_meta_value_by_key('tracking.analysis', $analysis->dbID));
        }
        print STDERR "Store meta key\n";
        $track_db->store_meta_key_value('tracking.analysis', $analysis->dbID);
        print STDERR "foreach...\n";
        my $module_name = $analysis->module;
        my $logic_name = $analysis->logic_name;
        if (! $module_name) {
            print STDERR 'It is not a perl module that run this analysis '.$logic_name, "\n";
            next;
        }
        foreach my $anal (@{$queue_config}) {
            if ($anal->{'logic_name'} eq $logic_name) {
                $runnable_path = $anal->{'runnabledb_path'}
                if (exists $anal->{'runnabledb_path'});
                last;
            }
        }
        my ($manme) = $module_name =~ /([^:]+)$/;
        require "$default_runnabledb_path/$module_name.pm";
        $module_name = $runnable_path.'::'.$module_name;
        my $module = $module_name->new(
                -db => $db,
                -analysis => $analysis,
                -input_id => 'input_id'
                );
        print STDERR "get databases...\n";
        my @a_input_dbs = qw(PAF_SOURCE_DB GENE_SOURCE_DB);
        my $ra_input_db_id = get_databases($module, \@a_input_dbs, $track_db);
        my @a_output_dbs = qw(OUTPUT_DB TARGET_DB);
        my $ra_output_db_id = get_databases($module, \@a_output_dbs, $track_db);
        print STDERR "Create AnalysisRun\n";
        my $current_analysis = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
                -analysis_id => $analysis->dbID,
                -input_db_id => join(':', @{$ra_input_db_id}),
                -output_db_id => join(':', @{$ra_output_db_id})
                );
        print STDERR Dumper($current_analysis);
        $analysisrun_adaptor->store($current_analysis);
    }
}

sub get_databases {
    my ($module, $ra_dbs, $track_db) = @_;

    my %h_ids;
    my $database_adaptor = $track_db->get_DatabaseAdaptor;
    foreach my $db (@{$ra_dbs}) {
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
            $h_ids{$edb->dbID} = $edb->dbID;
        }
    }
    my @a_ids = keys %h_ids;
    return \@a_ids;
}

=head2 

 Arg [1]    : 
 Example    : $;
 Description: 
 Returntype : 
 Exceptions : 


=cut

sub cleanup_meta_tracking {
    my $db = shift;

    my $dba = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
        -dbcon => $db->dbc
        );
    my @meta_keys = $dba->get_meta_values_by_key('tracking.analysis');
    while (my $meta_value = shift @meta_keys) {
        if(!@{$db->get_StateInfoContainer->list_input_ids_by_analysis($meta_value)}) {
            $dba->update_meta_key_by_value('tracking.analysis', 'tracking.done', $meta_value);
        }
    }
}

1;
