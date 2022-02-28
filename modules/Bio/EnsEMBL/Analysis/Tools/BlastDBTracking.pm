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

package Bio::EnsEMBL::Analysis::Tools::BlastDBTracking;

=head1 NAME - Bio::EnsEMBL::Analysis::Tools::BlastDBTracking

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::BlastDBTracking;

  my $bdbt    = Bio::EnsEMBL::Analysis::Tools::BlastDBTracking->new;
  my $entry   = $bdbt->by_filename('Pfam-A.hmm');
  my $version = $entry->version;

=head1 DESCRIPTION

A replacement for the orphaned BlastableVersion module, for querying
the blastdb tracking database.  For a description of the database, see
http://mediawiki.internal.sanger.ac.uk/wiki/index.php/Blastdbtrackingsystem.

=head1 AUTHOR

Michael Gray B<email> mg13@sanger.ac.uk

=cut

use warnings ;
use strict ;
use namespace::autoclean;
use Moose;
use MooseX::ClassAttribute;

use DBI;
use File::Basename qw(basename);
use File::Spec::Link;
use File::stat;
use POSIX qw(strftime);
use Readonly;
use Try::Tiny;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::Entry;

# Based on tjrc's BlastableVersion.pm, which at time of writing was
# installed in /software/perl-5.8.8/lib/site_perl/

Readonly my $DEFAULT_SYSTEM       => 'farm2';

Readonly my $DEFAULT_BLASTDB_HOST => 'cbi5d';
Readonly my $DEFAULT_BLASTDB_PORT => undef;
Readonly my $DEFAULT_BLASTDB_USER => 'blastdbro';
Readonly my $DEFAULT_BLASTDB_PASS => undef;
Readonly my $DEFAULT_BLASTDB_NAME => 'blastdb';

=head2 new

  All arguments are optional and have sensible defaults.

    system - system for which the blastdb is being queried (default: farm)

    db_host, db_user, db_pass, db_bame - blastdb connection parameters

    fallback_to_file_stat - if set and no blastdb entry, generate version from file mtime

=cut

# new() provided by Moose

has system  => ( is => 'ro', isa => 'Str', default => $DEFAULT_SYSTEM );

has db_host => ( is => 'ro', isa => 'Str',        default => $DEFAULT_BLASTDB_HOST );
has db_port => ( is => 'ro', isa => 'Maybe[Str]', default => $DEFAULT_BLASTDB_PORT );
has db_user => ( is => 'ro', isa => 'Str',        default => $DEFAULT_BLASTDB_USER );
has db_pass => ( is => 'ro', isa => 'Maybe[Str]', default => $DEFAULT_BLASTDB_PASS );
has db_name => ( is => 'ro', isa => 'Str',        default => $DEFAULT_BLASTDB_NAME );

has fallback_to_file_stat => ( is => 'ro', isa => 'Bool' );

class_has _Fake_version_for_otf   => ( is => 'rw', isa => 'Bool' );

has _dbh => (
    is       => 'ro',
    isa      => 'Object',
    builder  => '_connect',
    lazy     => 1,
    init_arg => undef,
    );

has _by_filename_sth => (
    is       => 'ro',
    isa      => 'Object',
    builder  => '_prepare_by_filename_sth',
    lazy     => 1,
    init_arg => undef,
    );

sub _connect {
    my $self = shift;
    my $dsn = sprintf("DBI:mysql:host=%s;database=%s", $self->db_host, $self->db_name);
    $dsn .= sprintf(";port=%d", $self->db_port) if $self->db_port;
    return DBI->connect($dsn, $self->db_user, $self->db_pass, { RaiseError => 1 });
}

sub _prepare_by_filename_sth {
    my $self = shift;
    return $self->_dbh->prepare(qq(
                                 SELECT f.filename,
                                        db.version,
                                        db.sanger_version,
                                        UNIX_TIMESTAMP(s.installation) AS installation,
                                        c.count,
                                        f.checksum
                                   FROM file      f
                                   JOIN component c ON c.com_id = f.com_id
                                   JOIN db          ON db.db_id = c.db_id
                                   JOIN system    s ON s.db_id  = db.db_id
                                  WHERE
                                        f.filename  = ?
                                    AND s.system    = ?
                                    AND s.available ='yes'
                                 ORDER BY s.installation DESC
                                 LIMIT 1
        ));
}

=head2 by_filename

  Search for entry by filename. 

=cut

sub by_filename {
    my ($self, $filename) = @_;

    # Resolve all symbolic links in the given filename
    $filename = File::Spec::Link->resolve_all($filename);

    my $sth = $self->_by_filename_sth;

    $sth->execute(basename($filename), $self->system);
    my $row = $sth->fetchrow_hashref;
    $sth->finish;

    if ($row) {
        return Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::Entry->new($row);
    }

    return unless $self->fallback_to_file_stat;

    return $self->_entry_from_file($filename);
}

sub _entry_from_file {
    my ($self, $filename) = @_;

    my $st = stat($filename);
    throw("Could not stat $filename: $!\n") unless $st;

    return Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::Entry->new({
        filename     => basename($filename),
        version      => strftime('%d-%b-%y', localtime($st->mtime)),
        installation => $st->mtime,
        count        => -1,
        from_file    => 1,
                                   });
}

=head2 get_db_version_mixin

  NOT a member function!

  This is designed to provide a minimal-weight replacement for get_db_version in various
  Bio::EnsEMBL::Analysis::Runnable[DB]::Finished modules.

  Example usage:

  sub get_db_version {
    my ($self, $db) = @_;
    return Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::get_db_version_mixin(
        $self, '_db_version_searched', $db);
  }

=cut

sub get_db_version_mixin {
    my ($caller, $attrib_key, $db_file ) = @_;

    if (__PACKAGE__->_Fake_version_for_otf) {
        return "1.0.0";
    }

    unless ( $caller->{$attrib_key} ) {
        if ($db_file) {

            my $entry;
            try {
                my $bdbt = Bio::EnsEMBL::Analysis::Tools::BlastDBTracking->new;
                $entry = $bdbt->by_filename($db_file);
            }
            catch {
                throw("Failed to get a BlastDBTracking entry for '$db_file': '$_'");
            };
            throw("BlastDBTracking entry not found for '$db_file'") unless $entry;

            $caller->{$attrib_key} = $entry->version;

        } else {
            throw("'$attrib_key' not set yet. I need to be called with a database filename first.");
        }
    }
    return $caller->{$attrib_key};
}

sub DEMOLISH {
    my $self = shift;
    $self->_dbh->disconnect if $self->_dbh;
    return;
}

1;
