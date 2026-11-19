
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs;

use warnings;
use strict;
use feature 'say';
use File::Spec;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
    my ($self) = @_;

    return {
        %{ $self->SUPER::param_defaults },
        decompress   => 0,
        create_faidx => 0,
    };
}

sub write_output {
    my ($self)       = @_;
    my $ftp_base_url = $self->param('ftp_base_url');
    my $fastq        = $self->param('iid');
    my $path         = $self->param('input_dir');
    my $srr;

    if ( -e $path . '/' . $fastq && $self->file_is_valid( $path, $fastq ) ) {
        if ( $self->param('decompress') ) {
            $fastq = $self->decompress( $path, $fastq );
        }

        if ( $self->param('create_faidx') ) {
            $self->create_faidx( $path, $fastq );
        }

        $self->complete_early('Input file already exists, will not download');
    }

    if ( $fastq =~ m/_/ ) {
        $srr = ( split /_/, $fastq )[0];
    }
    else {
        $srr = ( split /\./, $fastq )[0];
    }
    my $first  = substr $srr, 0, 6;
    my $second = '00' . ( substr $srr, -1, 1 );
    my $third  = '0' .  ( substr $srr, -2, 2 );

    my @urls = (
        "$ftp_base_url/$first/$second/$srr/$fastq",
        "$ftp_base_url/$first/$srr/$fastq",
        "$ftp_base_url/$first/$third/$srr/$fastq",
    );
    my $tmp = File::Spec->catfile( $path, ".${fastq}.download.$$." );
    my $res;

    # Never let wget write to the final filename.  The rename below is atomic,
    # so we can only see a complete download.
    for my $url (@urls) {
        unlink $tmp if -e $tmp;
        $res = $self->run_system_command(
            [
                'wget',          '-qq',
                '--tries=3',     '--timeout=60',
                '--waitretry=5', '-O',
                $tmp,            $url,
            ]
        );
        my $exit_code = $res >> 8;
        if ( !$res && -s $tmp && $self->file_is_valid( $path, $fastq, $tmp ) ) {
            last;
        }
        unlink $tmp if -e $tmp;
        $res = $exit_code << 8;
    }

    if ( !$res && -s $tmp && $self->file_is_valid( $path, $fastq, $tmp ) ) {

      # Remove an index belonging to an older version before replacing the file.
        unlink "$path/$fastq.fai" if -e "$path/$fastq.fai";
        rename( $tmp, "$path/$fastq" )
          or $self->throw(
"Could not move completed download into place: $tmp -> $path/$fastq: $!"
          );
    }
    else {
        unlink $tmp if -e $tmp;
        $res >>= 8  if $res;
        $self->throw(
            "Could not download a valid file $fastq; wget exit code is $res");
    }

    unless ( -e $path . '/' . $fastq ) {
        $self->throw(
                "Did not find the fastq file on the expected path. Path:\n"
              . $path . "/"
              . $fastq );
    }

    if ( $self->param('decompress') ) {
        $fastq = $self->decompress( $path, $fastq );
    }

    if ( $self->param('create_faidx') ) {
        $self->create_faidx( $path, $fastq );
    }
}

sub file_is_valid {
    my ( $self, $path, $fastq, $file ) = @_;
    $file ||= "$path/$fastq";

    return 0 unless -s $file;
    return 1 unless $file =~ /\.gz$/;

    # gzip -t catches truncated downloads which still have a final filename.
    my $res = $self->run_system_command( [ 'gzip', '-t', '--', $file ] );
    return !$res;
}

sub decompress {
    my ( $self, $path, $fastq ) = @_;
    my $cmd = 'gunzip -f ' . $path . '/' . $fastq;

    # Remove this in case indexing in the code block after this one
    if ( $fastq =~ s/\.gz$// ) {
        my $gunzip_res = system($cmd);
        if ($gunzip_res) {
            $self->throw( "Failed to decompress file. Command:\n" . $cmd );
        }
    }
    else {
        $self->warning(
"You selected decompress, but the file did not have a .gz extension, so will not try and decompress"
        );
    }

    # Update these in case the extension was removed
    $self->param( 'iid',        $fastq );
    $self->param( 'fastq_file', $fastq );

    return ($fastq);
}

sub create_faidx {
    my ( $self, $path, $fastq ) = @_;

    if ( -e $path . '/' . $fastq . '.fai' ) {
        $self->warning(
"You selected faidx, but the faidx exists, so will not try and create"
        );
        return;
    }

    my $cmd =
      $self->param_required('samtools_path') . ' faidx ' . $path . '/' . $fastq;
    my $faidx_res = system($cmd);
    if ($faidx_res) {
        $self->throw( "Failed to index file. Command:\n" . $cmd );
    }
}

1;
