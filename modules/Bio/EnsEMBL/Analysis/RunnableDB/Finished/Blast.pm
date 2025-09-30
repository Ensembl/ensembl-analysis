# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Finished::Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $blast   = Bio::EnsEMBL::Analysis::RunnableDB::Blast->new (
                                                    -db         => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
$blast->fetch_input();
$blast->run();
$blast->output();
$blast->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::Blast to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::Blast;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use base (
    'Bio::EnsEMBL::Analysis::RunnableDB::Blast',
    'Bio::EnsEMBL::Analysis::RunnableDB::Finished'
);

sub fetch_input {

    my ($self) = @_;
    my $slice =
      $self->fetch_sequence( $self->input_id, $self->db,
        $ANALYSIS_REPEAT_MASKING,
        $self->analysis->parameters =~ /lcmask/,    # Set flag for soft masking if lcmask is set
         );
    $self->query($slice);
    my %blast  = %{ $self->BLAST_PARAMS };
    my $parser = $self->make_parser;
    my $filter;
    if ( $self->BLAST_FILTER ) {
        $filter = $self->make_filter;
    }

    # Incremental updating of the embl blast db analysis
    # The embl blast dbs are made up of release files embl_*
    # and update files emnew_*. This block of code makes
    # sure that the analysis is only run against new version of either
    # of these files.

    my @files = split(",", $self->analysis->db_file);
    my @patches;
    if($files[-1] =~ /^embl_/){
        my $search_only_patch = 0;
        my $sic = $self->db->get_StateInfoContainer;
        my $db_version_saved = $sic->fetch_db_version($self->input_id, $self->analysis);
        my $db_version_current = $self->analysis->db_version;
        if($db_version_saved) {
            # split the embl blast db version "12-Mar-06 (85)" to
            # patch version "12-Mar-06" and release version "85"
            my ($patch_sv,$release_sv) = $db_version_saved =~ /^(\S+)\s+\((\d+)\)$/;
            my ($patch_cv,$release_cv) = $db_version_current =~ /^(\S+)\s+\((\d+)\)$/;
            if($release_sv eq $release_cv){
                $search_only_patch = 1;
                print STDOUT "blast db files [ @files ] version $release_sv already searched\n";
                # Just to make sure that nothing is going wrong with the incremental updating...
                throw("Problem with the embl blast db incremental updating, saved and current version identical !\n
                   saved [$db_version_saved] = current [$db_version_current]\n") unless($patch_sv ne $patch_cv)
            }
        }
        foreach my $file (@files) {
            my $patch_file = $file;
            $patch_file =~ s/^embl_/emnew_/g;
            $search_only_patch ? $file = $patch_file : push @patches,$patch_file;
        }
    }
    $self->analysis->db_file(join(",",@files,@patches));

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::Blast->new(
        -query    => $self->query,
        -program  => $self->analysis->program,
        -parser   => $parser,
        -filter   => $filter,
        -database => $self->analysis->db_file,
        -analysis => $self->analysis,
        %blast,
    );
    my $s = $self->runnable($runnable);
    return 1;

}

sub _createfiles {

    my ( $self, $dirname, $filenames ) = @_;
    my $unique = {};
    $unique = { map { $_, $unique->{$_}++ } @$filenames };
    my @files = ();
    $dirname ||= '/tmp';
    $dirname =~ s!(\S+)/$!$1!;
    foreach my $file (@$filenames) {
        if ( $unique->{$file} ) {

            #name not unique add random
            $file .= ".$$." . int( rand(200) );
            push( @files, "$dirname/$file" );
        }
        else {

            #name was unique just add it
            push( @files, "$dirname/$file.$$" );
        }
    }

    return @files;
}

=head2 run

    Title   :   write_output
    Usage   :   $self->write_output();
    Function:   Writes Features , hit_descriptions to database by calling parent write_output method
    Returns :   none
    Args    :   none

=cut

sub write_output {
    my ($self) = @_;
    $self->Bio::EnsEMBL::Analysis::RunnableDB::Finished::write_output();
    return 1;
}

=head2 db_version_searched

    Title   :  db_version_searched
               [ distinguished from Runnable::*::get_db_version() ]
    Useage  :  $self->db_version_searched('version string')
               $obj->db_version_searched()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished::Blast
               This is just a holding place for the string in this
               module
    Returns :  String or undef
    Args    :  String
    Caller  :  $self::run()
               Job::run_module()

=cut

sub db_version_searched {

    my ( $self, $arg ) = @_;
    $self->{'_db_version_searched'} = $arg if $arg;
    return $self->{'_db_version_searched'};

}

1;
