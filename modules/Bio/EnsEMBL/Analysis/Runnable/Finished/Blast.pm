package Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;

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

=head1 NAME - Bio::EnsEMBL::Analysis::Runnable::Finished::Blast

=head1 DESCRIPTION

Unlike Bio::EnsEMBL::Analysis::Runnable::Blast,
this module creates FeaturePairs from HSPs after
doing any depth filtering to save time and memory
when searching genomic sequences that generate
large numbers of blast matches.

=head2 usage of Bio::EnsEMBL::Analysis::Config::Blast with this module


BLAST_CONFIG =>
        {
            Uniprot =>
            {
            BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::Finished::BPliteWrapper',
            PARSER_PARAMS => {
                               -regex => '^\w+\s+(\w+)',
                               -query_type => 'pep',
                               -database_type => 'pep',
                               -threshold_type => 'PVALUE',
                               -threshold => 0.01,
                   -coverage => 10,
                   -discard_overlaps => 1,

                              },
            BLAST_FILTER =>'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
            FILTER_PARAMS => {
                              -min_score => 200,
                              -prune => 1,
                             },
            BLAST_PARAMS => {
                             -unknown_error_string => 'FAILED',
                             -type => 'wu',
                            },
            },
            DEFAULT =>
            {
             BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper',
             PARSER_PARAMS => {
                               -regex => '^(\w+)',
                               -query_type => undef,
                               -database_type => undef,
                              },
             BLAST_FILTER => undef,
             FILTER_PARAMS => {},
             BLAST_PARAMS => {
                              -unknown_error_string => 'FAILED',
                              -type => 'wu',
                             }
            },

           BLAST_AB_INITIO_LOGICNAME => 'Genscan'
       }


=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

modified by Sindhu K. Pillai B<email>sp1@sanger.ac.uk
=cut

use strict;
use warnings;
use Symbol;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Tools::BlastDBTracking;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use base ("Bio::EnsEMBL::Analysis::Runnable::Blast");

$ENV{BLASTMAT}      = $MAT_DIR;
$ENV{WUBLASTMAT}    = $MAT_DIR;
$ENV{BLASTFILTER}   = $BIN_DIR;
$ENV{WUBLASTFILTER} = $BIN_DIR;
$ENV{BLASTDB}       = $BLAST_DIR;

BEGIN {
    print STDERR "\nUSING " . __PACKAGE__ . "\n\n";
}

sub get_analysis {

    my ($self) = @_;

    my ($ana);
    unless ( $ana = $self->{'_analysis'} ) {
        my ($source) = $self->program =~ m{([^/]+)$}
          or throw(
            "Can't parse last element from path: '" . $self->program . "'" );
        $ana = $self->{'_analysis'} = Bio::EnsEMBL::Analysis->new(
            -db              => $self->database,
            -db_version      => 1,                 # ARUUGA!!!
            -program         => $source,
            -program_version => 1,
            -gff_source      => $source,
            -gff_feature     => 'similarity',
            -logic_name      => 'blast',
        );
    }
    return $ana;
}

sub run_analysis {
  my ($self) = @_;

  DB:foreach my $database (@{$self->databases}) {

    my $db = $database;
    $db =~ s/.*\///;
    $db =~ s/\W//g;
    #allow system call to adapt to using ncbi blastall.
    #defaults to WU blast
    my $command  = $self->program;
    my $blastype = "";
    my $filename = $self->queryfile;
    my $results_file = $self->create_filename($db, 'blast.out');
    $self->files_to_delete($results_file);
    $self->results_files($results_file);
    if ($self->type eq 'ncbi') {
      $command .= " -d $database -i $filename ";
    } else {
      $command .= " $database $filename gi ";
    }
    $command .= $self->options. ' 2>&1 > '.$results_file;

    info("Running blast ".$command);

    my $fh;
    unless (open($fh, "$command |")) {
      $self->delete_files;
      throw("Error opening Blast cmd <$command>." .
            " Returned error $? BLAST EXIT: '" .
            ($? >> 8) . "'," ." SIGNAL '" . ($? & 127) .
            "', There was " . ($? & 128 ? 'a' : 'no') .
            " core dump");
    }

    # this loop reads the STDERR from the blast command
    # checking for FATAL: messages (wublast) [what does ncbi blast say?]
    # N.B. using simple die() to make it easier for RunnableDB to parse.
    while(<$fh>){
      if(/FATAL:(.+)/){
        my $match = $1;
        next DB if $match =~ /There is nothing in the requested database to search/;
        # clean up before dying
        $self->delete_files;
        if($match =~ /no valid contexts/){
          die qq{"VOID"\n}; # hack instead
        }elsif($match =~ /Bus Error signal received/){
          die qq{"BUS_ERROR"\n}; # can we work out which host?
        }elsif($match =~ /Segmentation Violation signal received./){
          die qq{"SEGMENTATION_FAULT"\n}; # can we work out which host?
        }elsif($match =~ /Out of memory;(.+)/){
          # (.+) will be something like "1050704 bytes were last
          #requested."
          die qq{"OUT_OF_MEMORY"\n};
          # resenD to big mem machine by rulemanager
        }elsif($match =~ /the query sequence is shorter than the word length/){
          #no valid context
          die qq{"VOID"\n}; # hack instead
        }else{
          warning("Something FATAL happened to BLAST we've not ".
                  "seen before, please add it to Package: "
                  . __PACKAGE__ . ", File: " . __FILE__."\n[$match]\n");
          die ($self->unknown_error_string."\n");
          # send appropriate string
          #as standard this will be failed so job can be retried
          #when in pipeline
        }
      }elsif(/WARNING:(.+)/){
        # ONLY a warning usually something like hspmax=xxxx was exceeded
        # skip ...
      }elsif(/^\s{10}(.+)/){ # ten spaces
        # Continuation of a WARNING: message
        # Hope this doesn't catch more than these.
        # skip ...
      }
    }
    unless(close $fh){
      # checking for failures when closing.
      # we should't get here but if we do then $? is translated
      #below see man perlvar
      $self->delete_files;
      warning("Error running Blast cmd <$command>. Returned ".
              "error $? BLAST EXIT: '" . ($? >> 8) .
              "', SIGNAL '" . ($? & 127) . "', There was " .
              ($? & 128 ? 'a' : 'no') . " core dump");
      die ($self->unknown_error_string."\n");
    }
  }
}

sub parse_results {

    my ($self)           = @_;
    my $results          = $self->results_files;
    my $bplites          = $self->parser->parse_files($results);
    my $threshold_type   = $self->parser->threshold_type;
    my $threshold        = $self->parser->threshold;
    my $discard_overlaps = $self->parser->discard_overlaps;
    my $coverage         = $self->parser->coverage;
    my $hits             =
      $self->parser->get_best_hits( $bplites, $threshold_type, $threshold );
    my $query_length = $self->query->length
      or throw("Couldn't get query length");
    my $output =
      $self->parser->_apply_coverage_filter( $query_length, $hits,
        $threshold_type, $threshold, $coverage, $discard_overlaps );
    $self->output($output);
    return $output;
}

sub clean_databases {
    my ( $self) = @_;
    $self->{'_databases'} = [];
}

sub clean_results_files {
  my ($self) = @_;
  $self->{'results_files'} = [];
}

sub databases {
    my ($self, $val) = @_;

    if ($val) {

        # Trim any leading and trailing space
        $val =~ s/(^\s+|\s+$)//g;
        if ($val =~ /\s/) {

            # We have a list of databases to pass to blast as a virtual db
            my @db_list = map { m{^/} ? $_ : "$ENV{BLASTDB}/$_" } split /\s+/, $val;
            $self->{'_databases'} = [qq{'@db_list'}];
            $self->get_db_version($db_list[0]);     ### Should check that all dbs are at same version!
        }
        else {
            $self->{'_databases'} = [];
            my @databases;
            foreach my $dbname (split(/,/, $val)) {

                # allows the use of a comma separated list in $self->database
                # prepend the environment variable $BLASTDB if
                # database name is not an absoloute path
                unless ($dbname =~ m!^/!) {
                    $dbname = $ENV{BLASTDB} . "/" . $dbname;
                }

                # If the expanded database name exists put this in
                # the database array.
                #
                # If it doesn't exist then see if $database-1,$database-2 exist
                # and put them in the database array
                if (-f $dbname) {
                    push(@databases, $dbname);
                }
                else {
                    my $count = 1;
                    my $db_filename;
                    while (-f ($db_filename = "${dbname}-${count}")) {
                        push(@databases, $db_filename);
                        $count++;
                    }
                    $! = undef;    # to stop pollution as it will be "No such file or directory" after while loop above.
                }
            }
            if (scalar(@databases) == 0) {
                throw("No databases exist for " . $val);
            }
            else {
                foreach my $db_name (@databases) {
                    $self->get_db_version($db_name) if $db_name =~ /emnew_/;
                }
                $self->get_db_version($databases[0]);
                push @{ $self->{'_databases'} }, @databases;
            }
        }
    }

    return $self->{'_databases'};
}


=head2 get_db_version

    Title   :  get_db_version
               [ distinguished from RunnableDB::*::db_version_searched() ]
    Useage  :  $self->get_db_version('/data/base/path')
               $obj->get_db_version()
    Function:  Set a blast database version from the supplied path
               Get a blast database version from previously supplied path
    Returns :  String
    Args    :  String (should be a full database path)
    Caller  :  $self::fetch_databases()
               RunnableDB::Finished_EST::db_version_searched()

=cut

sub get_db_version {
    my ($self, $db) = @_;
    my $ver = Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::get_db_version_mixin($self, '_db_version_searched', $db,);

    printf STDERR "B:E:A:Runnable::Finished::Blast::get_db_version '%s' => '%s'\n", $db || '<read_only>', $ver;

    return $ver;
}

sub DESTROY
{
    my ( $self ) = @_;
    # just do the cleanup
    $self->delete_files;
}

1;

__END__

