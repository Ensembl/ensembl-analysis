package Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;

use strict;
use warnings;
use BlastableVersion;
use Symbol;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use base ("Bio::EnsEMBL::Analysis::Runnable::Blast");

$ENV{BLASTMAT} 		= $MAT_DIR;
$ENV{WUBLASTMAT}	= $MAT_DIR;
$ENV{BLASTFILTER} 	= $BIN_DIR;
$ENV{WUBLASTFILTER} = $BIN_DIR;
$ENV{BLASTDB} 		= $BLAST_DIR;

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
      $command .= " $database $filename -gi ";
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
        }elsif($match =~ /the query sequence is shorter
               than the word length/){
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
	$self->{databases} = [];
}

sub clean_results_files {
  my ($self) = @_;
  $self->{'results_files'} = [];
}

sub databases {
	my ( $self, @vals ) = @_;

	if ( not exists $self->{databases} ) {
		$self->{databases} = [];
	}

	foreach my $val (@vals) {
		my $db_names = $val;
		my @databases;
		$db_names =~ s/\s//g;

		foreach my $dbname ( split( ",", $db_names ) )
		{    # allows the use of a comma separated list in $self->database
			    # prepend the environment variable $BLASTDB if
			    # database name is not an absoloute path
			unless ( $dbname =~ m!^/! ) {
				$dbname = $ENV{BLASTDB} . "/" . $dbname;
			}

			# If the expanded database name exists put this in
			# the database array.
			#
			# If it doesn't exist then see if $database-1,$database-2 exist
			# and put them in the database array
			if ( -f $dbname ) {
				push( @databases, $dbname );
			}
			else {
				my $count = 1;
				my $db_filename;
				while ( -f ( $db_filename = "${dbname}-${count}" ) ) {
					push( @databases, $db_filename );
					$count++;
				}
				$! = undef
				  ; # to stop pollution as it will be "No such file or directory" after while loop above.
			}
		}
		if ( scalar(@databases) == 0 ) {
			throw( "No databases exist for " . $db_names );
		} else {
			foreach my $db_name (@databases){
				$self->get_db_version($db_name) if $db_name =~ /emnew_/;
			}
			$self->get_db_version($databases[0]);
			push @{$self->{databases}}, @databases;
		}
	}

	return $self->{databases};

}

=head2 get_db_version

    Title   :  get_db_version
               [ distinguished from RunnableDB::*::db_version_searched() ]
    Useage  :  $self->get_db_version('/data/base/path')
               $obj->get_db_version()
    Function:  Set a blast database version from the supplied path
               Get a blast database version from previously supplied path
               Uses tjrc''s BlastableVersion module.
    Returns :  String
    Args    :  String (should be a full database path)
    Caller  :  $self::fetch_databases()
               RunnableDB::Finished_EST::db_version_searched()

=cut

sub get_db_version {
	my ( $self, $db ) = @_;
	my $debug_this = 0;    # this just shows debug info.
	my $force_dbi  = 0;    # this will force a dbi call SLOW!!!!!!
	unless ( $self->{'_db_version_searched'} ) {
		if ($db) {
			$BlastableVersion::debug = $debug_this;
			warning
			  "BlastableVersion is cvs revision ".$BlastableVersion::revision." \n"
			  if $debug_this;

			my $ver = eval {
				my $blast_ver = BlastableVersion->new();
				$blast_ver->set_hostname('farm2');
				$blast_ver->force_dbi($force_dbi);    # if set will be SLOW.
				$blast_ver->get_version($db);
				$blast_ver;
			};
			throw("I failed to get a BlastableVersion for $db") if $@;

			my $dbv  = $ver->version();
			my $sgv  = $ver->sanger_version();
			my $name = $ver->name();
			my $date = $ver->date();
			unless ($dbv) {
				throw(  "I know nothing about $db I tried to find out:\n"
					  . " - name <". $name . ">\n"
					  . " - date <". $date . ">\n"
					  . " - version <". $dbv . ">\n"
					  . " - sanger_version <". $sgv. ">\n" );
			}
			$self->{'_db_version_searched'} = $dbv;
		}
		else {
			throw(  "You've asked about what I searched, but I don't know."
				  . " It's not set. I need to be called with a database filename first"
			);

			# The code probably got here because of a problem with the MLDBM
			# cache file on the machine this was running on.
			# the cache file is stored @ /var/tmp/blast_versions
			# try <rm -f /var/tmp/blast_versions>
		}
	}
	return $self->{'_db_version_searched'};
}

sub DESTROY
{
	my ( $self ) = @_;
	# just do the cleanup
	$self->delete_files;
}

1;

__END__

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

