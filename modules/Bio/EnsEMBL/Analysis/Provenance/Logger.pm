package Bio::EnsEMBL::Analysis::Provenance::Logger;
use strict;
use warnings;
use JSON;
use File::Path qw(make_path);
use Time::HiRes qw(gettimeofday);
use POSIX qw(strftime);

sub new {
  my ($class, %args) = @_;
  
  # Generate run ID based on current datetime
  my $run_id = $args{run_id} || $class->_generate_run_id();
  
  my $self = bless {
    log_dir           => $args{log_dir} || 'logs',
    run_id            => $run_id,
    json              => JSON->new->utf8->canonical->pretty(0)->allow_blessed(1)->convert_blessed(1),
  }, $class;

  $self->{log_dir} =~ s/\/$//; # remove trailing slash if present
  
  # Create run-specific directory structure
  $self->{run_dir} = "$self->{log_dir}/runs/$run_id";
  
  $self->_initialize_directories();

  return $self;
}

sub log {
    my ($self, $stage, $data) = @_;
    
    # Add timestamp if not provided
    $data->{timestamp} ||= $self->_timestamp();
    
    # Add stage and run_id to data
    $data->{stage} = $stage;
    $data->{run_id} = $self->{run_id};
    
    # Convert to JSON
    my $json_string;
    eval {
        $json_string = $self->{json}->encode($data);
    };
    if ($@) {
        # Handle serialization errors gracefully
        warn "JSON encoding error: $@";
        # Create a safe version of the data with problematic objects removed or simplified
        $data = $self->_sanitize_for_json($data);
        $json_string = $self->{json}->encode($data);
    }
    
    # Write to both run-specific and stage-specific log files
    $self->_write_to_log_file($stage, $json_string);
    
    return 1;
}

# Get the current run directory path
sub get_run_dir {
    my ($self) = @_;
    return $self->{run_dir};
}

# Get the current run ID
sub get_run_id {
    my ($self) = @_;
    return $self->{run_id};
}

# List all previous runs
sub list_runs {
    my ($self) = @_;
    my $runs_dir = "$self->{log_dir}/runs";
    return [] unless -d $runs_dir;
    
    opendir(my $dh, $runs_dir) or die "Could not open $runs_dir: $!";
    my @runs = grep { -d "$runs_dir/$_" && $_ !~ /^\./ } readdir($dh);
    closedir($dh);
    
    return [sort @runs];
}

# ## Example Call 
# $logger->log(
#     'MiniMap',
#     {
#         feature_id => 'ENSG00000139618',
#         feature_type => 'gene',
#         result => 'hit',
#         details => {
#             identity => 98.7,
#             coverage => 99.2,
#             next_stage => 'ProteinValidation'
#         },
#         message => 'MiniMap alignment successful'
#     }
# );

# ## Example with custom run_id:
# my $logger = Bio::EnsEMBL::Analysis::Provenance::Logger->new(
#     log_dir => '/path/to/logs',
#     run_id => 'manual_run_001'
# );

###########################
##### Private methods #####
###########################

sub _generate_run_id {
    my ($class) = @_;
    
    # Create run ID with format: YYYYMMDD_HHMMSS_microseconds
    my ($seconds, $microseconds) = gettimeofday();
    my $datetime = strftime("%Y%m%d_%H%M%S", localtime($seconds));
    my $micro_suffix = sprintf("%06d", $microseconds);
    
    return "${datetime}_${micro_suffix}";
}

sub _initialize_directories {
    my ($self) = @_;
    
    my $base_dir = $self->{log_dir};
    my $run_dir = $self->{run_dir};
    
    # Create main directories
    foreach my $dir (
        "$base_dir/runs",
        "$run_dir/by_stage",
        "$run_dir/indexes",
        "$base_dir/current" # symlink target for current run
    ) {
        eval { make_path($dir) };
        if ($@) {
            die "Could not create log directory $dir: $@";
        }
    }
    
    # Create/update symlink to current run
    my $current_link = "$base_dir/current";
    if (-l $current_link) {
        unlink($current_link) or warn "Could not remove existing symlink: $!";
    }
    
    # Create relative symlink to current run
    my $relative_path = "runs/" . $self->{run_id};
    symlink($relative_path, $current_link) or warn "Could not create symlink to current run: $!";
    
    # Write run metadata
    $self->_write_run_metadata();
}

sub _write_run_metadata {
    my ($self) = @_;
    
    my $metadata = {
        run_id => $self->{run_id},
        start_time => $self->_timestamp(),
        log_dir => $self->{log_dir},
        run_dir => $self->{run_dir},
        perl_version => $^V ? $^V->stringify : $],
        hostname => $ENV{HOSTNAME} || `hostname 2>/dev/null` || 'unknown',
        user => $ENV{USER} || $ENV{USERNAME} || 'unknown',
        working_directory => `pwd 2>/dev/null` || 'unknown'
    };
    
    # Clean up newlines
    chomp($metadata->{hostname});
    chomp($metadata->{working_directory});
    
    my $metadata_file = "$self->{run_dir}/run_metadata.json";
    my $json_string = $self->{json}->encode($metadata);
    
    open(my $fh, ">", $metadata_file) or die "Could not create metadata file $metadata_file: $!";
    print $fh $json_string . "\n";
    close($fh) or die "Could not close metadata file: $!";
}

sub _write_to_log_file {
    my ($self, $stage, $json_string) = @_;
    
    my $file = "$self->{run_dir}/by_stage/$stage.log";
    
    # Append to file (each line is a complete JSON object)
    open(my $fh, ">>", $file) or die "Could not open $file for appending: $!";
    print $fh $json_string . "\n";
    close($fh) or die "Could not close $file: $!";
}

sub _timestamp {
    my ($self) = @_;
    
    my ($seconds, $microseconds) = gettimeofday();
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = gmtime($seconds);
    
    return sprintf("%04d-%02d-%02dT%02d:%02d:%02d.%03dZ",
        $year + 1900, $mon + 1, $mday, $hour, $min, $sec, int($microseconds / 1000));
}

# Helper method to sanitize objects that can't be serialized to JSON
sub _sanitize_for_json {
    my ($self, $data) = @_;
    
    if (ref($data) eq 'HASH') {
        my %clean_data;
        foreach my $key (keys %$data) {
            if (defined $data->{$key}) {
                if (ref($data->{$key})) {
                    $clean_data{$key} = $self->_sanitize_for_json($data->{$key});
                } else {
                    $clean_data{$key} = $data->{$key};
                }
            } else {
                $clean_data{$key} = undef;
            }
        }
        return \%clean_data;
    }
    elsif (ref($data) eq 'ARRAY') {
        my @clean_data;
        foreach my $item (@$data) {
            if (defined $item) {
                if (ref($item)) {
                    push @clean_data, $self->_sanitize_for_json($item);
                } else {
                    push @clean_data, $item;
                }
            } else {
                push @clean_data, undef;
            }
        }
        return \@clean_data;
    }
    elsif (ref($data)) {
        # Handle blessed objects or other references
        return "$data"; # Convert object to string representation
    }
    
    return $data;
}

1;