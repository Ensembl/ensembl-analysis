package Bio::EnsEMBL::Analysis::Provenance::Logger;
use strict;
use warnings;
use JSON;
use File::Path qw(make_path);
use Time::HiRes qw(gettimeofday);


sub new {
  my ($class, %args) = @_;
  my $self = bless {
    log_dir           => $args{log_dir} || 'logs', # default log directory - should typically be analysis directory
    json              => JSON->new->utf8->canonical->pretty(0)->allow_blessed(1)->convert_blessed(1),
  }, $class;

  $self->{log_dir} =~ s/\/$//; # remove trailing slash if present
  $self->_initialize_directories();

  return $self;
}

sub log {
    my ($self, $stage, $data) = @_;
    
    # Add timestamp if not provided
    $data->{timestamp} ||= $self->_timestamp();
    
    # Add stage to data
    $data->{stage} = $stage;
    
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
    
    # Create stage directory if needed
    my $dir = "$self->{log_dir}/by_stage";
    make_path($dir) unless -d $dir;
    
    # Write to log file
    $self->_write_to_log_file($stage, $json_string);
    
    return 1;
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


###########################
##### Private methods #####
###########################

sub _initialize_directories {
    my ($self) = @_;
    
    my $base_dir = $self->{log_dir};
    
    # Create main directories
    foreach my $dir (
        "$base_dir/by_stage",
        "$base_dir/indexes"
    ) {
        eval { make_path($dir) };
        if ($@) {
            die "Could not create log directory $dir: $@";
        }
    }
}

sub _write_to_log_file {
    my ($self, $stage, $json_string) = @_;
    
    my $file = "$self->{log_dir}/by_stage/$stage.log";
    
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