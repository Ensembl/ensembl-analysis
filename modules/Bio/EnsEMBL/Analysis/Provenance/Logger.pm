

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
    json              => JSON->new->utf8->canonical->pretty(0),
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
    my $json_string = $self->{json}->encode($data);
    
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