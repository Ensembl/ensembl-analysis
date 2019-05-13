=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::ConfigWriter

=head1 SYNOPSIS

  my $bq_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
      -modulename => 'Bio::EnsEMBL::Analysis::Config::Databases',
      -is_example => 1,
     );
  $bq_cfg->key_by_parent('REFERENCE_DB', '-dbname', 'homo_sapiens_core_68_37');
  $bq_cfg->write_config(1);

=head1 DESCRIPTION

  This module update or create config files from a template (.example) or from
  an existing config file.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Tools::ConfigWriter;

use warnings ;
use strict;
use Data::Dumper;
use File::Copy;
use File::Path;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA) ;
@ISA = qw();

=head2 new

  Arg [-MODULENAME] : String, module name
  Arg [-IS_EXAMPLE] : Int, 1 or 0, use the example file like Bio/EnsEMBL/Analysis/Config/Databases.pm.example
  Arg [-MODULEDIR]  : String, path to the root of the library, ie '/path/ensembl-analysis/modules'.
                        Look in @INC by default
  Arg [-BACKUPDIR]  : String, path to the directory which will store the backup. It's the module directory by
                        default.
                       NOTE: It will write the module like: /my/backup/path/Databases.pm.123333332
  Arg [-VERBOSE]    : Int
  Function          : Create the object and parse the config file
  Returntype        : Bio::EnsEMBL::Analysis::Tools::ConfigWriter;
  Exceptions        : Throw if no module name is given
  Example           : $cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(-modulename => 'Bio::EnsEMBL::Analysis::Config::Databases', -is_example => 1);

=cut

sub new {
  my ($class,@args) = @_;
  my $self = bless {},$class;

  my ($config_module, $is_example, $moddir, $backup_dir, $verbose) = rearrange (['MODULENAME', 'IS_EXAMPLE', 'MODULEDIR', 'BACKUPDIR', 'VERBOSE'], @args);

  throw("No module name provided...\n") unless $config_module;
  verbose('EXCEPTION');
  verbose($verbose) if $verbose;
  $self->modulename($config_module);
  $self->moduledir($moddir);
  $self->is_example($is_example);
  $self->parse($config_module);
  if ($backup_dir) {
      $self->backupdir($backup_dir);
  }
  else {
      $self->backupdir($self->moduledir);
  }
  return $self;
}

=head2 parse

  Function  : Parse the config file
  Returntype: 
  Exceptions: Throw if it can't open the file
  Example   : $cfg->parse();

=cut

sub parse {
    my $self = shift;

    my $path = $self->modulename;
    $path .= '.example' if ($self->is_example);
    if ($self->moduledir and -e $self->moduledir.'/'.$path) {
        $path = $self->moduledir.'/'.$path;

    }
    else {
        foreach my $tmppath (@INC) {
            if (-e $tmppath.'/'.$path) {
               $path = $tmppath.'/'.$path;
               $self->moduledir($tmppath) unless $self->moduledir;
               last;
            }
        }
    }

    my $head;
    my $config_hash;
    my $open_config = 0;
    my $import;
    my $header = 1;
    my $hash = 0;
    my %Config;

    open(FH, $path) || throw("Could not open $path\n");
    while(<FH>) {
        if (/%Config\s*=\s*\(/) {
            $header = 0;
            $hash = 1;
        }

        elsif(/%Config\s*=/) {
            $header = 0;
            $open_config = 1;
        }

       # The following two conditionals are for cases when the first opening
       # bracket is not on the same line as Config\s*=
       elsif($open_config && $_ =~ /\s*\(/) {
            $open_config = 0;  
            $hash = 1;
        }

        elsif($open_config) {
          # Discard these lines 
        }

        elsif (/\s*\)\s*;/ and $hash) {
            $hash = 0;
        }

        elsif ($header) {
            $head .= $_;
        }

        elsif ($hash) {
            $config_hash .= $_;
        }

        else {
            $import .= $_;
        }
    }
    close(FH);

    unless($config_hash) {
        throw("Could not parse config hash out of module correctly, check config structure!\n");
    }

    unless($import) {
        throw("Could not parse import block out of module correctly, check config structure!\n");
    }


    $self->header($head);
    $self->tail($import);
    eval '%Config = ('.$config_hash.');';

    if($@) {
        throw("Problem with eval of config hash, could be a problem with config structure!".
              "\n\nError from eval: ".$@."\n");
    }
  
    $self->config(\%Config);
}

=head2 write_config

  Function  : Write the configuration file
  Returntype: 
  Exceptions: Throw if it can't write the file
  Example   : $cfg->write_config;

=cut

sub write_config {
    my ($self, $backup_if_exists) = @_;

    local $Data::Dumper::Indent = 1;
    local $Data::Dumper::Quotekeys = 0;
    local $Data::Dumper::Deepcopy = 1;
    local $Data::Dumper::Sortkeys = \&_sort_keys;
    my $config_hash = Dumper($self->config);
    $config_hash =~ s/^\$VAR1 = \{(.*)\};/%Config = ($1);/s;
    $config_hash =~ s/'([^']+)' =>/$1 =>/gs;
    $config_hash =~ s/'(\d+)'/$1/gs;
    $config_hash =~ s/'undef'/undef/gs;
    $config_hash =~ s/\\'//gs;
    $config_hash =~ s/\\\\([^\\])/\\$1/gs;
    $config_hash =~ s/ '"([^"']*)"'/ '$1'/gs;
    my $backup;
    $backup = $self->backup if ($backup_if_exists and -e $self->get_module_path);
    open(FW, '>'.$self->get_module_path) || throw('Could not open '.$self->get_module_path."\n");
    print FW $self->header, $config_hash, $self->tail;
    close(FW);
    return $backup;
}

=head2 backup

  Function  : Back up the configuration file by adding the Unix time as an extension
  Returntype: String, path to back up file
  Exceptions: Throw if it can't copy the file
  Example   : $cfg->backup;

=cut

sub backup {
    my ($self) = @_;

    my $modulename = $self->get_module_path;
    print STDERR $modulename, "\n";
    my $backup_name = $modulename;
    if ($self->backupdir and $self->backupdir ne $self->moduledir) {
        $self->modulename =~ /([^\/]+)$/;
        $backup_name = $self->backupdir.'/'.$1;
    }
    $backup_name .= '.'.time;
    copy($modulename, $backup_name) || throw('Could not back up '.$modulename."\n");
    return $backup_name;
}

=head2 create_path

  Function  : Create the path to the module, uses $self->moduledir
  Returntype: 
  Exceptions: Throw if it can't create the path
  Example   : $cfg->create_path;

=cut

sub create_path {
    my $self = shift;

    my ($path) = $self->get_module_path =~ '(.*)/';
    eval {
        mkpath($path);
    };
    throw('Could not create '.$path) if ($@);
}

=head2 key_by_parent

  Arg[1]    : String, name of the father
  Arg[2]    : String, name of the hash key
  Arg[3]    : (optional) Value, the value to assign to the key, NO check is performed!!!
  Function  : Get the value for a key knowing it's parent, to avoid ambiguity
  Returntype: Whatever value is there
  Exceptions:
  Example   : $cfg->key('REFERENCE_DB', '-dbname');

=cut

sub key_by_parent {
    my ($self, $parent, $key, $value) = @_;

    my $father = _key_exists($self->config, $parent);
    if (exists $$father->{$key}) {
        $$father->{$key} = $value if ($value);
        return $$father->{$key};
    }
    return;
}

=head2 delete_key_by_parent

  Arg[1]    : String, name of the father
  Arg[2]    : String, name of the hash key
  Function  : Delete the value for a key knowing it's parent, to avoid ambiguity
  Returntype: The deleted key
  Exceptions:
  Example   : $cfg->delete_key_by_parent('REFERENCE_DB', '-pass');

=cut

sub delete_key_by_parent {
    my ($self, $parent, $key) = @_;

    my $father = _key_exists($self->config, $parent);
    return delete $$father->{$key};
}

=head2 root_value

  Arg[1]    : String, name of the hash key
  Arg[2]    : (optional) Value, Whatever data you want to add, NO check is done
  Function  : Add or retrieve a value to a key for an analysis. Runnable config file
  Returntype: Whatever data is there, it will be a reference for an array or a hash
  Exceptions:
  Example   : $cfg->root_value('DNA_DBNAME', 'REFERENCE_DB');

=cut

sub root_value {
    my ($self, $key, $value) = @_;

    $self->config->{$key} = $value if ($value);
    return $self->config->{$key} if (exists $self->config->{$key});
}

=head2 delete_root_key

  Arg[1]    : String, name of the hash key
  Function  : Delete the value for a root key.
  Returntype: The deleted key
  Exceptions:
  Example   : $cfg->delete_root_key('DNA_DBNAME');

=cut

sub delete_root_key {
    my ($self, $key) = @_;

    return delete $self->config->{$key};
}

=head2 value_by_logic_name

  Arg[1]    : String, logic name
  Arg[2]    : String, key name
  Arg[3]    : (optional) Value, Whatever data you want to add, NO check is done
  Function  : Add or retrieve a value to a key for an analysis. Runnable config file
  Returntype: Whatever data is there, it will be a reference for an array or a hash
  Exceptions:
  Example   : $cfg->value_by_logic_name('cdna_update', 'OUTDB', 'EXONERATE_DB');

=cut

sub value_by_logic_name {
    my ($self, $logic_name, $key, $value) = @_;

    foreach my $hkey  (keys %{$self->config}) {
        next unless $hkey =~ /CONFIG_BY_LOGIC/;
        next unless (exists $self->config->{$hkey}->{$logic_name});
        my $tmp = _key_exists($self->config->{$hkey}->{$logic_name}, $key);
        if ($tmp) {
            $$tmp = $value if ($value);
            return $$tmp;
        }
        else {
            $self->config->{$hkey}->{$logic_name}->{$key} = $value;
            return $self->config->{$hkey}->{$logic_name}->{$key};
        }
    }
    return;
}

=head2 default_value

  Arg[1]    : String, name of the hash key
  Arg[2]    : (optional) Value, Whatever data you want to add, NO check is done
  Function  : Add or retrieve a value to a key for an analysis. Runnable config file
  Returntype: Whatever data is there, it will be a reference for an array or a hash
  Exceptions:
  Example   : $cfg->default_value('OUTDB', 'GENEBUILD_DB');

=cut

sub default_value {
    my ($self, $key, $value) = @_;

    return $self->value_by_logic_name('DEFAULT', $key, $value);
}

=head2 copy_analysis_from_config

  Arg[1]    : String, logic name from an existing analysis
  Arg[2]    : String, logic name of the new analysis
  Function  : Copy an analysis, useful when you want to add several (almost) identical
                analysis.
  Returntype: None
  Exceptions: Throw if the "old" analysis does not exists or if the "new" one exists
  Example   : $cfg->copy_analysis_from_config('cdna_update', 'cdna_update_2');

=cut

sub copy_analysis_from_config {
    my ($self, $logic_name, $new_analysis) = @_;

    throw( "You need to provide a new logic_name!\n") unless ($new_analysis);

    foreach my $hkey (%{$self->config}) {
 
        next unless $hkey =~ /CONFIG_BY_LOGIC/;
        my $tmp = _key_exists($self->config->{$hkey}, $logic_name);
        if (!$tmp) {
            throw ("The analysis $logic_name does not exists!\n");
        }
        else {
            my $tmp_hash = _clone($$tmp);
            $self->add_analysis_to_config($new_analysis, $tmp_hash);
        }
    }
}

=head2 add_analysis_to_config

  Arg[1]    : String, logic name of the analysis to add
  Arg[2]    : (optional) Value, Whatever data you want to add, NO check is done
  Function  : Add an analysis to the config file for a RunnableDB module
  Returntype: None
  Exceptions: Throw if the analysis exists
  Example   : my $new_hash = _clone(\%hash);

=cut

sub add_analysis_to_config {
    my ($self, $logic_name, $key) = @_;

    foreach my $hkey (%{$self->config}) {
        next unless $hkey =~ /CONFIG_BY_LOGIC/;
        my $tmp = _key_exists($self->config->{$hkey}, $logic_name);
        if ($tmp) {
            throw ("The analysis $logic_name already exists!\n");
        }
        else {
            if ($key) {
                throw("Should be a hash for $logic_name\n") unless (ref($key) eq 'HASH');
                $self->config->{$hkey}->{$logic_name} = $key;
            }
            else {
                $self->config->{$hkey}->{$logic_name} = {};
            }
        }
    }
}

=head2 delete_analysis

  Arg[1]    : String, logic name of the analysis
  Function  : Delete an analysis from a config file
  Returntype: The deleted key
  Exceptions:
  Example   : $cfg->delete_analysis('cdna_update');

=cut

sub delete_analysis {
    my ($self, $logic_name) = @_;

    foreach my $hkey (%{$self->config}) {
        next unless $hkey =~ /CONFIG_BY_LOGIC/;
        return delete $self->config->{$hkey}->{$logic_name};
    }
}

=head2 analysis_from_batchqueue

  Arg[1]    : String, logic name of the analysis
  Arg[2]    : String, key to retrieve
  Arg[3]    : (optional) Value, Whatever data you want to add, NO check is done
  Function  : Add or retrieve a value to an analysis to/from the batch queue file
  Returntype: None
  Exceptions:
  Example   : $self->analysis_from_batchqueue('cdna_update');

=cut

sub analysis_from_batchqueue {
    my ($self, $logic_name, $key, $value) = @_;

    throw("Need a key to set $value\n") if ($value and !$key);
    foreach my $qc_hash (@{$self->queue_config}) {
        next unless ($qc_hash->{logic_name} eq $logic_name);
        if ($value) {
            $qc_hash->{$key} = $value;
        }
        elsif ($key) {
            return $qc_hash->{$key};
        }
        return $qc_hash;
    }
}

=head2 add_analysis_to_batchqueue

  Arg[1]    : String, logic name of the analysis
  Arg[2]    : HashRef, hash reference
  Function  : Add an analysis to the batch queue file
  Returntype: None
  Exceptions:
  Example   : $self->add_analysis_to_batchqueue('cdna_update');

=cut

sub add_analysis_to_batchqueue {
    my ($self, $logic_name, $value) = @_;

    throw("You must provide a logic_name\n") unless $logic_name;
    foreach my $qc_hash (@{$self->queue_config}) {
        throw("You already have a $logic_name analysis in your BatchQueue.pm!\n") if ($qc_hash->{logic_name} eq $logic_name);
    }
    if ($value) {
        throw('You should pass a hash reference instead of a '.ref($value)) if (ref($value) ne 'HASH');
        $value->{logic_name} = $logic_name;
        push(@{$self->queue_config}, $value);
    }
    else {
        push(@{$self->queue_config}, {logic_name => $logic_name});
    }
}

=head2 copy_analysis_from_batchqueue

  Arg[1]    : String, logic name from an existing analysis
  Arg[2]    : String, logic name of the new analysis
  Function  : Copy an analysis, useful when you want to add several (almost) identical
                analysis.
  Returntype: None
  Exceptions: Throw if the "old" analysis does not exists or if the "new" one exists
  Example   : $cfg->copy_analysis_from_batchqueue('cdna_update', 'cdna_update_2');

=cut

sub copy_analysis_from_batchqueue {
    my ($self, $logic_name, $new_logic_name) = @_;

    my $analysis;
    foreach my $qc_hash (@{$self->queue_config}) {
        throw("You already have a $new_logic_name analysis in your BatchQueue.pm!\n") if ($qc_hash->{logic_name} eq $new_logic_name);
        $analysis = $qc_hash if ($qc_hash->{logic_name} eq $logic_name);
    }
    throw("Could not find analysis with logic_name $logic_name") unless $analysis;
    $self->add_analysis_to_batchqueue($new_logic_name,  _clone($analysis));
}

=head2 delete_analysis_key_from_batchqueue

  Arg[1]    : String, name of the father
  Arg[2]    : String, name of the hash key
  Function  : Delete the value for a key knowing it's parent, to avoid ambiguity
  Returntype: The deleted key
  Exceptions:
  Example   : $cfg->delete_analysis_key_from_batchqueue('cdna_update', 'OUTDB');

=cut

sub delete_analysis_key_from_batchqueue {
    my ($self, $logic_name, $key) = @_;

    my $analysis;
    foreach my $qc_hash (@{$self->queue_config}) {
        $analysis = $qc_hash if ($qc_hash->{logic_name} eq $logic_name);
    }
    throw("Could not find analysis with logic_name $logic_name") unless $analysis;
    return delete $analysis->{$key};
}

=head2 empty_queue_config

  Arg[1]    : (optional) ArrayRef, list of analysis to keep
  Function  : Delete all the analysis from the batch queue file except
                the ones passed in argument
  Returntype: None
  Exceptions:
  Example   : $cfg->empty_queue_config(['cdna_update', 'cdna_update_2']);

=cut

sub empty_queue_config {
    my ($self, $config_to_keep) = @_;

    if ($config_to_keep) {
        my %hash = map { $_ => $_ } @$config_to_keep;
        my @config_array = grep { exists $hash{$_->{logic_name}} } @{$self->queue_config};
        $self->queue_config(\@config_array);
    }
    else {
        $self->queue_config([]);
    }
}

=head2 delete_databases

  Arg[1]    : (optional) Array, list of database names to keep
  Function  : Delete from the database hash all the databases reference
  Returntype: None
  Exceptions:
  Example   : $cfg->delete_databases('REFERENCE_DB', 'EXONERATE_DB');

=cut

sub delete_databases {
    my ($self, $arrayref) = @_;
    my %to_keep = map { $_ => $_ } @$arrayref;

    foreach my $key (keys %{$self->get_databases}) {
        delete $self->get_databases->{$key} unless (exists $to_keep{$key});
    }
}


###################
# Getters/Setters #
###################

=head2 get_databases

  Function  : Get the DATABASES hash
  Returntype: Hashref
  Exceptions:
  Example   : $cfg->get_databases;

=cut

sub get_databases {
    my $self = shift;

    return $self->config->{DATABASES};
}

=head2 queue_config

  Function  : Get the QUEUE_CONFIG hash
  Returntype: Hashref
  Exceptions:
  Example   : $cfg->queue_config;

=cut

sub queue_config {
    my ($self, $value) = @_;

    if ($value) {
        $self->config->{QUEUE_CONFIG} = $value;
    }
    return $self->config->{QUEUE_CONFIG};
}

sub get_module_path {
    my $self = shift;

    return $self->moduledir.'/'.$self->modulename;
}

sub header {
    my ($self, $arg) = @_;

    if ($arg) {
        $self->{'_config_header'} = $arg;
    }
    return $self->{'_config_header'};
}

sub config {
    my ($self, $h_arg) = @_;

    if ($h_arg) {
        throw('Not a hash ref element!') unless (ref($h_arg) eq 'HASH');
        $self->{'_config_hash'} = $h_arg;
    }
    return $self->{'_config_hash'};
}

sub tail {
    my ($self, $arg) = @_;

    if ($arg) {
        $self->{'_config_import'} = $arg;
    }
    return $self->{'_config_import'};
}

sub moduledir {
    my ($self, $arg) = @_;

    if ($arg) {
        $self->{'_config_moddir'} = $arg;
    }
    return $self->{'_config_moddir'};
}

sub is_example {
    my ($self, $arg) = @_;

    if ($arg) {
        $self->{'_config_is_example'} = $arg;
    }
    return $self->{'_config_is_example'};
}

sub modulename {
    my ($self, $arg) = @_;

    if ($arg) {
        $arg =~ s'::'/'g;
        $self->{'_config_modname'} = $arg.'.pm';
    }
    return $self->{'_config_modname'};
}

sub backupdir {
    my ($self, $arg) = @_;

    if ($arg) {
        $self->{'_config_backupdir'} = $arg;
    }
    return $self->{'_config_backupdir'};
}


###################
# Private methods #
###################

=head2 _key_exists

  Arg[1]    : Hashref, reference to hash
  Arg[2]    : String, name of the hash key
  Function  : Check if this key exists, return the first one found
  Returntype: Boolean
  Exceptions:
  Example   : _key_exists($hash, 'DATABASES');

=cut

sub _key_exists {
    my ($hash, $key) = @_;

    foreach my $hkey (keys %$hash) {
        return \$hash->{$hkey} if ($key eq $hkey);
        if (ref($hash->{$hkey}) eq 'HASH') {
            my $res = _key_exists($hash->{$hkey}, $key);
            return $res if ($res);
        }
    }
    return;
}

=head2 _clone

  Arg[1]    : Reference, reference to a hash to be deep cloned
  Function  : Deep cloning
  Returntype: HashRef
  Exceptions: 
  Example   : my $new_hash = _clone(\%hash);

=cut

sub _clone {
    my $hashref = shift;

    my %hash;
    foreach my $key (keys %$hashref) {
        if (ref($hashref->{$key}) eq 'ARRAY') {
            my @tmp = @{$hashref->{$key}};
            $hash{$key} = \@tmp;
        }
        elsif (ref($hashref->{$key}) eq 'HASH') {
            $hash{$key} = _clone($hashref->{$key});
        }
        else {
            $hash{$key} = $hashref->{$key};
        }
    }
    return \%hash;
}

sub _sort_keys {
    my ($hash) = @_;

    my @res = sort _sort_hash_key keys %$hash;
    return \@res;
}

sub _sort_hash_key {

    return -1 if ($a eq 'DEFAULT'
               or $a eq 'logic_name'
               or $a eq 'REFERENCE_DB'
               or $b eq 'QUEUE_CONFIG'
               or $b eq 'KILL_LIST_DB');
    return 1 if ($b eq 'DEFAULT'
              or $b eq 'logic_name'
              or $b eq 'REFERENCE_DB'
              or $a eq 'QUEUE_CONFIG'
              or $a eq 'KILL_LIST_DB');
    return $a cmp $b;
}

1;
