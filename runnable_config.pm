@RUNNABLE_CONFIG = ();
$ANALYSIS_WORK_DIR = '/nfs/users/nfs_t/th3/enscode/git_test/ensembl-analysis';
$ANALYSIS_INPUT_DIR = '/nfs/users/nfs_t/th3/enscode/git_test/ensembl-analysis';
$ANALYSIS_TARGET_DIR = '/nfs/users/nfs_t/th3/enscode/git_test/ensembl-analysis';sub import { my ($callpack) = caller(0); my $pack = shift; my @vars = @_ ? @_ : keys(%Config); return unless @vars; eval "package $callpack; use vars qw(" . join(" ", map { "\$".$_ } @vars) . ")"; die $@ if $@; foreach (@vars) { if (defined $Config{ $_ }) { no strict "refs"; *{"${callpack}::$_"} = \$Config{ $_ }; } else { die "Error: Config: $_ not known\n"; } } }
1;
