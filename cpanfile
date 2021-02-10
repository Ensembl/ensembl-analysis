requires 'Bio::DB::HTS';
requires 'Path::Tiny';
requires 'LWP::UserAgent';
requires 'Proc::ProcessTable';

on 'test' => sub {
  requires 'Test::More', '>= 0.96, < 2.0';
  requires 'Test::Most', '>= 0.37';
};
