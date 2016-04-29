use strict;
use warnings;

use Test::More;
use FindBin '$Bin';

printf("Running in $Bin\n");

eval 'use Test::Code::TidyAll';
plan skip_all => "Test::Code::TidyAll required to check if the code is clean." if $@;

tidyall_ok( conf_file => 'tidyallrc' );

