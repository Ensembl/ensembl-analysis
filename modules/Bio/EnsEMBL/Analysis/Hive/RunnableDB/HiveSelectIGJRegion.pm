=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectIGJRegion

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectIGJRegion;

use strict;
use warnings;

use Bio::EnsEMBL::IO::Parser::Fasta;

use parent 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB';


sub fetch_input {
  my ($self) = @_;

  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($self->param_required('inputfile'));
  $self->param('fasta_parser', $parser);
}


sub run {
  my ($self) = @_;

  my $count = 0;
  my $parser = $self->param('fasta_parser');
  while ($parser->next) {
    my $seq = $parser->getSequence;
    my $name = $parser->getHeader;
    if ($seq =~ /[GT][GT][GT][TC][TC][TAGC][TG][GC][TGC]\w{21,24}[CGT][AC][CATG][TC][G][TG][GC]\w{3,10}TT[TC]GG\w{4}GG\w/) {
      $self->say_with_header("IGHJ $name");
      if ($seq =~ /[C][A][C][A][GA][TAC][G]\w{21,24}[TAG][C][ACT][GAC][AG][AG][A][C][CG]/) {
        $self->say_with_header("IGHV $name");
        $self->output([$name]);
        ++$count;
      }
    }
    elsif ($seq =~ /[GA][GA][T][T][T][TC][T][G][T]\w{21,24}[C][A][CT][T][G][T][G]\w{3,10}TT[TC]GG\w{4}GG\w/) {
      $self->say_with_header("IGKJ $name ");
      if ($seq =~ /[C][A][C][AT][G][T][G]\w{11,13}[A][C][A][ATGC][A][A][A][CT][C]/) {
        $self->say_with_header("IGKV $name");
        $self->output([$name]);
        ++$count;
      }
    }
    elsif ($seq =~ /[GA][G][TG][T][T][TG][TG][G][T]\w{11,13}[C][A][C][ATC][G][TC][GA]\w{3,10}TT[TC]GG\w{4}GG\w/) {
      $self->say_with_header("IGLJ $name");
      if ($seq =~ /[C][AT][C][AG][G][TA][G]\w{21,24}[AC][CT][ACT][ACGT][AG][A][ACT][CT][CTA]/) {
        $self->say_with_header("IGLV $name");
        $self->output([$name]);
        ++$count;
      }
    }
  }
  $self->say_with_header("$count sequences matches");
}

sub write_output {
  my ($self) = @_;

  my @iids;
  foreach my $header (@{$self->output}) {
    push(@iids, {iid => $header});
  }
  if (@iids) {
    $self->dataflow_output_id(\@iids, $self->param('_branch_to_flow_to'));
  }
  else {
    $self->input_job->autoflow($self->param('_auto_flow'));
  }
}

1;

#if ($seq =~ /[GT][GT][GT][TC][TC][TAGC][TG][GC][TGC]\w{21,24}[CGT][AC][CATG][TC][G][TG][GC]\w{3,10}TT[TC]GG\w{4}GG\w/) {
#elsif ($seq =~ /[GA][GA][T][T][T][TC][T][G][T]\w{21,24}[C][A][CT][T][G][T][G]\w{3,10}TT[TC]GG\w{4}GG\w/) {
#elsif ($seq =~ /[GA][G][TG][T][T][TG][TG][G][T]\w{11,13}[C][A][C][ATC][G][TC][GA]\w{3,10}TT[TC]GG\w{4}GG\w/) {
# salmon
#ssa01
#ssa02
#ssa03
#ssa04
#ssa05
#ssa06
#ssa07
#ssa08
#ssa09
#ssa10
#ssa11
#ssa12
#ssa13
#ssa14
#ssa15
#ssa16
#ssa17
#ssa18
#ssa19
#ssa20
#ssa21
#ssa22
#ssa23
#ssa24
#ssa25
#ssa26
#ssa27
#ssa28
#ssa29
#AGKD04011014.1
#AGKD04012080.1
#AGKD04012485.1
#AGKD04014939.1
#AGKD04015363.1
#AGKD04016225.1
#AGKD04018123.1
#AGKD04018253.1
#AGKD04019579.1
#AGKD04023920.1
#AGKD04023923.1
#AGKD04024231.1
#AGKD04024391.1
#AGKD04024479.1
#AGKD04024749.1
#AGKD04029839.1
#AGKD04032958.1
#AGKD04034868.1
#AGKD04035730.1
#AGKD04035779.1
#AGKD04036018.1
#AGKD04043507.1
#AGKD04044888.1
#AGKD04046147.1
#AGKD04051894.1
#AGKD04059584.1
#AGKD04066888.1
#AGKD04082585.1
#AGKD04107488.1
#AGKD04114817.1
#AGKD04120812.1
#AGKD04138386.1
#AGKD04139384.1
#AGKD04164950.1
#AGKD04195220.1
#AGKD04210049.1
#AGKD04227737.1
#AGKD04228884.1
#AGKD04229233.1
#AGKD04318915.1
#AGKD04323652.1
#AGKD04326048.1
#AGKD04332399.1
#AGKD04369470.1
#AGKD04377089.1
#AGKD04397534.1
#AGKD04404353.1
#AGKD04464111.1
#AGKD04492084.1
#AGKD04517319.1
#AGKD04573689.1
#AGKD04593735.1
#AGKD04661153.1
#AGKD04669115.1
#AGKD04673533.1
#AGKD04705521.1
#AGKD04707035.1
#AGKD04713183.1
#AGKD04727288.1
#AGKD04731014.1
#AGKD04771704.1
#AGKD04819777.1
#
# Lamprey
# PIZI01000001.1
# PIZI01000002.1
# PIZI01000003.1
# PIZI01000004.1
# PIZI01000005.1
# PIZI01000006.1
# PIZI01000007.1
# PIZI01000008.1
# PIZI01000009.1
# PIZI01000010.1
# PIZI01000011.1
# PIZI01000012.1
# PIZI01000013.1
# PIZI01000014.1
# PIZI01000015.1
# PIZI01000016.1
# PIZI01000017.1
# PIZI01000018.1
# PIZI01000019.1
# PIZI01000020.1
# PIZI01000021.1
# PIZI01000022.1
# PIZI01000023.1
# PIZI01000024.1
# PIZI01000025.1
# PIZI01000026.1
# PIZI01000027.1
# PIZI01000028.1
# PIZI01000029.1
# PIZI01000030.1
# PIZI01000031.1
# PIZI01000032.1
# PIZI01000033.1
# PIZI01000034.1
# PIZI01000035.1
# PIZI01000036.1
# PIZI01000037.1
# PIZI01000038.1
# PIZI01000039.1
# PIZI01000040.1
# PIZI01000041.1
# PIZI01000042.1
# PIZI01000043.1
# PIZI01000044.1
# PIZI01000045.1
# PIZI01000046.1
# PIZI01000047.1
# PIZI01000048.1
# PIZI01000049.1
# PIZI01000050.1
# PIZI01000051.1
# PIZI01000052.1
# PIZI01000053.1
# PIZI01000054.1
# PIZI01000055.1
# PIZI01000056.1
# PIZI01000057.1
# PIZI01000058.1
# PIZI01000059.1
# PIZI01000060.1
# PIZI01000061.1
# PIZI01000062.1
# PIZI01000063.1
# PIZI01000064.1
# PIZI01000065.1
# PIZI01000066.1
# PIZI01000070.1
# PIZI01000071.1
# PIZI01000074.1
# PIZI01000080.1
# PIZI01000082.1
# PIZI01000084.1
# PIZI01000086.1
# PIZI01000090.1
# PIZI01000091.1
# PIZI01000096.1
# PIZI01000100.1
# PIZI01000104.1
# PIZI01000107.1
# PIZI01000108.1
# PIZI01000110.1
# PIZI01000121.1
# PIZI01000123.1
# PIZI01000125.1
# PIZI01000140.1
# PIZI01000141.1
# PIZI01000163.1
# PIZI01000172.1
# PIZI01000182.1
# PIZI01000183.1
# PIZI01000187.1
# PIZI01000188.1
# PIZI01000194.1
# PIZI01000227.1
# PIZI01000234.1
# PIZI01000236.1
# PIZI01000239.1
# PIZI01000254.1
# PIZI01000327.1
# PIZI01000340.1
# PIZI01000374.1
# PIZI01000391.1
# PIZI01000463.1
# PIZI01000482.1
# PIZI01000524.1
# PIZI01000558.1
# PIZI01000662.1
# PIZI01000670.1
# PIZI01000709.1
# PIZI01000876.1
# PIZI01000966.1
# PIZI01001137.1
# PIZI01001271.1
# PIZI01001429.1
# PIZI01001443.1
# PIZI01001612.1
# PIZI01001631.1
# PIZI01001685.1
# PIZI01001952.1
# PIZI01002231.1
# PIZI01002528.1
# PIZI01002870.1
# PIZI01003367.1
# PIZI01003560.1
# PIZI01003810.1
# PIZI01004071.1
# PIZI01004144.1
# PIZI01004647.1
# PIZI01004851.1
# PIZI01005394.1
# PIZI01005483.1
# PIZI01005497.1
# PIZI01005859.1
# PIZI01005881.1
# PIZI01005961.1
# PIZI01006183.1
# PIZI01006255.1
# PIZI01006806.1
# PIZI01006855.1
# PIZI01007292.1
# PIZI01008040.1
# PIZI01008186.1
# PIZI01008998.1
# PIZI01010312.1
# PIZI01010867.1
# PIZI01010946.1
# PIZI01011795.1

# if ($seq =~ /GGTTTTTGT\w{21,24}CACTGTG\w{3,10}TT[TC]GG\w{4}GG\w/) {
# Salmon
# ssa01
# ssa03
# ssa06
# ssa07
# ssa12
# ssa17
# AGKD04012080.1
# AGKD04016225.1
# AGKD04023920.1
# AGKD04023923.1
# AGKD04032958.1
# AGKD04035730.1
# AGKD04036018.1
# AGKD04043507.1
# AGKD04727288.1
# AGKD04819777.1
# Lamprey
# 0
#
#
#
#
