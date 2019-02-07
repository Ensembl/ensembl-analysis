=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::IGDFinder

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::IGDFinder;

use strict;
use warnings;

use parent 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB';

sub fetch_input {
  my ($self) = @_;

  my $dna_db = $self->get_database_by_name('dna_db');
  my $db = $self->get_database_by_name('target_db', $dna_db);
  $self->query($self->fetch_sequence($self->input_id, $db));
}

sub run {
  my ($self) = @_;

  my $seq = $self->query->seq;
  while($seq =~ /([GAC][GC][TA]T[TA][TGCA][TGA][GA][ATG])\w{10,14}([CGT]AC[TAG]GT[GC])/g) {
    my $pos = pos($seq);
    if ($1 eq 'GGTTTTTGA' and $2 eq 'CACTGTG') {
      $self->say_with_header('Strict 5 '.__LINE__." $pos");
    }
    else {
      my $nonamer = $1;
      my $heptamer = $2;
      if ($nonamer =~ /GGTTT[TC]TG[ATG]/ and $heptamer =~ /CAC[TA]TGTG/) {
        $self->say_with_header('Relax 5 '.__LINE__." $pos");
      }
      elsif ($nonamer =~ /GG[TA]T[TA][TC][TG]G[ATG]/ and $heptamer =~ /[CG]AC[TA]TGTG/) {
        $self->say_with_header('Relax 5 '.__LINE__." $pos");
      }
      elsif ($nonamer =~ /GG[TA]T[TA][TCG][TGA]G[ATG]/ and $heptamer =~ /[CGT]AC[TAG]TGT[GC]/) {
        $self->say_with_header('Relax 5 '.__LINE__." $pos");
      }
      else {
        $self->say_with_header('Relax 5 '.__LINE__." $pos $nonamer   $heptamer");
      }
    }
  }
}

sub write_output {
  my ($self) = @_;
}
1;

