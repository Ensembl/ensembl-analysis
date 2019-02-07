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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT;

use strict;
use warnings;

use Bio::EnsEMBL::IO::Parser::EMBL;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_seq_length => 10,
    column_names => ['iid'],
  }
}


sub fetch_input {
  my $self = shift;

  my $files = $self->param_required('iid');
  my $input_seq_count = 0;
  my $below_min_length_count = 0;
  my $contains_stop = 0;

  my $min_seq_length = $self->param('min_seq_length');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));
  if (ref($files) ne 'ARRAY') {
    $files = [$files];
  }
  my @iids;
  foreach my $file_path (@$files) {
    $self->throw("The input id doesn't exist, offending path:\n$file_path")
      unless(-e $file_path);
    my $parser = Bio::EnsEMBL::IO::Parser::EMBL->open($file_path);

    my @features;
    my %ful_gene_stats;
    while($parser->next()) {
      ++$input_seq_count;
      my @header = split(';', $parser->get_raw_id->[0]);
      next if ($header[4] eq ' PAT' or $header[5] eq ' INV' or $header[5] eq ' SYN' or $header[5] eq ' UNC');
      next if ($parser->get_raw_description->[0] =~ /uPA/);
      my $open_feature = 0;
      my $seq_feature;
      my $seq_name;
      my $open_seq = 0;
      my $type;
      my $type_index = 1;
      foreach my $feature (@{$parser->get_raw_features}) {
        if ($feature =~ /^FT {3}(\S+)/) {
          my $segment = $1;
#          print STDOUT __LINE__, ' ', $segment, "\n";
          if ($open_feature) {
            if ($seq_feature) {
              my $biotype = $type;
              my $seq;
              if ($seq_name) {
#            print STDOUT __LINE__, ' ', $seq_name, "\n";
                if ($seq_name =~ /(IG|TR)\w*[VJDCM]\d?/) {
                  $biotype = $1.'_'.$type;
                }
                $seq_name =~ s/\*/./g;
                $seq_name .= '.'.$type_index++;
              }
              else {
                $biotype = "IG_$type";
                $seq_name = $type_index++;
              }
              my $accession = substr($header[0], 5).'.'.substr($header[1], 4).".$seq_name";
              $seq = $seq_feature;
              if ($type eq 'full_gene') {
                if ($seq =~ /C\w{10,18}W\w{40,48}[IVLFCMA]\w{3,5}C\w{10,14}[FW]G\wG/) {
                  $biotype = 'cwhcj_strict';
                  ++$ful_gene_stats{cwhcj_strict};
                }
                elsif ($seq =~ /C\w{7,20}W\w{40,50}[IVLFCMA]\w{3,8}C\w{8,19}[FW]G\wG/) {
                  $biotype = 'cwhcj_relax';
                  ++$ful_gene_stats{cwhcj_relax};
                }
                elsif ($seq =~ /W\w{40,48}[IVLFCMA]\w{3,5}C\w{10,14}[FW]G\wG/) {
                  $biotype = 'whcj_strict';
                  ++$ful_gene_stats{whcj_strict};
                }
                elsif ($seq =~ /W\w{40,50}[IVLFCMA]\w{3,8}C\w{10,19}[FW]G\wG/) {
                  $biotype = 'whcj_relax';
                  ++$ful_gene_stats{whcj_relax};
                }
                elsif ($seq =~ /[IVLFCMA]\w{3,5}C\w{10,14}[FW]G\wG/) {
                  $biotype = 'hcj_strict';
                  ++$ful_gene_stats{hcj_strict};
                }
                elsif ($seq =~ /[IVLFCMA]\w{3,8}C\w{10,19}[FW]G\wG/) {
                  $biotype = 'hcj_relax';
                  ++$ful_gene_stats{hcj_relax};
                }
                elsif ($seq =~ /C\w{10,14}[FW]G\wG/) {
                  $biotype = 'cj_strict';
                  ++$ful_gene_stats{cj_strict};
                }
                elsif ($seq =~ /C\w{10,19}[FW]G\wG/) {
                  $biotype = 'cj_relax';
                  ++$ful_gene_stats{cj_relax};
                }
                elsif ($seq =~ /[FW]G\wG/) {
                  $biotype = 'cj_strict';
                  ++$ful_gene_stats{cj_strict};
                }
                elsif ($seq =~ /[FW]G\wG/) {
                  $biotype = 'cj_relax';
                  ++$ful_gene_stats{cj_relax};
                }
                elsif ($seq =~ /C\w*W\w*[IVLFCMA]*C*[FW]G\wG/) {
                  $biotype = 'cwhcj';
                  ++$ful_gene_stats{cwhcj};
                }
                elsif ($seq =~ /C\w{10,18}W\w{40,48}[IVLFCMA]\w{3,5}C/) {
                  $biotype = 'cwhc_strict';
                  ++$ful_gene_stats{cwhc_strict};
                }
                elsif ($seq =~ /C\w{10,20}W\w{40,50}[IVLFCMA]\w{3,8}C/) {
                  $biotype = 'cwhc_relax';
                  ++$ful_gene_stats{cwhc_relax};
                }
                elsif ($seq =~ /W\w{40,48}[IVLFCMA]\w{3,5}C/) {
                  $biotype = 'whc_strict';
                  ++$ful_gene_stats{whc_strict};
                }
                elsif ($seq =~ /W\w{40,50}[IVLFCMA]\w{3,8}C/) {
                  $biotype = 'whc_relax';
                  ++$ful_gene_stats{whc_relax};
                }
                elsif ($seq =~ /[IVLFCMA]\w{3,5}C/) {
                  $biotype = 'hc_strict';
                  ++$ful_gene_stats{hc_strict};
                }
                elsif ($seq =~ /[IVLFCMA]\w{3,8}C/) {
                  $biotype = 'hc_relax';
                  ++$ful_gene_stats{hc_relax};
                }
                elsif ($seq =~ /C\w{55,63}C/) {
                  $biotype = 'g_strict';
                  ++$ful_gene_stats{g_strict};
                }
                elsif ($seq =~ /C\w{55,66}C/) {
                  $biotype = 'g_relax';
                  ++$ful_gene_stats{g_relax};
                }
                else {
                  $biotype = 'other';
                  ++$ful_gene_stats{other};
                  $self->warning("OTHER $accession")
                }
              }
              $self->say_with_header("$accession $biotype $seq");
              push(@features, {accession => $accession, biotype => $biotype, seq => $seq, source_db => 'imgt', pe_level => 2, group_name => 'imgt'});
              undef $type;
              undef $seq_feature;
              undef $seq_name;
            }
            else {
              $self->warning('Something went wrong for '.$header[0]." $seq_name $feature");
            }
          }
          if ($segment =~ /^([CDJVM])-REGION/) {
            ++$open_feature;
            $type = $1.'_gene';
          }
          elsif ($segment eq 'CDS') {
            ++$open_feature;
            $type = 'full_gene';
          }
          else {
            $open_feature = 0;
          }
        }
        elsif ($open_feature) {
          if ($feature =~ '/translation="(\w+)("?)') {
            $seq_feature = $1;
            ++$open_seq unless ($2);
          }
          elsif ($feature =~ '/(IMGT_)?allele="([^" ]+)') {
            $seq_name = $2;
          }
          elsif ($feature =~ '/(IMGT_)?gene="([^" ]+)') {
            $seq_name = $2 unless ($seq_name);
          }
          elsif ($feature =~ '/protein_id="([^"]+)"') {
            if ($type eq 'full_gene' or !$seq_name) {
              $seq_name = $1;
            }
          }
          elsif ($open_seq) {
            $feature =~ /^FT {4,}(\w*)("?)/;
            $seq_feature .= $1;
            $open_seq = 0 if ($2);
          }
        }
      }
      if ($type) {
        if ($seq_feature) {
          my $accession;
          my $biotype = $type;
          my $seq;
          if ($seq_name) {
            if ($seq_name =~ /(IG|TR)\w*[VJDCM]\d?/) {
              $biotype = $1.'_'.$type;
            }
            $seq_name =~ s/\*/./g;
            $seq_name .= '.'.$type_index++;
          }
          else {
            $biotype = "IG_$type";
            $seq_name = $type_index++;
          }
          if ($seq and $type eq 'full_gene') {
            if ($seq =~ /C\w{10,18}W\w{40,48}[IVLFCMA]\w{3,5}C\w{10,14}[FW]G\wG/) {
              $biotype = 'cwhcj_strict';
              ++$ful_gene_stats{cwhcj_strict};
            }
            elsif ($seq =~ /C\w{7,20}W\w{40,50}[IVLFCMA]\w{3,8}C\w{8,19}[FW]G\wG/) {
              $biotype = 'cwhcj_relax';
              ++$ful_gene_stats{cwhcj_relax};
            }
            elsif ($seq =~ /W\w{40,48}[IVLFCMA]\w{3,5}C\w{10,14}[FW]G\wG/) {
              $biotype = 'whcj_strict';
              ++$ful_gene_stats{whcj_strict};
            }
            elsif ($seq =~ /W\w{40,50}[IVLFCMA]\w{3,8}C\w{10,19}[FW]G\wG/) {
              $biotype = 'whcj_relax';
              ++$ful_gene_stats{whcj_relax};
            }
            elsif ($seq =~ /[IVLFCMA]\w{3,5}C\w{10,14}[FW]G\wG/) {
              $biotype = 'hcj_strict';
              ++$ful_gene_stats{hcj_strict};
            }
            elsif ($seq =~ /[IVLFCMA]\w{3,8}C\w{10,19}[FW]G\wG/) {
              $biotype = 'hcj_relax';
              ++$ful_gene_stats{hcj_relax};
            }
            elsif ($seq =~ /C\w{10,14}[FW]G\wG/) {
              $biotype = 'cj_strict';
              ++$ful_gene_stats{cj_strict};
            }
            elsif ($seq =~ /C\w{10,19}[FW]G\wG/) {
              $biotype = 'cj_relax';
              ++$ful_gene_stats{cj_relax};
            }
            elsif ($seq =~ /[FW]G\wG/) {
              $biotype = 'cj_strict';
              ++$ful_gene_stats{cj_strict};
            }
            elsif ($seq =~ /[FW]G\wG/) {
              $biotype = 'cj_relax';
              ++$ful_gene_stats{cj_relax};
            }
            elsif ($seq =~ /C\w*W\w*[IVLFCMA]*C*[FW]G\wG/) {
              $biotype = 'cwhcj';
              ++$ful_gene_stats{cwhcj};
            }
            elsif ($seq =~ /C\w{10,18}W\w{40,48}[IVLFCMA]\w{3,5}C/) {
              $biotype = 'cwhc_strict';
              ++$ful_gene_stats{cwhc_strict};
            }
            elsif ($seq =~ /C\w{10,20}W\w{40,50}[IVLFCMA]\w{3,8}C/) {
              $biotype = 'cwhc_relax';
              ++$ful_gene_stats{cwhc_relax};
            }
            elsif ($seq =~ /W\w{40,48}[IVLFCMA]\w{3,5}C/) {
              $biotype = 'whc_strict';
              ++$ful_gene_stats{whc_strict};
            }
            elsif ($seq =~ /W\w{40,50}[IVLFCMA]\w{3,8}C/) {
              $biotype = 'whc_relax';
              ++$ful_gene_stats{whc_relax};
            }
            elsif ($seq =~ /[IVLFCMA]\w{3,5}C/) {
              $biotype = 'hc_strict';
              ++$ful_gene_stats{hc_strict};
            }
            elsif ($seq =~ /[IVLFCMA]\w{3,8}C/) {
              $biotype = 'hc_relax';
              ++$ful_gene_stats{hc_relax};
            }
            elsif ($seq =~ /C\w{55,63}C/) {
              $biotype = 'g_strict';
              ++$ful_gene_stats{g_strict};
            }
            elsif ($seq =~ /C\w{55,66}C/) {
              $biotype = 'g_relax';
              ++$ful_gene_stats{g_relax};
            }
            else {
              $biotype = 'other';
              ++$ful_gene_stats{other};
              $self->warning("OTHER $accession")
            }
          }
          $accession = substr($header[0], 5).'.'.substr($header[1], 4).".$seq_name";
          $accession =~ tr/ /_/;
          $seq = $seq_feature;
          $self->say_with_header("$accession $biotype $seq");
          push(@features, {accession => $accession, biotype => $biotype, seq => $seq, source_db => 'imgt', pe_level => 2, group_name => 'imgt'});
        }
        else {
          $self->warning('Something went wrong for '.$header[0]." $seq_name");
        }
      }
    }
    foreach my $data (@features) {
      $self->say_with_header($data->{accession}.' '.$data->{biotype}.' '.$data->{seq});
      push(@iids, $data->{accession});
    }
    $table_adaptor->store(\@features);
    foreach my $key (keys %ful_gene_stats) {
      $self->warning("$key\t".$ful_gene_stats{$key});
    }
  }

  $self->warning("Total sequences: $input_seq_count\nSequences below $min_seq_length: $below_min_length_count\nSequences with stops: $contains_stop\nStored sequences: ".scalar(@iids)."\n");
  $self->param('inputlist', \@iids);
}

1;
