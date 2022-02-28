=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Infernal - 

=head1 SYNOPSIS

    my $runnableDB = Bio::EnsEMBL::Analysis::RunnableDB::Infernal->new(
            	-db       => $db_adaptor,
		-input_id => $flag_id,		
		-analysis => $analysis,
	        );
    $runnabledb->fetch_input();
    $runnabledb->run();
    $runnabledb->write_output();


=head1 DESCRIPTION

RunnableDB to wrap cmsearch - part of the Infernal suite of programs by Sean Eddy.
Uses RFAM blast hits to identify regions where cmsearch should be run. Blast hits
are run using the RfamBlast module and then are filtered on a familly by familly basis.
Blast hits that look promising are flagged by predict_ncRNAs.pl and the flags are used
as input_ids.
Uses the Bio::EnsEMBL::Analysis::Config::Databases config file. Writes non coding genes to 
the $GB_FINALDB database and gets DNA from $GB_ database. The blast hits from BlastRfam are
stored in the pipeline database;
Creates single exon non-coding genes with gene descriptions obtained from Rfam.descriptions
file. The dna align feature representing the initial blast hit is added as a supporting feature
and the RNA secondary structure predicted by Infernal is added as a transcript attribute
in a length encoded string form to take up less room in the db.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterncRNAs;

use strict;
use warnings;
use feature 'say';

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

my $runnable;

=head2 fetch_input

  Title   :   fetch_input
  Usage   :   $self->fetch_input
  Function:   Opens connections to the final gene build database to store the genes in.
          :   Also opens connection to the genebuild DNA database to fetch DNA from.
          :   Parses flags use as input_ids and fetches dna_align_features to run cmsearch
          :   on.
  Returns :   none
  Args    :   none

=cut

sub fetch_input{
  my ($self)=@_;

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                           );
  $self->analysis($analysis);


  # The output db should be the one that the dafs to check have been written to
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));


  if($dna_dba) {
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $daf_adaptor = $output_dba->get_DnaAlignFeatureAdaptor;


  my $mirna_filtered_ids = $self->filter_mirna_hits($output_dba,$daf_adaptor);


  say "Filtered miRNA count: ".scalar(@{$mirna_filtered_ids});

  $self->final_mirna_db_ids($mirna_filtered_ids);

  return 1;
}

=head2 fetch_input

  Title   :   fetch_input_deprecated // removed 12/2018
  Usage   :   $self->fetch_input
  Function:   Opens connections to the final gene build database to store the genes in.
          :   Also opens connection to the genebuild DNA database to fetch DNA from.
          :   Parses flags use as input_ids and fetches dna_align_features to run cmsearch
          :   on.
  Returns :   none
  Args    :   none

=cut

sub fetch_input_deprecated{
  my ($self)=@_;

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                           );
  $self->analysis($analysis);


  # The output db should be the one that the dafs to check have been written to
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));


  if($dna_dba) {
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $daf_adaptor = $output_dba->get_DnaAlignFeatureAdaptor;


  my $rfam_filtered_ids = $self->filter_rfam_hits($output_dba,$daf_adaptor);
  my $mirna_filtered_ids = $self->filter_mirna_hits($output_dba,$daf_adaptor);


  say "Filtered miRNA count: ".scalar(@{$mirna_filtered_ids});
  say "Filtered Rfam count: ".scalar(@{$rfam_filtered_ids});

  $self->final_mirna_db_ids($mirna_filtered_ids);
  $self->final_rfam_db_ids($rfam_filtered_ids);

  return 1;
}


=head2 run

  Args       : none
  Description: Runs the runnable.
  Exceptions : Throws if the the runnable module is not set
  Returntype : scalar

=cut

sub run{
  my ($self) = @_;

  return 1;
}

=head2 write_output

  Args       : none
  Description: Writes the single exon ncRNA genes into the final genebuild databse, 
             : also stores the structure attribute associated with the transcript
  Exceptions : Throws if gene or transcript attribute fail to write to the database
  Returntype : scalar

=cut
sub write_output{
  my ($self) = @_;

  my $mirna_db_ids = $self->final_mirna_db_ids();

  my $mirna_batch = [];
  my $mirna_batch_size = 1000;
  my $batch_count = 0;
  my $output_hash = {};
  my $mirna_branch_code = 2;

  # Should put these into a single sub that takes in the ids, the branch code and the batch size at some point
  foreach my $output_id (@{$mirna_db_ids}) {
    if($batch_count == $mirna_batch_size) {
      $output_hash = {};
      $output_hash->{'iid'} = $mirna_batch;
      $self->dataflow_output_id($output_hash,$mirna_branch_code);
      $mirna_batch = [];
      push(@{$mirna_batch},$output_id);
      $batch_count = 1;
    } else {
      push(@{$mirna_batch},$output_id);
      $batch_count++;
    }
  }

  # Send out last batch
  if(scalar(@{$mirna_batch})) {
    $output_hash = {};
    $output_hash->{'iid'} = $mirna_batch;
    $self->dataflow_output_id($output_hash,$mirna_branch_code);
  }


  return 1;
}

sub write_output_deprecated{
  my ($self) = @_;

  my $mirna_db_ids = $self->final_mirna_db_ids();
  my $rfam_db_ids = $self->final_rfam_db_ids();

  my $mirna_batch = [];
  my $mirna_batch_size = 1000;
  my $batch_count = 0;
  my $output_hash = {};
  my $mirna_branch_code = 2;
  my $rfam_branch_code = 3;

  # Should put these into a single sub that takes in the ids, the branch code and the batch size at some point
  foreach my $output_id (@{$mirna_db_ids}) {
    if($batch_count == $mirna_batch_size) {
      $output_hash = {};
      $output_hash->{'iid'} = $mirna_batch;
      $self->dataflow_output_id($output_hash,$mirna_branch_code);
      $mirna_batch = [];
      push(@{$mirna_batch},$output_id);
      $batch_count = 1;
    } else {
      push(@{$mirna_batch},$output_id);
      $batch_count++;
    }
  }

  # Send out last batch
  if(scalar(@{$mirna_batch})) {
    $output_hash = {};
    $output_hash->{'iid'} = $mirna_batch;
    $self->dataflow_output_id($output_hash,$mirna_branch_code);
  }

  my $rfam_batch = [];
  my $rfam_batch_size = 25;
  $batch_count = 0;
  $output_hash = {};
  foreach my $output_id (@{$rfam_db_ids}) {
    if($batch_count == $rfam_batch_size) {
      $output_hash = {};
      $output_hash->{'iid'} = $rfam_batch;
      $self->dataflow_output_id($output_hash,$rfam_branch_code);
      $rfam_batch = [];
      push(@{$rfam_batch},$output_id);
      $batch_count = 1;
    } else {
     push(@{$rfam_batch},$output_id);
     $batch_count++;
   }
  }

  # Send out last batch
  if(scalar(@{$rfam_batch})) {
    $output_hash = {};
    $output_hash->{'iid'} = $rfam_batch;
    $self->dataflow_output_id($output_hash,$rfam_branch_code);
  }

  return 1;
}


sub filter_rfam_hits {
  my ($self,$output_dba,$daf_adaptor) = @_;

  my %rfam_blasts;
  my %rfam_threshold;
  my @rfam_filtered_ids = ();

  my $rfam_logic_name = 'rfamblast';
  my $rfam_analysis_id = ${sql("select analysis_id from analysis where logic_name='".$rfam_logic_name."';",$output_dba)}[0]->[0];

  my %thr;
  my @domcount =  @{sql("select count(*) , LEFT(hit_name,7) from dna_align_feature where  analysis_id = ".$rfam_analysis_id.
                        " group by LEFT(hit_name,7) ;",$output_dba)};

  foreach my $dom (@domcount){
    $thr{$dom->[1]} = $dom->[0];
  }

  my $total = scalar(keys (%thr));
  my $complete;
  my $last = 0;
  my $count = 0;
  print "\n0_________10________20________30________40________50________60________70________80________90_______100\n";
  foreach my $domain (sort keys %thr){
    $count++;
    my @dafs;
    if($thr{$domain} > 2000) {
      my @top_scores = @{sql("SELECT evalue FROM dna_align_feature WHERE analysis_id = ".$rfam_analysis_id.
                             " AND LEFT(hit_name,7) = '".$domain."' ORDER BY  evalue ASC limit 2000;",$output_dba)};
      my $cutoff = pop (@top_scores)->[0];
      @dafs = @{$daf_adaptor->generic_fetch("left(hit_name,7) = \"".$domain."\" AND evalue <= ".$cutoff)};
    } else {
      @dafs = @{$daf_adaptor->generic_fetch("left(hit_name,7) = \"".$domain."\"")};
    }
    $complete = int($count/$total*100);
    if ($complete > $last){
      my $num = $complete -  $last;
      foreach (my $i = 0; $i < $num; $i++){
        print "=";
      }
    }
    $last = $complete;
    @dafs = sort {$a->p_value <=> $b->p_value} @dafs if (scalar(@dafs) > 2000 );
    foreach my $daf(@dafs) {
      if ($rfam_blasts{$domain} && scalar(@{$rfam_blasts{$domain}}) >= 2000 ) {
        last;
      }

      if ($daf->p_value > 0.01) {
        next;
      }
      push @{$rfam_blasts{$domain}},$daf;
    }

  }

  foreach my $domain (keys %rfam_blasts){
    my @hits = @{$rfam_blasts{$domain}};
    foreach my $hit (@hits){
      push @rfam_filtered_ids,$hit->dbID;
    }
  }

  return(\@rfam_filtered_ids);
}


sub filter_mirna_hits {
  my ($self,$output_dba,$daf_adaptor) = @_;

  my @mirna_filtered_ids = ();
  my $mirna_logic_name = 'blastmirna';
  my $mirna_analysis_id = ${sql("select analysis_id from analysis where logic_name='".$mirna_logic_name."';",$output_dba)}[0]->[0];
  my %mirbase_blasts;

  print "Filtering miRNAs";
  print "\n0_________10________20________30________40________50________60________70________80________90_______100\n";

  my %mi_types;
  my @domcount =  @{sql("select count(*) , hit_name from dna_align_feature where  analysis_id = ".$mirna_analysis_id.
                        " group by hit_name;",$output_dba)};
  foreach my $dom (@domcount){
    $mi_types{$dom->[1]} = $dom->[0];
  }
  my $count = 0;
  my $total = scalar(keys %mi_types);
  my $last = 0;
  foreach my $type (keys %mi_types){
    $count++;
    my @dafs;
    if($mi_types{$type} > 50) {
      my @top_scores = @{sql("SELECT evalue FROM dna_align_feature WHERE analysis_id = ".$mirna_analysis_id.
                             " AND hit_name = '".$type."' ORDER BY  evalue ASC limit 50;",$output_dba)};
        my $cutoff = pop (@top_scores)->[0];
        @dafs = @{$daf_adaptor->generic_fetch("hit_name = \"".$type."\" AND evalue <= ".$cutoff)};
    } else {
      @dafs = @{$daf_adaptor->generic_fetch("hit_name = \"".$type."\" ")};
    }

    my $complete = int($count/$total*100);
    if($complete > $last) {
      my $num = $complete -  $last;
      foreach (my $i = 0; $i < $num; $i++){
        print "=";
      }
    }
    $last = $complete;
    @dafs = sort {$a->p_value <=> $b->p_value} @dafs if (scalar(@dafs) > 50 );

    foreach my $daf(@dafs){
      if ($mirbase_blasts{$type} &&  scalar(@{$mirbase_blasts{$type}}) >= 50 ) {
        last;
      }
      push @{$mirbase_blasts{$type}},$daf;
    }
  }

  print "\nGenerating Flags\n";

  foreach my $domain(keys %mirbase_blasts){
    my @hits = @{$mirbase_blasts{$domain}};
    foreach my $hit (@hits){
      push @mirna_filtered_ids,$hit->dbID;
    }
  }

  return(\@mirna_filtered_ids);
}

sub final_mirna_db_ids {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_final_mirna_db_ids',$val);
  }

  return($self->param('_final_mirna_db_ids'));
}


sub final_rfam_db_ids {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_final_rfam_db_ids',$val);
  }

  return($self->param('_final_rfam_db_ids'));
}

sub sql {
  my ($query,$db) = @_;
  my @array;

  print "FM2 QUERY: ".$query;

  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  while ( my @row = $sth->fetchrow_array ) {
    push @array , \@row;
  }
  return \@array;
}


# I wrote this as a conventional way of getting the ids using the API but it seems to be too slow
# when there are lots of hits
sub fetch_input_slow {
  my ($self)=@_;

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                           );
  $self->analysis($analysis);

  # The output db should be the one that the dafs to check have been written to
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));


  if($dna_dba) {
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $daf_adaptor = $output_dba->get_DnaAlignFeatureAdaptor;
  my $slice_adaptor = $output_dba->get_SliceAdaptor;

  my $mirna_dafs = $daf_adaptor->fetch_all_by_logic_name('mirna_blast');
  say "Unfiltered miRNA hit count: ".scalar(@{$mirna_dafs});

  my $mirna_types;
  my $mirna_filtered_ids = [];
  my $rfam_types;
  my $rfam_filtered_ids = [];

  my $slices = $slice_adaptor->fetch_all('toplevel');

  foreach my $slice (@{$slices}) {
    my $mirna_dafs = $daf_adaptor->fetch_all_by_Slice($slice,'mirna_blast');
    say "Working on slice: ".$slice->name;
    # First group all the db ids for each hit name in the hash and store the evalue
    foreach my $daf (@{$mirna_dafs}) {
      my $evalue = $daf->p_value;
      my $hit_name = $daf->hseqname;
      my $db_id =$daf->dbID;
      $mirna_types->{$hit_name}->{$db_id} = $evalue;
    }
    undef($mirna_dafs);

    my $rfam_dafs = $daf_adaptor->fetch_all_by_Slice($slice,'rfam_blast');
    foreach my $daf (@{$rfam_dafs}) {
      my $evalue = $daf->p_value;
      my $hit_name = $daf->hseqname;
      my $db_id =$daf->dbID;
      $rfam_types->{$hit_name}->{$db_id} = $evalue;
    }
    undef($rfam_dafs);
  }

  # Now loop through each hit name, sort the underlying db ids based on evalue, loop through them until
  # the evalue it higher than the evalue of the first key
  foreach my $mirna_group (keys(%$mirna_types)) {
    my $db_id_hash = $mirna_types->{$mirna_group};
    my @sorted_db_ids = (sort {$db_id_hash->{$a} <=> $db_id_hash->{$b} } keys %$db_id_hash);
    my $best_evalue = $db_id_hash->{$sorted_db_ids[0]};
    foreach my $db_id (@sorted_db_ids) {
      my $current_evalue = $db_id_hash->{$db_id};
      if($current_evalue > $best_evalue) {
        last;
      } else {
        push(@{$mirna_filtered_ids},$db_id);
      }
    }
  }


  # Now loop through each hit name, sort the underlying db ids based on evalue, loop through them until
  # the evalue it higher than the evalue of the first key
  foreach my $rfam_group (keys(%$rfam_types)) {
    my $hit_count = 0;
    my $db_id_hash = $rfam_types->{$rfam_group};
    my @sorted_db_ids = (sort {$db_id_hash->{$a} <=> $db_id_hash->{$b} } keys %$db_id_hash);
    my $best_evalue = $db_id_hash->{$sorted_db_ids[0]};
    foreach my $db_id (@sorted_db_ids) {
      my $current_evalue = $db_id_hash->{$db_id};
      if($current_evalue > $best_evalue || $hit_count > 2000 || $current_evalue > 0.01) {
        last;
      } else {
        push(@{$rfam_filtered_ids},$db_id);
        $hit_count++;
      }
    }
  }

  say "Filtered miRNA count: ".scalar(@{$mirna_filtered_ids});
  say "Filtered Rfam count: ".scalar(@{$rfam_filtered_ids});

  $self->final_mirna_db_ids($mirna_filtered_ids);
  $self->final_rfam_db_ids($rfam_filtered_ids);

  return 1;
}

1;
