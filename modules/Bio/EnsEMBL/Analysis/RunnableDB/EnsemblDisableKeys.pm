=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::EnsemblDisableKeys - 

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::EnsemblDisableKeys->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->run();  #writes to DB

=head1 DESCRIPTION

This module has 2 functions:
 - disable/enable multiples keys in a database by giving an input id like :
  DBKEY:table_name:enable
 - dump, sort and reload a table from a database: 
  DBKEY:table_name:sort_and_dump:NEW_DB_KEY
   if the input id does not have a value for the new database, it backup the table
   and stores the sorted table in the "normal" table

  
=head1 METHODS


=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::EnsemblDisableKeys;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
#use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::General qw(BIN_DIR ANALYSIS_WORK_DIR);

use vars qw(@ISA);

#This class inhertis from BaseGeneBuild and BGB inhertis from RunnableDB.pm  
@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild); 

## Constants
# The suffix of the backup table if a new database has not been added in the input_id
my $suffix_sorted_table = '_nonsorted';
# The suffix of the sorted dumped file
my $suffix_sorted = '.sorted';


sub new {
  my ( $class, @args ) = @_; 
 
  my $self = $class->SUPER::new(@args); 

  throw("Your input id is not DBKEY:table_name:action : ".$self->input_id()."!") unless ($self->input_id() =~ /^\S+:\S+:\S+/);
  # split input id 
  ($self->{DB_HKEY}, $self->{TABLE_NAME}, $self->{MODIFIER}, $self->{NEW_DB_KEY}) = split(':', $self->input_id()) ; 

  $self->db->disconnect_when_inactive(1); 
  return $self;
}


sub fetch_input{  
  my ($self) = @_;   
  return ; 
} 

sub run {
  my ($self) = @_; 
  throw ("Can't run - no runnable objects") unless ( $self->runnable );
  print STDERR @{$self->runnable}."\n";

  my $ref_db = $self->get_dbadaptor($self->DB_HKEY);   

  if ( !defined $ref_db ) { 
     throw ("db " . $self->DB_HKEY . " does not exist !!!\n");  
  }  
  if ($self->MODIFIER eq 'disable' or $self->MODIFIER eq 'enable') {
      $self->toggle_indexation_multkeys($ref_db);
  }
  elsif ($self->MODIFIER eq 'sort_and_dump') { 
      my $temp_dir = $ANALYSIS_WORK_DIR;
      $self->sort_and_dump_table($ref_db, $temp_dir);
  }
  else {
      throw("We do not know this input id: ".$self->input_id()."\n");
  }
}


=head2 DB_HKEY

    Getter for the database hash name

=cut

sub DB_HKEY {
    my $self = shift;
    return $self->{DB_HKEY};
}


=head2 TABLE_NAME

    Getter for the name of the table

=cut

sub TABLE_NAME {
    my $self = shift;
    return $self->{TABLE_NAME};
}


=head2 MODIFIER

    Getter for the action to do:
     - enable
     - disable
     - sort_and_dump

=cut

sub MODIFIER {
    my $self = shift;
    return $self->{MODIFIER};
}


=head2 NEW_DB_KEY

    Getter for the database hash name of the new database when 
    we sort and dump a table

=cut

sub NEW_DB_KEY {
    my $self = shift;
    return $self->{NEW_DB_KEY};
}


# FUNCTIONS:

=head2 toggle_indexation_multkeys

    Enable/disable the indexation of multiple keys for a specific table
    The name of the table and the name of the DB are in the input id

    $self->toggle_indexation_multkeys($ref_db)

=cut

sub toggle_indexation_multkeys {
    my ($self, $ref_db) = @_;

    my $sth = $ref_db->prepare('alter table '.$self->TABLE_NAME.' '.$self->MODIFIER.' keys') ;
    $sth->execute();

    # do check 
    $sth = $ref_db->prepare('show index from '.$self->TABLE_NAME);
    $sth->execute();
    foreach my $row (@{$sth->fetchall_arrayref}) {
      next unless ($row->[1] == 1);
      if (($self->MODIFIER eq 'disable' and $row->[11] ne 'disabled')
       or ($self->MODIFIER eq 'enable' and $row->[11] ne '')) {
          throw ("Keys were not ".$self->MODIFIER."d for ".$self->TABLE_NAME."\n");
      }
    }
    $sth->finish();
}


=head2 sort_and_dump_table

    Dump the table, create a new database based on a refenrence or backup the table
    and load the sorted table in the new DB or the table

    $self->sort_and_dump_table(ref_db, temp_dir)

=cut

sub sort_and_dump_table {
    my ($self, $ref_db, $temp_dir) = @_;

    my $file_name = $self->dump_table($temp_dir);
    $self->sort_table($file_name);
    if (defined $self->NEW_DB_KEY) {
        $self->create_new_database($self->main_reference_db, $self->NEW_DB_KEY );
        $self->load_sorted_table($self->get_dbadaptor($self->NEW_DB_KEY), $file_name.$suffix_sorted);
        $self->check_new_table($ref_db, $self->get_dbadaptor($self->NEW_DB_KEY));
    }
    else {
        $self->backup_table($ref_db); 
        $self->load_sorted_table($ref_db, $file_name.$suffix_sorted);
        $self->check_new_table($ref_db);
    }
}


=head2 dump_table

    Dump the dna align feature table of a database to disk 
    Returns the path to the file

    $filename = $self->dump_table($temp_dir)

=cut

sub dump_table { 
    my ($self, $temp_dir) = @_;  

    # get db parameters

    my $dbname = $self->get_dbadaptor($self->DB_HKEY)->dbname();  
    my $host = $self->get_dbadaptor($self->DB_HKEY)->host();
    my $user = $self->get_dbadaptor($self->DB_HKEY)->username();
    my $pass= $self->get_dbadaptor($self->DB_HKEY)->password();
    my $port = $self->get_dbadaptor($self->DB_HKEY)->port();


    my $cmd = -x $BIN_DIR."/mysql" ? $BIN_DIR.'/mysql' : 'mysql'; 
    $cmd .= " -NB -quick -u$user -h$host -P$port " ;  

    if ( $pass ne '' ) { 
        $cmd.=" -p$pass ";
    } 
    my $dump_file_name = $temp_dir."/".$self->DB_HKEY."_".$dbname."_".$self->TABLE_NAME."_".$host."_".$port."_".$$.".dump"; 

    print "Dumping data from ".$self->DB_HKEY." : $dbname\@$host\n" ; 

    $cmd .= " -D$dbname  -e\"select * from ".$self->TABLE_NAME;
    # Works only with dna_align_feature table
    $cmd .= "\" | cut -f2- > $dump_file_name " ;   
    info("cmd : $cmd\n");
    system($cmd) == 0 or throw("ERROR - can't continue. command failed :\" $cmd \"\n");
    print $self->TABLE_NAME." dumped to file :  $dump_file_name\n" ; 
    return $dump_file_name; 
}


=head2 sort_table

    Sort the file which contains the table previously dumped

    $self->sort_table($file_name)

=cut

sub sort_table {   
    my ($self, $file) = @_;  

    my $cmd ="sort -u $file | sort -n -k1 -n -k2 -n -k3 |sed 's/^/\\N\t/' > ".$file.$suffix_sorted;
    info($cmd . "\n");   
    system($cmd) == 0 or throw("ERROR - can't continue. command failed :\" $cmd \"\n");
}


=head2 backup_table

    Backup the table and rename the table before we load the sorted data

    $self->backup_table($ref_db)

=cut

sub backup_table {
    my ($self, $ref_db) = @_;

    #check if backup table already exists
    my $sth = $ref_db->prepare('show tables like "'.$self->TABLE_NAME.$suffix_sorted_table.'"');
    $sth->execute();
    throw("The table ".$self->TABLE_NAME.$suffix_sorted_table." already exists!") if (scalar(@{$sth->fetchall_arrayref}));

    $sth = $ref_db->prepare('create table '.$self->TABLE_NAME.$suffix_sorted_table.' select * from '.$self->TABLE_NAME);
    $sth->execute();

    $sth->finish();
}


=head2 load_sorted_table

    Load the sorted table back in the DB

    $self->load_sorted_table($ref_db, $filename)

=cut

sub load_sorted_table {
    my ($self, $ref_db, $file) = @_;

    my $sth = $ref_db->prepare("truncate ".$self->TABLE_NAME);
    $sth->execute();

    $sth = $ref_db->prepare("load data local infile '$file' into table ".$self->TABLE_NAME);
    $sth->execute();

    unlink $file;
    $file =~ s/$suffix_sorted//;
    unlink $file;
}


=head2 check_new_table

    Check if the new table has the same number of rows as the original

    $self->check_new_table($old_db, $new_db)

=cut

sub check_new_table {
    my ($self, $old_db, $new_db) = @_;

    my $sth = $old_db->prepare('select count(*) from '.$self->TABLE_NAME);
    $sth->execute();
    my ($rows_count) = @{$sth->fetchall_arrayref->[0]};
    my $sorted_count;
    if (defined $new_db) {
        $sth = $new_db->prepare('select count(*) from '.$self->TABLE_NAME);
        $sth->execute();
        ($sorted_count) = @{$sth->fetchall_arrayref->[0]};
    }
    else {
        $sth = $old_db->prepare('select count(*) from '.$self->TABLE_NAME.$suffix_sorted_table);
        $sth->execute();
        ($sorted_count) = @{$sth->fetchall_arrayref->[0]};
    }
    $sth->finish();
    warning("The new table ".$self->TABLE_NAME."has not the same number of rows as its original! $rows_count => $sorted_count") if ($sorted_count != $rows_count);
}


1; 

