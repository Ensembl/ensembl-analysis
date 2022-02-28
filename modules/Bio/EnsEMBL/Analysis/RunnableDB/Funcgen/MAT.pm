# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT->new
     (
         -db       => $db,
         -input_id => 'chromosome::20:1:100000:1',
         -analysis => $analysis,
     );
  $runnable->fetch_input;
  $runnable->run;
  $runnable->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl functional genomics 
database and the Runnable MAT which wraps the program MAT (Model-based 
Analysis of Tiling-array, see http://chip.dfci.harvard.edu/~wli/MAT).

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::MAT;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::MAT;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Funcgen::Importer;

use vars qw(@ISA); 
@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);

=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new MAT runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT object
  Exceptions  : 
  Example     : 

=cut

sub new {

    print "Analysis::RunnableDB::Funcgen::MAT::new\n";
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->read_and_check_config($CONFIG);

    return $self;
	
}

=head2 query

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : integer
  Function  : container for chip number 
  Returntype: integer
  Exceptions: throws if passed the incorrect value
  Example   : 

=cut


sub query {
    my ($self, $chip) = @_;
    if($chip){
        throw("Must pass RunnableDB::Funcgen::MAT::query an integer ".
              "specifying the chip to process not ".$chip) 
            unless($chip =~ m/^\d+$/);
        $self->{'chip'} = $chip;
    }
    return $self->{'chip'};
}

=head2 check_Sets

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : fetch and set ResultSets of interest
  Returntype: 1
  Exceptions: none
  Example   : 

=cut


sub check_Sets
{

    warn("\nNEED TO IMPLEMENT SETTING OF RESULT/DATA/FEATURE SET!!!\n\n");

}

sub fetch_input {

    my ($self) = @_;
    print "Analysis::RunnableDB::Funcgen::MAT::fetch_input\n";

    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'
        .$self->analysis->module;

    warn("No input to fetch since we are going to work directly on the CEL files.");

    $runnable = $runnable->new
        (
         -query => $self->query,
         -program => $self->analysis->program_file,
         -analysis => $self->analysis,
         );
    
    $self->runnable($runnable);

    return 1;

}

sub write_output{
    
    print "RunnableDB::Funcgen::MAT::write_output\n";
    my ($self) = @_;

    my @result_files = map { $_->resultsfile } @{$self->runnable};
    
    print Dumper @result_files;

    my $Importer = Bio::EnsEMBL::Funcgen::Importer->new
        (
         -name        => $ENV{EXPERIMENT},
         -vendor      => $ENV{VENDOR},
         -format      => $ENV{FORMAT},
         -db          => $self->db,
         -host        => $self->db->host,
         -port        => $self->db->port,
         -user        => $self->db->username,
         -pass        => $self->db->password,
         -dbname      => $self->db->dbname,
         -species     => $ENV{SPECIES},
         -data_version => $ENV{DATA_VERSION},
         -input_dir   => $ENV{ANALYSIS_WORK_DIR},
         -output_dir  => $ENV{ANALYSIS_WORK_DIR},
         -result_files => \@result_files,
         -experimental_set_name => $ENV{EXP_SET},
         -exp_date     => $ENV{EXP_DATE},
         -feature_analysis => $self->analysis->logic_name,
         -feature_type_name => $ENV{FEATURE_TYPE},
         -cell_type_name => $ENV{CELL_TYPE},
         -location    => $ENV{LOCATION},
         -contact     => $ENV{CONTACT},
         -group       => $ENV{EGROUP},
         -recover     => 1,
         -verbose     => 1,
#         -ssh         =>  $ssh,
#         -array_set   => $array_set,
#         -array_name  => $array_name,
#         -result_set_name => $rset_name, #not implemented yet
#         -write_mage    => $write_mage,
#         -update_xml => $update_xml,
#         -no_mage => $no_mage,
#         -data_root   => $data_dir,
#         -dump_fasta  => $fasta,
#         -norm_method => $nmethod,
#         -farm => $farm,
#         -input_dir   => $input_dir,
#         -exp_date     => $exp_date,
#         -result_files => \@result_files,
#         -old_dvd_format => $old_dvd_format,
         #Exp does not build input dir, but could
         #This allows input dir to be somewhere 
         #other than efg dir structure
         );

    $Importer->register_experiment();
    
}


1;
