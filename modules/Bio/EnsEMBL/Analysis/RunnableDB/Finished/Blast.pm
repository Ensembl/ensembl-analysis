# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Finished::Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $blast   = Bio::EnsEMBL::Pipeline::RunnableDB::Blast->new ( 
                                                    -db         => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
$blast->fetch_input();
$blast->run();
$blast->output();
$blast->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::Blast to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::Blast;

use strict;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::Config::General;
use base 'Bio::EnsEMBL::Analysis::RunnableDB::Blast';


sub fetch_input{

  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db, 
                                    $ANALYSIS_REPEAT_MASKING);
  $self->query($slice);
  my %blast = %{$self->BLAST_PARAMS};
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::Blast->new
    (
     -query => $self->query,
     -program => $self->analysis->program,
     -parser => $parser,
     -filter => $filter,
     -database => $self->analysis->db_file,
     -analysis => $self->analysis,
     %blast,
    );
  my $s=$self->runnable($runnable);
  return 1;

}


sub _createfiles {

    my ($self, $dirname, $filenames) = @_;
    my $unique = {};
    $unique    = { map { $_, $unique->{$_}++ } @$filenames };
    my @files  = ();
    $dirname ||= '/tmp';
    $dirname   =~ s!(\S+)/$!$1!;
    foreach my $file(@$filenames){
        if($unique->{$file}){
            #name not unique add random
            $file .= ".$$.".int(rand(200));
            push(@files, "$dirname/$file");
        }else{
            #name was unique just add it
            push(@files, "$dirname/$file.$$");
        }
    }

    return @files;
}

=head2 run

    Title   :   write_output
    Usage   :   $self->write_output();
    Function:   Writes Features , hit_descriptions to database by calling parent write_output method
    Returns :   none
    Args    :   none

=cut

sub write_output{

    my ($self) = @_;
    $self->SUPER::write_output();
    foreach my $runnable (@{$self->runnable}) {
      my $db_version = $runnable->get_db_version if $runnable->can('get_db_version');
      $self->db_version_searched($db_version); # make sure we set this here
      if ( my @output = @{$runnable->output} ) {
	my $dbobj      = $self->db;
        my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
        my %ids        = map { $_->hseqname, 1 } @output;
        $seqfetcher->write_descriptions( $dbobj, keys(%ids) );
      }

   }
   return 1;

}



=head2 db_version_searched

    Title   :  db_version_searched
               [ distinguished from Runnable::*::get_db_version() ]
    Useage  :  $self->db_version_searched('version string')
               $obj->db_version_searched()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished::Blast
               This is just a holding place for the string in this
               module
    Returns :  String or undef
    Args    :  String
    Caller  :  $self::run()
               Job::run_module()

=cut

sub db_version_searched{

    my ($self, $arg) = @_;
    $self->{'_db_version_searched'} = $arg if $arg;
    return $self->{'_db_version_searched'};

}

1;
