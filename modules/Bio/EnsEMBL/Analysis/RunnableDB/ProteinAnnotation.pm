
=head1 NAME

ProteinAnnotation.pm - DESCRIPTION of Object

=head1 SYNOPSIS

this is the base class for the ProteinAnnotation runnabledbs

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;
use Bio::EnsEMBL::Analysis::Config::ProteinAnnotation;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


################################
sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  throw("Analysis object required") unless ($self->analysis);

  $self->read_and_check_config;

  return $self;
}


################################
sub fetch_input {
  my ($self) = @_;  
  my $input_id;

  my $db;
  if ($self->GENEDB) {
    $db = $self->get_dbadaptor($self->GENEDB);
  } else {
    $db = $self->db;
  }
  if ( $self->DNA_DB ) {  
     my $m_dna_db = $self->get_dbadaptor($self->DNA_DB); 
     $m_dna_db->disconnect_when_inactive(1); 
     $db->dnadb($m_dna_db);
     $db->dnadb($m_dna_db);
  } 
  my $error = "For logic " . $self->analysis->logic_name . ", Input id must be:\n";
  $error .= "(a) a valid translation dbID; or\n";
  $error .= "(b) the fully qualified directory/name of a fasta peptide file; or\n";
  $error .= "(c) the name only of a peptide fasta file (with BASE_DIR defined)\n";
  $error .= $self->input_id . " is none of these\n";

  if ($self->input_id =~ /^(\d+)$/)  {
    my $prot;
    eval {
      $prot = $db->get_TranslationAdaptor->fetch_by_dbID($self->input_id);
    };
    if($@ or not defined $prot) {
      throw($error);
    }
    
    $input_id  =  Bio::PrimarySeq->new(-seq         => $prot->seq,
				       -id          => $self->input_id,
				       -accession   => $self->input_id,
				       -moltype     => 'protein');
  } elsif ($self->input_id =~ /^\// and -e $self->input_id) {
    # assume fully-qualified file name
    $input_id = $self->input_id;
  } elsif (defined $self->BASE_DIR and 
           $self->BASE_DIR and
           -e $self->BASE_DIR . "/" . $self->input_id) {    
    $input_id = $self->BASE_DIR . "/" . $self->input_id;
  } else {
    throw($error);
  }
  
  $self->query($input_id);
}


##################################
sub write_output {
  my($self) = @_;
  
  my @features = @{$self->output()};
    
  my $db;
  if ($self->GENEDB) {
    $db = $self->get_dbadaptor($self->GENEDB);
  } else {
    $db =$self->db;
  }

  my $adap = $db->get_ProteinFeatureAdaptor;
  
  foreach my $feat(@features) {
    $adap->store($feat, $feat->seqname);
  }
  
  return 1;
}



##################################
sub run {
  my ($self,$dir) = @_;

  throw("Runnable module not set") unless ($self->runnable());

  my $db;
  if ($self->GENEDB) {
    $db = $self->get_dbadaptor($self->GENEDB);
  } else {
    $db =$self->db;
  }

  $db->dbc->disconnect_when_inactive(1);
  my ($run) = @{$self->runnable};
  $run->run($dir);
  $db->dbc->disconnect_when_inactive(0);

  $self->output($run->output);
}


###################################
sub query{
  my ($self, $query) = @_;

  if(defined $query){
    if (ref($query)) {      
      if (not $query->isa('Bio::PrimarySeqI')) {
        throw("Must pass RunnableDB:query a Bio::PrimarySeqI " 
              . "not a ".$query);
      } 
    } elsif (not -e $query) {
      throw("Must pass RunnableDB::query a filename that exists " . ref($query));
    }
    $self->{_query} = $query;        
  }
  return $self->{_query};
}

#####################################
sub output {
  my ($self, $output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($output){
    if(ref($output) ne 'ARRAY'){
      throw('Must pass RunnableDB:output an array ref not a '.$output);
    }
    $self->{'output'} = $output;
  }
  return $self->{'output'};
}


####################################
#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($PROTEINANNOTATION_CONFIG_BY_LOGIC);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values

  if (defined $self->BASE_DIR and 
      $self->BASE_DIR and
      not -d $self->BASE_DIR) {
    throw("BASE_DIR " . $self->BASE_DIR . " could not be found")
  }
}


sub BASE_DIR {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_base_dir} = $val;
  }

  if (not exists $self->{_base_dir}) {
    return undef;
  } else {
    return $self->{_base_dir};
  }

}



sub GENEDB {
  my ($self, $val) = @_;


  if (defined $val) {
    $self->{_gene_db} = $val;
  }

  if (not exists $self->{_gene_db}) {
    return undef;
  } else {
    return $self->{_gene_db};
  }
} 


sub DNA_DB {
  my ($self, $val) = @_;


  if (defined $val) {
    $self->{_dna_db} = $val;
  }

  if (not exists $self->{_dna_db}) {
    return undef;
  } else {
    return $self->{_dna_db};
  }
}


1;
