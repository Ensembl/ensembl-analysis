
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature




=head1 SYNOPSIS

my $clone = 
  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$clone->fetch_input();
$clone->run();
$clone->write_output(); #writes to DB

=head1 DESCRIPTION

This object maps clone sequences to a genome,
and writing the results as Features. 

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::Config::ExonerateAlignFeature;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

############################################################
sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  
  #Because we dont know whether the sort of adapter (compara, hive, core)
  #we'll be passed, it's better to just remake the core dbadaptor and
  #keep it on ourselves as the dna db - this has to be the intent of the
  #dataadaptor input through the rulemanager / beekeeper / whatever

  $self->read_and_check_config($EXONERATE_ALIGNFEAT_CONFIG);
  return $self;
}

sub fetch_input {
  my ($self) = @_;

  my $logic = $self->analysis->logic_name;

  ##########################################
  # set up the target (genome)
  ##########################################

  my $target = $self->GENOMICSEQS;
  my $query =  $self->QUERYSEQS;

  if ( -e $target ){
    if(-d $target ) {
      warn ("Target $target is a directory of files\n");
    }elsif (-s $target){
      warn ("Target $target is a whole-genome file\n");
    }else{
      throw("'$target' isn't a file or a directory?");
    }
  } else {
    throw("'$target' could not be found");
  }
  
  my ($chunk_number, $chunk_total);

  if (-e $query and -d $query) {
    # query seqs is a directory; input id will be a file in that directory
    $query = "$query/" . $self->input_id;
    if (not -e $query) {
      throw( "Query file '$query' does not exist'\n");
    }
  }
  elsif (-e $query and -s $query) {
    # query seqs is a single file; input id will correspond to a chunk number
    my $iid_regexp = $self->IIDREGEXP;
    
    throw("When your input ids are not filenames, you must define ".
          "IIDREGEXP in config to enable inference of chunk number and total")
        if not defined $iid_regexp;
    
    ($chunk_number, $chunk_total) = $self->input_id =~ /$iid_regexp/;
  } else {
    throw("'$query' refers to something that could not be made sense of\n");
  }

  ##########################################
  # setup the runnable
  ##########################################

  my %parameters = %{ $self->parameters_hash };
  
  if (not exists( $parameters{-options} )
      and defined $self->OPTIONS) {
    $parameters{-options} = $self->OPTIONS;
  }
  
  #print STDERR "PROGRAM FILE: ".$self->analysis->program_file."\n";
  my $program = $self->PROGRAM;
  $program = $self->analysis->program_file if not defined $program;
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature->new(
    -analysis           => $self->analysis,
    -program            => $program,
    -query_type         => $self->QUERYTYPE,
    -query_file         => $query,
    -query_chunk_number => $chunk_number ? $chunk_number : undef,
    -query_chunk_total  => $chunk_total ? $chunk_total : undef,
    -target_file        => $target,
    %parameters,
  );
  
  $self->runnable($runnable);

}

############################################################

sub run {
  my ($self) = @_;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my ($runnable) = @{$self->runnable};
  
  $runnable->run;
  my $features = $runnable->output;

  if ($self->filter) {
    $features = $self->filter->filter_results($features);
  }
  
  my $output_db = $self->create_output_db;
  $self->output_db($output_db) ;  

  $self->process_features($features);
  $self->output($features);
}

############################################################

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb = $self->output_db; 
  my $fa = $outdb->get_DnaAlignFeatureAdaptor;

  foreach my $f (@{$self->output}){  
    eval{ 
      $fa->store($f);
    };
    if ($@) {
      $self->throw("Unable to store clone feature!\n $@");
    }
  }
}

############################################################

sub process_features {
  my ( $self, $flist,$output_db  ) = @_;

  my $slice_adaptor = $self->output_db->get_SliceAdaptor;

  my %genome_slices;

  foreach my $f (@$flist) {

    $f->analysis( $self->analysis );

    # get the slice based on the seqname stamped on in the runnable
    my $slice_id = $f->seqname;

     if ( not exists $genome_slices{$slice_id} ) {
      # assumes genome seqs were named in the Ensembl API Slice naming
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
   }
    my $slice = $genome_slices{$slice_id};
 
    $f->slice($slice);  
  }

}

sub create_output_db {
  my ($self) = @_;

  my $outdb;
  my $dnadb;

  if ( $self->OUTDB) {
    $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{ $self->OUTDB }); 
  } else {
    $outdb = $self->db;
  } 
  $self->db->disconnect_when_inactive(1); 
  return $outdb;
}


sub output_db {
 my ( $self,$arg ) = @_ ; 
  
 if ( defined $arg && $arg->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")) {  
   $self->{'output_db'} = $arg ; 
 } 
 return $self->{'output_db'} ; 
} 



sub filter {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_feature_filter} = $val;
  }
  return $self->{_feature_filter};
}


#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (
    qw(
      QUERYSEQS
      QUERYTYPE
      GENOMICSEQS      
    )
  ){
    if ( not defined $self->$config_var ){
      throw("You must define $config_var in config for logic '$logic'");
    }
  }

  # output db does not have to be defined, but if it is, it should be a hash
  if ($self->OUTDB && ref( $self->OUTDB ) ne "HASH") {
    throw("OUTDB in config for '$logic' must be a hash ref of db connection pars.");
  }

  if ($self->FILTER) {
    if (not ref($self->FILTER) eq "HASH" or
        not exists($self->FILTER->{OBJECT}) or
        not exists($self->FILTER->{PARAMETERS})) {
          
      throw("FILTER in config fo '$logic' must be a hash ref with elements:\n" . 
            "  OBJECT : qualified name of the filter module;\n" .
            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
    } else {
      my $module = $self->FILTER->{OBJECT};
      my $pars   = $self->FILTER->{PARAMETERS};
      
      (my $class = $module) =~ s/::/\//g;
      eval{
        require "$class.pm";
      };
      throw("Couldn't require ".$class." ExonerateAlignFeature:require_module $@") if($@);
    
      $self->filter($module->new(%{$pars}));
    }
  }

}


sub QUERYSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYSEQS'} ) ) {
    return $self->{'_CONFIG_QUERYSEQS'};
  } else {
    return undef;
  }
}


sub QUERYTYPE {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYTYPE'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYTYPE'} ) ) {
    return $self->{'_CONFIG_QUERYTYPE'};
  } else {
    return 'dna';
  }
}


sub GENOMICSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_GENOMICSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_GENOMICSEQS'} ) ) {
    return $self->{'_CONFIG_GENOMICSEQS'};
  } else {
    return undef;
  }
}

sub IIDREGEXP {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_IIDREGEXP'} = $value;
  }

  if ( exists( $self->{'_CONFIG_IIDREGEXP'} ) ) {
    return $self->{'_CONFIG_IIDREGEXP'};
  } else {
    return undef;
  }
}

sub FILTER {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_FILTER'} = $value;
  }
  
  if (exists($self->{'_CONFIG_FILTER'})) {
    return $self->{'_CONFIG_FILTER'};
  } else {
    return undef;
  }
}

sub OUTDB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OUTDB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OUTDB'} ) ) {
    return $self->{'_CONFIG_OUTDB'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OPTIONS'} ) ) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}


sub PROGRAM {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_PROGRAM'} = $value;
  }

  if ( exists( $self->{'_CONFIG_PROGRAM'} ) ) {
    return $self->{'_CONFIG_PROGRAM'};
  } else {
    return undef;
  }
}



###############################################
###     end of config
###############################################

1;
