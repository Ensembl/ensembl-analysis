=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM->new(
					     -genomic => $seq,
					     -features => \@features,
								);
    $obj->run

    my @newfeatures = $obj->output;

or something similar????

=head1 DESCRIPTION

Finds which pfam domains a particular swissprot hit matches, finds the related hmms and runs genewiseHmm

=head1 CONTACT

rds@sanger.ac.uk
refactored by Sindhu K. Pillai B<sp1@sanger.ac.uk>
=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM;


use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Tools::BlastDBTracking;
use Bio::EnsEMBL::Analysis::Tools::FeatureFilter;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Bio::Tools::BPlite;
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Runnable::Finished::GenewiseHmm;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
#use DB_File;
use Fcntl;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

    Arg      : Bio:Seq obj, reference to array of uniprot features, location of hmmfetch program, location of hmmdb, location of dbm index
    Function : Make a new HalfwiseHMM object defining the above variables
    Returntype: Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM
    Exception: Will throw an exception if no genomic sequence is defined or no features are passed
    Caller   :
    Example  : $halfwise = Bio::EnsEMBL::Analysis::Runnable::HalfwiseHMM->new(genomic => $seq
                                                                              features => \@features);

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    $self->{'_query'} = undef;
    $self->{'_input_features'} = [];
    $self->{'_output_features'} = [];
    $self->{'_hmmfetch'} = undef;
    $self->{'_hmmdb'} = undef;
    $self->{'_hmmfilename'} = undef; #file to store protein profile hmms
    $self->{'_errorfile'} = undef; #file to store any errors hmmfetch throw
    $self->{'_dbm_file'} = undef;
    #print STDERR "args = @args\n";
    my ($genomic, $features, $hmmfetch,$hmmconvert, $hmmdb, $pfamdb, $memory, $options, $analysis,$program) = rearrange([qw(QUERY
											   FEATURES
											   HMMFETCH
											   HMMCONVERT
											   HMMDB
											   PFAMDB
											   MEMORY
											   OPTIONS
											   ANALYSIS
											   PROGRAM)], @args);
    throw("No genomic sequence input") unless defined($genomic);
    $self->query($genomic);
    throw("No features input") unless defined($features);
    $self->add_input_features($features);
    throw("No pfam database") unless $pfamdb;
    $self->pfamDB($pfamdb);
    $self->options($options);
    $self->program($program);
    throw("No analysis") unless defined($analysis);
    $self->analysis($analysis);
    if ($hmmfetch){
	$self->hmmfetch($hmmfetch);
    } else {
	$self->hmmfetch('/software/pubseq/bin/hmmer/hmmfetch');
    }
    if ($hmmconvert){
	$self->hmmconvert($hmmconvert);
    } else {
	$self->hmmconvert('/software/pubseq/bin/hmmer/hmmconvert');
    }
    if($hmmdb){
	$self->hmmdb($hmmdb);
    } else {
	$self->hmmdb('/data/blastdb/Finished/Pfam-A.hmm');

    }
    $self->memory ($memory)    if (defined($memory));
    return $self; # success - we hope!
}




#################
#GET/SET METHODS#
#################

sub analysis {
    my ( $self, $analysis ) = @_;
    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}


sub pfamDB{
    my ($self, $dbobj) = @_;
    $self->{'_pfamDB'} = $dbobj if $dbobj;
    throw("Not a Bio::EnsEMBL::DBSQL::DBAdaptor")
        unless $self->{'_pfamDB'}->isa("Bio::EnsEMBL::DBSQL::DBConnection");
    return $self->{'_pfamDB'};
}

sub pfam_db_version{
    my ($self) = @_;
    unless($self->{'_pfam_db_version'}){
        my $db = $self->pfamDB();
        $self->{'_pfam_db_version'} = $db->get_meta_value_by_key('version');
    }
    return $self->{'_pfam_db_version'};
}

sub pfam_ls_version{
    my ($self) = @_;
    unless($self->{'_pfam_ls_version'}){
	my $ver = Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::get_db_version_mixin(
            $self, '_pfam_ls_version', $self->hmmdb
            );
    }
    return $self->{'_pfam_ls_version'};
}

sub program{
    my ($self, $program) = @_;
    $self->{'_program'} = $program if $program;
    return $self->{'_program'};
}
=head2 query

  Arg      : Bio:Seq object
  Function : get/set method for query
  Exception: throws an exception if not passed a Bio:PrimarySeq
  Caller   :

=cut

sub query {
    my( $self, $value ) = @_;
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_query'} = $value;
    }
    return $self->{'_query'};
}

=head2 add_input_features

    Arg      : reference to array of features
    Function : get/set method for input features
    Exception: throws an exception if not passed an array ref
    Caller   :

=cut

sub add_input_features{
    my ($self, $features) = @_;
    if (ref($features) eq "ARRAY") {
      push(@{$self->{'_input_features'}},@$features);
    } else {
      throw("[$features] is not an array ref.");
    }
}

=head2 all_input_features

    Arg      : none
    Function : returns all input features
    Exception: none
    Caller   :

=cut


sub all_input_features{

  my ($self) = @_;

  return @{$self->{'_input_features'}};

}


#methods that probably won't be needed


=head2 hmmfetch

    Arg      : location of hmmfetch program
    Function : gets/sets location of hmmfetch program
    Exception: none
    Caller   :

=cut

sub hmmfetch {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmfetch'} = $args;
    }
    return $self->{'_hmmfetch'};
}

=head2 hmmconvert

    Arg      : location of hmmconvert program
    Function : gets/sets location of hmmfconvert program
    Exception: none
    Caller   :

=cut

sub hmmconvert {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmconvert'} = $args;
    }
    return $self->{'_hmmconvert'};
}

=head2 hmmdb

    Arg      : location of hmmdb
    Function : gets/sets location of hmmdb
    Exception: none
    Caller   :

=cut

sub hmmdb {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmdb'} = $args;
    }
    return $self->{'_hmmdb'};
}


=head2 hmmfilename

    Arg      : hmm filename
    Function : gets/sets hmmfilename
    Exception: none
    Caller   :

=cut

sub hmmfilename {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_hmmfilename'} = $args;
    }

    return $self->{'_hmmfilename'};

}

sub hmm3filename {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_hmm3filename'} = $args;
    }

    return $self->{'_hmm3filename'};

}

sub hmm2filename {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_hmm2filename'} = $args;
    }

    return $self->{'_hmm2filename'};

}
=head2 dbm_file

    Arg      : dbm filename
    Function : gets/sets dbm filename
    Exception: none
    Caller   :

=cut


sub dbm_file {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_dbm_file'} = $args;
    }

    return $self->{'_dbm_file'};

}

=head2 memory

    Arg      : value memory to be set to, this is an option of how much memory genewise can use when it is run
    Function : gets/sets memory
    Exception: none
    Caller   :

=cut

sub memory {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_memory'} = $arg;
    }

    return $self->{'_memory'} || 1000000;
}

##########
#ANALYSIS#
##########

=head2 run

    Arg      : directory if to be set to anything other than tmp
    Function : runs functions which run genewisehmm
    Exception: none
    Caller   :

=cut

sub run {

    my ($self, $dir) = @_;

    my $swissprot_ids = $self->get_swissprot_ids();
    my $pfam_ids = $self->get_pfam_hmm($swissprot_ids,$dir);

    #$self->create_genewisehmm($pfam_ids, $dir);
    #$self->create_genewisehmm_individually($pfam_ids, $dir);
    $self->create_genewisehmm_complete($pfam_ids);
}

=head2 get_swissprot_ids

    Arg      : none
    Function : gets swissprot ids and strands from input features and returns a hash of these
    Exception: none
    Caller   :

=cut

### 1, -1 and 0 are used to repesent strand, 1 and -1 have the same meaning as in teh ensembl databases 1 is the forward strand and -1 is the reverse
### 0 is used if a particular id has hit on both strands

sub get_swissprot_ids{
    my ($self) = @_;
    my @features = $self->all_input_features();
    my $swissprot_ids; #{}
    foreach my $f(@features){
        my $hname = $f->hseqname();
        my $strand = $f->strand();
        #print "swissprot id = ".$hname."\n";
        if(!defined($swissprot_ids->{$hname})){
            $swissprot_ids->{$hname} = $strand;
        }else{
            $swissprot_ids->{$hname} = ($strand eq $swissprot_ids->{$hname} ? $strand : 0);
        }
    }
    return $swissprot_ids;
}

sub get_pfam_hmm {
	my ($self, $uniprot_ids, $dir) = @_;
    my $pfam_accs;
    my $pfam_lookup;
    my $db = $self->pfamDB();

    $self->workdir('/tmp') unless($self->workdir($dir));
    $self->checkdir();

	my $sql = "SELECT
	    			DISTINCT(CONCAT(a.pfamA_acc,'.',a.version)),
	    			a.pfamA_id,
	    			a.description
				FROM
				    pfamA_reg_full_significant f,
				    pfamseq p,
				    pfamA a
				WHERE
				    f.auto_pfamseq = p.auto_pfamseq AND
				    f.in_full = 1 AND
				    a.auto_pfamA = f.auto_pfamA AND
				    p.pfamseq_acc IN ";

    my @accessions = keys(%$uniprot_ids);
    map ( s/(\w*(-\d+)?)(\.\d+)?/'$1'/,@accessions);

    # deal with empty @accessions
    next unless defined @accessions;

    $sql .= "(".join(',',@accessions).")";

	#print STDOUT $sql."\n";

	my $sth = $db->prepare($sql);
    $sth->execute();
    my ($pfam_acc, $pfam_id, $description);
    $sth->bind_columns(\($pfam_acc, $pfam_id, $description));

    while(my $row = $sth->fetchrow_arrayref()){
    	$pfam_lookup->{$pfam_id}	= [$pfam_acc, $description];
    	$pfam_accs->{$pfam_acc}		= 1;
    }
    $sth->finish();

	$self->get_hmmdb($pfam_lookup);
	$self->pfam_lookup($pfam_lookup);

    return $pfam_accs;
}

=head2 get_pfam_ids

    Arg      : reference to a hash of swissprotids and stands
    Function : gets all pfam ids for a particular swissprot id and puts them in a hash along with the strand the swissprot id is found on and returns its
    Exception: warns if swissprot id has no pfam domains
    Caller   :

=cut

sub get_pfam_ids{

    my ($self, $swissprot_ids) = @_;
    my $pfam_accs;
    my $pfam_lookup;
    # CREATE TEMPORARY TABLE
    my $tbl_name = "pipeline_tmp_$$";
    my $create_table = qq{CREATE TEMPORARY TABLE $tbl_name(
                            pfamseq_id varchar(12) NOT NULL PRIMARY KEY,
                            strand enum('1','0','-1') DEFAULT '0'
                            )TYPE = HEAP
			}; # should this be HEAP?
                            # There's never gonna be that many matches
                            # to exceed tmp_table_size = 10048576 ???
    my $db = $self->pfamDB();
    my $sth = $db->prepare($create_table);
    $sth->execute();
    $sth->finish();
    # INSERT
    my (@binds, @values);
    my $sql = qq{INSERT IGNORE INTO $tbl_name (pfamseq_id, strand) VALUES };
    while (my ($swiss_id, $strand) = each(%$swissprot_ids)){
        #print STDERR "$swiss_id : $strand\n";
    # strip off sequence version from uniprot ID (introduced in Uniprot release 7.0)
    $swiss_id =~ s/\.\d+//;
	push(@binds, $swiss_id, $strand);
	push(@values, qq{ (?, ?)});
    }
    if(scalar(@values)){
	$sql .= join(", ", @values);
	#warning $sql;
	$sth = $db->prepare($sql);
	$sth->execute(@binds);
	$sth->finish();
    }
    # SELECT
#     my $select = qq{SELECT a.pfamA_acc, t.strand, t.pfamseq_id, a.pfamA_id, a.description
#                         FROM pfamA a, pfamA_reg_full f, pfamseq p, $tbl_name t
#                         WHERE p.auto_pfamseq = f.auto_pfamseq
#                         && a.auto_pfamA = f.auto_pfamA
#                         && p.pfamseq_id = t.pfamseq_id
#                         && f.in_full = 1};
    my $select = qq{SELECT CONCAT(a.pfamA_acc,".",a.version), t.strand, t.pfamseq_id, a.pfamA_id, a.description
			FROM pfamA_reg_full f, pfamseq p, $tbl_name t,  pfamA a
			WHERE f.auto_pfamseq = p.auto_pfamseq
			&& p.pfamseq_acc     = t.pfamseq_id
			&& f.in_full         = 1
			&& a.auto_pfamA      = f.auto_pfamA;};
    $sth = $db->prepare($select);
    $sth->execute();
    my ($pfam_acc, $strand, $swall, $pfam_id, $description);
    $sth->bind_columns(\($pfam_acc, $strand, $swall, $pfam_id, $description));

    while(my $row = $sth->fetchrow_arrayref()){
        print STDERR "$swall : $strand : $pfam_acc : $pfam_id\n";
        # making pfam_lookup for later ...
        $pfam_lookup->{$pfam_id} = [$pfam_acc, $description];
        if(defined($pfam_accs->{$pfam_acc})){
            $pfam_accs->{$pfam_acc} = ($strand eq $pfam_accs->{$pfam_acc} ? $strand : 0);
        }else{
            $pfam_accs->{$pfam_acc} = $strand;
        }
    print "\n pfam_acc: $pfam_acc, strand: $strand, uniprot: $swall,pfamid: $pfam_id,description: $description";

    }
    $self->pfam_lookup($pfam_lookup); # store the lookup somewhere useful.
    return $pfam_accs;
}

# get/set for lookup multi-valued hash { pfam_id => [pfam_acc, pfam_desc], ... }
# can append multiple to the lookup (see { %{$self->{'_pfam_lookup'}}, %{$hash_ref} })
sub pfam_lookup{
    my ($self, $hash_ref) = @_;
    if(ref($hash_ref) eq 'HASH'){
        $self->{'_pfam_lookup'} ||= {};
        $self->{'_pfam_lookup'}   = { %{$self->{'_pfam_lookup'}}, %{$hash_ref} };
    }
    return $self->{'_pfam_lookup'};
}

=head2 create_genewisehmm_individually

    Arg      : reference to pfam_id hash and directory if it is to be set to anything other than temp
    Function : runs through each pfam id runs, get_hmm and creates and runs the GenewiseHmm runnable
    Exception: throws an exception if anything other than 1, -1 or 0 is found for strand
    Caller   :

=cut

sub create_genewisehmm_individually{
    my ($self, $pfam_ids, $dir) = @_;
    print STDERR "getting hmms\n"; ##########
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
    my $memory = $self->memory;

    print STDERR "there are ".scalar(keys(%$pfam_ids))." pfam ids to run\n"; ##########

    while(my ($pfam_id, $strand) = each(%$pfam_ids)){
        print STDERR "doing the hmm for id: $pfam_id\n"; ##########
        $self->get_hmm($pfam_id);
        if (-z $self->hmmfilename){
            warning("hmm file not created :$!");
            next;
        }
        my @strands = ();
        if($strand == 1){
            push(@strands, $strand);
        }elsif($strand == -1){
            push(@strands, $strand);
        }elsif($strand == 0){
            push(@strands, -1);
            push(@strands, 1);
        }else{
            throw("strand = ".$strand." something funnies going on : $!\n");
        }
        foreach my $s (@strands){
            $self->run_genewisehmm($s);
            #print STDERR "running on strand: $s\n"; ##########
            #my $genewisehmm = $self->get_GenewiseHMM($s, $memory);
            #$genewisehmm->run();
            #my @features = $genewisehmm->output();
            #$self->display_genes(\@features); ##########
            #print STDERR "adding ".scalar(@features)." to output\n"; ##########
            #$self->add_output_features(\@features);
        }
        print STDERR "removing hmmfile: ".$self->hmmfilename()."\n"; ##########
        unlink $self->hmmfilename();
    }
    throw("finished this bit");
    return undef;
}

=head2 create_genewisehmm_complete

    Arg      : reference to pfam_id hash and directory if it is to be set to anything other than temp
    Function : creates a hmm databases of pfam ids, then creates and runs the GenewiseHmm runnable
    Exception:
    Caller   :

=cut

sub create_genewisehmm_complete{
    my ($self, $pfam_ids) = @_;
    print STDERR "getting hmms\n"; ##########
    print STDERR "\nthere are ".scalar(keys(%$pfam_ids))." pfam ids in database\n"; ##########

    return unless scalar(keys(%$pfam_ids));

    print STDERR "doing the hmm for ids: ". join(" ", keys(%$pfam_ids)) . "\n"; ##########

    $self->run_genewisehmm(0);
    print STDERR "removing hmm files:\n"; ##########
    unlink $self->hmm3filename();
    unlink $self->hmm2filename();

    return undef;
}

=head2 run_genewisehmm

    Arg       :
    Function  :
    Exception :
    Caller    : $self->create_genewisehmm_complete, $self->create_genewisehmm_individually
=cut

sub run_genewisehmm{
    my ($self, $strand) = @_;
    my $memory = $self->memory;
    print STDERR "running on strand: $strand\n"; ##########
    my $genewisehmm = $self->get_GenewiseHMM($strand, $memory);
    $genewisehmm->run();
    my @features = $genewisehmm->output();
    $self->display_genes(\@features); ##########
    print STDERR "adding ".scalar(@features)." to output\n"; ##########
    $self->add_output_features(\@features);
}

=head2 get_GenewiseHMM

    Arg      : the strand of the hit and the memory to be used by genewisehmm
    Function : runs the new method of GenewiseHmm and returns the runnable
    Exception: none
    Caller   :

=cut

sub get_GenewiseHMM{
  my ($self, $strand, $memory) = @_;
  my $reverse = ($strand == -1 ? 1 : undef);
  print STDERR "creating genewisehmm strand $strand reverse $reverse\n"; ##########
  print STDERR "OPTIONS To Genewise: ".$self->options()."\n"; ##########
#  $genewisehmm->set_environment("/usr/local/ensembl/data/wisecfg/");
  $ENV{WISECONFIGDIR} = "/usr/local/ensembl/data/wisecfg/";

  my $genewisehmm =
    Bio::EnsEMBL::Analysis::Runnable::Finished::GenewiseHmm->new(
    							'-query'    => $self->query(),
                                '-memory'   => $memory,
                                '-hmmfile'  => $self->hmm2filename(),
                                '-reverse'  => $reverse,
                                '-genewise' => $self->program(),
                                '-options'  => $self->options(),
								'-analysis' => $self->analysis(),
    );
  return $genewisehmm;

}


=head2 get_hmm

    Arg      : id of hmm to be fetched normally a pfam id
    Function : runs hmmfetch on given id
    Exception: thows if no hmmfile is created
    Caller   :

=cut

sub get_hmm{
  my ($self, $id) = @_;
  #print "getting hmms\n";
  $self->hmmfilename($id.".".$$.".hmm");
  #print "getting hmm for id ".$id."\n";
  my $command =  $self->hmmfetch." ".$self->hmmdb." ".$id." > ".$self->hmmfilename;

  #print "command = ".$command."\n";
  #system ('pwd');
  eval{
    system($command);
  };
  if($@){
    warning("hmmfetch threw error : $@\n");
  }
  if (-z $self->hmmfilename){
    warning("hmm file not created :$!");
  }elsif(-e $self->hmmfilename){
    open(ERROR, $self->hmmfilename) or die "couldn't open error file : $!";
    while(<ERROR>){
      if(/no such hmm/i){
	print STDERR "error message : ".$_."\n";
      }
    }
    close(ERROR) or die "couldn't close error file : $!";
  }
}

sub get_hmmdb{

    my ($self, $pfam_ids) = @_;
    my @pfamIds = keys(%$pfam_ids);
    print STDERR "getting the hmms for ids: ". join(" ", @pfamIds) . "\n"; ##########
    $self->hmm3filename("halfwise.$$.hmmdb3");
    $self->hmm2filename("halfwise.$$.hmmdb2");

    foreach my $id(@pfamIds){
        print "getting hmm for id ".$id."\n";
        my $command =  $self->hmmfetch." ".$self->hmmdb." ".$id." >> ".$self->hmm3filename;
        eval{
          system($command);
        };
        if($@){
          warning("hmmfetch threw error : $@\n");
        }
        if (-z $self->hmm3filename){
          warning("hmm file not created :$!");
        }elsif(-e $self->hmm3filename){
          open(ERROR, $self->hmm3filename) or die "couldn't open error file : $!";
          while(<ERROR>){
            if(/no such hmm/i){
              print STDERR "error message : ".$_."\n";
            }
          }
          close(ERROR) or die "couldn't close error file : $!";
        }
    }
    # convert the hmm file from HMMER3 to HMMER2 format for genewise
    my $command = $self->hmmconvert." -2 ".$self->hmm3filename." > ".$self->hmm2filename;
 	print STDERR "$command\n";
    eval{
          system($command);
    };
    if($@){
          warning("hmmconvert threw error : $@\n");
    }



}
=head2 add_output_features

    Arg      : array ref to array of output features
    Function : adds the array of features to the output feature array
    Exception: throws an exception if not passed an array ref
    Caller   :

=cut

sub add_output_features{
    my ($self, $features) = @_;
    #print "adding ".scalar(@$features)." to output\n";
    if (ref($features) eq "ARRAY") {
      push(@{$self->{'_output_features'}},@$features);
    } else {
      throw("[$features] is not an array ref.");
    }
    #print "total feature no = ".scalar(@{$self->{'_output_features'}})."\n";

}




################
#OUTPUT METHODS#
################


=head2 output

    Arg      : none
    Function : returns the array of output features
    Exception: none
    Caller   :

=cut


sub output{
  my ($self) = @_;
  #print "outputing data\n";
  #print "returning ".scalar(@{$self->{'_output_features'}})." to output\n";
  return @{$self->{'_output_features'}};

}

sub display_genes {
  my ($self, $result) = @_;
  #Display output
  my @results = @$result;
  foreach my $obj (@results)
    {
      print STDERR ("gene:  ".$obj->gffstring. "\n");
      if ($obj->sub_SeqFeature)
	{
	  foreach my $exon ($obj->sub_SeqFeature)
	    {
	      print STDERR "Exon:  ".$exon->gffstring."\n";
	      if ($exon->sub_SeqFeature)
		{
		  foreach my $sub ($exon->sub_SeqFeature){
		    print STDERR "supporting features:  ".$sub->gffstring."\n";
		  }
		}
	    }
	}
    }
}


1;
