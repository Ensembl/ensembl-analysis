package Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;

use strict;
use warnings;
use FileHandle;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Analysis::Tools::BPlite;

use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($filename, $regex, $query, $target) = rearrange(['FILENAME', 'REGEX',
                                                       'QUERY_TYPE',
                                                       'DATABASE_TYPE',
                                                       ], @args);
  
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->regex('(\w+)\s+');
  ######################

  $self->filename($filename) if($filename);
  $self->regex($regex) if(defined $regex);
  $self->query_type($query) if($query);
  $self->database_type($target) if($target);
  return $self;
}


sub filename{
  my $self = shift;
  $self->{'filename'} = shift if(@_);
  return $self->{'filename'};
}


sub regex{
  my $self = shift;
  $self->{'regex'} = shift if(@_);
  return $self->{'regex'};
}


sub query_type{
  my ($self, $qtype) = @_; 
  if($qtype){
    throw("Query type must be either dna or pep not ".$qtype)
      unless(($qtype =~ /dna/i) || ($qtype = /pep/i));
    $self->{'query_type'} = $qtype;
  }
  return $self->{'query_type'};
}

sub database_type{
  my ($self, $dtype) = @_; 
  if($dtype){
    throw("Database type must be either dna or pep not ".$dtype)
      unless(($dtype =~ /dna/i) || ($dtype = /pep/i));
    $self->{'database_typxe'} = $dtype;
  }
  return $self->{'database_type'};
}


sub output{
  my ($self, $output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($output){
    throw("Must pass and arrayref not a ".$output." BPliteWrapper:output")
      unless(ref($output) eq 'ARRAY');
    push(@{$self->{'output'}}, @$output);
  }
  return $self->{'output'};
}

sub feature_factory{
  my ($self, $feature_factory) = @_;
  if($feature_factory){
    $self->{'feature_factory'} = $feature_factory;
  }
  if(!$self->{'feature_factory'}){
    $self->{'feature_factory'} = Bio::EnsEMBL::Analysis::Tools::FeatureFactory
      ->new();
  }
  return $self->{'feature_factory'};
}

sub parse_file{
  my ($self, $file) = @_;
  if(!$file){
    $file = $self->filename;
  }
  $self->filename($file);
  throw("File ".$file." must exist to be parsed BPliteWrapper:parse_file ")
    unless(-e $file);
 
  my $bplite = $self->get_parser($file);
  $self->get_hsps($bplite);
  return $self->output;
}


sub get_parser{
  my ($self, $file) = @_;
  my $fh = new FileHandle;
  $fh->open("<". $file);
  my $bplite = Bio::EnsEMBL::Analysis::Tools::BPlite->new
    (
     -fh => $fh,
    );
  return $bplite;
}

sub get_hsps{
  my ($self, $parser) = @_;
  my $regex = $self->regex;
  my @output;
 NAME: while(my $sbjct = $parser->nextSbjct){
    my ($name) = $sbjct->name =~ /$regex/;
    throw("Error parsing name from ".$sbjct->name." check your ".
          "blast setup and blast headers") unless($name);
  HSP: while (my $hsp = $sbjct->nextHSP) {
      push(@output, $self->split_hsp($hsp, $name));
    }
  }
  $self->output(\@output);
}



sub split_hsp {
    my ($self,$hsp,$name) = @_;
    my $qstrand = $hsp->query->strand;
    my $hstrand = $hsp->subject->strand;
    my ($qinc,   $hinc)    = $self->find_increments($qstrand,$hstrand);
    
    my @qchars = split(//,$hsp->querySeq);  # split alignment into array of
                                            # chars
    my @hchars = split(//,$hsp->sbjctSeq);  # ditto for hit sequence
    
    my $qstart = $hsp->query->start(); # Start off the feature pair start
    my $hstart = $hsp->subject->start(); # ditto
    my $qend   = $hsp->query->start(); # Set the feature pair end also
    my $hend   = $hsp->subject->start(); # ditto
   if ($qstrand == -1) {
      $qstart = $hsp->query->end;
      $qend   = $hsp->query->end;
    }
    if ($hstrand == -1) {
      $hstart = $hsp->subject->end;
      $hend   = $hsp->subject->end;
    }
    
    my $count = 0; # counter for the bases in the alignment
    my $found = 0; # flag saying whether we have a feature pair
    

    my @tmpf;

    while ($count <= $#qchars) {
      # We have hit an ungapped region.  Increase the query and hit 
      #counters and flag that we have a feature pair.
      
      if ($qchars[$count] ne '-' &&
          $hchars[$count] ne '-') {
        
        $qend += $qinc;
        $hend += $hinc;
        
        $found = 1;
      } else {
        
        # We have hit a gapped region.  If the feature pair flag is set 
        # ($found) then make a feature pair, store it and reset the start 
        # and end variables.
        
        if ($found == 1) {
          my $fp = $self->convert_to_featurepair($qstart, $qend, $qstrand, 
                                                 $qinc, $hstart, $hend, 
                                                 $hstrand, $hinc, $name,
                                                 $hsp->score, 
                                                 $hsp->percent, $hsp->P,
                                                 $hsp->positive, 
                                                 $hsp->match);
          push(@tmpf,$fp);
        }
        
        # We're in a gapped region.  We need to increment the sequence that
        # doesn't have the gap in it to keep the coordinates correct.
        # We also need to reset the current end coordinates.
        
        if ($qchars[$count] ne '-') {
          $qstart = $qend   + $qinc;
        } else {
          $qstart = $qend;
        }
        if ($hchars[$count] ne '-') {
          $hstart = $hend   + $hinc;
        } else {
          $hstart = $hend;
        }
        
        $qend = $qstart;
        $hend = $hstart;
        
        $found = 0;
      }
      $count++;
    }
    # Remember the last feature
    if ($found == 1) {
      my $fp = $self->convert_to_featurepair($qstart, $qend, $qstrand, 
                                             $qinc, $hstart, $hend, 
                                             $hstrand, $hinc, $name,
                                             $hsp->score, 
                                             $hsp->percent, $hsp->P,
                                             $hsp->positive, 
                                             $hsp->match);
      push(@tmpf,$fp);
    }
    my $fp;
    
    
    $qinc = abs( $qinc );
    $hinc = abs( $hinc );
    
    if( $qinc == 3 && $hinc == 1 ) {
      $fp = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tmpf);
    } elsif( $qinc == 1 && $hinc == 3 ) {
      $fp = Bio::EnsEMBL::PepDnaAlignFeature->new(-features => \@tmpf);
    } elsif( $qinc == 1 && $hinc == 1 ) {
      $fp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tmpf);
    } else {
      $self->throw( "Hardcoded values wrong?? " );
    }
    
    # helps debugging subsequent steps
    $fp->{'qseq'} = $hsp->querySeq();
    $fp->{'sseq'} = $hsp->sbjctSeq();
    
    # for compara
    $fp->positive_matches($hsp->positive);
    $fp->identical_matches($hsp->match);
    return $fp;
  }
#sub split_hsp{
#  my ($self, $hsp, $name) = @_;
#  print "\nSpliting hsp ".$name."\n";
#  my @qchars = split(//,$hsp->querySeq);  # split alignment into array of 
#                                          #characters
#  my @hchars = split(//,$hsp->sbjctSeq);  # ditto for hit sequence
#  my $qstart = $hsp->query->start();
#  my $qend = $hsp->query->start;
#  my $qstrand = $hsp->query->strand;
#  my $hstart = $hsp->subject->start;
#  my $hend = $hsp->subject->start;
#  my $hstrand = $hsp->subject->strand;
#  if($qstrand == -1){
#    $qstart = $hsp->query->end;
#    $qend = $hsp->query->start();
#  }
#  if($hstrand == -1){
#    $hstart = $hsp->subject->end;
#    $hend = $hsp->subject->start;
#  }
#  my ($qinc, $hinc) = $self->find_increments($qstrand, $hstrand);
#  my $count = 0; #counter for how may basepairs though the hit you are
#  my $found = 0; #marker to indicate a feature pair has been found
#  my @tmpfs;
#  print "There are ".@qchars." query characters and ".@hchars.
#    " hit characters\n";
#  print "query start = ".$qstart." query end ".$qend." length ".
#    ($qend - $qstart + 1)."\n";
#  print "hit start = ".$hstart." hit end ".$hend."length ".
#    ($hend - $hstart + 1)."\n";
#  while($count <= $#qchars){
#    if ($qchars[$count] ne '-' && $hchars[$count] ne '-') {
#      $qend += $qinc;
#      $hend += $hinc;
#      $found = 1;
#    }else{
#      if($found){
#        my $fp = $self->convert_to_featurepair($qstart, $qend, $qstrand,
#                                               $qinc, $hstart, $hend, 
#                                               $hstrand, $hinc, $name,
#                                               $hsp->score, $hsp->percent,
#                                               $hsp->P, $hsp->positive,
#                                               $hsp->match);
#        push(@tmpfs, $fp);
#      }
#      if ($qchars[$count] ne '-') {
#        $qstart = $qend   + $qinc;
#      } else {
#        $qstart = $qend;
#      }
#      if ($hchars[$count] ne '-') {
#        $hstart = $hend   + $hinc;
#      } else {
#        $hstart = $hend;
#      }
      
#      $qend = $qstart;
#      $hend = $hstart;
      
#      $found = 0;
#    }
#    $count++;
#  }
#  if($found){
#    my $fp = $self->convert_to_featurepair($qstart, $qend, $qstrand,
#                                           $qinc, $hstart, $hend, 
#                                           $hstrand, $hinc, $name,
#                                           $hsp->score, $hsp->percent,
#                                           $hsp->P, $hsp->positive,
#                                           $hsp->match);
#    push(@tmpfs, $fp);
#  }
#  foreach my $f(@tmpfs){
#    print "query ".$f->start." ".$f->end." ".$f->strand." length ".
#      ($f->end - $f->start + 1)." hit ".$f->hstart." ".
#      $f->hend." ".$f->strand." length ". ($f->hend - $f->hstart + 1)."\n";
#  }
#  $qinc = abs( $qinc );
#  $hinc = abs( $hinc );
#  my $fp;
#  if( $qinc == 3 && $hinc == 1 ) {
#    $fp = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tmpfs);
#  } elsif( $qinc == 1 && $hinc == 3 ) {
#    $fp = Bio::EnsEMBL::PepDnaAlignFeature->new(-features => \@tmpfs);
#  } elsif( $qinc == 1 && $hinc == 1 ) {
#    $fp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tmpfs);
#  } else {
#    throw( "Hardcoded values wrong?? " );
#  }
  
#  # helps debugging subsequent steps
#  $fp->{'qseq'} = $hsp->querySeq();
#  $fp->{'sseq'} = $hsp->sbjctSeq();
#  #for compara
#  $fp->positive_matches($hsp->positive);
#  $fp->identical_matches($hsp->match); 	 
#  return $fp;
#}

sub find_increments{
  my ($self, $qstrand, $hstrand) = @_;
  my $qinc   = 1 * $qstrand;
  my $hinc   = 1 * $hstrand;

  my $qtype = lc($self->query_type);
  my $htype = lc($self->database_type);

  if ($qtype eq 'dna' && $htype eq 'pep') {
    $qinc = 3 * $qinc;
  } 
  if ($qtype eq 'pep' && $htype eq 'dna') {
    $hinc = 3 * $hinc;
  }
  
  return ($qinc,$hinc);
}



sub convert_to_featurepair{
  my ($self, $qstart, $qend, $qstrand, $qinc, $hstart, $hend, $hstrand, 
      $hinc, $name, $score, $percent, $pvalue, $positive, $matches) = @_;
  my $tmpqend = $qend; $tmpqend -= $qinc;
  my $tmphend = $hend; $tmphend -= $hinc;
    
  my $tmpqstart = $qstart;
  my $tmphstart = $hstart;
  
  # This is for dna-pep alignments.  The actual end base
  # will be +- 2 bases further on.
  if (abs($qinc) > 1) {
    $tmpqend += $qstrand * 2;
  }
  if (abs($hinc) > 1) {
    $tmphend += $hstrand * 2;
  }
  
  # Make sure start is always < end
  if ($tmpqstart > $tmpqend) {
    my $tmp    = $tmpqstart;
    $tmpqstart = $tmpqend;
    $tmpqend   = $tmp;
  }
  if ($tmphstart > $tmphend) {
    my $tmp    = $tmphstart;
    $tmphstart = $tmphend;
    $tmphend   = $tmp;
  }
  my $fp = $self->feature_factory->create_feature_pair($tmpqstart, 
                                                       $tmpqend, $qstrand,
                                                       $score, $tmphstart,
                                                       $tmphend, $hstrand,
                                                       $name, $percent, 
                                                       $pvalue, undef, undef,
                                                       undef, $positive,
                                                       $matches);
  return $fp;
}
