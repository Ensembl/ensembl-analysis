

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut


package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;


@ISA = qw (Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub new {
  my ($class, @args) = @_;
 
  my $self = $class->SUPER::new(@args);
 
  my ($confirms) =  rearrange(['CONFIRM'], @args);

  if (defined $confirms) {
    $self->confirm_file($confirms);
  }

  if (not defined $self->database) {
    throw("You must supply a databaser to search");
  } 

  return $self;
}



sub run {
  my ($self, @args) = @_;

  my $scanning_code = 
      $self->create_scanning_code_from_patterns($self->database,
                                                $self->confirm_file);
  

  my @fps;

  if (-s $self->query) {
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -file => '<'.$self->query);
    while(my $seq = $seqio->next_seq) {
      push @fps, @{$self->scan_sequence($scanning_code, $seq)};
    }
    $seqio->close;
  } elsif (ref($self->query) and
           $self->query->isa("Bio::PrimarySeqI")) {
    push @fps, @{$self->scan_sequence($scanning_code, $self->query)};
  }

  $self->output(\@fps);
}


sub scan_sequence {
  my ($self, $code, $seq) = @_;

  my (@RESULTS);

  my $SEQID = $seq->display_id;
  my $SEQ = $seq->seq;

  eval($code);
  warn($@) if $@;

  my @features;
  foreach my $res (@RESULTS) {
    my $fp = $self->create_protein_feature($res->{start},
                                           $res->{end},
                                           $res->{score},
                                           $res->{seqid},
                                           0, 0,
                                           $res->{acc},
                                           $self->analysis,
                                           0, 0);
    push @features, $fp;
  }

  return \@features;
}

sub create_scanning_code_from_patterns {
  my ($self, $pattern_file, $confirm_file) = @_;

  my $confirm_hash = {}; 
  if (defined $confirm_file) {
    $self->_read_confirm_patterns($confirm_file,
                                  $confirm_hash);
  }

  my $scan_code = "\n";       # Perl-code to be constructed
  
  open (PAT,"<$pattern_file")
      or throw("Cannot open pattern file $pattern_file");
  while (<PAT>){
    my ($acc,$pattern,$name,$taxonrange) = split(/\s+/);
    $taxonrange =~ s/\?//g; 
    
    $scan_code .= "while(\$SEQ =~ /$pattern/g){\n";
    $scan_code .= "  my (\$match,\$start,\$end,\$confirmed) = (\$&,length(\$\`)+1,pos(\$SEQ),0);\n";
    $scan_code .= "  my \$result = {};\n";
    $scan_code .= "  \$result->{seqid}  = \$SEQID;\n";
    $scan_code .= "  \$result->{acc}  = \"$acc\";\n";
    $scan_code .= "  \$result->{name} = \"$name\";\n";
    $scan_code .= "  \$result->{start} = \$start;\n";
    $scan_code .= "  \$result->{end}   = \$end;\n";
    $scan_code .= "  \$result->{score} = 0;\n";

    if (exists($confirm_hash->{$acc})){    	
      foreach my $con_pat (@{$confirm_hash->{$acc}}) {
        $scan_code .= "  \$result->{score} = 1 if \$match =~ /$con_pat/;\n";
      }
    } 
    $scan_code .= "  push \@RESULTS, \$result;\n}\n";
  }
  close (PAT);

  #print $scan_code;

  return $scan_code;
}


sub _read_confirm_patterns {
  my ($self, $confirm_file, $confirms_hash) = @_;
  
  open (PAT,"<$confirm_file") 
      or throw("Cannot open file of confirms '$confirm_file'");
  while (<PAT>){
    my ($acc,$pattern)=split(/\s+/);
    push @{$confirms_hash->{$acc}}, $pattern;
  }
  close (PAT);
}


sub confirm_file {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_confirm_file} = $val;
  }
  if (exists $self->{_confirm_file}) {
    return $self->{_confirm_file};
  } else {
    return undef;
  }
}

1;
