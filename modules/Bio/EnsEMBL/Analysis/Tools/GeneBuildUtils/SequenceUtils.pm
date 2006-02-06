package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::SequenceUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id);
use vars qw (@ISA  @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw(
             create_filename
             write_seq_file
            );



=head2 create_filename

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string, stem of filename
  Arg [3]   : string, extension of filename
  Arg [4]   : directory file should live in
  Function  : create a filename containing the PID and a random number
  with the specified directory, stem and extension
  Returntype: string, filename
  Exceptions: throw if directory specifed doesnt exist
  Example   : my $queryfile = $self->create_filename('seq', 'fa');

=cut



sub create_filename{
  my ($stem, $ext, $dir) = @_;
  if(!$dir){
    $dir = '/tmp';
  }
  $stem = '' if(!$stem);
  $ext = '' if(!$ext);
  throw($dir." doesn't exist SequenceUtils::create_filename") 
    unless(-d $dir);
  my $num = int(rand(100000));
  my $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  while(-e $file){
    $num = int(rand(100000));
    $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  }
  return $file;
}


sub write_seq_file{
  my ($seq, $filename) = @_;
  throw("Need a Bio::Seq object not a ".$seq)
    if(!$seq || !$seq->isa('Bio::Seq'));
  $filename = create_filename('seq', 'fa', '/tmp') 
    if(!$filename);
  my $seqout = Bio::SeqIO->new(
                               -file => ">".$filename,
                               -format => 'Fasta',
                              );
  eval{
    $seqout->write_seq($seq);
  };
  if($@){
    throw("FAILED to write $seq to $filename Runnable:write_seq_file $@");
  }
  return $filename;
}

1;
