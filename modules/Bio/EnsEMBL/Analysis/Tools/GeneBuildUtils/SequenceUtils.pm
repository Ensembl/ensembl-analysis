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
             create_file_name
             write_seqfile
            );



=head2 create_file_name

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



sub create_file_name{
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



=head2 write_seq_file

  Arg [1]   : Bio::Seq
  Arg [2]   : string, filename
  Function  : This uses Bio::SeqIO to dump a sequence to a fasta file
  Returntype: string, filename
  Exceptions: throw if failed to write sequence
  Example   : 

=cut

sub write_seqfile{
  my ($seq, $filename) = @_;
  throw("Need a Bio::Seq object not a ".$seq)
    if(!$seq || !$seq->isa('Bio::Seq'));
  $filename = create_file_name('seq', 'fa', '/tmp') 
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
