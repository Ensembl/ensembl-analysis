from Sequence import Sequence

if __name__ == '__main__':

  sequence = Sequence(0,20,'-',"seq1","AAAA")
#  t = sequence.get_sequence
  print(sequence.sequence)
  print(sequence.start)
#  sequence.sequence = "TEST";
  sequence.get_sequence('test_fasta.fa')
  print(sequence.start)
  print(sequence.sequence)
#  print(t)
