from Exon import Exon

if __name__ == '__main__':

  exon = Exon(1,100,1,"Y")
  print(exon.start)
  print(exon.end)
  print(exon.strand)
  print(exon.location_name)

