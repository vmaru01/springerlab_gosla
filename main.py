from Bio import SeqIO
from Bio.Blast import NCBIWWW
NCBIWWW.email = "maruyinks@gmail.com"

filename = "file/FAT81153_pass_944e9f34_212.fastq"
line_count = 5
records = list((SeqIO.parse(filename, "fastq")))

for i in range(line_count):
    print(records[i].id, records[i].seq[:10], records[i].seq[-10:])

# define function to blast records;
