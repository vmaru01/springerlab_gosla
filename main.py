# Dependencies

from Bio import SeqIO  # file parser
from Bio.Blast import NCBIWWW  # Blast to the web
from Bio.Blast import NCBIXML  # and back
NCBIWWW.email = "maruyinks@gmail.com"

filename = "file/FAT81153_pass_944e9f34_212.fastq"

# Blast sequence and write result to a xml file
# OUT: Parse file, index sequences and blast using index
records = list(SeqIO.parse(filename, "fastq"))
result_handle = NCBIWWW.qblast("blastn", "nt", records[2].seq)  # get sequence read using index location

# IN: save result to xml format
with open("test_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())  # use str result as write object
result_handle.close()

result_handle = open("test_blast.xml")
blast_record = NCBIXML.read(result_handle)

e_val_thresh = 0.04
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < e_val_thresh:
            print("**Alignment**")
            print(f"Sequence: {alignment.title}")
            print(f"Length: {alignment.length}")
            print(f"e value: {hsp.expect}")
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
        else:
            print("BLAST Failed :(")
