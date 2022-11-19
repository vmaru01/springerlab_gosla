from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
NCBIWWW.email = "maruyinks@gmail.com"

filename = "file/FAT81153_pass_944e9f34_212.fastq"


# define function to blast records;
# def blastseq_request(line_count):
#     records = list((SeqIO.parse(filename, "fastq")))
#     for i in range(line_count):
#         print(f"ID: {records[i].id} | Sequence: {records[i].seq[:20]}...{records[i].seq[-10:]}")


# blastseq_request(7)
# return sequence id for blast-ing

# def nth_line_seq(line_number):
#     with open(filename, "r") as handle:
#         records = list(SeqIO.parse(handle, "fastq"))
#         result = records[line_number].seq[:20]
#         return result


# print(nth_line_seq(12))
records = list(SeqIO.parse(filename, "fastq"))
result_handle = NCBIWWW.qblast("blastn", "nt", records[2].seq)

with open("test_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
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
