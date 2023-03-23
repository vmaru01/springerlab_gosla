# FIND RELATED SPECIES; IDENTIFY POLYMORPHISMS

import os
import gzip
from Bio import Entrez
from Bio import SeqIO  # file parser
from Bio.Blast import NCBIWWW  # Blast to the web and back
from Bio.Blast import NCBIXML  # sort and filter blast output

Entrez.email = 'vmaru01@outlook.com'
NCBIWWW.email = 'vmaru01@outlook.com'
fastq_dir = '/Users/springerlab/Desktop/Nanopore Test Data/20220817-sla-luc/no_sample/20220817_1750_MC-113190_FAT81153_9a4aacb0/fastq_pass'

# To BLAST and Back
i = 0
for item in os.listdir(fastq_dir):
    while i == 0:
        try:
            if item.lower().endswith(".gz"):
                fullpath = os.path.join(fastq_dir, item)
                print(fullpath)
                with gzip.open(fullpath, 'rt') as file:
                    records = SeqIO.parse(file, "fastq")
                    for index, line in enumerate(records):
                        if index == 0:
                            print(f'Sequence Length: {len(line.seq)}\nSequence: {line.seq}')
                            blast_handle = NCBIWWW.qblast("blastn", "nt", line.seq)
                    with open("test_blast.xml", "w") as blast_xml:
                        blast_xml.write(blast_handle.read())
                        blast_handle.close()
        except TypeError:
            print("Error: Wrong file type")
            pass
        i += 1

# Filter
with open("test_blast.xml") as xml_read:
    blast_result = NCBIXML.parse(xml_read)
    eval_threshold = float(input("Enter the e-value to filter by (ex. 2.3e-4): "))
    for alignment in blast_result.alignments:
        for hsp in alignment.hsps:
            print(hsp.id)
            print(hsp.score)
            print(hsp.identities)
            if hsp.evalue < eval_threshold:
                print(hsp)   # revise

# entrez_handle = Entrez.esearch(db='Nucleotide')

# Filter
# with open("test_blast.xml", "r") as xml_read:
#     blast_result = SearchIO.parse("test_blast.xml", "blast-xml")
#     eval_threshold = float(input("Enter the e-value to filter by (ex. 2.3e-4): "))
#     for alignment in blast_result:
#         print(alignment)
#         for hsp in alignment.hsps:
#             if hsp.evalue < eval_threshold:
#                 print(f"HSP: {hsp}")
#                 print("---------------------"*4)

# # unzip folder for BLAST
# if zipfile.is_zipfile(fastqzip_dir):  # zipped file check
#     with zipfile.ZipFile(fastqzip_dir, "r") as zipped_file:
#         for filename in zipped_file.namelist():  # list zipped file content
#             print(f".namelist(): {zipped_file.namelist()} \nCurrent filename: {filename}")
#             if filename.lower().endswith((".fastq.gz", "fasta")):  # fastq, fasta file check
#                 print("Checking file format...file matches")
#                 file = filename.read()
#                 # OUT: BLAST - Parse file, index sequences and blast using sequence index
#                 #
#                 #
#             else:
#                 print("Not a file")
# else:
#     print("Not a zipped file")
#
#                 print("---------------------"*4)
#
# # unzip folder for BLAST
# if zipfile.is_zipfile(fastqzip_dir):  # zipped file check
#     with zipfile.ZipFile(fastqzip_dir, "r") as zipped_file:
#         for filename in zipped_file.namelist():  # list zipped file content
#             print(f".namelist(): {zipped_file.namelist()} \nCurrent filename: {filename}")
#             if filename.lower().endswith((".fastq.gz", "fasta")):  # fastq, fasta file check
#                 print("Checking file format...file matches")
#                 file = filename.read()
#                 # OUT: BLAST - Parse file, index sequences and blast using sequence index
#                 #
#                 #
#             else:
#                 print("Not a file")
# else:
#     print("Not a zipped file")
