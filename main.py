# FIND RELATED SPECIES; IDENTIFY POLYMORPHISMS

import os
import gzip
from Bio import Entrez
from Bio import SearchIO
from Bio import SeqIO  # file parser
from Bio.Blast import NCBIWWW  # Blast to the web and back
from Bio.Blast import NCBIXML  # sort and filter blast output

Entrez.email = 'vmaru01@outlook.com'
NCBIWWW.email = 'vmaru01@outlook.com'
fastq_dir = "C:/Users/steam/Downloads/Nanopore Data"

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
xml_read = open("test_blast.xml")
blast_result = SearchIO.parse(xml_read, "blast-xml")
eval_threshold = float(input("Enter the e-value to filter by (ex. 2.3e-4): "))
for item in blast_result:
    print(item)
    for hsp in item.hsps:
        if hsp.evalue < eval_threshold:
            print(hsp)
            print(f"BITSCORE: {hsp.bitscore}")
            print(f"HIT ID: {hsp.hit_id}")
            print(f"EVALUE: {hsp.evalue}")
            print("**********" * 4)
          # revise
xml_read.close()

#ENTREZ
handle = Entrez.esearch(db="nucleotide", term="txid63221")
record = Entrez.read(handle)
for result in record:
    print(record)
