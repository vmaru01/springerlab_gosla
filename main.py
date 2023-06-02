# FIND RELATED SPECIES; IDENTIFY POLYMORPHISMS

import os
import gzip
from Bio import Entrez
from Bio import SearchIO
from Bio import SeqIO  # file parser
from Bio.Blast import NCBIWWW  # Blast to the web and back
# from Bio.Blast import NCBIXML  # sort and filter blast output

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
            # taxids = item.accession.split('.')[0]
            # print(f"TAXID: {taxids}")
            print(hsp)
            print(f"BITSCORE: {hsp.bitscore}")
            print(f"HIT ID: {hsp.hit_id}")
            print(f"EVALUE: {hsp.evalue}")
            print("**********" * 6)
        # revise
xml_read.close()

# ENTREZ
search_handle = Entrez.esearch(db="taxonomy", term="63221")
search_record = Entrez.read(search_handle)
search_handle.close()
print(search_record)

link_handle = Entrez.elink(db="taxonomy", id=search_record['IdList'][0])
link_record = Entrez.read(link_handle)
search_handle.close()
print(link_record)

post_handle = Entrez.esummary(db="taxonomy", id=link_record[0]['LinkSetDb'][0]['Link'][0])
post_record = Entrez.parse(post_handle)
search_handle.close()
print(post_record)
