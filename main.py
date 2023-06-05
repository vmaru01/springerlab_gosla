# FIND RELATED SPECIES; IDENTIFY POLYMORPHISMS

import os
import gzip
from Bio import Entrez
from Bio import SearchIO
from Bio import SeqIO  # file parser
from Bio.Blast import NCBIWWW  # Blast to the web and back

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
eval_threshold = float(input("Enter the e-value to filter by (ex. 1.3e-1); Lower e-values will give more "
                             "significant results: "))
gi_list = []
for item in blast_result:
    print(f"BLAST SUMMARY: \n {item}")
    for hsp in item.hsps:
        if hsp.evalue < eval_threshold:
            gi_list.append(hsp.hit_id.split("|")[1])
xml_read.close()

# ENTREZ
for gi_number in range(len(gi_list)):
    search_handle = Entrez.esearch(db="nuccore", term=gi_list[gi_number])
    search_record = Entrez.read(search_handle)
    search_handle.close()

    for search_id in search_record["IdList"]:
        fetch_handle = Entrez.efetch(db="nuccore", id=search_id, rettype='gb', retmode="xml")
        fetch_record = Entrez.parse(fetch_handle)
        search_handle.close()
        for record in fetch_record:
            if "Viridiplantae" in record['GBSeq_taxonomy']:
                print(record["GBSeq_organism"])
