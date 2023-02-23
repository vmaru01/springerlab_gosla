# FIND RELATED SPECIES; IDENTIFY POLYMORPHISMS
# Parse fasta/fastq file
# Blast output to NCBI
# Write result to xml file
# Sort alignment scores; filter at e = 1e-7 (approx. 80% bit score)

import os
import zipfile
from Bio import SeqIO  # file parser
from Bio.Blast import NCBIWWW  # Blast to the web and back
from Bio import SearchIO  # sort and filter blast output

NCBIWWW.email = "maruyinks@gmail.com"
fastq_dir = "file/FAT81153_pass_944e9f34_212.fastq"
fastqzip_dir = "file/test_fastq.gz.zip"

# unzip folder for BLAST
if zipfile.is_zipfile(fastqzip_dir):  # verify that file is zipped
    with zipfile.ZipFile(fastqzip_dir, "r") as zipped_file:
        for filename in zipped_file.namelist():  # list zipped file content
            print(f"Printing content in zipped_file...{filename}")
            if filename.endswith(".fastq.gz"):  # verify that you're working with a file
                print("Checking file format...Correct match!")
                file = zipped_file.read(filename)
                print("")
                # OUT: BLAST - Parse file, index sequences and blast using sequence index
                records = list(SeqIO.parse(file, "fastq"))
                blast_handle = NCBIWWW.qblast("blastn", "nt", records[5].seq)
            else:
                print("Not a file")
else:
    print("Not a zipped file")

# IN: XML - write BLAST result to xml file
with open("test_blast.xml", "w") as blast_xml:
    blast_xml.write(blast_handle.read())  # use str result as write object
blast_handle.close()

# Filter
with open("test_blast.xml", "r") as xml_read:
    blast_result = SearchIO.parse("test_blast.xml", "blast-xml")
    eval_threshold = float(input("Enter the e-value to filter by (ex. 2.3e-4): "))
    for alignment in blast_result:
        print(alignment)
        for hsp in alignment.hsps:
            if hsp.evalue < eval_threshold:
                print(f"HSP: {hsp}")
                print("---------------------"*4)
