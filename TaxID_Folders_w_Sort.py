import os
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import shutil
from io import StringIO

fasta_dir = "/Users/Ryan/Desktop/FastATest2"
output_dir = "/Users/Ryan/Desktop/Test_Results"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Set your email address for Entrez
Entrez.email = "srcheverie@upei.ca"

for filename in os.listdir(fasta_dir):
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        fasta_file = os.path.join(fasta_dir, filename)
        with open(fasta_file) as file:
            sequence_query = file.read()

        blast_program = "blastn"
        database = "nt"
        word_size = 28
        expect = 0.001
        gapcosts = "5 2"

        result_handle = NCBIWWW.qblast(
            blast_program,
            database,
            sequence_query,
            expect=expect,
            word_size=word_size,
            gapcosts=gapcosts
        )
        
        # Save the result content
        result_content = result_handle.read()

        output_file = os.path.join(output_dir, filename.replace(".fasta", ".xml"))
        with open(output_file, "w") as out:
            out.write(result_content)

        # Create a new file-like object for parsing
        result_content_handle = StringIO(result_content)
        blast_records = NCBIXML.parse(result_content_handle)

        try:
            # Process and analyze the top hit
            blast_record = next(blast_records)
            alignment = blast_record.alignments[0]
            print("Alignment description:", alignment.title)
            # Extract the accession number
            accession = alignment.accession

            # Retrieve tax ID using Entrez efetch
            handle = Entrez.efetch(db="nuccore", id=accession, retmode="xml")
            record = Entrez.read(handle)
            tax_id = record[0]['GBSeq_taxonomy']
            print("Tax ID:", tax_id)

            # Create a folder based on the tax ID if it doesn't exist
            tax_folder = os.path.join(output_dir, tax_id)
            if not os.path.exists(tax_folder):
                os.makedirs(tax_folder)

            # Move the FASTA file to the tax folder
            shutil.move(fasta_file, os.path.join(tax_folder, filename))

            hsp = alignment.hsps[0]
            print("Score:", hsp.score)
            print("E-value:", hsp.expect)
            print("Query sequence:", hsp.query)
            print("Alignment sequence:", hsp.match)
            print("Subject sequence:", hsp.sbjct)
            print()

        except IndexError:
            print("No hits found for:", fasta_file)
            continue