import os
import glob
import gzip
from Bio import SeqIO

# Change the working directory to the folder containing the FASTQ files
os.chdir('/Users/Ryan/Desktop/Raw_FastQ')

# Get the list of all .fastq.gz files in the folder
fastq_files = glob.glob("*.fastq.gz")

# Iterate over the files and convert them to FASTA
combined_sequences = []
for fastq_file in fastq_files:
    # Extract the file name without the extension
    file_name = fastq_file.replace(".fastq.gz", "")
    
    # Read the fastq.gz file and convert sequences to FASTA format
    with gzip.open(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            combined_sequences.append(record)

# Write the combined sequences to a single FASTA file
output_file = "/Users/Ryan/Desktop/FastA_Merged/merged.fasta"
with open(output_file, "w") as handle:
    SeqIO.write(combined_sequences, handle, "fasta")
