#!/usr/bin/env python3

import subprocess
import os

# Set the path to the folder containing FASTA files
folder_path = "/Users/Ryan/Desktop/Contigs/"

# Get a list of all the FASTA files in the folder
fasta_files = [file for file in os.listdir(folder_path) if file.endswith(".fasta")]

# Iterate over each FASTA file
for fasta_file in fasta_files:
    # Create the full path to the input FASTA file
    fasta_file_path = os.path.join(folder_path, fasta_file)

    # Set the output file path
    output_file = fasta_file_path + '.blast.xml'

    # Construct the command to run blastn
    cmd = f'blastn -query "{fasta_file_path}" -db nt -out "{output_file}" -outfmt 5'

    # Execute the blastn command
    subprocess.run(cmd, shell=True)

    print("BLAST search completed for:", fasta_file_path)
    print("Results saved to:", output_file)