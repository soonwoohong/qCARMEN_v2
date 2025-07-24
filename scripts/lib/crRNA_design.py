import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# first get the primer


class crRNA_Design:
    def __init__(self,
                 output_dir: str
                 ):
        self.output_dir = output_dir
        self.crRNA_dir = os.path.join(output_dir, "crRNA")
        os.makedirs(self.crRNA_dir, exist_ok=True)

    def read_valid_primers(self,
                           gene_name,
                           valid_primers):
        amplicons = valid_primers['aplicon_seq']

        # Save amplicons as a temp FASTA file for running badgers
        # these amplicons will be the input for badgers
        temp_fasta = os.path.join(self.crRNA_dir, "temp.fasta")

        # Individual FASTA files
        fasta_records = []
        for gb_file in os.listdir(gene_dir):
            fasta_records.append(SeqIO.read(os.path.join(gene_dir, gb_file), "genbank"))

        filename = os.path.join(fasta_dir, f"{gene_name}.fasta")
        with open(filename, "w") as f:
            SeqIO.write(fasta_records, f, "fasta")

