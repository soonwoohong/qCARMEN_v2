"""

Specificity check based on BLAST

Original

Run the following code in command line

- MAC / linux
* install blast
conda install bioconda::blast
* install database (ex. mus musculus or mouse genome)
update_blastdb.pl --decompress mouse_genome


* check this cookbook in more detail:
https://www.ncbi.nlm.nih.gov/books/NBK279690/pdf/Bookshelf_NBK279690.pdf

- SpecificityChecking: narrowing down primers that are specific to the target
                       against species.


- NOTE that the direction of sequences is from 5' to 3' in most cases.

"""
import os
import subprocess
import sys
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from typing import List, Dict, Tuple, Optional
import logging

# from local
from .primer_design import PrimerPair

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PrimerBlast:
    def __init__(self, output_dir: str, organism: str = "Mus musculus"):
        """
               Initialize PrimerBLAST

               Args:

               """
        if organism == "Mus musculus":
            self.filter = 10090 # organism filter based on NCBI taxonomy id
        self.organism = organism
        self.output_dir = output_dir
        self.genbank_dir = os.path.join(output_dir, "genbank")
        self.fasta_dir = os.path.join(output_dir, "fasta")
        self.primer_dir = os.path.join(output_dir, "primer")
        self.align_dir = os.path.join(output_dir, "aligned_fasta")
        self.specificity_dir = os.path.join(output_dir, "specificity")
        os.makedirs(self.specificity_dir, exist_ok=True)
    def run_blast_locally(self, gene, primers):

        # unique_names = [0_fwd, 1_fwd ... 0_rev, 1_rev ....]
        unique_primers, unique_names = self._load_primer(primers)

        # make a temp fsata file for blast
        temp_file = os.path.join(self.specificity_dir, "temp.fasta")
        temp_records = []
        for unique_idx in range(len(unique_primers)):
            temp_records.append(SeqRecord(Seq(unique_primers[unique_idx]), id=gene+"_"+unique_names[unique_idx]))

        with open(temp_file, "w") as f:
            SeqIO.write(temp_records, f, "fasta")

        output_file = os.path.join(self.specificity_dir, f"{gene}_specificity.xml")
        cmd = self._make_blast_cmd(temp_file, output_file)
        logger.info(f"Checking primers specificity for {gene}...")
        subprocess.run(cmd, shell=True)

        return

    def _load_primer(self, primers):
        """
        args:
        primers: Dataframe of PrimerPair
                                gene: str
                                forward_seq: str
                                reverse_seq: str
                                forward_tm: float
                                reverse_tm: float
                                forward_gc: float
                                reverse_gc: float
                                amplicon_seq: str
                                amplicon_len: int
                                primer_class: int
        :return:
        """
        # load primers
        fwd_primers = primers['forward_seq'].to_list()
        fwd_name = [f"{i}_fwd" for i in range(len(fwd_primers))]
        rev_primers = primers['reverse_seq'].to_list()
        rev_name = [f"{i}_rev" for i in range(len(rev_primers))]
        total_primers = fwd_primers + rev_primers
        total_names = fwd_name + rev_name
        # merge them together and remove duplicates while preserving order
        unique_primers = []
        unique_names = []
        seen = set()

        for primer, name in zip(total_primers, total_names):
            if primer not in seen:
                seen.add(primer)
                unique_primers.append(primer)
                unique_names.append(name)

        return unique_primers, unique_names

    def _make_blast_cmd(self, temp_file, output_file):

        cmd = [
            "blastn",
            "-db", self.db,
            "-query", temp_file,
            '-outfmt', '5',  # XML output
            '-out', output_file,
            '-word_size', '7',
            '-evalue', '1000',
            '-num_alignments', '100',
            '-task', 'blastn-short'  # Optimized for short queries
        ]

        return cmd