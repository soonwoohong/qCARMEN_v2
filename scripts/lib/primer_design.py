"""

Design primers based on NCBI primer-BLAST

- CommonPrimerDesign: designing common primers across isoforms that are specific to the target
                      against other off-target genes and species.

"""

import numpy as np
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import primer3
import logging
import os

from typing import List, Dict, Tuple, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class CommonPrimerDesign:
    def __init__(self,
                 genbank_dir: str,
                 fasta_dir: str,
                 primer_dir: str,
                 organism: str = "Mus musculus",
                 num_primers: int = 10
                 ):
        """
        Arg:
            organism: Target organism for specificity checking
            primer_dir: primer directory for primers
            num_primers: Number of primer pairs to return

        """
        self.organism = organism
        self.genbank_dir = genbank_dir
        self.fasta_dir = fasta_dir
        self.primer_dir = primer_dir
        os.makedirs(self.primer_dir, exist_ok=True)
        self.gene_list = [gb_file for gb_file in os.listdir(genbank_dir) if os.path.isdir(os.path.join(genbank_dir, gb_file))]

        self.default_params = {
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_PRODUCT_MIN': 70,
            'PRIMER_PRODUCT_MAX': 1000,
            'PRIMER_NUM_RETURN': num_primers,
            'PRIMER_MAX_POLY_X': 5,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_SALT_DIVALENT': 1.5,
            'PRIMER_DNTP_CONC': 0.2,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_SELF_ANY': 8.0,
            'PRIMER_MAX_SELF_END': 3.0,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,
            'PRIMER_PAIR_MAX_COMPL_END': 3.0,
        }
        def find_conserved_regions(self,
                                   fasta_dir: str,
                                   min_length: int = 100,
                                   min_identity: float = 0.95):
            """
            Find conserved regions across all isoforms

            """
            sequences = []
            for gb_file in self.gene_list:
                gb_path = os.path.join(self.genbank_dir, gb_file)
                for record in SeqIO.parse(gb_path, "genbank"):
                    sequences.append(record.seq)