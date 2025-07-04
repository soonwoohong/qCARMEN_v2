"""

Design primers based on NCBI primer-BLAST

- CommonPrimerDesign: desigining common primers across isoforms that are specific to the target
                      against other off-target genes and species.

"""

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import primer3
import logging
import os

from typing import List, Dict, Tuple, Optional

class CommonPrimerDesign:
    def __init__(self,
                 genbank_dir: str,
                 output_dir: str,
                 organism: str = "Mus musculus"
                 ):
        self.organism = organism
        self.genbank_dir = genbank_dir
        self.output_dir = output_dir
        self.gene_list = [gb_file for gb_file in os.listdir(genbank_dir)]
        print(self.gene_list)

