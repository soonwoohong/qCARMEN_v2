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
import subprocess
from typing import List, Dict, Tuple, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class CommonPrimerDesign:
    def __init__(self,
                 output_dir: str,
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
        self.output_dir = output_dir
        self.genbank_dir = os.path.join(output_dir, "genbank")
        self.fasta_dir = os.path.join(output_dir, "fasta")
        self.primer_dir = os.path.join(output_dir, "primer")
        self.align_dir = os.path.join(output_dir, "aligned_fasta")
        os.makedirs(self.primer_dir, exist_ok=True)
        os.makedirs(self.align_dir, exist_ok=True)
        self.gene_list = [gb_file for gb_file in os.listdir(self.genbank_dir) if os.path.isdir(os.path.join(self.genbank_dir, gb_file))]

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
                               min_length: int = 100,
                               min_identity: float = 0.95)-> Dict[str, tuple]:
        """
        Find conserved regions across all isoforms of a gene

        """
        conserved_regions = {}

        # align the sequences

        for gene in self.gene_list:
            fasta_file = os.path.join(self.fasta_dir, gene+".fasta")
            align_file = os.path.join(self.align_dir, f"aligned_{gene}.fasta")
            # count how many sequences in a fasta file
            fasta_len = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
            # if a fasta file has one sequence, it doesn't require alignment.
            if fasta_len == 1:
                with open(align_file, "w") as f:
                    SeqIO.write(SeqIO.read(fasta_file, "fasta"), f, "fasta")

            else:
                # alignment using mafft
                subprocess.run(['mafft', '--auto', '--quiet', fasta_file],
                               stdout=open(align_file, 'w'))
            logger.info(f"Saved aligned fasta files for {gene} ({fasta_len} sequence(s) aligned)")
            # Find conserved regions
            alignment = AlignIO.read(align_file, "fasta")
            conserved_regions.update({gene: self._identify_conserved_regions(
                alignment, min_length, min_identity
            )})

        # conserved_regions = {gene name: [(start_pos, end_pos, identity)]}
        # positions here are based on the first sequences in the alignment fasta.
        return conserved_regions

    def _identify_conserved_regions(self, alignment: MultipleSeqAlignment,
                                    min_length: int,
                                    min_identity: float) -> List[Tuple[int, int, float]]:
        """Identify conserved regions in the alignment"""
        conserved_regions = []
        seq_length = alignment.get_alignment_length()

        # Scan through alignment to find conserved stretches
        i = 0
        while i < seq_length:
            # Check if position is conserved
            column = alignment[:, i]
            if self._calculate_identity(column) >= min_identity:
                # Start of conserved region
                start = i
                while i < seq_length and self._calculate_identity(alignment[:, i]) >= min_identity:
                    i += 1
                end = i

                # Check if region is long enough
                if end - start >= min_length:
                    # Calculate average identity
                    total_identity = sum(
                        self._calculate_identity(alignment[:, j])
                        for j in range(start, end)
                    )
                    avg_identity = total_identity / (end - start)

                    # Convert to ungapped coordinates
                    ungapped_start = self._gapped_to_ungapped(alignment[0], start)
                    ungapped_end = self._gapped_to_ungapped(alignment[0], end)

                    conserved_regions.append((ungapped_start, ungapped_end, avg_identity))
            else:
                i += 1

        return conserved_regions

    def _calculate_identity(self, column: str) -> float:
        """Calculate identity for an alignment column"""
        # Remove gaps
        non_gap_chars = [c for c in column if c != '-']
        if not non_gap_chars:
            return 0.0

        # Count most common character
        char_counts = {}
        for char in non_gap_chars:
            char_counts[char] = char_counts.get(char, 0) + 1

        max_count = max(char_counts.values())
        return max_count / len(non_gap_chars)

    def _gapped_to_ungapped(self, seq_record, gapped_pos: int) -> int:
        """Convert gapped alignment position to ungapped sequence position"""
        ungapped_pos = 0
        for i in range(gapped_pos):
            if seq_record.seq[i] != '-':
                ungapped_pos += 1
        return ungapped_pos
