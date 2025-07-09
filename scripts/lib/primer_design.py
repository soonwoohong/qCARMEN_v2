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
from dataclasses import dataclass

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

@dataclass
class PrimerPair:
    """Store primer pair information from Primer-BLAST"""
    forward_seq: str
    reverse_seq: str
    forward_tm: float
    reverse_tm: float
    forward_gc: float
    reverse_gc: float
    product_size: int
    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int
    penalty: float
    target_coverage: List[str]  # Which isoforms this primer pair covers
    off_targets: int
    specificity_info: Dict

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
                               min_identity: float = 0.95)-> (Dict[str, tuple], Dict[str, int]):
        """
        Find conserved regions across all isoforms of a gene

        conserved_regions = {gene name: [(start_pos, end_pos, identity)]}
        positions here are based on the first sequences in the alignment fasta.

        num_fasta = gene name: fasta_length

        """
        conserved_regions = {}
        num_fasta = {}

        # align the sequences

        for gene in self.gene_list:
            fasta_file = os.path.join(self.fasta_dir, gene+".fasta")
            align_file = os.path.join(self.align_dir, f"aligned_{gene}.fasta")
            # count how many sequences in a fasta file
            fasta_len = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
            num_fasta[gene] = fasta_len
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
            conserved_regions[gene] = self._identify_conserved_regions(
                alignment, min_length, min_identity)


        return conserved_regions, num_fasta

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

    def design_primers(self, conserved_regions):
        """
        This is the main method to use for designing common primers
        :return:
        """
        all_primers_dict = {}
        for gene in self.gene_list:
            # 1. load conserved regions
            temp_conserved = conserved_regions[gene]
            logger.info(f"Found {len(temp_conserved)} conserved regions")
            fasta_file = os.path.join(self.fasta_dir, gene+".fasta")
            # 2: create base sequence
            records = list(SeqIO.parse(fasta_file, "fasta"))
            base_seq = str(records[0].seq)
            base_id = str(records[0].id)

            # 3: design primers in conserved regions
            all_primers =  []
            parent_path = os.path.join(self.genbank_dir,gene)
            genbank_files = [os.path.join(parent_path, genfile) for genfile in os.listdir(parent_path) if ".gb" in genfile]
            for region_start, region_end, identity in temp_conserved:
                primers = self._design_with_primer3(base_seq, region_start, region_end)

                validated_primers = self._validate_primers_across_isoforms(primers, genbank_files)

                all_primers.extend(validated_primers)

            all_primers_dict[gene] = all_primers
        return all_primers_dict

    def _design_with_primer3(self, base_seq: str, region_start: int, region_end: int):
        # Define target region (see below for information in detail)
        # TODO: (070825 - target on exon-exon junction
        # at least one of primers should be on exon-exon junction.

        target_len = min(200, region_end - region_start)
        target_start = region_start + (region_end - region_start - target_len) // 2
        """
                Conserved region: [========================================] (500bp)
                Position:         100                                      600

                target_len = min(200, 500) = 200bp

                Extra space = 500 - 200 = 300bp
                Half of extra = 300 // 2 = 150bp

                target_start = 100 + 150 = 250

                Result:
                [===========|████████████████|===========]
                100         250            450         600
                            └─── 200bp ────┘
                            (Target Region)
        """

        seq_args = {
            'SEQUENCE_ID': 'consensus',
            'SEQUENCE_TEMPLATE': base_seq,
            'SEQUENCE_TARGET': [target_start, target_len],
        }
        # Design primers
        results = primer3.bindings.design_primers(seq_args, self.default_params)

        # Extract primer pairs
        primers = []
        num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)

        for i in range(num_returned):
            primer_data = {
                'forward_seq': results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'reverse_seq': results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'forward_tm': results[f'PRIMER_LEFT_{i}_TM'],
                'reverse_tm': results[f'PRIMER_RIGHT_{i}_TM'],
                'forward_gc': results[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                'reverse_gc': results[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
                'product_size': results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                'penalty': results[f'PRIMER_PAIR_{i}_PENALTY'],
                'forward_pos': results[f'PRIMER_LEFT_{i}'][0],
                'reverse_pos': results[f'PRIMER_RIGHT_{i}'][0],
            }
            primers.append(primer_data)

        return primers


    # no flanking allowed
    # TODO: 070825 - flanking allowed
    def _validate_primers_across_isoforms(self,
                                          primers: List[Dict],
                                          genbank_files: List[str]) -> List[PrimerPair]:
        """Validate that primers work across all isoforms"""
        validated = []

        for primer_data in primers:
            forward_seq = primer_data['forward_seq']
            reverse_seq = primer_data['reverse_seq']

            # Check each isoform
            coverage = []
            works_for_all = True

            for gb_file in genbank_files:
                record = SeqIO.read(gb_file, "genbank")
                seq_str = str(record.seq).upper()

                # Check if both primers are present
                fwd_pos = seq_str.find(forward_seq.upper())
                rev_rc = str(Seq(reverse_seq).reverse_complement()).upper()
                rev_pos = seq_str.find(rev_rc)

                if fwd_pos >= 0 and rev_pos >= 0 and rev_pos > fwd_pos:
                    product_size = rev_pos + len(reverse_seq) - fwd_pos
                    if 50 <= product_size <= 2000:  # Reasonable product size
                        coverage.append(record.id)
                    else:
                        works_for_all = False
                        break
                else:
                    works_for_all = False
                    break

            if works_for_all:
                primer_pair = PrimerPair(
                    forward_seq=forward_seq,
                    reverse_seq=reverse_seq,
                    forward_tm=primer_data['forward_tm'],
                    reverse_tm=primer_data['reverse_tm'],
                    forward_gc=primer_data['forward_gc'],
                    reverse_gc=primer_data['reverse_gc'],
                    product_size=primer_data['product_size'],
                    forward_start=primer_data['forward_pos'],
                    forward_end=primer_data['forward_pos'] + len(forward_seq),
                    reverse_start=primer_data['reverse_pos'],
                    reverse_end=primer_data['reverse_pos'] + len(reverse_seq),
                    penalty=primer_data['penalty'],
                    target_coverage=coverage,
                    off_targets=0,  # Will be updated by specificity check
                    specificity_info={}
                )
                validated.append(primer_pair)

        logger.info(f"Validated {len(validated)} primer pairs across all isoforms")
        return validated






