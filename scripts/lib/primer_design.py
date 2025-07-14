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
class ExonJunction:
    """Information about an exon-exon junction"""
    junction_seq: str  # Sequence around junction
    present_in_isoforms: List[str]
    junction_pos_in_seq: int  # Position of junction in junction_seq

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
    primer_class: int # 1: primers on exon-exon junction; 2: amplicon spanning exon-exon junction 3. primers on conserved region
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
            'PRIMER_MAX_TM_DIFF': 2,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_PRODUCT_MIN': 70,
            'PRIMER_PRODUCT_MAX': 1000,
            'DESIRED_AMPLICON_LENGTH': 300,
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
            'PRIMER_EXON_5_MIN_MATCH': 7,
            'PRIMER_EXON_3_MIN_MATCH': 4,
            'PRIMER_EXON_3_MAX_MATCH': 8
        }
    def design_primers(self, genbank_files: List[str]) -> List[PrimerPair]:
        """
        Main method to design junction-targeting primers

        Returns primers in priority order:
        1. ON junction: primer_class = 1
        2. SPANNING junction: primer_class = 2
        3. Conserved region : primer_class = 3
            3A: no primer available on junction and spanning junction
            3B: no junction found because only one exon exists

        save all the possible primers for designing crRNAs later.
        top n primers are saved in the primer directory. (default: {output_dir}/primers)
        """

        logger.info(f"Designing junction primers for {len(genbank_files)} files")

        #extract sequences and exon structures
        isoforms = self._load_isoforms(genbank_files)

        # find all exon junctions
        #junctions = self._find_all_junctions(isoforms)
    def _load_isoforms(self, genbank_files: List[str]) -> Dict[str, Dict]:
        """
        Load isoforms from genbank files and extract exon information.
        """
        isoforms = {}
        for gb_file in genbank_files:
            record = SeqIO.read(gb_file, "genbank")

            # extract exons
            exons = []


            for feature in record.features:
                if feature.type == "exon":
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    exons.append((start, end))

            # if no exons annotated (ex. Ifi44l, NM_031367.1), consider the entire gene as one exon.
            if not exons:
                exons.append((0, len(record.seq) - 1))

            # sort exons by position (x[0] = start position)
            exons.sort(key=lambda x: x[0])

            # create mRNA sequence (this is normally not necessary, because the entire sequence (here, record.seq) represents mRNA)
            mRNA_seq = ""
            for start, end in exons:
                mRNA_seq += str(record.seq[start:end])

            isoforms[record.id] = {
                'mRNA_seq': mRNA_seq,
                'exons': exons,
            }

        return isoforms

    def _find_all_junctions(self, isoforms: Dict[str, Dict]) -> List[ExonJunction]:
        """Find all unique exon-exon junctions across isoforms"""
        junctions = {}
        junction_count = 0
        duplicate_check = []
        # isoform_id = isoform id
        # isoform data = mRNA_seq / exons
        for isoform_id, isoform_data in isoforms.items():
            exons = isoform_data['exons']
            mRNA_seq = isoform_data['mRNA_seq']

            # current position in mRNA sequence
            mRNA_pos = 0

            for i in range(len(exons) - 1):
                exon_length = exons[i][1] - exons[i][0]
                junction_pos = mRNA_pos + exon_length

                # extract sequence around junction (30 bases each side)
                # max and min to avoid junction sequences beyond mRNA sequence
                start = max(0, junction_pos - 30)
                end = min(len(mRNA_seq), junction_pos + 30)
                junction_seq = mRNA_seq[start:end]
                junction_key = str(junction_count)
                if junction_seq not in duplicate_check:
                    junctions[junction_key] = ExonJunction(
                        junction_seq=junction_seq,
                        present_in_isoforms=[isoform_id],
                        junction_pos_in_seq=junction_pos - start
                    )
                    junction_count += 1
                    duplicate_check.append(junction_seq)
                else:
                    update_key = duplicate_check.index(junction_seq)
                    junctions[str(update_key)].present_in_isoforms.append(isoform_id)

                mRNA_pos += exon_length
        return list(junctions.values())

    def _design_on_junction_primers(self, gene: str, junctions: List[ExonJunction],
                                    isoforms: Dict[str, Dict]) -> List[PrimerPair]:
        """Design primers on exon-exon junctions"""
        primers =[]
        junction_count = 0
        valid_junctions = []

        min_len = self.default_params['PRIMER_MIN_SIZE']
        max_len = self.default_params['PRIMER_MAX_SIZE']
        exon_5_min_match = self.default_params['PRIMER_EXON_5_MIN_MATCH']


        for junction in junctions:
            # only design if junction is present in all isoforms
            if len(junction.present_in_isoforms) != len(isoforms):
                continue

            junction_count += 1
            valid_junctions.append(junction)
        logger.info(
            f'{gene} has {len(isoforms)} isoform(s) and common {junction_count} junction(s) across all isoform(s).')

        # design primers that overlap the junction
        # at least one primer (fwd or rev) must cross the junction

        # Case 1: num_junctions >= 1
        # Try different primer positions around junction
        for val_junction in valid_junctions:
            # vary position relative to junction
            # check the forward primer overlapping >{PRIMER_EXON_5_MIN_MATCH} bases on the exon
            # EX
            # [-----exon1----]--[-----exon2-----]
            #        {-----fwd-----}
            #        <overlap>
            conserved_down_seq = self._find_conserved_regions_around_junction(val_junction, isoforms)
            for overlap in range(exon_5_min_match, exon_5_min_match+10):
                fw_start = val_junction.junction_pos_in_seq - overlap



                for primer_len in range(-min_len, max_len+1):
                    fw_end = fw_start + primer_len
                    if fw_start >=0 and fw_end < len(val_junction.junction_seq):
                        fw_seq = val_junction.junction_seq[fw_start:fw_end]

                        # find reverse primer targeting conserved downstream region.

    def _find_conserved_regions_around_junction(self, junction: ExonJunction,
                                                isoforms: Dict[str, Dict]):
        """Find conserved regions around a junction in isoforms"""
        junction_seq = junction.junction_seq
        junction_pos_in_seq = junction.junction_pos_in_seq
        conserved_len = round(self.default_params['DESIRED_AMPLICON_LENGTH'] * 1.5)

        conserved_regions = {}

        for isoform_id, isoform_data in isoforms.items():
            mRNA_seq = isoform_data['mRNA_seq']
            junction_pos_in_mRNA = mRNA_seq.find(junction_seq)+junction_pos_in_seq
            # [exon1] --- junction --- [exon2]
            #                 ------------> downstream
            # <--------------- upstream
            # checking length
            downstream_seq = mRNA_seq[junction_pos_in_mRNA:min(junction_pos_in_mRNA+conserved_len,len(mRNA_seq))]
            upstream_seq = mRNA_seq[max(0, junction_pos_in_mRNA-conserved_len):junction_pos_in_mRNA]
            conserved_regions[isoform_id] = {"isoform_id": isoform_id,
                                             "junction_seq": junction_seq,
                                             "junction_pos_in_seq": junction_pos_in_seq,
                                            "junction_pos_in_mRNA": junction_pos_in_mRNA,
                                             "upstream_seq": upstream_seq,
                                             "downstream_seq": downstream_seq}


        # validate conserved regions across isoforms
        self._validate_conserved_regions_around_junction(conserved_regions, 0.95, self.default_params['PRIMER_PRODUCT_MIN'])

    def _validate_conserved_regions_around_junction(self, conserved_regions: Dict,
                                                    min_identity: float,
                                                    min_length: int):
        # validate conserved regions around junction across isoforms
        valid_conserved_regions = {}
        downstream_seq = []
        upstream_seq = []
        for isoform_id, isoform_data in conserved_regions.items():
            downstream_seq.append(isoform_data["downstream_seq"])
            upstream_seq.append(isoform_data["upstream_seq"])
        downstream_min_length = min(map(len, downstream_seq))
        upstream_min_length = min(map(len, upstream_seq))
        alignment_down =MultipleSeqAlignment([SeqRecord(Seq(down_seq[:downstream_min_length])) for down_seq in downstream_seq])


        # Scan through donwstream sequnces to find valid conserved downstream seq
        i = 0
        while i < downstream_min_length:
            # Check if position is conserved
            column = alignment_down[:, i]
            if self._calculate_identity(column) >= min_identity:
                # Start of conserved region
                start = i
                while i < downstream_min_length and self._calculate_identity(alignment_down[:, i]) >= min_identity:
                    i += 1
                end = i

                # Check if region is long enough
                if end - start >= min_length:
                    # Calculate average identity
                    total_identity = sum(
                        self._calculate_identity(alignment_down[:, j])
                        for j in range(start, end)
                    )
                    avg_identity = total_identity / (end - start)

                    # Convert to  coordinates
                    valid_conserved_regions['downstream'] = (downstream_seq[0], start, end, avg_identity)

            else:
                i += 1
        # Scan through upstrea sequnces to find valid conserved downstream seq
        # upstream sequence is a little tricky because they need to be checked backward.
        # [exon1] --- junction --- [exon2]
        # <--------------- upstream
        #                ^
        #      check the identitiy from here
        alignment_up = MultipleSeqAlignment([SeqRecord(Seq(up_seq[-1:-1-upstream_min_length:-1])) for up_seq in upstream_seq])
        i = 0
        while i < upstream_min_length:
            # Check if position is conserved
            column = alignment_up[:, i]
            if self._calculate_identity(column) >= min_identity:
                # Start of conserved region
                start = i
                while i < upstream_min_length and self._calculate_identity(alignment_up[:, i]) >= min_identity:
                    i += 1
                end = i

                # Check if region is long enough
                if end - start >= min_length:
                    # Calculate average identity
                    total_identity = sum(
                        self._calculate_identity(alignment_up[:, j])
                        for j in range(start, end)
                    )
                    avg_identity = total_identity / (end - start)

                    # Convert to coordinates
                    valid_conserved_regions['upstream'] = (upstream_seq[0], -1-start, -end, avg_identity)
                    # actual sequence = upstream_seq[-end:-1-start]
            else:
                i += 1

        return valid_conserved_regions

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

    def design_primers_old(self, conserved_regions):
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






