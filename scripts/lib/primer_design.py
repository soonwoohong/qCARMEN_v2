"""

Design primers based on NCBI primer

- CommonPrimerDesign: designing common primers across isoforms that are specific to the target


- NOTE that the direction of sequences is from 5' to 3' in most cases.

"""

import numpy as np
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.SeqUtils import gc_fraction
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass, asdict

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
    gene: str
    forward_seq: str
    reverse_seq: str
    forward_tm: float
    reverse_tm: float
    forward_gc: float
    reverse_gc: float
    amplicon_seq: str
    amplicon_len: int
    primer_class: int # 1: fwd primer on junction; 2: rev primer on junction 3. primers on conserved region


class CommonPrimerDesign:
    def __init__(self,
                 output_dir: str,
                 num_primers: int = 10
                 ):
        """
        Arg:
            primer_dir: primer directory for primers
            num_primers: Number of primer pairs to return

        """
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
            'PRIMER_MIN_TM': 58.0,
            'PRIMER_MAX_TM': 64.0,
            'PRIMER_OPT_TM': 61,
            'PRIMER_MAX_TM_DIFF': 1.5,
            'PRIMER_MIN_GC': 0.4,
            'PRIMER_MAX_GC': 0.6,
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
    def design_primers(self, gene_name: str, genbank_files: List[str]):
        """
        Main method to design junction-targeting primers

        Returns primers in priority order:
        1. ON junction: primer_class = 1
        2. SPANNING junction: primer_class = 2
        3. Conserved region : primer_class = 3

        save all the possible primers for designing crRNAs later.
        top n primers are saved in the primer directory. (default: {output_dir}/primers)
        """

        logger.info(f"Designing primers for {gene_name}")
        #extract sequences and exon structures
        isoforms = self._load_isoforms(genbank_files)
        junctions = self._find_all_junctions(isoforms)
        # primersA: fwd primer on junction
        # primersB: rev primer on junction
        primers_A, primers_B = self._design_on_junction_primers(gene_name, junctions,isoforms)
        # primersC: primers on conserved regions

        primers_C = self._design_on_conserved_primers(gene_name)

        primers = pd.DataFrame([asdict(obj) for lst in [primers_A, primers_B, primers_C] for obj in lst])
        primers = primers.drop_duplicates(subset=['amplicon_seq'], ignore_index=True)

        return primers


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
                'mRNA_seq': mRNA_seq.upper(),
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
                                    isoforms: Dict[str, Dict]) -> Tuple[List[PrimerPair], List[PrimerPair]]:
        """Design primers on exon-exon junctions"""
        primers_A = [] # forward primer on junction
        primers_B = [] # reverse primer on junction
        junction_count = 0
        valid_junctions = []

        min_len = self.default_params['PRIMER_MIN_SIZE']
        max_len = self.default_params['PRIMER_MAX_SIZE']
        min_GC = self.default_params['PRIMER_MIN_GC']
        max_GC = self.default_params['PRIMER_MAX_GC']
        exon_5_min_match = self.default_params['PRIMER_EXON_5_MIN_MATCH']
        exon_3_min_match = self.default_params['PRIMER_EXON_3_MIN_MATCH']
        exon_3_max_match = self.default_params['PRIMER_EXON_3_MAX_MATCH']

        Tm_min = self.default_params['PRIMER_MIN_TM']
        Tm_max = self.default_params['PRIMER_MAX_TM']
        Tm_max_diff = self.default_params['PRIMER_MAX_TM_DIFF']

        desired_amplicon_length = self.default_params['DESIRED_AMPLICON_LENGTH']


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

        # Case 1: num_junctions >= 1, fwd primer on junction / rev primer on conserved regions
        # Try different primer positions around junction

        for val_junction in valid_junctions:
            # vary position relative to junction
            # check the forward primer overlapping >{PRIMER_EXON_5_MIN_MATCH} bases on the exon
            # EX
            # [-----exon1----]--[-----exon2-----]
            #        {-----fwd-----}
            #        <overlap>
            try:
                conserved_down_seq, down_start, down_end, down_identity = self._find_conserved_regions_around_junction(val_junction, isoforms)['downstream']
            except:
                logger.info(f"The following junction has no valid downstream conserved sequence: {val_junction}")
                continue

            for overlap in range(exon_5_min_match, exon_5_min_match+10):
                fwd_start = val_junction.junction_pos_in_seq - overlap
                # 5' exon match = overlap
                for fwd_primer_len in range(min_len, max_len+1):
                    # 3' exon match = fwd_primer_len - overlap
                    if exon_3_min_match<= fwd_primer_len-overlap <= exon_3_max_match:
                        fwd_end = fwd_start + fwd_primer_len
                        if fwd_start >=0 and fwd_end < len(val_junction.junction_seq):
                            fwd_seq = val_junction.junction_seq[fwd_start:fwd_end]

                            # check melting temperature and gc ratio of forward primer
                            # note that here gc_fraction is floating number (0<gc_fraction<1)
                            fwd_primer_tm = primer3.calc_tm(fwd_seq)
                            fwd_gc = gc_fraction(fwd_seq)
                            if Tm_min <= fwd_primer_tm <= Tm_max and min_GC <= fwd_gc <= max_GC:
                                # find reverse primer targeting conserved downstream region.
                                # conserved_down_seq starts from the junction
                                # amplicon length = overlap + rev_end_pos (on the conserved region)
                                # [exon1] - junction - [exon2]
                                #   <----------------------> +- 30 bases: junction_seq
                                #              ------------------------------->: from junction downstream_seq
                                # ----------> fwd
                                #  (overlap)               rev <--------
                                #                                       ^
                                #                                   (rev_end_pos)
                                # <------------------------------------> amplicon length

                                amplicon_length = min(desired_amplicon_length, down_end+overlap)
                                for rev_primer_len in range(min_len, max_len + 1):
                                    rev_end = amplicon_length - overlap
                                    rev_start = rev_end - rev_primer_len
                                    if rev_start >= 0 and rev_end < down_end:
                                        rev_seq = conserved_down_seq[rev_start:rev_end]
                                        rev_primer_tm = primer3.calc_tm(rev_seq)
                                        rev_gc = gc_fraction(rev_seq)
                                        tm_diff = abs(rev_primer_tm - fwd_primer_tm)
                                        if Tm_min <= rev_primer_tm <= Tm_max and min_GC <= rev_gc <= max_GC and tm_diff <= Tm_max_diff:
                                            rev_seq = Seq(rev_seq).reverse_complement()
                                            amplicon_seq = fwd_seq[:overlap]+conserved_down_seq[:rev_end]
                                            primers_A.append(
                                                PrimerPair(
                                                    gene,
                                                    fwd_seq,
                                                    str(rev_seq).upper(),
                                                    fwd_primer_tm,
                                                    rev_primer_tm,
                                                    fwd_gc*100,
                                                    rev_gc*100,
                                                    amplicon_seq,
                                                    amplicon_length,
                                                    primer_class=1))

        # Case 2: num_junctions >= 1, rev primer on junction / fwd primer on conserved regions

        for val_junction in valid_junctions:
            # vary position relative to junction
            # check the rev primer overlapping >{PRIMER_EXON_3_MIN_MATCH} bases on the exon
            # EX
            # [-----exon1----]--[-----exon2-----]
            #        {-----rev------------}
            #                    <overlap>
            try:
                conserved_up_seq, up_start, up_end, up_identity = self._find_conserved_regions_around_junction(val_junction, isoforms)['upstream']
            except:
                logger.info(f"The following junction has no valid upstream conserved sequence: {val_junction}")
                continue
            for overlap in range(exon_5_min_match, exon_5_min_match+10):
                rev_end = val_junction.junction_pos_in_seq + overlap
                for rev_primer_len in range(min_len, max_len+1):
                    # 3' exon match = fwd_primer_len - overlap
                    if exon_3_min_match <= rev_primer_len - overlap <= exon_3_max_match:
                        rev_start = rev_end - rev_primer_len
                        if rev_start >=0 and rev_end < len(val_junction.junction_seq):
                            rev_seq = val_junction.junction_seq[rev_start:rev_end]
                            rev_primer_tm = primer3.calc_tm(rev_seq)
                            rev_gc = gc_fraction(rev_seq)
                            if Tm_min <= rev_primer_tm <= Tm_max and min_GC <= rev_gc <= max_GC:
                                overlap_seq = rev_seq[-overlap:]
                                rev_seq = Seq(rev_seq).reverse_complement()
                                """ 
                                # find forward primer targeting conserved upstream region.
                                # conserved_up_seq starts from the junction
                                # amplicon length = overlap + fwd_end_pos (on the conserved region)
                                # [exon1] - junction - [exon2]
                                #   <----------------------> +- 30 bases: junction_seq
                    ----------------------------->: from junction upstream_seq
                    ----------> fwd              (overlap)
                                             rev <----------
                   ^ (fwd_start_pos)
                    <--------------------------------------> amplicon length
                                """
                                amplicon_length = min(desired_amplicon_length, up_end+overlap)
                                for fwd_primer_len in range(min_len, max_len+1):
                                    # fwd_primer will be upstream_seq[-fwd_start:-fwd_end]
                                    # position on upstream_seq (backward)
                                    fwd_start = amplicon_length - overlap
                                    fwd_end = fwd_start - fwd_primer_len
                                    if fwd_start >= 0 and fwd_end < up_end:
                                        fwd_seq = conserved_up_seq[-fwd_start:-fwd_end]
                                        fwd_primer_tm = primer3.calc_tm(fwd_seq)
                                        fwd_gc = gc_fraction(fwd_seq)
                                        tm_diff = abs(fwd_primer_tm - rev_primer_tm)
                                        if Tm_min <= fwd_primer_tm <= Tm_max and min_GC <= fwd_gc <= max_GC and tm_diff <= Tm_max_diff:
                                            amplicon_seq = conserved_up_seq[-fwd_start:]+overlap_seq
                                            primers_B.append(
                                                PrimerPair(gene,
                                                    fwd_seq,
                                                    str(rev_seq).upper(),
                                                    fwd_primer_tm,
                                                    rev_primer_tm,
                                                    fwd_gc*100,
                                                    rev_gc*100,
                                                    amplicon_seq,
                                                    amplicon_length,
                                                    primer_class=2))


                                #            upstream seq
                                #  ----------------------------> <junction>
                                #  ^                          ^
                                # (end)                  (start)
                                # actual sequence = upstream_seq[-end:]

        return primers_A, primers_B

    def _find_conserved_regions_around_junction(self, junction: ExonJunction,
                                                isoforms: Dict[str, Dict]):
        """Find conserved regions around a junction in isoforms"""
        junction_seq = junction.junction_seq.upper()
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
                                             #"junction_seq": junction_seq,
                                             #"junction_pos_in_seq": junction_pos_in_seq,
                                            #"junction_pos_in_mRNA": junction_pos_in_mRNA,
                                             "upstream_seq": upstream_seq,
                                             "downstream_seq": downstream_seq}


        # validate conserved regions across isoforms
        return self._validate_conserved_regions_around_junction(conserved_regions, 0.95, self.default_params['PRIMER_PRODUCT_MIN'])

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
                valid_conserved_regions['downstream'] = (downstream_seq[0].upper(), start, end, avg_identity)

        # Scan through upstrea sequnces to find valid conserved downstream seq
        # upstream sequence is a little tricky because they need to be checked backward.
        # [exon1] --- junction --- [exon2]
        # <--------------- upstream
        #                ^
        #      check the identitiy from here
        alignment_up = MultipleSeqAlignment([SeqRecord(Seq(up_seq[-1:-1-upstream_min_length:-1])) for up_seq in upstream_seq])
        i = 0

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
                valid_conserved_regions['upstream'] = (upstream_seq[0].upper(), start, end, avg_identity)
                # actual sequence = upstream_seq[-end:]
                # start position doesn't matter; it's always 0

        return valid_conserved_regions

    def _design_on_conserved_primers(self, gene: str)->List[PrimerPair]:
        conserved_regions = self._find_conserved_regions(gene,
                                                         min_length = round(self.default_params['DESIRED_AMPLICON_LENGTH']/2))
        # both primers are not on junction but on conserved regions
        buffer = 20 # avoid primers on 5'end or 3'end of the conserved region.
        primers_C = []
        logger.info(f"{gene} has {len(conserved_regions)} conserved regions across isoforms.")

        for conserved_seq, identity in conserved_regions:
            target_start = buffer
            target_end = len(conserved_seq)-buffer
            target_len = target_end - target_start
            target_seq = conserved_seq[target_start:target_end]
            opt_size = min(target_len, self.default_params['DESIRED_AMPLICON_LENGTH'])

            seq_args = {
                'SEQUECE_ID': gene,
                'SEQUENCE_TEMPLATE': target_seq
            #    'SEQUENCE_TARGET': [target_start, target_end]
            }

            global_args = {
                'PRIMER_OPT_SIZE': self.default_params['PRIMER_OPT_SIZE'],
                'PRIMER_MIN_SIZE': self.default_params['PRIMER_MIN_SIZE'],
                'PRIMER_MAX_SIZE': self.default_params['PRIMER_MAX_SIZE'],
                'PRIMER_OPT_TM': self.default_params['PRIMER_OPT_TM'],
                'PRIMER_MIN_TM': self.default_params['PRIMER_MIN_TM'],
                'PRIMER_MAX_TM': self.default_params['PRIMER_MAX_TM'],
                'PRIMER_MIN_GC': self.default_params['PRIMER_MIN_GC']*100,
                'PRIMER_MAX_GC': self.default_params['PRIMER_MAX_GC']*100,
                'PRIMER_PRODUCT_SIZE_RANGE': [[int(opt_size*0.85), int(opt_size)*1.15]],
                'PRIMER_NUM_RETURN': self.default_params['PRIMER_NUM_RETURN']*2,

            }

            # design primers
            results = primer3.design_primers(seq_args, global_args)

            # Extract primer pairs
            num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)

            for i in range(num_returned):
                left_pos, _ = results[f'PRIMER_LEFT_{i}']
                right_pos, _ = results[f'PRIMER_RIGHT_{i}']
                primers_C.append(
                    PrimerPair(gene,
                               results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                               results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                               results[f'PRIMER_LEFT_{i}_TM'],
                               results[f'PRIMER_RIGHT_{i}_TM'],
                               results[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                               results[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
                               str(target_seq[left_pos:right_pos+1]).upper(),
                               results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                               primer_class = 3
                               )


                    )


        return primers_C

    def _find_conserved_regions(self,
                                gene,
                                min_length: int = 100,
                                min_identity: float = 0.95)-> List[Tuple[str, float]]:
        """
        Find conserved regions across all isoforms of a gene

        conserved_regions = {gene name: [(start_pos, end_pos, identity)]}
        positions here are based on the first sequences in the alignment fasta.

        skip if there is only one isoform

        """

        # align the sequences
        fasta_file = os.path.join(self.fasta_dir, gene+".fasta")
        align_file = os.path.join(self.align_dir, f"aligned_{gene}.fasta")
        # count how many sequences in a fasta file
        fasta_len = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        num_fasta = fasta_len
        # if a fasta file has one sequence, it doesn't require alignment.
        if fasta_len == 1:
            with open(align_file, "w") as f:
                temp_record = SeqIO.read(fasta_file, "fasta")
                SeqIO.write(temp_record, f, "fasta")
                temp_seq = temp_record.seq
            return [(temp_seq.upper(),1.0)]

        else:
            # alignment using mafft
            subprocess.run(['mafft', '--auto', '--quiet', fasta_file],
                           stdout=open(align_file, 'w'))
        #logger.info(f"Saved aligned fasta files for {gene} ({fasta_len} sequence(s) aligned)")
        # Find conserved regions
        alignment = AlignIO.read(align_file, "fasta")
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
                    #ungapped_start = self._gapped_to_ungapped(alignment[0], start)
                    #ungapped_end = self._gapped_to_ungapped(alignment[0], end)
                    #base_id = alignment[0].id
                    gapped_conserved_sequence = alignment[0].seq
                    upgapped_conserved_sequence = gapped_conserved_sequence.upper()[start:end].replace("-","")

                    conserved_regions.append((upgapped_conserved_sequence, avg_identity))
            else:
                i += 1

        return conserved_regions


    def _calculate_identity(self, column: str) -> float:
        """Calculate identity for an alignment column"""
        # Remove gaps

        if '-' in column:
            return 0.0
        else:
            non_gap_chars = [c for c in column]

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


