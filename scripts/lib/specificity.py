"""

Specificity check based on BLAST

Originally, I was planning to run BLAST locally, but I realized the database is a bit larger than I thought. It requires > 600 GB...
This version uses Biopyton library, Blast.qblast to run blast remotely.

* check this cookbook in more detail:
https://www.ncbi.nlm.nih.gov/books/NBK279690/pdf/Bookshelf_NBK279690.pdf

- SpecificityChecking: narrowing down primers that are specific to the target
                       against species.
- TODO: against off-target (not necessary in most case but required for checking specificity across different species)

- NOTE that the direction of sequences is from 5' to 3' in most cases.

"""
import os
from Bio import SeqIO, Blast
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import qblast
from dataclasses import dataclass
import pandas as pd
from typing import List, Dict, Tuple, Optional
import logging

# from local
from .primer_design import PrimerPair

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PrimerBlast:
    def __init__(self,
                 output_dir: str,
                 organism: str = "Mus musculus",
                 database: str = "core_nt",
                 identity_threshold: float = 0.8,
                 coverage_threshold: float = 0.8,
                 max_hits: int = 100,
                 delay_between_queries: float = 1.0
                 ):
        """
               Initialize PrimerBLAST

        Args:

        organism: Target organism for specificity (e.g., "Mus musculus")
        database: NCBI database to search against ("nt", "refseq_rna", etc.)
        identity_threshold: Maximum allowed identity with off-targets (0.0-1.0)
        coverage_threshold: Minimum coverage for considering a hit (0.0-1.0)
        max_hits: Maximum number of BLAST hits to retrieve
        delay_between_queries: Delay between BLAST queries (seconds) to be nice to NCBI

        """
        self.organism = organism
        self.database = database
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.max_hits = max_hits
        self.delay = delay_between_queries
        # Map organism names to NCBI taxonomy IDs for better filtering

        self.organism_taxids = {
            "Mus musculus": "10090",
            "Homo sapiens": "9606",
            "Rattus norvegicus": "10116",
            "Danio rerio": "7955",
            "Drosophila melanogaster": "7227",
            "Caenorhabditis elegans": "6239",
            "Arabidopsis thaliana": "3702",
            "Saccharomyces cerevisiae": "4932"
        }
        self.taxid = self.organism_taxids.get(organism, "")
        if not self.taxid:
            logger.warning(f"No taxonomy ID found for {organism}. Using organism name filter.")

        self.output_dir = output_dir
        self.genbank_dir = os.path.join(output_dir, "genbank")
        self.fasta_dir = os.path.join(output_dir, "fasta")
        self.primer_dir = os.path.join(output_dir, "primer")
        self.align_dir = os.path.join(output_dir, "aligned_fasta")
        self.specificity_dir = os.path.join(output_dir, "specificity")
        os.makedirs(self.specificity_dir, exist_ok=True)

    def check_primer_specificity(self,
                                 gene,
                                 primers):
        """
        Args:
        primers: dataframe of PrimerPairs
        gene: name of target gene

        """
        results =[]
        logger.info(f"Checking specificity for {len(primers)} primer pairs from {gene}")

        # merge primers and exclude any duplicate primers
        unique_primers, unique_names = self._load_primer(primers)

        # create multi-sequence fasta
        fasta_content = ""
        for i in range(len(unique_primers)):
            primer =unique_primers[i]
            name = unique_names[i]
            fasta_content += f">{name}\n{primer}\n"

        self.output_file = os.path.join(self.specificity_dir, f"{gene}_specificity.xml")

        results = self._blast_primer_remote(fasta_content, gene)

    def _blast_primer_remote(self,
                             primer_seq: str,
                             query_id: str):

        # Prepare BLAST parameters
        blast_params = {
            'program': 'blastn',
            'short_query': True,
            'database': self.database,
            'sequence': primer_seq,
            'expect': 1000,  # Higher E-value for short sequences
            'word_size': 7,  # Smaller word size for short primers
            'filter': 'F',  # No low complexity filtering
            'hitlist_size': self.max_hits
        }

        # add organism filter
        if self.taxid:
            blast_params['entrez_query'] = f'(txid{self.taxid}[ORGN])'

        results = qblast(**blast_params)
        with open(self.output_file, "wb") as f:
                f.write(results.read())

        #with Blast.parse(self.output_file) as blast_records:
        #    for blast_record in blast_records:
        #        print(blast_record)





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

