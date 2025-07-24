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
import time
from Bio import SeqIO, Blast
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import qblast
from dataclasses import dataclass, asdict
import pandas as pd
from typing import List, Dict, Tuple, Optional
import logging

# from local
from .primer_design import PrimerPair

@dataclass
class BlastHit:
    """Store BLAST hit information"""
    query_id: str
    subject_id: str
    subject_def: str
    evalue: float
    identity: str
    midline: str
    query_start: int
    query_end: int

@dataclass
class SpecificPrimerPair(PrimerPair):
    """Store primer pair information from Primer-BLAST"""
    specificity_class: int # 1:both primers don't have any unintended target hit; 2:each primer hits a different target


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PrimerBlast:
    def __init__(self,
                 output_dir: str,
                 organism: str = "Mus musculus",
                 database: str = "refseq_rna",
                 max_hits: int = 400,
                 ):
        """
               Initialize PrimerBLAST

        Args:

        organism: Target organism for specificity (e.g., "Mus musculus")
        database: NCBI database to search against ("nt", "refseq_rna", etc.)
        identity_threshold: Maximum allowed identity with off-targets (0.0-1.0)
        coverage_threshold: Minimum coverage for considering a hit (0.0-1.0)
        max_hits: Maximum number of BLAST hits to retrieve

        """
        self.organism = organism
        self.database = database
        self.max_hits = max_hits
        # if there are homologous genes in the target... the homologous gene should be excluded for blast
        self.homologus_gene = {"Tgtp2": ["Tgtp1"]}
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
        self.blast_dir = os.path.join(self.specificity_dir, "raw_blast_xml")
        self.blast_hit_dir = os.path.join(self.specificity_dir, "blast_hits")
        os.makedirs(self.specificity_dir, exist_ok=True)
        os.makedirs(self.blast_dir, exist_ok=True)
        os.makedirs(self.blast_hit_dir, exist_ok=True)


    def check_primer_specificity(self,
                                 gene,
                                 primers):
        """
        Args:
        primers: dataframe of PrimerPairs
        gene: name of target gene

        """

        logger.info(f"Checking specificity for {len(primers)} primer pairs from {gene}")

        # merge primers and exclude any duplicate primers
        unique_primers, unique_names, alias_name = self._load_primer(primers)


        # create multi-sequence fasta
        fasta_content = ""
        for i in range(len(unique_primers)):
            primer =unique_primers[i]
            name = unique_names[i]
            fasta_content += f">{name}\n{primer}\n"

        self.output_file = os.path.join(self.blast_dir, f"{gene}_blast_results.xml")
        # blast the primers
        blast_results = self._blast_primer_remote(fasta_content, gene, alias_name)
        blast_results_df = pd.DataFrame(blast_results)
        blast_results_df.to_csv(os.path.join(self.blast_hit_dir, f"{gene}_blast_hits.csv"), index=False)

        """
        now categorize and analyze the blast hits 
        
        Case1: both primers don't have any unintended target hit
        -> should be included in the primer list and considered as top candidates
        Case2: each primer hits a different target
        -> should be included in the primer list but marked for further review
        Case3: both primers hit the same target
        -> should be excluded from the primer list
        """

        valid_primers = []
        for primer_idx in range(len(primers)):
            fwd_id = f"{primer_idx}_fwd"
            rev_id = f"{primer_idx}_rev"
            fwd_blast_hits = blast_results_df[blast_results_df['query_id'] == fwd_id]['subject_id'].to_list()
            rev_blast_hits = blast_results_df[blast_results_df['query_id'] == rev_id]['subject_id'].to_list()
            common_hits = list(set(fwd_blast_hits) & set(rev_blast_hits))
            if len(common_hits) > 0:
                specificity_class = 3
            else:
                if len(fwd_blast_hits) + len(rev_blast_hits) != 0:
                    specificity_class = 2
                else:
                    specificity_class = 1
                primer = primers.iloc[primer_idx]
                valid_primers.append(SpecificPrimerPair(**primer.to_dict(),
                                                        specificity_class=specificity_class
                                                        ))
        return valid_primers


    def _blast_primer_remote(self,
                             primer_seq: str,
                             gene: str,
                             alias_name: Dict[str, List[str]],
                             min_mismatch_total: int = 2,
                             min_mismatch_3end: int = 2,
                             hotspot_3end: int = 5):
        """
        primer_seq = fasta sequence
        ex. >fwd_0 \n AGCTGT... \n >rev_0 \n AGCTGT... \n
        fasta id will be query id in blast results
        because duplicate has been already removed, it's necessary to assign duplicate ids to different sequences using alilas list.

        From NCBI Blast
        This requires at least one primer (for a given primer pair) to have the specified number of mismatches to unintended targets.
        The larger the mismatches (especially those toward 3' end) are between primers and the unintended targets,
        the more specific the primer pair is to your template (i.e., it will be more difficult to anneal to unintended targets).
        However, specifying a larger mismatch value may make it more difficult to find such specific primers.
        Try to lower the mismatch value in such case.

        return: blast hit that should be later excluded in the filtered primer list.


        :return:
        """
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

        # if the blast result already exists, skip this part
        # the following lines sometimes fail. I think it's because of timeout from running blast remotely.
        # to solve this issue, I added the retry part after 2 minutes rest.
        if not os.path.exists(self.output_file):
            try:
                results = qblast(**blast_params)
            except:
                time.sleep(120)
                results = qblast(**blast_params)

            with open(self.output_file, "wb") as f:
                f.write(results.read())

        hits = []
        # the following also make an error sometimes...
        # the only way is to rerun the code
        try:
            blast_records = Blast.parse(self.output_file)
        except:
            results = qblast(**blast_params)
            with open(self.output_file, "wb") as f:
                f.write(results.read())
            blast_records = Blast.parse(self.output_file)


        for blast_record in blast_records:
            for hsps in blast_record:
                for hsp in hsps:
                    # HSP stands for "High-scoring Segment Pair"
                    # it represents a local alignment between your query sequence (primer) and a subject sequence (database hit).
                    identity_num = hsp.annotations['identity'] # Number of identical positions
                    # checking identity at the hotspot from the 3' end
                    primer_len = len(hsp.query.seq)
                    hotspot_start =  primer_len-hotspot_3end # position from 5'end
                    identity_num_3end = hsp.annotations['midline'][hotspot_start:].count("|")

                    total_mismatch = primer_len - identity_num
                    end_mismatch = hotspot_3end - identity_num_3end
                    target_description = hsp.target.description

                    query_id = hsp.query.description
                    query_alias = alias_name[query_id] # a list of duplicate primer

                    # Primer specificity stringency

                    # if there is a homologous gene that needs to be excluded for blast search
                    homologous_genes = []
                    try:
                        homologous_genes.extend(self.homologus_gene[gene])
                        homologous_genes.append(gene)
                    except:
                        homologous_genes.append(gene)

                    if not any(gene_name.upper() in target_description.upper() for gene_name in homologous_genes):
                        if total_mismatch < min_mismatch_total and end_mismatch < min_mismatch_3end:
                            for alias in query_alias:
                                hit = BlastHit(
                                    query_id=alias,
                                    subject_id=hsp.target.id,
                                    subject_def=hsp.target.description,
                                    evalue=hsp.annotations['evalue'],
                                    identity=f"{identity_num} / {hsp.length}",
                                    midline=hsp.annotations['midline'],
                                    query_start=hsp.coordinates[1][0],
                                    query_end=hsp.coordinates[1][1]
                                )

                                hits.append(hit)



        return hits

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
        alias_names = {}
        seen = set()

        for primer, name in zip(total_primers, total_names):
            if primer not in seen:
                seen.add(primer)
                unique_primers.append(primer)
                unique_names.append(name)
                alias_names[name] = [name]
            else:
                alias_names[unique_names[unique_primers.index(primer)]].append(name)

        return unique_primers, unique_names, alias_names

