"""

Fetch FASTA and GenBank files from NCBI based on gene names and (optional) isoforms.
Note that it cannot be used on Della because it doesn't allow for internet access.

"""
from typing import List, Tuple, Dict
import logging
from Bio import Entrez, SeqIO
import pandas as pd
import time
import os
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
# Output: 2024-01-10 14:30:45 - INFO - Processing gene BRCA1

logger = logging.getLogger(__name__)


class NCBIGeneFetcher:
    """Fetch gene sequences from NCBI"""

    def __init__(self, email: str = None, api_key: str = None, organism: str = "Mus musculus"):
        """
        Initialize the fetcher with NCBI credentials

        Args:
            email: User email address for NCBI (default is not set)
            api_key: Optional NCBI API key for higher rate limits
                     (only 3 queries per second are allowed without an API key, 10 queries/sec with a key)
            organism: Default organism (species) name
        """
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.organism = organism
        self.delay = 0.34 if not api_key else 0.11  # Rate limiting

    def search_gene(self, gene_name: str, organism: str = None) -> List[str]:
        """
        Search for gene records in NCBI

        Returns a list of accession IDs
        """
        if organism is None:
            organism = self.organism

        # Build search query - look for mRNA sequences
        query = (
            f'({gene_name}[Gene Name] NOT PREDICTED[Title] AND srcdb_refseq[PROP]) '
            f'AND "{organism}"[porgn] AND biomol_mrna[PROP]'
        )

        #logger.info(f"Searching for {gene_name} in {organism}")

        try:
            # Search in nucleotide database
            with Entrez.esearch(db="nuccore", term=query, retmax=100) as handle:
                results = Entrez.read(handle)

            id_list = results.get("IdList", [])
            logger.info(f"Found {len(id_list)} sequences for {gene_name}")

            time.sleep(self.delay)  # Rate limiting
            return id_list

        except Exception as e:
            logger.error(f"Error searching for {gene_name}: {str(e)}")
            return []

    def fetch_sequences(self, id_list: List[str], gene_name: str) -> List:
        """
        Fetch sequences in both GenBank and FASTA format

        Returns list of genbank_records
        """
        if not id_list:
            return []

        genbank_records = []


        # Fetch in batches to avoid timeouts
        batch_size = 20
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]

            try:
                # Fetch GenBank format
                #logger.info(f"Fetching GenBank records {i + 1}-{min(i + batch_size, len(id_list))} for {gene_name}")

                with Entrez.efetch(
                        db="nuccore",
                        id=",".join(batch_ids),
                        rettype="gbwithparts",
                        retmode="text"
                ) as handle:
                    for record in SeqIO.parse(handle, "genbank"):
                        # Filter out predicted sequences
                        if "PREDICTED" not in record.description.upper():
                            genbank_records.append(record)

                time.sleep(self.delay)


            except Exception as e:
                logger.error(f"Error fetching sequences for {gene_name}: {str(e)}")
                continue

        return genbank_records

    def filter_isoforms(self, records: List, isoform_ids: List[str]) -> List:
        """
        Filter records to only include specified isoforms
        """

        isoform_set = set(isoform_ids)
        return [
            record for record in records
            if record.id in isoform_set or record.name in isoform_set
        ]


    def save_sequences(self,
                       gene_name: str,
                       genbank_records: List,
                       output_dir: str):
        """
        Save sequences to files
        """
        # Create gene directory
        gb_dir = os.path.join(output_dir, "genbank")
        os.makedirs(gb_dir, exist_ok=True)
        gene_dir = os.path.join(gb_dir, gene_name)
        os.makedirs(gene_dir, exist_ok=True)
        # Save GenBank files

        for record in genbank_records:
            filename = os.path.join(gene_dir, f"{record.id}.gb")
            with open(filename, "w") as f:
                SeqIO.write(record, f, "genbank")

        # Save FASTA files
        fasta_dir = os.path.join(output_dir, "fasta")
        os.makedirs(fasta_dir, exist_ok=True)

        # Individual FASTA files
        fasta_records = []
        for gb_file in os.listdir(gene_dir):
            fasta_records.append(SeqIO.read(os.path.join(gene_dir, gb_file), "genbank"))

        filename = os.path.join(fasta_dir, f"{gene_name}.fasta")
        with open(filename, "w") as f:
            SeqIO.write(fasta_records, f, "fasta")


        logger.info(f"Saved {len(genbank_records)} GenBank files for {gene_name}")

        return gene_dir
    def csv_loader(self, csv_file: str):
        """
        process a csv file with gene information

        Args:
            csv_file: path to csv file
                      requires columns: gene, (optional) isoforms
        return
        """

        gene_data = pd.read_csv(csv_file)
        logger.info(f"Processing {len(gene_data)} genes from {csv_file}")
        gene_name = gene_data['gene'].tolist()
        isoform_dict = {}
        for i in range(len(gene_name)):
            temp_isoform = gene_data['isoforms'].iloc[i]
            if not pd.isna(temp_isoform):
                isoform_dict.update = {gene_name[i]: temp_isoform}

        return (gene_name, isoform_dict)

