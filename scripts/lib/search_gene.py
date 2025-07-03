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
                logger.info(f"Fetching GenBank records {i + 1}-{min(i + batch_size, len(id_list))} for {gene_name}")

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
        if not isoform_ids:
            return records

        filtered = []
        for record in records:
            # Check if record ID matches any of the specified isoforms
            for isoform in isoform_ids:
                if isoform in record.id or isoform in record.name:
                    filtered.append(record)
                    break

        return filtered

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
        fasta_dir = os.path.join(gene_dir, "fasta")
        os.makedirs(fasta_dir, exist_ok=True)

        # Individual FASTA files
        fasta_records = []
        for gb_file in os.listdir(gene_dir):
            fasta_records.append(SeqIO.read(os.path.join(gene_dir, gb_file), "genbank"))

        filename = os.path.join(fasta_dir, f"{gene_name}.fasta")
        with open(filename, "w") as f:
            SeqIO.write(fasta_records, f, "fasta")


        logger.info(f"Saved {len(genbank_records)} GenBank and {len(fasta_records)} FASTA files for {gene_name}")

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


        '''
        # Parse isoforms
        if isoforms_str and isoforms_str.strip():
            isoform_list = [iso.strip() for iso in isoforms_str.split(',') if iso.strip()]
        else:
            isoform_list = []

        logger.info(f"\n[{i}/{len(genes_data)}] Processing {gene_name} (identifier: {identifier})")
        if isoform_list:
            logger.info(f"  Specific isoforms requested: {', '.join(isoform_list)}")
        else:
            logger.info(f"  Fetching all isoforms")

        # Search for gene
        id_list = fetcher.search_gene(gene_name)

        if not id_list:
            logger.warning(f"No sequences found for {gene_name}")
            summary_data.append({
                'gene': gene_name,
                'identifier': identifier,
                'status': 'not_found',
                'num_genbank': 0,
                'num_fasta': 0,
                'isoforms_found': ''
            })
            continue

        # Fetch sequences
        genbank_records, fasta_records = fetcher.fetch_sequences(id_list, gene_name)

        # Filter isoforms if specified
        if isoform_list:
            logger.info(f"  Filtering for specific isoforms: {', '.join(isoform_list)}")
            genbank_filtered = fetcher.filter_isoforms(genbank_records, isoform_list)
            fasta_filtered = fetcher.filter_isoforms(fasta_records, isoform_list)

            if not genbank_filtered:
                logger.warning(f"  No matching isoforms found for {gene_name}")
                # If no specific isoforms found, use all
                genbank_filtered = genbank_records
                fasta_filtered = fasta_records
        else:
            genbank_filtered = genbank_records
            fasta_filtered = fasta_records

        # Save sequences
        if genbank_filtered or fasta_filtered:
            gene_dir = fetcher.save_sequences(
                gene_name,
                genbank_filtered,
                fasta_filtered,
                output_dir,
                identifier
            )

            # Get list of isoforms found
            isoforms_found = list(set([r.id for r in genbank_filtered]))


        return gene_data

def process_gene_list(csv_file: str, output_dir: str, email: str,
                      api_key: str = None, organism: str = "Mus musculus"):
    """
    Process a CSV file with gene information
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize fetcher
    fetcher = NCBIGeneFetcher(email, api_key, organism)

    # Read CSV file
    genes_data = []
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes_data.append(row)

    logger.info(f"Processing {len(genes_data)} genes from {csv_file}")

    # Create summary file
    summary_file = os.path.join(output_dir, "fetch_summary.csv")
    summary_data = []

    # Process each gene
    for i, gene_info in enumerate(genes_data, 1):
        gene_name = gene_info['gene']
        identifier = gene_info['ident']
        isoforms_str = gene_info.get('isoforms', '')

        # Parse isoforms
        if isoforms_str and isoforms_str.strip():
            isoform_list = [iso.strip() for iso in isoforms_str.split(',') if iso.strip()]
        else:
            isoform_list = []

        logger.info(f"\n[{i}/{len(genes_data)}] Processing {gene_name} (identifier: {identifier})")
        if isoform_list:
            logger.info(f"  Specific isoforms requested: {', '.join(isoform_list)}")
        else:
            logger.info(f"  Fetching all isoforms")

        # Search for gene
        id_list = fetcher.search_gene(gene_name)

        if not id_list:
            logger.warning(f"No sequences found for {gene_name}")
            summary_data.append({
                'gene': gene_name,
                'identifier': identifier,
                'status': 'not_found',
                'num_genbank': 0,
                'num_fasta': 0,
                'isoforms_found': ''
            })
            continue

        # Fetch sequences
        genbank_records, fasta_records = fetcher.fetch_sequences(id_list, gene_name)

        # Filter isoforms if specified
        if isoform_list:
            logger.info(f"  Filtering for specific isoforms: {', '.join(isoform_list)}")
            genbank_filtered = fetcher.filter_isoforms(genbank_records, isoform_list)
            fasta_filtered = fetcher.filter_isoforms(fasta_records, isoform_list)

            if not genbank_filtered:
                logger.warning(f"  No matching isoforms found for {gene_name}")
                # If no specific isoforms found, use all
                genbank_filtered = genbank_records
                fasta_filtered = fasta_records
        else:
            genbank_filtered = genbank_records
            fasta_filtered = fasta_records

        # Save sequences
        if genbank_filtered or fasta_filtered:
            gene_dir = fetcher.save_sequences(
                gene_name,
                genbank_filtered,
                fasta_filtered,
                output_dir,
                identifier
            )

            # Get list of isoforms found
            isoforms_found = list(set([r.id for r in genbank_filtered]))

            summary_data.append({
                'gene': gene_name,
                'identifier': identifier,
                'status': 'success',
                'num_genbank': len(genbank_filtered),
                'num_fasta': len(fasta_filtered),
                'isoforms_found': ';'.join(isoforms_found),
                'output_dir': gene_dir
            })
        else:
            summary_data.append({
                'gene': gene_name,
                'identifier': identifier,
                'status': 'no_valid_sequences',
                'num_genbank': 0,
                'num_fasta': 0,
                'isoforms_found': ''
            })

    # Write summary
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"\nSummary saved to {summary_file}")

    # Print summary statistics
    total_success = len(summary_df[summary_df['status'] == 'success'])
    total_not_found = len(summary_df[summary_df['status'] == 'not_found'])
    total_no_valid = len(summary_df[summary_df['status'] == 'no_valid_sequences'])

    print("\n" + "=" * 60)
    print("FETCH SUMMARY")
    print("=" * 60)
    print(f"Total genes processed: {len(genes_data)}")
    print(f"Successfully fetched: {total_success}")
    print(f"Not found in NCBI: {total_not_found}")
    print(f"No valid sequences: {total_no_valid}")
    print(f"\nOutput directory: {output_dir}")
    print(f"Summary file: {summary_file}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch gene sequences from NCBI based on CSV input"
    )
    parser.add_argument(
        "csv_file",
        help="CSV file with gene information (columns: gene, ident, isoforms)"
    )
    parser.add_argument(
        "-o", "--output",
        default="gene_sequences",
        help="Output directory (default: gene_sequences)"
    )
    parser.add_argument(
        "-e", "--email",
        required=True,
        help="Email address for NCBI"
    )
    parser.add_argument(
        "-k", "--api-key",
        help="NCBI API key (optional, but recommended for faster downloads)"
    )
    parser.add_argument(
        "-s", "--species",
        default="Mus musculus",
        help="Species/organism name (default: Mus musculus)"
    )

    args = parser.parse_args()

    # Add timestamp to output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"{args.output}_{timestamp}"

    # Process the gene list
    process_gene_list(
        csv_file=args.csv_file,
        output_dir=output_dir,
        email=args.email,
        api_key=args.api_key,
        organism=args.species
    )


'''