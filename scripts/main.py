#!/usr/bin/env python3

import os
import argparse
import logging

# from local
from lib import NCBIGeneFetcher


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
# Output: 2024-01-10 14:30:45 - INFO - Processing gene BRCA1
logger = logging.getLogger(__name__)


def main():
    # command-line interface
    parser = argparse.ArgumentParser()
    parser.add_argument("input_type", help="csv or genebank folder", type=str)
    parser.add_argument("-o", "--output_dir", help="enter the output directory", type=str, required=True)

    parser.add_argument("-g", "--genbank_dir", help="enter the genbank directory", type=str)
    parser.add_argument("-t", "--target", help="input file", type=str)
    parser.add_argument("-n", "--ncbi", help="API key for NCBI", type=str)
    parser.add_argument("-u", "--user", help="user email for NCBI", type=str)
    #parser.add_argument("-p", "--project_name", help="project_name", type=str, required=True)
    #parser.add_argument("-t", "--target", help="input file", type=str, required=True)

    #parser.add_argument("-org", "--organism", help="organism species", type=str, required=False)


    args = parser.parse_args()
    input_type = args.input_type
    target_file = args.target
    user_email = args.user
    ncbi_key = args.ncbi
    output_dir = args.output_dir

    if input_type == 'csv':
        if not ncbi_key:
            print("if NCBI key is not provided, only 3 queries per second are allowed.")
        else:
            logger.info(f"Searching for GenBank files from NCBI")

        # load csv file and search gene information from NCBI
        # It will automatically save fasta and genbank files of a gene.
        # If you specify isoforms, the output only includes the filtered fasta/genbank otherwise all isoforms will be exported.

        fetcher = NCBIGeneFetcher()
        csv_data = fetcher.csv_loader(target_file)
        gene_list = csv_data[0]
        isoform_dict = csv_data[1]
        for gene in gene_list:
            try:
                isoform_list = isoform_dict[gene]
            except:
                isoform_list = []
            id_list = fetcher.search_gene(gene)
            genbank_records = fetcher.fetch_sequences(id_list, gene)

            # filter if isoforms are specified.
            if len(isoform_list)>0:
                genbank_records = fetcher.filter_isoforms(genbank_records, isoform_list)

            # save fasta, genbank files under the output directory.
            fetcher.save_sequences(gene, genbank_records, output_dir)

        genbank_dir = os.path.join(output_dir, "genbank")

    elif input_type == "genbank":
        try:
            genbank_dir = args.genbank_dir
        except:
            raise ValueError("genbank_dir is required when input_type is genbank")

    else:
        raise ValueError(("input_type must be csv, fasta, or genbank, but you entered unknown input type '%s'") % input_type)

    #project_name = args.project_name
    #target_file = args.target

    #organism_species = args.organism
    #num_cpu = args.num_cpu



if __name__ == "__main__":
    main()