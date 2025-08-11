#!/usr/bin/env python3

import os
import argparse
import logging
from Bio import SeqIO
import pandas as pd
from dataclasses import dataclass, asdict
import subprocess

# from local
from lib import NCBIGeneFetcher
from lib import CommonPrimerDesign
from lib import PrimerBlast
from lib import crRNA_Design

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
    parser.add_argument("-pad", "--padding", help="padding for desiging crRNA on amplicon", type=int, required=True, default=20)
    parser.add_argument("-g", "--genbank_dir", help="enter the genbank directory", type=str)
    parser.add_argument("-t", "--target", help="input file", type=str)
    parser.add_argument("-n", "--ncbi", help="API key for NCBI", type=str)
    parser.add_argument("-u", "--user", help="user email for NCBI", type=str)
    #parser.add_argument("-p", "--project_name", help="project_name", type=str, required=True)
    #parser.add_argument("-t", "--target", help="input file", type=str, required=True)

    parser.add_argument("-org", "--organism", help="organism species", type=str, required=False)


    args = parser.parse_args()
    input_type = args.input_type
    target_file = args.target
    user_email = args.user
    ncbi_key = args.ncbi
    output_dir = args.output_dir
    organism = args.organism
    padding = args.padding

    if input_type == 'csv':
        if not ncbi_key:
            print("if NCBI key is not provided, only 3 queries per second are allowed.")
        else:
            logger.info(f"Searching for GenBank files from NCBI")

        # load csv file and search gene information from NCBI
        # It will automatically save fasta and genbank files of a gene.
        # If you specify isoforms, the output only includes the filtered fasta/genbank otherwise all isoforms will be exported.

        fetcher = NCBIGeneFetcher(organism=organism)
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
        fasta_dir = os.path.join(output_dir, "fasta")
    elif input_type == "genbank":
        try:
            genbank_dir = args.genbank_dir
        except:
            raise ValueError("genbank_dir is required when input_type is genbank")
        gene_list = [gb_file for gb_file in os.listdir(genbank_dir) if os.path.isdir(os.path.join(genbank_dir, gb_file))]

        # Save FASTA files
        fasta_dir = os.path.join(output_dir, "fasta")
        os.makedirs(fasta_dir, exist_ok=True)

        # Individual FASTA
        for gene_name in gene_list:
            fasta_records = []
            gene_dir = os.path.join(genbank_dir, gene_name)
            for gb_file in os.listdir(gene_dir):
                if gb_file.endswith(".gb"):
                    fasta_records.append(SeqIO.read(os.path.join(gene_dir, gb_file), "genbank"))

            filename = os.path.join(fasta_dir, f"{gene_name}.fasta")
            with open(filename, "w") as f:
                SeqIO.write(fasta_records, f, "fasta")

    else:
        raise ValueError(("input_type must be csv, fasta, or genbank, but you entered unknown input type '%s'") % input_type)


    primer_dir = os.path.join(output_dir, "primer")
    valid_primer_dir = os.path.join(output_dir, "valid_primer")
    os.makedirs(valid_primer_dir, exist_ok=True)
    primer_design = CommonPrimerDesign(output_dir)
    specificity_checker = PrimerBlast(output_dir, organism)
    crRNA_design = crRNA_Design(output_dir, padding=padding)


    for gene_name in gene_list:
        parent_path = os.path.join(genbank_dir, gene_name)
        genbank_files = [os.path.join(parent_path, genfile) for genfile in os.listdir(parent_path) if ".gb" in genfile]

        # designing primers
        primers = primer_design.design_primers(gene_name, genbank_files)
        primers.to_csv(os.path.join(primer_dir, gene_name+"_primers.csv"))

        # specificity checking
        valid_primers = specificity_checker.check_primer_specificity(gene_name, primers)
        valid_primers_df = pd.DataFrame(valid_primers)
        valid_primers_df.to_csv(os.path.join(valid_primer_dir, gene_name+"_valid_primers.csv"))

        # design crRNA
        crRNA_design.design_crRNA(gene_name, valid_primers_df)




    #conserved_regions, num_fasta = primer_design.find_conserved_regions()

    #print(gene_isoform_dict)

    #project_name = args.project_name
    #target_file = args.target

    #organism_species = args.organism
    #num_cpu = args.num_cpu



if __name__ == "__main__":
    main()