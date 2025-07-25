import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from local
from .primer_design import CommonPrimerDesign

# first get the primer


class crRNA_Design:
    def __init__(self,
                 output_dir: str,
                 max_amplicon_num: int = 10):
        self.output_dir = output_dir
        self.crRNA_dir = os.path.join(output_dir, "crRNA")
        os.makedirs(self.crRNA_dir, exist_ok=True)
        self.max_amplicon_num = max_amplicon_num
        _load_params = CommonPrimerDesign(output_dir)
        self.default_params = _load_params.default_params

    def read_valid_primers(self,
                           gene_name,
                           valid_primers):
        # rank them based on primer/specificity class
        # 1. primer_class <= 2 and specificity_class ==1
        # : either fwd or rev primers are on a exon-exon junction; they both don't have any blast hits
        selected_primers = valid_primers.query('primer_class<=2 and specificity_class==1')
        if len(selected_primers) < self.max_amplicon_num:
            additional_primers = valid_primers.query('primer_class<=2 and specificity_class==2')
            rank_df = additional_primers.copy()

            #additional_primers = additional_primers.sample(self.max_amplicon_num-len(selected_primers))

            selected_primers = pd.concat([selected_primers, additional_primers])

        amplicons = valid_primers['amplicon_seq']


        # Save amplicons as a temp FASTA file for running badgers
        # these amplicons will be the input for badgers
        temp_fasta_file = os.path.join(self.crRNA_dir, "temp.fasta")
        fasta_records = []
        for i in range(len(amplicons)):
            amplicon_seq = amplicons[i]
            amplicon_id = f"{gene_name}_amplicon_{i}"
            fasta_records.append(SeqRecord(Seq(amplicon_seq), id=amplicon_id, description=""))

        with open(temp_fasta_file, "w") as f:
            SeqIO.write(fasta_records, f, "fasta")

        print(fasta_records)

