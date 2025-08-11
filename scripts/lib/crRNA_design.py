import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

# from local
from .primer_design import CommonPrimerDesign





class crRNA_Design:
    def __init__(self,
                 output_dir: str,
                 padding: int = 20, # to avoid crRNA located on the edge of amplicon
                 max_amplicon_num: int = 5):
        self.output_dir = output_dir
        self.crRNA_dir = os.path.join(output_dir, "crRNA")
        os.makedirs(self.crRNA_dir, exist_ok=True)
        self.max_amplicon_num = max_amplicon_num
        self.padding = padding
        _load_params = CommonPrimerDesign(output_dir)

        self.default_params = _load_params.default_params
    def design_crRNA(self,
                     gene_name,
                     valid_primer):
        filtered_primers = self._read_valid_primers(gene_name,
                                             valid_primer)
        primer_idx = filtered_primers.index.to_list()
        primer_amplicon = filtered_primers['amplicon_seq']


        for i in range(len(filtered_primers)):
            primer_id = f"{gene_name}_{primer_idx[i]}"
            BADGERS_dir = os.path.join(self.crRNA_dir, gene_name, primer_id)
            os.makedirs(BADGERS_dir, exist_ok=True)
            # making a temp file for BADGERS
            temp_fasta_file = os.path.join(BADGERS_dir, "temp.fasta")
            temp_range_file = os.path.join(BADGERS_dir, "range.tsv")
            if os.path.exists(os.path.join(BADGERS_dir, "final_results.tsv")):
                continue

            temp_amplicon = primer_amplicon.iloc[i]
            if len(temp_amplicon) < self.padding*2:
                range_start, range_end = 0, len(temp_amplicon)
            else:
                range_start, range_end = self.padding, len(temp_amplicon)-self.padding
            with open(temp_fasta_file, "w") as f:
                SeqIO.write(SeqRecord(Seq(temp_amplicon), id=primer_id, description=""), f, "fasta")
            with open(temp_range_file, "w") as f:
                f.write(f"{range_start}\t{range_end}\n")

            BADGERS_script = f'python3 scripts/badgers-cas13/design_guides.py multi both {temp_fasta_file} {BADGERS_dir} --use_range {temp_range_file}'
            BADGERS_script = BADGERS_script + f' --n_top_guides 10'
            BADGERS_script = BADGERS_script + f' --n_top_guides_per_site 3'
            #print(BADGERS_script)
            subprocess.run(BADGERS_script, shell=True)
        # Save amplicons as a temp FASTA file for running badgers
        # these amplicons will be the input for badgers
        """
        temp_fasta_file = os.path.join(self.crRNA_dir, "temp.fasta")
        fasta_records = []
        for i in range(len(amplicons)):
            amplicon_seq = amplicons[i]
            amplicon_id = f"{gene_name}_amplicon_{i}"
            fasta_records.append(SeqRecord(Seq(amplicon_seq), id=amplicon_id, description=""))

        with open(temp_fasta_file, "w") as f:
            SeqIO.write(fasta_records, f, "fasta")
        """

    def _read_valid_primers(self,
                           gene_name,
                           valid_primers):
        desired_amplicon_length = self.default_params["DESIRED_AMPLICON_LENGTH"]
        desired_melt_temp = self.default_params["PRIMER_OPT_TM"]

        # rank them based on primer/specificity class
        # 1. primer_class <= 2 and specificity_class ==1
        # : either fwd or rev primers are on a exon-exon junction; they both don't have any blast hits
        # 2. primer_class <=2 and specificity_class == 2
        # : either fwd or rev primers are on a exon-exon junction, but one of them has blast hits but they don't hit the same target.
        #selected_primers = valid_primers.query('primer_class<=2 and specificity_class==1')
        valid_primers_temp = valid_primers.assign(
            diff_amp=abs(valid_primers["amplicon_len"] - desired_amplicon_length),
            diff_melt=abs(valid_primers["forward_tm"] - desired_melt_temp) + abs(valid_primers["reverse_tm"] - desired_melt_temp),
            avg_gc=(valid_primers["forward_gc"] + valid_primers["reverse_gc"]) / 2)
        selected_primers = valid_primers_temp.query('primer_class<=2')
        # sort them by
        # (0) specificity_class
        # (1) amplicon length (close to the desired amplicon length),
        # (2) melt temp (close to the desired melt temp),
        # (3) low GC ratio
        selected_primers = selected_primers.sort_values(by=['specificity_class', 'diff_amp', 'diff_melt', 'avg_gc'],
                                       ascending=[True, True, True, True]
                                       )

        primer_count = len(selected_primers)
        # if the number of valid primers satisfying this criteria < max_amplicon_num
        # 3. then consider primer_class = 3 (primers on conserved regions)

        if primer_count < self.max_amplicon_num:
            req_num = self.max_amplicon_num - primer_count
            additional_primers = valid_primers_temp.query('primer_class==3')
            additional_primers.sort_values(by=['specificity_class','diff_amp', 'diff_melt', 'avg_gc'],
                                           ascending=[True, True, True, True],
                                           inplace=True)
            selected_primers = pd.concat([selected_primers, additional_primers.iloc[:min(len(additional_primers), req_num)]])
        else:
            selected_primers = selected_primers.iloc[:self.max_amplicon_num]

        selected_primers.drop(columns=['diff_amp', 'diff_melt', 'avg_gc'], inplace=True)

        return selected_primers



