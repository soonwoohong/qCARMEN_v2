# this will add T7 promoter, blocker, and convert usual primers to rhPrimer GEN1
# example: DDDDDDDDrDDDDMx
# check the idt website about rhprimer: https://www.idtdna.com/pages/products/qpcr-and-pcr/custom-primers/rhpcr-primers
# the proto-spacer of crRNA into the actual crRNA with DR.

import pandas as pd
import os


def final_mod(crRNA_df):
    T7_promoter = "TAATACGACTCACTATAggg"
    blocker = "/iSpC3/"
    direct_repeat = "rCrUrArArArUrCrUrGrArUrGrGrGrGrUrUrUrUrUrGrCrUrUrCrCrCrCrUrGrArUrUrUrUrG"
    fwd = []
    rev =[]
    crRNA = []
    mismatch = {"A": "G",
                "T": "C",
                "C": "A",
                "G": "T"
                }
    base_pair = {"A": "T",
                 "T": "A",
                 "C": "G",
                 "G": "C"
                 }
    base_pair_RNA = {"A": "rU",
                     "T": "rA",
                     "C": "rG",
                     "G": "rC"
                     }

    for primer_idx in range(len(crRNA_df)):
        temp_fwd = crRNA_df['forward_seq'][primer_idx]
        temp_rev = crRNA_df['reverse_seq'][primer_idx]
        temp_amplicon = crRNA_df['amplicon_seq'][primer_idx]
        temp_guide = crRNA_df['guide_sequence'][primer_idx]
        len_fwd = len(temp_fwd)
        len_rev = len(temp_rev)

        fwd_add_on = "r"+temp_amplicon[len_fwd:len_fwd+5] # convert the first DNA base to RNA base
        fwd_add_on = fwd_add_on[:-1]+mismatch[fwd_add_on[-1].upper()] # make the intentional mismatch at the final base
        final_fwd = T7_promoter + temp_fwd + fwd_add_on + blocker

        # for reverse primer, it needs to get the reverse complement of amplicon
        rev_comp_amplicon = temp_amplicon[::-1].upper().translate(str.maketrans(base_pair))
        rev_add_on = "r"+rev_comp_amplicon[len_rev:len_rev+5]
        rev_add_on = rev_add_on[:-1]+mismatch[rev_add_on[-1].upper()]
        final_rev = temp_rev + rev_add_on + blocker

        # for final crRNA
        temp_guide_rev_com = temp_guide[::-1].translate(str.maketrans(base_pair_RNA))
        final_crRNA = direct_repeat + temp_guide_rev_com

        fwd.append(final_fwd)
        rev.append(final_rev)
        crRNA.append(final_crRNA)

    final_primers = pd.DataFrame({'final_fwd_primer': fwd, 'final_rev_primer': rev, 'final_crRNA_w/_DR': crRNA})
    final_output = pd.concat([final_primers, crRNA_df], axis=1)
    final_output.insert(0,'gene',final_output.pop('gene'))

    return final_output