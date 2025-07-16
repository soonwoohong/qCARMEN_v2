"""

Specificity check based on local BLAST
install blast first on command line
conda install bioconda::blast


- SpecificityChecking: narrowing down primers that are specific to the target
                       against other off-target genes and species.


- NOTE that the direction of sequences is from 5' to 3' in most cases.

"""
import os
import subprocess
import sys

class PrimerBlast:


