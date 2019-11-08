import re

import numpy as np


def parse_molcas_couplings(text):
    regex = re.compile("Total derivative coupling(.+?)\-+\s+norm", re.DOTALL)
    blocks = regex.findall(text)
    points = len(blocks)
    nac_lines = [
        [line.split() for line in block.split("\n")[8:-1]] for block in blocks
    ]
    nac_lines = np.array(nac_lines)
    nacs = nac_lines[:,:,1:].astype(float)
    # atoms = nac_lines[:,:,0]

    return nacs
