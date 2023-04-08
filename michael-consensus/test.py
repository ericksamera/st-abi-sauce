from Bio import AlignIO, SeqIO
import numpy as np
from collections import Counter

aln_vals = list(AlignIO.parse('michael-consensus/16S-with-gaps.fasta','fasta'))

consensus = ''
for pos in range(aln_vals[0].get_alignment_length()):
    nuc_at_pos = [seq.seq[pos] for seq in aln_vals[0]]
    counts_at_pos = Counter(nuc_at_pos).most_common(5)
    primary_nuc, secondary_nuc = counts_at_pos[0][0], counts_at_pos[0][1]
    
    nuc_to_add = counts_at_pos[0][0][0]
    if counts_at_pos[0][0][0] == '-': nuc_to_add = '-'
    consensus += nuc_to_add
print(consensus)