# This file is part of the csbioscripts package
#
# Copyright (c) 2025 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

from Bio import Align

def scoresequencealignment(seqA, seqB):
    aligner = Align.PairwiseAligner()
    score = aligner.score(seqA, seqB)
    return score

def alignsequences(refseq, queryseq, gapopenpenalty=-0.5, gapextendpenalty=-0.1):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = gapopenpenalty
    aligner.extend_gap_score = gapextendpenalty
    alignments = aligner.align(refseq, queryseq)
    return alignments
