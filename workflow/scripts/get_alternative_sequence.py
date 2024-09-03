import os, sys, re, argparse, glob
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

sequences = SeqIO.to_dict(SeqIO.parse(snakemake.input.fasta, "fasta"))

for seq_id, sequence in sequences.items():
    seq = str(sequence.seq)
    nseq = seq[:snakemake.params.window_length] + snakemake.params.locus['alt']+ seq[snakemake.params.window_length+len(snakemake.params.locus['ref']):]
    sequence.seq = Seq(nseq)

SeqIO.write(sequences.values(), snakemake.output.fasta, "fasta")