import sys
import os, glob
import pandas as pd

# logging
sys.stderr = open(snakemake.log[0], "w")

complement_nucleotide = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def reverse_complement(sequence):
    return ''.join(complement_nucleotide.get(base, base) for base in reversed(sequence.upper()))
def reverse_complement_fasta(input_file):
    # Read the input FASTA file
    sequence = read_fasta(input_file)
    return reverse_complement(sequence)
def read_fasta(file_path):
    sequences = {}
    current_id = None
    current_seq = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id:
        sequences[current_id] = ''.join(current_seq)

    return sequences[current_id]
def read_primer3_out(path):
    # read table 
    if os.path.exists(path):
        table = pd.read_csv(path, delim_whitespace=True, skiprows=3, header=None)
        columns = ['ID','sequence','pi','length','N','GC','Tm','self_any_th','self_end_th','harpin','quality']
        table.columns = columns
        return table

primers = read_primer3_out(snakemake.input.primer3_out)

for row_id, primer in primers.iterrows():
    name = "{base_name}_{primer_id}".format(base_name = snakemake.params.name, primer_id =  primer.ID)

    if snakemake.params.orientation == '35':
        seq = reverse_complement_fasta(snakemake.input.fasta)
        name += '_35'
        
        complement_allele = ''.join([complement_nucleotide [nucleotide] for nucleotide in snakemake.params.locus['alt']])
        primer_seq = primer.sequence[:-len(snakemake.params.locus['ref'])]
        primer_seq += complement_allele[::-1]
        
    else:
        seq = read_fasta(snakemake.input.fasta)
        name += '_53'
        primer_seq = primer.sequence[:-len(snakemake.params.locus['ref'])]
        primer_seq += snakemake.params.locus['alt']


    if not os.path.exists(snakemake.output.alt_input):
        os.makedirs(snakemake.output.alt_input)


    input_name = '{variant_id}_{allele}_{primer_id}.txt'.format(variant_id = snakemake.params.name.split('_')[0],
                                                        allele = "alt",
                                                        primer_id = primer.ID)
    outpath = os.path.join(snakemake.output.alt_input, input_name)

    with open(outpath, 'w') as file:
        file.write("SEQUENCE_ID={name}\n".format(name=name))
        file.write("SEQUENCE_TEMPLATE={seq}\n".format(seq=seq))
        file.write("SEQUENCE_PRIMER={primer_seq}\n".format(primer_seq=primer_seq))
        file.write("P3_FILE_FLAG=1\n")
        file.write("PRIMER_PICK_RIGHT_PRIMER=0\n")
        file.write("PRIMER_PICK_LEFT_PRIMER=1\n=")