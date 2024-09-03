import sys

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

if snakemake.params.orientation == '35':
    seq = reverse_complement_fasta(snakemake.input.fasta)
    position = len(seq) - snakemake.params.target_pos
else:
    seq = read_fasta(snakemake.input.fasta)
    position = snakemake.params.target_pos + len(snakemake.params.allele)
    
with open(snakemake.output.primer3_in, 'w') as file:
    file.write("SEQUENCE_ID={name}\n".format(name=snakemake.params.name))
    file.write("SEQUENCE_TEMPLATE={seq}\n".format(seq=seq))
    file.write("SEQUENCE_TARGET={pos},{sep}\n".format(
        pos=position, sep=snakemake.params.sep))
    file.write("P3_FILE_FLAG=1\n=")


        