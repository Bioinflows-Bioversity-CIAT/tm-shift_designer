from snakemake.utils import validate
import pandas as pd
import glob as glob
import os


# Snakemake configuration
configfile: "config/config.yaml"

def read_variants_file(file_path):
    if os.path.exists(file_path):
        print(f"The file '{file_path}' exists.")
        variants = pd.read_csv(file_path, sep='\t')
        variants.set_index('variant_id', drop=False, inplace=True)
        
        print("{markers} variants are detected for primer design".format(
            markers = str(variants.shape[0])))
        
        return variants
    else:
        print(f"The file '{file_path}' does not exist.")

# Read assemblies
assemblies = pd.read_csv(config['assemblies_units_path'], sep = '\t')
validate(assemblies, schema = "../schemas/assembly_unit.yaml")

# Read target variants
variants = read_variants_file(config['target_variants_path'])
validate(variants, schema = "../schemas/target_variant.yaml")

# Set indices 
assemblies.set_index("assembly_id", inplace = True, drop=False)
variants.set_index("variant_id", inplace = True,drop=False)

# UTIL FUNCTIONS



def get_variant_info(wildcards):
    info = variants.loc[wildcards.variant_id]
    locus_info = {
        'chrom': info.chrom,
        'pos': info.pos,
        'pi': info.pos - config['window_length'],
        'pf': info.pos + config['window_length'],
        'variant_id': info.variant_id,
        'ref': info.ref,
        'alt': info.alt,
        'assembly_id': info.assembly_id}
    return locus_info

def get_assembly(wildcards):
    info = get_variant_info(wildcards)
    return assemblies.loc[info['assembly_id']].path

def read_fasta(file_path):
    sequences = {}
    current_sequence = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line

    return sequences[current_sequence]

def get_all_primers():
    print("test")
    paths = list()
    for n, row in variants.iterrows():
        directory = 'results/{variant_id}/summary'.format(variant_id = n)
        paths.append(directory)
    return paths

def get_primer_summaries():
    paths = list()
    for n, row in variants.iterrows():
        path = 'results/{variant_id}/selected/{variant_id}_summary.csv'.format(variant_id = n)
        paths.append(path)
    return paths

def get_all_primers_outs_per_qtl(wildcards):
    
    ref_common = expand('results/{variant_id}/primer3/ref_allele/{variant_id}_{allele}_{sep}_{orient}.rev',
                        variant_id = wildcards.variant_id,
                        allele = 'ref',
                        sep = range(1,25),
                        orient = ['53','35'])
    ref_allele_specific = expand('results/{variant_id}/primer3/ref_allele/{variant_id}_{allele}_1_{orient}.for',
                        variant_id = wildcards.variant_id,
                        allele = 'ref',
                        sep = range(1,25),
                        orient = ['53','35'])
    alt_allele_specific = expand('results/{variant_id}/primer3/alt_allele/{variant_id}_{allele}_{sep}_{orient}_{primer_type}_out',
                        variant_id = wildcards.variant_id,
                        allele = 'alt',
                        sep = '1',
                        orient = ['53','35'],
                        primer_type = 'for')
    return ref_common + ref_allele_specific + alt_allele_specific 

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

def reverse_complement(sequence):
    return ''.join(complement_nucleotide.get(base, base) for base in reversed(sequence.upper()))


def read_primer3_out(path):
    # read table 
    if os.path.exists(path):
        table = pd.read_csv(path, delim_whitespace=True, skiprows=3, header=None)
        columns = ['ID','sequence','pi','length','N','GC','Tm','self_any_th','self_end_th','harpin','quality']
        table.columns = columns
        return table
    
def get_alternative_input(wildcards):
    out_dir = checkpoints.build_primer3_input_primercheck.get(**wildcards).output.alt_input
    input_files = glob.glob(out_dir + '/*.txt')
    return input_files


# UTIL DICT
complement_nucleotide = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

wildcard_constraints:
    variant_id = "|".join(variants['variant_id'].unique()),
    assembly_id = "|".join(assemblies['assembly_id']),
    sep = "|".join([str(i) for i in range(1,25)]),
    orient = "|".join(['53','35']),
    allele = "|".join(['ref', 'alt']),
    primer_type = "|".join(['for', 'rev'])