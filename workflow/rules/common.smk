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
                        sep = range(1,config['max_dist']),
                        orient = ['53','35'])
    ref_allele_specific = expand('results/{variant_id}/primer3/ref_allele/{variant_id}_{allele}_1_{orient}.for',
                        variant_id = wildcards.variant_id,
                        allele = 'ref',
                        sep = range(1,config['max_dist']),
                        orient = ['53','35'])
    alt_allele_specific = expand('results/{variant_id}/primer3/alt_allele/{variant_id}_{allele}_{sep}_{orient}_{primer_type}_out',
                        variant_id = wildcards.variant_id,
                        allele = 'alt',
                        sep = '1',
                        orient = ['53','35'],
                        primer_type = 'for')
    return ref_common + ref_allele_specific + alt_allele_specific 

def get_alternative_input(wildcards):
    out_dir = checkpoints.build_primer3_input_primercheck.get(**wildcards).output.alt_input
    input_files = [os.path.join(out_dir,file) for file in os.listdir(out_dir) if file[-3:] == 'txt']
    return input_files

wildcard_constraints:
    variant_id = "|".join(variants['variant_id'].unique()),
    assembly_id = "|".join(assemblies['assembly_id']),
    sep = "|".join([str(i) for i in range(1,25)]),
    orient = "|".join(['53','35']),
    allele = "|".join(['ref', 'alt']),
    primer_type = "|".join(['for', 'rev'])