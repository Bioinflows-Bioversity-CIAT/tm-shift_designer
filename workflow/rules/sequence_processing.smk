rule extract_roi:
    input:
        ref = get_assembly,
        qtns_table = config['target_variants_path']
    params:
        locus = lambda wildcards: get_variant_info(wildcards)
    output:
        fasta = 'results/{variant_id}/{variant_id}_ref.fasta'
    log:
        'results/log/samtools_faidx_{variant_id}'
    conda:
        '../envs/primer3.yaml'
    shell:
        """
        samtools faidx {input.ref} {params[locus][chrom]}:{params[locus][pi]}-{params[locus][pf]} > {output.fasta} 2> {log}
        """
        
rule get_alternative_roi:
    input:
        fasta = 'results/{variant_id}/{variant_id}_ref.fasta'
    output:
        fasta = 'results/{variant_id}/{variant_id}_alt.fasta'
    params:
        locus = lambda wildcards: get_variant_info(wildcards),
        window_length = config['window_length']
    log:
        'results/log/get_alternative_roi_{variant_id}'
    conda:
        '../envs/biopython.yaml'
    script:
        '../scripts/get_alternative_sequence.py'