rule extract_roi:
    input:
        ref = get_assembly,
        qtns_table = config['target_variants_path']
    params:
        locus = lambda wildcards: get_variant_info(wildcards)
    output:
        fasta = 'results/{variant_id}/{variant_id}_ref.fasta'
    conda:
        '../envs/primer3.yaml'
    shell:
        """
        samtools faidx {input.ref} {params[locus][chrom]}:{params[locus][pi]}-{params[locus][pf]} > {output.fasta}
        """
        
rule get_alternative_roi:
    input:
        fasta = 'results/{variant_id}/{variant_id}_ref.fasta'
    output:
        fasta = 'results/{variant_id}/{variant_id}_alt.fasta'
    params:
        locus = lambda wildcards: get_variant_info(wildcards)
    conda:
        '../envs/biopython.yaml'
    shell:
        """
        python workflow/scripts/utils.py get_alt_sequence {input.fasta} \
        {params[locus][ref]} \
        {params[locus][alt]} \
        {config[window_length]} \
        {output.fasta} 
        """