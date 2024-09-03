rule build_primer3_input:
    input:
        fasta = 'results/{variant_id}/{variant_id}_{allele}.fasta'
    output:
        primer3_in = 'results/{variant_id}/primer3/ref_allele/config/{variant_id}_{allele}_{sep}_{orient}_primer3.txt'
    params:
        sep = lambda wildcards: wildcards.sep,
        orientation = lambda wildcards: wildcards.orient,
        allele = lambda wildcards: get_variant_info(wildcards)[wildcards.allele],
        name = lambda wildcards: "{vid}_{al}_{sep}_{ori}".format(
            vid = wildcards.variant_id,
            al = wildcards.allele,
            sep = wildcards.sep,
            ori = wildcards.orient),
        target_pos = config['window_length'],
    log:
        'results/log/build_primer3_input_{variant_id}_{allele}_{sep}_{orient}.log'
    conda:
        '../envs/biopython.yaml'
    shell:
        """
        python workflow/scripts/utils.py build_primer3_input \
        {params[name]} \
        {input.fasta} \
        {params[target_pos]} \
        {params[allele]} \
        {params[orientation]} \
        {params[sep]} \
        {output.primer3_in} 2> {log}
        """
            
checkpoint build_primer3_input_primercheck:
    input:
        primer3_out = 'results/{variant_id}/primer3/ref_allele/{variant_id}_ref_{sep}_{orient}.{primer_type}',
        fasta = 'results/{variant_id}/{variant_id}_alt.fasta'
    output:
        alt_input = directory('results/{variant_id}/primer3/alt_allele/config/{variant_id}_{allele}_{sep}_{orient}_{primer_type}')
    params:
        sep = lambda wildcards: wildcards.sep,
        orientation = lambda wildcards: wildcards.orient,
        locus = lambda wildcards: get_variant_info(wildcards),
        name = lambda wildcards: "{vid}_{al}_{sep}_{ptype}".format(
            vid = wildcards.variant_id,
            al = wildcards.allele,
            sep = wildcards.sep,
            ptype = wildcards.primer_type),
        target_pos = config['window_length'],
    log:
        'results/log/build_primer3_check_input_{variant_id}_{allele}_{sep}_{orient}_{primer_type}.log'
    conda:
        '../envs/biopython.yaml'
    shell:
        """
        python workflow/scripts/utils.py build_primer3_check_primer \
        {params[name]} \
        {input.primer3_out} \
        {input.fasta} \
        {params[orientation]} \
        {params[locus][ref]} \
        {params[locus][alt]} \
        {params[sep]} \
        {output.alt_input} 2> {log}
        """
           
            
rule primer3_search_alt:
    input:
        alt_input = get_alternative_input
    output:
        out_dir = directory('results/{variant_id}/primer3/alt_allele/{variant_id}_{allele}_{sep}_{orient}_{primer_type}_out'),
    conda:
        '../envs/primer3.yaml'
    log:
        'results/log/primer3_check_input_{variant_id}_{allele}_{sep}_{orient}_{primer_type}.log'
    shell:
        """
        if [ ! -d  {output.out_dir} ]; then
            mkdir {output.out_dir}
        fi
        for primer3_in in {input.alt_input}; do
            echo "${{primer3_in}}";
            primer3_core \
            --p3_settings_file resources/primer_check.txt \
            ${{primer3_in}} 2>> {log} && \
            mv {wildcards.variant_id}_{wildcards.allele}_{wildcards.sep}_{wildcards.primer_type}*_{wildcards.orient}.for \
       {output.out_dir}
        done
        """

rule primer3_search:
    input:
        primer3_in = 'results/{variant_id}/primer3/ref_allele/config/{variant_id}_{allele}_{sep}_{orient}_primer3.txt'
    output:
        out = 'results/{variant_id}/primer3/ref_allele/{variant_id}_{allele}_{sep}_{orient}_out.txt',
        forward = 'results/{variant_id}/primer3/ref_allele/{variant_id}_{allele}_{sep}_{orient}.for',
        rev = 'results/{variant_id}/primer3/ref_allele/{variant_id}_{allele}_{sep}_{orient}.rev'
    conda:
        '../envs/primer3.yaml'
    log:
        'results/log/primer3_input_{variant_id}_{allele}_{sep}_{orient}.log'
    shell:
        """
        primer3_core --output={output.out} \
        --p3_settings_file resources/primer_discriminative.txt \
        {input.primer3_in} 2> {log} && \
        mv {wildcards.variant_id}_{wildcards.allele}_{wildcards.sep}_{wildcards.orient}.for {output.forward}
        mv {wildcards.variant_id}_{wildcards.allele}_{wildcards.sep}_{wildcards.orient}.rev {output.rev}
        """