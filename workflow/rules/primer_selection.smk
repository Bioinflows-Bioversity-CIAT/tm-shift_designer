rule select_primers:
    input:
        primer3_outs = get_all_primers_outs_per_qtl
    output:
        summary = directory('results/{variant_id}/summary'),
    log:
        'results/log/select_primers_{variant_id}.log'
    conda:
        '../envs/biopython.yaml'
    shell:
        """
        python workflow/scripts/select_primers_tm-shift_summarizer.py summarize {wildcards.variant_id} \
        results/{wildcards.variant_id}/primer3 \
        {output.summary} 2> {log}
        """ 