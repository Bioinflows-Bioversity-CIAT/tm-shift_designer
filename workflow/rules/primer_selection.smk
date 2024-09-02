rule select_primers:
    input:
        primer3_outs = get_all_primers_outs_per_qtl
    output:
        summary = directory('results/{variant_id}/summary'),
    shell:
        """
        python workflow/scripts/select_primers_tm-shift_summarizer.py summarize {wildcards.variant_id} \
        results/{wildcards.variant_id}/primer3 \
        {output.summary}
        """ 