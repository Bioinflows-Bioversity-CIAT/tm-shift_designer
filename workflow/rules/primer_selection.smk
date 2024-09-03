rule select_primers:
    input:
        primer3_outs = get_all_primers_outs_per_qtl
    output:
        summary = directory('results/{variant_id}/summary'),
    log:
        'results/log/select_primers_{variant_id}.log'
    params:
        basedir = lambda wildcards: "results/{vid}/primer3".format(vid = wildcards.variant_id)
    conda:
        '../envs/biopython.yaml'
    script:
        '../scripts/select_primers_tm-shift_summarizer.py'