rule build_primer3_input:
    input:
        fasta = 'results/{variant_id}/{variant_id}_{allele}.fasta'
    output:
        primer3_in = 'results/{variant_id}/primer3/ref_allele/config/{variant_id}_{allele}_{sep}_{orient}_primer3.txt'
        
    params:
        sep = lambda wildcards: wildcards.sep,
        target_pos = config['window_length']
    run:
        name = wildcards.variant_id+'_'+ wildcards.allele + '_{sep}'.format(sep=str(wildcards.sep))
        
        variant_info = variants.loc[wildcards.variant_id]
        
        
        if wildcards.orient == '35':
            seq = reverse_complement_fasta(input.fasta)
            name += '_35'
            position = len(seq) - params.target_pos 
        else:
            seq = read_fasta(input.fasta)
            name += '_53'
            position = params.target_pos 
        
            if wildcards.allele == 'ref':
                position += len(variant_info.ref)
            else:
                position += len(variant_info.alt)
            
        with open(output.primer3_in, 'w') as file:
            file.write("SEQUENCE_ID={name}\n".format(name=name))
            file.write("SEQUENCE_TEMPLATE={seq}\n".format(seq=seq))
            file.write("SEQUENCE_TARGET={pos},{sep}\n".format(
                pos=position, sep=wildcards.sep))
            file.write("P3_FILE_FLAG=1\n=".format(name=name))
            
            
checkpoint build_primer3_input_primercheck:
    input:
        primer3_out = 'results/{variant_id}/primer3/ref_allele/{variant_id}_ref_{sep}_{orient}.{primer_type}',
        fasta = 'results/{variant_id}/{variant_id}_alt.fasta'
    output:
        alt_input = directory('results/{variant_id}/primer3/alt_allele/config/{variant_id}_{allele}_{sep}_{orient}_{primer_type}')
    params:
        sep = lambda wildcards: wildcards.sep,
        target_pos = config['window_length']
    run:
        import pandas as pd
        primers = read_primer3_out(input.primer3_out)
        
        variant_info = variants.loc[wildcards.variant_id]
        
        for row_id, primer in primers.iterrows():
            name = wildcards.variant_id+'_'+ wildcards.allele 
            name += '_{sep}_{primer_type}_{primer_id}'.format(sep=str(wildcards.sep),
                                                              primer_type = wildcards.primer_type, primer_id=primer.ID)
            if wildcards.orient == '35':
                seq = reverse_complement_fasta(input.fasta)
                name += '_35'
                
                complement_allele = ''.join([complement_nucleotide [nucleotide] for nucleotide in variant_info.alt])
                primer_seq = primer.sequence[:-len(variant_info.ref)]
                primer_seq += complement_allele[::-1]
                
            else:
                seq = read_fasta(input.fasta)
                name += '_53'
                primer_seq = primer.sequence[:-len(variant_info.ref)]
                primer_seq += variant_info.alt

            
            if not os.path.exists(output.alt_input):
                os.makedirs(output.alt_input)
            
            
            input_name = '/{variant_id}_{allele}_{primer_id}.txt'.format(variant_id = wildcards.variant_id,
                                                                  allele = wildcards.allele,
                                                                  primer_id = primer.ID)
            
            with open(output.alt_input + input_name, 'w') as file:
                file.write("SEQUENCE_ID={name}\n".format(name=name))
                file.write("SEQUENCE_TEMPLATE={seq}\n".format(seq=seq))
                file.write("SEQUENCE_PRIMER={primer_seq}\n".format(primer_seq=primer_seq))
                file.write("P3_FILE_FLAG=1\n")
                file.write("PRIMER_PICK_RIGHT_PRIMER=0\n")
                file.write("PRIMER_PICK_LEFT_PRIMER=1\n=")
            
            
            
rule primer3_search_alt:
    input:
        alt_input = get_alternative_input
    output:
        out_dir = directory('results/{variant_id}/primer3/alt_allele/{variant_id}_{allele}_{sep}_{orient}_{primer_type}_out'),
    conda:
        '../envs/primer3.yaml'
    shell:
        """
        if [ ! -d  {output.out_dir} ]; then
            mkdir {output.out_dir}
        fi
        for primer3_in in {input.alt_input}; do
            echo "${{primer3_in}}";
            primer3_core \
            --p3_settings_file resources/primer_check.txt \
            ${{primer3_in}} && \
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
    shell:
        """
        primer3_core --output={output.out} \
        --p3_settings_file resources/primer_discriminative.txt \
        {input.primer3_in} && \
        mv {wildcards.variant_id}_{wildcards.allele}_{wildcards.sep}_{wildcards.orient}.for {output.forward}
        mv {wildcards.variant_id}_{wildcards.allele}_{wildcards.sep}_{wildcards.orient}.rev {output.rev}
        """