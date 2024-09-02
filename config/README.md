# General settings

To configure the workflow, modify `config/config.yaml` according to yout needs, following the explanations provided in the file.

# Assembly units sheet

Add assemblies to file `config/assembly_units.txt`. This file is tab-separated file with two columns: `assembly_id` and `path`. The former is the alias of the assembly and the latter the path of .fasta assembly file. Both columns are mandatory and cannot be duplicated. 


# Target variants sheet

Add target variants where you want to design the Tm-shit primers in the file `config/target_variants.txt`. This file is tab-sepparated. In the column `variant_id` specify a unique name to identify the variant. The `assembly_id` uses the `assembly_id` specified in `config/assembly_units.txt` and the columns `chrom` and `pos` define the physical position on the given assembly. In columns `ref` and `alt` is stored the reference and alternative allele that the primers going to genotype.