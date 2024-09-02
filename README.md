# Snakemake workflow: Tm-shift Designer

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/Bioinflows-Bioversity-CIAT/tm-shift_designer/workflows/Tests/badge.svg?branch=main)](https://github.com/Bioinflows-Bioversity-CIAT/tm-shift_designer/actions?query=branch%3Amain+workflow%3ATests)


This Snakemake workflow assist in designing of Tm (melting temperature) shift markers for SNPs (Single Nucleotide Polymorphisms) genotyping. A Tm-shift markers consists of two allele-specific primers  (each of wich contains a 3'-terminal base that corresponds to one of the two SNP allelic variants) and a common reverse primer that aplifies both alleles. For allele-specific primers are attached GC-rich tails of different lengths such that SNP alleles in genomic DNA samples can be discriminated by the Tm of the PCR products (Wang *et al*., 2005). 

## Description

This worflow takes as input one or multiple assemblies (config/assembly_units.txt) and a list of genomic variants referencing the given assemblies (config/target_variants.txt). For each variant a region of interest (roi) is extracted from the given reference. This roi is centered in the target variant and extend a distance of  **window** (50 bp) upstream and downstream, representing the search space. This file is located at `results/{variant_id}/{variant_id}_ref.fasta`. The file `results/{variant_id}/{variant_id}_alt.fasta` only differ with the former at the position **window+1** where the nucleotide is the alterinative specified in config/target_variants.txt

![roi_definition](resources/assets/img/roi_def.svg)

Once the roi is produced its necesary to build the input file for [primer3](https://github.com/primer3-org/primer3). The input files are located in `results/{variant_id}/ref_allele/config/{variant_id}_ref_{distance}_{orientation}_primer3.txt`. In this input file is specified the searching space (roi) in 5'3' and 3'5' orientation, aditionally is used the primer3 paramenter `SEQUENCE_TARGET` that restrict the search of primers for only those pairs that flank the target region. The target region is expanded with the intention of seek more common primer alternatives. This parameter could range between 1 - `max_dist` (20 bp, defalult), specified in the config file.

![primer3_search](resources/assets/img/primer3_input.svg)

For allele-specific and common primers search, primer3 use a preset configuration file in `resources/primer_discriminative.txt` where the following parameters are setted:
- Primer size (pb): 14(min)/21(opt)/35(max) 
- Primer melting temperature (°C): 58(min)/60(opt)/75(max) 
- Amplicon size (pb): 40-100

Each primer3 call produce the forward and reverse file `results/{variant_id}/primer3/ref_allele/{variant_id}_ref_{dist}_{orientation}.[for|rev]` that specifies the sequences and quality statistics for each potential primer. 

Now for each orientation (5'3',3'5') primers are selected with the script `workflow/scripts/select_primers_tm-shift_summarizer.py` acording to the Tm criteria stated by (Wang *et al*,. 2005):


| Primer type  | Min  | Optimal  | Max  |
|---|---|---|---|
| Allele-specific  | 58.3°C | 59-62°C  |  63.5°C |
| Common  | 62°C | 63-70°C  |  75°C |

The primers that meet those critera are listed in the file `results/{variant_id}/summary/{variant_id}_primers.csv`. Here an example:

|  pr_name |  sequence |  length |  tm | GC_content  |
|---|---|---|---|---|
|  Pv21_49613701_**01_Fa** | gcgggcTCGCCCGCCGACGT  | 14  | 60.852  | 78.57  |
| Pv21_49613701_**01_Fb**  | gcgggcagggcggcTCGCCCGCCGACGG  |  14 |  61.802 |  80.0 |
|  Pv21_49613701_**R01** |  GGGCAGAGACGGCGTGG |  17 | 62.928  | 76.47  |
|  Pv21_49613701_**01_Ra** |  gcgggcGGTCTGGGCGTGGCCA |  16 | 61.827  | 62.686  |
|  Pv21_49613701_**01_Rb** |  gcgggcagggcggcGGTCTGGGCGTGGCCC |  16 | 62.686  | 81.25  |
|  Pv21_49613701_**F01** |  GCCCTCGACGAGGCCG |  16 | 62.438  | 81.25  |

The allele-specific primers in 5'3' orientation ends with the suffix **x_Fa** and **x_Fb** and the common primers ends with **_Rx** (_R01, _R02, ...) whereas the allele-specific primers in 3'5' orientations ends with the suffix **x_Ra** and **x_Rb** and the common primers with **_Fx** (_F01, _F02, ...). To the sequences of allele-specific primers are added a GC-rich tails of different length. The long 14 bp GC tail has the sequence `gcgggcagggcggc` and the short 6 bp GC tail has the sequence `gcgggc` (Wang *et al*,. 2005). The long tail is attach to the allele-specific primer with higher Tm. 

The file `results/{variant_id}/summary/{variant_id}_combinations.csv` takes all the primers for each orientation and list all the trio combinations (two allele-specific and one common). Each row represent a trio and is specified for each primer:
- id
- Tm (witout considering the attached GC tail)
- GC content (witout considering the attached GC tail)
- Self any score (tendency of a primer to bind to itself)
- Self end score (tries to bind the 3'-END to a identical primer and scores the best binding it can find)
- Hairpin (tendency of a primer to fold to itself)
- Penalty score

Additionally was calculated for the trio:
- The distance from allele-specific to the begining of common primer (`distance_common`)
- Amplicon size (whithout GC tails, `amplicon_size_[ref|alt]`)
- Tm difference between the reference allele-specific primer and the alternative, `Tm_diff`. 
