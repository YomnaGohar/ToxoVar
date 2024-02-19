# ToxoVar
This repository contains a Snakemake pipeline for the comparative analysis of three isolates from three clones of T. gondii ME49, which have been obtained and maintained under different culture conditions. The pipeline is designed to analyze genomic variations and perform pairwise comparisons.

Key features and processes of the pipeline include:

    Calling SNPs and indels using Medaka.
    Detecting structural variations (SVs) with Sniffles.
    Integrating called variations for all isolates into a variation graph constructed by VG.
    Aligning reads to the graph for each sample independently.
    Outputting a set of high-quality variants for each sample.

The goal is to provide a comprehensive framework for understanding the genetic diversity and evolutionary dynamics of T. gondii isolates, leveraging the power of graph genomes and long-read sequencing data to uncover detailed insights into genomic variations.

# Dependencies
    Medaka 1.11.3: This needs to be installed by the user.
    Sniffles 2.2:  This will be automatically installed by the Snakemake pipeline.
    Ensemble_vep:  This needs to be installed by the user. 
    SnpEff:        This needs to be installed by the user.
    GraphAligner:  This will be automatically installed by the Snakemake pipeline.
    
Here's a proofread and refined version of the "Usage" section for better clarity and correctness:

---

### Usage
To utilize the pipeline, users need to modify the `config.yaml` file located in the config directory. This modification includes specifying the directory path to the FASTQ files in `samples:`. It is crucial that the sample name in the `config.yaml` file matches other entries in `pileup:` and `models:`. Additionally, the paths to the reference and GFF files should also be included in the config file.

#### To Call SNPs and Indels
Execute the following command:
```bash
snakemake --profile profiles/ call_small_variations
```

#### To Call Structural Variations (SVs)
Use this command:
```bash
snakemake --profile profiles/ call_small_variations call_SVs
```

    

