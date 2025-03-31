# ToxoVar
This repository contains a Snakemake pipeline for the comparative analysis of three isolates from three clones of T. gondii ME49, which have been obtained and maintained under different culture conditions. The pipeline is designed to analyze genomic variations and perform pairwise comparisons.

Key features and processes of the pipeline include:

    Calling SNPs and indels using Medaka.
    Detecting structural variations (SVs) with Sniffles.
    Integrating called variations for all isolates into a variation graph constructed by VG.
    Aligning reads to the graph for each sample independently.
    Outputting a set of high-quality variants for each sample.

The goal is to provide a comprehensive framework for understanding the genetic diversity and evolutionary dynamics of T. gondii isolates, leveraging the power of graph genomes and long-read sequencing data to uncover detailed insights into genomic variations.
## Resources

The `resources/` directory contains essential reference files required by the pipeline:

- `2015T_assembly.fa`: A high-quality genome assembly of the *T. gondii* ME49 2015T isolate. (variants are called against it.)
- `annotation.gff3.gz`: Genome annotation file generated using **Companion**
- `homo.bed`: Genomic positions of **homopolymeric regions**
- `TR.bed`: Genomic positions of **tandem repeats (TRs)**
- `numts.bed`: Genomic positions of **nuclear mitochondrial DNA segments (NUMTs)**
- `positions_to_mask_variants_in.bed`: Genomic positions that expected to contain false positive variant calls

These BED files are used for variant filtering and quality control steps within the workflow.

# Dependencies
    Medaka 1.11.3:
    Sniffles 2.2:  
    Ensemble_vep: 
    SnpEff:       
    GraphAligner:  
    
---

### Usage
To utilize the pipeline, users need to modify the `config.yaml` file located in the config directory. This modification includes specifying the directory path to the FASTQ files in `samples:`. It is crucial that the sample name in the `config.yaml` file matches other entries in `pileup:` and `models:`. Additionally, the paths to the reference, GFF and BED files should also be included in the config file. For information about each entry in the config.yaml file, a sample file (config_sample.yaml) is provided in the config/ directory

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

    

