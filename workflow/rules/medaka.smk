rule call_small_variantions:
  input:
    expand("../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vep", sample=config["samples"]),
    expand("../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.assignedID.vcf",sample=config["samples"])
rule medaka:
  input:
    ref=config["Files"]["ref"],
    concat=lambda wildcards: config["samples"][wildcards.sample]
  output:
    dir=directory("../analysis/medaka/medaka_{sample}/"),
    vcf="../analysis/medaka/medaka_{sample}/medaka.annotated.vcf", 
    bam="../analysis/medaka/medaka_{sample}/calls_to_ref.bam" # I had to add this line because SV.smk is dependent on this output
  params: 
    model=lambda wildcards: config["models"][wildcards.sample],
    num_threads=4,
  log: 
    "logs/medaka_{sample}_log/medaka_{sample}.log"
  shell:
    """
    ulimit -n 8000
    medaka_haploid_variant -r {input.ref} -i {input.concat} -m {params.model} -t {params.num_threads} -o {output.dir} 2> {log}
    """
rule SNPSift_VarType:
  #run SNP_Sift vartype on the medaka.annotated.vcfs
  input:
    medaka_annotated_infile="../analysis/medaka/medaka_{sample}/medaka.annotated.vcf"
  output:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf"
  log:
    "logs/medaka_{sample}_log/snpSift_VarType.log"
  params:
    snpEff_path=config["tools"]["snpEff_path"]
  shell:
    """
    java -jar {params.snpEff_path} varType {input.medaka_annotated_infile} > {output} 2> {log}
    """
rule assign_id_medaka_all:
  input: 
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf" 
  output:
    corrected_head="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.vcf",
    uncompressed="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.assignedID.vcf"
  params:
    col_name="Toxo_ME49_{sample}"
  shell: 
    """
    awk 'BEGIN {{OFS="\t"}} /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/ {{$NF = "{params.col_name}"}} 1' {input} > {output.corrected_head}
    cat {output.corrected_head} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
    """  
rule SNPSift_filter:
  #run SNPSift_filter for the medaka.annotated_with_Vartype.vcf
  input:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf"
  output: 
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1.vcf"
  params:
    criteria='"(QUAL >= 1)"',
    snpEff_path=config["tools"]["snpEff_path"]
  log: 
    "logs/medaka_{sample}_log/snpSift_filter.log" 
  shell: 
    """
    java -jar {params.snpEff_path} filter {params.criteria} {input} > {output} 2> {log}
    """
   
rule counting_variant_types_before_filtering:
  input:
   All="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf",
   qual_1="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1.vcf"
  output:
    "../analysis/medaka/medaka_{sample}/variant.statistics.txt"
  shell:
    """
    python3 scripts/2023.12.11_variants_counting_from_vcf_file.py {input.All} > {output}
    python3 scripts/2023.12.11_variants_counting_from_vcf_file.py {input.qual_1} >> {output}
    """
rule pileup:
  input:
   bam="../analysis/medaka/medaka_{sample}/calls_to_ref.bam",
   ref=config["Files"]["ref"]
  output:
    "../analysis/pileup/{sample}_mpileup.txt"
  params:
    samtools_bin=config["tools"]["samtools"]   
  shell:  
   """
   perl ../scripts/extractPileUpFrequencies.pl --samtools_bin {params.samtools_bin} --BAM {input.bam} --outputFile {output} --referenceGenome {input.ref}     
   """
rule filtering_vcf_for_valid_variants:
  input:
   vcf="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1.vcf",
   pileup= "../analysis/pileup/{sample}_mpileup.txt", #lambda wildcards: config["pileup"][wildcards.sample],
   stats="../analysis/medaka/medaka_{sample}/variant.statistics.txt"
  output:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.vcf"
  shell:
    """
    python3 scripts/extract_valid_SNPs_Indels_basic_for_snakemake.py {input.vcf} {input.pileup} {output} >> {input.stats}
    """
rule correct_header_medaka:
  input:
   "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.vcf"
  output:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.vcf"
  params:
    col_name="Toxo_ME49_{sample}"
  shell:
    """
    awk 'BEGIN {{OFS="\t"}} /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/ {{$NF = "{params.col_name}"}} 1' {input} > {output}
    """    
rule sort_vcf_medaka:
  input: 
   "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.vcf"
  output:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf"
  shell:
    "bcftools sort {input} > {output}"     
rule assign_id_medaka:
  input:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf" 
  output:
    uncompressed="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.assignedID.vcf", 
    compressed="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.assignedID.vcf.gz" 
  shell: 
    """
    cat {input} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
    bgzip -c {output.uncompressed} > {output.compressed}
    tabix -p vcf {output.compressed}
    """  
rule variant_annotation:
  input:
    val="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf",
    gff=config["Files"]["gff"],
    ref=config["Files"]["ref"]
  output:
    vep="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vep" 
  params:
    vep_path=config["tools"]["Ensemble_vep"] 
  shell:  
   """
   {params.vep_path} -i {input.val} -gff {input.gff} --fasta {input.ref} -o {output.vep}
   """
rule bgzip_concat_vcf_medaka:
  input: 
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf"
  output:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf.gz"
  shell:
    "bgzip {input}"

rule index_concat_vcf_medaka:
  input:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf.gz"
  output:
    "../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.vcf.gz.tbi" 
  shell:
    "tabix -p vcf {input}"