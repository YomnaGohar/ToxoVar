rule call_small_variantions:
  input:
     expand("{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf", out=config["output_dir"],sample=config["samples"]),
     expand("{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.assignedID.vcf", out=config["output_dir"],sample=config["samples"]),
     expand("{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.vcf", out=config["output_dir"],sample=config["samples"]),
     expand("{out}/analysis/pileup/{sample}_mpileup.txt", out=config["output_dir"],sample=config["samples"] ),
     expand("{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf.gz.tbi",out=config["output_dir"],sample=config["samples"])
   
rule medaka:
  input:
    ref=config["Files"]["ref"],
    concat=lambda wildcards: config["samples"][wildcards.sample]
  output:
    dir=directory("{out}/analysis/medaka/medaka_{sample}/"),
    vcf="{out}/analysis/medaka/medaka_{sample}/medaka.annotated.vcf",
    bam="{out}/analysis/medaka/medaka_{sample}/calls_to_ref.bam"
  params: 
    model=lambda wildcards: config["models"][wildcards.sample],
    num_threads=4
  log: 
    "{out}/analysis/logs/medaka_{sample}_log/medaka_{sample}.log"
  threads: 10
  shell:
    """
    ulimit -n 8000
    medaka_variant -r {input.ref} -i {input.concat} -m {params.model} -t {params.num_threads} -o {output.dir} -b 50 2> {log}
    """

rule exclude_chromosomes:
    input:
        vcf="{out}/analysis/medaka/medaka_{sample}/medaka.annotated.vcf"
    output:
        vcf="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_no_api.vcf"
    params:
        exclude_string=lambda wildcards: " || ".join([f'CHROM="{seq}"' for seq in config["Files"]["sequences_to_exclude"].split(",")])
    shell:
        """
        bcftools view --exclude '{params.exclude_string}' {input.vcf} -o {output.vcf}
        """
rule SNPSift_VarType:
  input:
    medaka_annotated_infile="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_no_api.vcf"
  output:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf"
  log:
    "{out}/analysis/logs/medaka_{sample}_log/snpSift_VarType.log"
  params:
    snpEff_path=config["tools"]["snpEff_path"]
  threads: 4   
  shell:
    """
    java -jar {params.snpEff_path} varType {input.medaka_annotated_infile} > {output} 2> {log}
    """
rule assign_id_medaka_all:
  input: 
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf" 
  output:
    corrected_head="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.vcf",
    uncompressed="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.assignedID.vcf"
  params:
    col_name="Toxo_ME49_{sample}"
  threads: 4   
  shell: 
    """
    awk 'BEGIN {{OFS="\t"}} /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/ {{$NF = "{params.col_name}"}} 1' {input} > {output.corrected_head}
    cat {output.corrected_head} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
    """  
rule pileup:
  input:
   bam="{out}/analysis/medaka/medaka_{sample}/calls_to_ref.bam",
   ref=config["Files"]["ref"]
  output:
    "{out}/analysis/pileup/{sample}_mpileup.txt"
  params:
    samtools_bin=config["tools"]["samtools"]   
  threads: 10 
  shell:  
   """
   perl scripts/extractPileUpFrequencies_v2.pl --samtools_bin {params.samtools_bin} --BAM {input.bam} --outputFile {output} --minMappingQuality 0 --minBaseQuality 0 --referenceGenome {input.ref}    
   """
rule filtering_vcf_for_valid_variants:
  input:
   vcf="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType.vcf",
   pileup= "{out}/analysis/pileup/{sample}_mpileup.txt", #lambda wildcards: config["pileup"][wildcards.sample],
  output:
    valid="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.vcf",
    stat="{out}/analysis/medaka/medaka_{sample}/variant.statistics.txt"
  threads: 4  
  params:
     depth=config["filter"]["depth"],
     gq=config["filter"]["gq"],
     snp_del_overlap=config["filter"]["snp_del_overlap"]
  shell:
    """
    python3 scripts/2023.12.11_variants_counting_from_vcf_file.py {input.vcf} > {output.stat}
    python3 scripts/medaka_filtering.py {input.vcf} {input.pileup} {params.depth} {params.gq} {params.snp_del_overlap}
    python3 scripts/2023.12.11_variants_counting_from_vcf_file.py {output.valid} >> {output.stat}
    """               
rule correct_header_medaka:
  input:
   "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.vcf"
  output:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.vcf"
  params:
    col_name="Toxo_ME49_{sample}"
  threads: 4   
  shell:
    """
    awk 'BEGIN {{OFS="\t"}} /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/ {{$NF = "{params.col_name}"}} 1' {input} > {output}
    """    
rule sort_vcf_medaka:
  input: 
   "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.vcf"
  output:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf"
  threads: 4   
  shell:
    "bcftools sort {input} > {output}"     
rule assign_id_medaka:
  input:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf" 
  output:
    uncompressed="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID.vcf", 
    compressed="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID.vcf.gz" 
  threads: 4    
  shell: 
    """
    cat {input} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
    bgzip -c {output.uncompressed} > {output.compressed}
    tabix -p vcf {output.compressed}
    """  
rule bgzip_concat_vcf_medaka:
  input: 
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf"
  output:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf.gz"
  threads: 4   
  shell:
    "bgzip {input}"

rule index_concat_vcf_medaka:
  input:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf.gz"
  output:
    "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.vcf.gz.tbi" 
  threads: 4     
  shell:
    "tabix -p vcf {input}"
    

        
