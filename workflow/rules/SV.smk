rule call_SVs:
  input:
    expand("../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf", sample=config["samples"])
rule sniffles:
    input:
        bam="../analysis/medaka/medaka_{sample}/calls_to_ref.bam",
        ref=config["Files"]["ref"]
    output:
        "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference.vcf"
    log:
        "logs/sniffles/log_{sample}_sniffles"
    conda: 
        config["conda_env"]["sniffles"]
    shell:
        "sniffles --input {input.bam} --reference {input.ref} --vcf {output} > {log} 2>&1"

rule filter_vcf:
    input:
        vcf="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference.vcf"
    output:
        valid="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid.vcf",
        statistics="../analysis/Sniffles/{sample}/statistics.csv"
    shell:
        "python scripts/extract_valid_sniffles_SV_for_snakemake.py {input.vcf} {output.valid} {output.statistics}"

rule correct_vcf:
    input:
        vcf="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid.vcf",
        ref=config["Files"]["ref"]
    output:
        corrected="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.vcf"
    shell:
        "python scripts/correct_sniffles_vcf.py {input.vcf} {input.ref} {output.corrected}"
rule correct_header_SV:
    input:
     "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.vcf"
    output:
     "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.vcf"
    params:
     col_name="Toxo_ME49_{sample}"
    shell:
     """
     awk 'BEGIN {{OFS="\t"}} /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/ {{$NF = "{params.col_name}"}} 1' {input} > {output}
     """
rule sort_vcf_SV:
    input: 
      "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.vcf"
    output:
     "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf"
    shell:
     "bcftools sort {input} > {output}"     
rule bgzip_concat_vcf_sniffles:
  input: 
    "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf"
  output:
    "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf.gz"
  shell:
    "bgzip {input}"

rule index_concat_vcf_sniffles:
  input:
    "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf.gz"
  output:
    "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf.gz.tbi" 
  shell:
    "tabix -p vcf {input}"     