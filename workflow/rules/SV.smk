rule call_SVs:
  input:
    expand("{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID_sorted.vcf.gz.tbi",out=config["output_dir"], sample=config["samples"]) ,
    expand("{out}/analysis/Sniffles/{sample}/statistics.csv",out=config["output_dir"], sample=config["samples"]) ,
rule map:
     input:
        ref=config["Files"]["ref"],
        concat=lambda wildcards: config["samples"][wildcards.sample]
     output:
        sam="{out}/analysis/Sniffles/{sample}/calls_to_ref.sam",
        bam="{out}/analysis/Sniffles/{sample}/calls_to_ref.bam",
        sort="{out}/analysis/Sniffles/{sample}/calls_to_ref_sorted.bam",
        index="{out}/analysis/Sniffles/{sample}/calls_to_ref_sorted.bam.bai"
     shell:
        """
        minimap2 -ax map-ont  {input.ref} {input.concat} > {output.sam}
        samtools view -bS {output.sam} > {output.bam}
        samtools sort {output.bam} > {output.sort}
        samtools index {output.sort}
        """                            
rule sniffles:
    input:
        bam="{out}/analysis/Sniffles/{sample}/calls_to_ref_sorted.bam",
        ref=config["Files"]["ref"]
    output:
        "{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference.vcf"
    log:
        "{out}/analysis/logs/sniffles/log_{sample}_sniffles"
    shell:
        "sniffles --input {input.bam} --reference {input.ref} --vcf {output} > {log} 2>&1"
       
rule exclude_chromosomes_sniffles:
    input:
        vcf="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference.vcf"
    output:
        vcf="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_no_api.vcf"
    params:
        exclude_string=lambda wildcards: " || ".join([f'CHROM="{seq}"' for seq in config["Files"]["sequences_to_exclude"].split(",")])
    shell:
        """
        bcftools view --exclude '{params.exclude_string}' {input.vcf} -o {output.vcf}
        """        
rule correct_vcf:
    input:
        vcf="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_no_api.vcf",
        ref=config["Files"]["ref"]
    output:
        corrected="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.vcf"
    shell:
        "python scripts/correct_sniffles_vcf.py {input.vcf} {input.ref} {output.corrected}"
rule filter_vcf:
    input:
        vcf="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.vcf"
    output:
        statistics="{out}/analysis/Sniffles/{sample}/statistics.csv"
    shell:
        "python scripts/sniffles_count.py {input.vcf} {output.statistics}"        

rule assign_id_SV_all:
  input:
      "{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.vcf"
  output:
      corrected="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected_newHead.vcf",
      uncompressed="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID.vcf",
      sorted = "{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID_sorted.vcf",
      compressed="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID_sorted.vcf.gz",
      idex="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID_sorted.vcf.gz.tbi"
  params:
     col_name="Toxo_ME49_{sample}"    
  shell:
      """
      awk 'BEGIN {{OFS="\t"}} /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/ {{$NF = "{params.col_name}"}} 1' {input} > {output.corrected}
      cat {output.corrected} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
      bcftools sort {output.uncompressed} > {output.sorted}
      bgzip -c {output.sorted} > {output.compressed}
      tabix {output.compressed}
      """ 
