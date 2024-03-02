ruleorder: statistics_from_vep > vep > zip_sorted_vcf > sort_separated_vcf > separate_vcf
ruleorder: getting_statistics_from_vep > overlapping_vep > overlapping_var
rule generate_figures:
    input:
        expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{var}.{length}.variant_consequences_in_coding_areas.txt", var=["SNP","INS","DEL"],sample=config["samples"],length=["less_than_50","greater_than_50"]),
        expand("../analysis/Graph/graph_construction/{var}.{length}.variant_consequences.all.pdf", var=["SNP","INS","DEL"],length=["less_than_50","greater_than_50"]),
        expand("../analysis/Graph/graph_construction/{var}.{length}.variant_consequences_in_coding_areas.pdf", var=["SNP","INS","DEL"],length=["less_than_50","greater_than_50"])
rule separate_vcf:
    input:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter.vcf"
    output:
        var="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.vcf"
    params: 
        vartype=lambda wildcards: wildcards.var.upper(),  
        length=lambda wildcards: wildcards.length
    run:
        comparison_operator = "<" if params.length == "less_than_50" else ">="
        shell(
            """
            bcftools view -i 'INFO/VARTYPE="{params.vartype}" && (strlen(ALT)-strlen(REF)){comparison_operator}50' {input} -o {output.var}
            """
        )
rule sort_separated_vcf:
    input:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.vcf"
    output:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vcf"
    shell:
        """
        bcftools sort {input} > {output}
        """
rule zip_sorted_vcf:
    input:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vcf"
    output:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vcf.gz"
    shell:
        """
        bgzip -c {input} > {output}
        tabix {output}
        """
rule vep:
    input:
      var="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vcf",
      gff=config["Files"]["gff"],
      ref=config["Files"]["ref"]
    output:
      vep="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vep",
      stat="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vep_summary.txt"
    params:
      vep_path=config["tools"]["Ensemble_vep"]    
    shell:
     """
     {params.vep_path} -i {input.var} -gff {input.gff} --fasta {input.ref} --stats_text -o {output.vep}
     """
rule statistics_from_vep:
    input:
     stat="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.{length}.sorted.vep_summary.txt"
    output:
     all_cons="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{var}.{length}.variant_consequences.all.txt",
     coding="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{var}.{length}.variant_consequences_in_coding_areas.txt"
    shell:
     """
     python3 scripts/vep.statistics.py {input.stat} {output.all_cons} {output.coding}
     """        
rule figure:
    input:
     f=lambda wildcards: [f"../analysis/Graph/graph_construction/{sample}_graph_Alignment/{wildcards.var}.{wildcards.length}.variant_consequences.all.txt" for sample in config["samples"]]
    output:
     o= "../analysis/Graph/graph_construction/{var}.{length}.variant_consequences.all.pdf"
    shell:
     """
     python3 scripts/variant_barplot.py {output.o} {input.f} 
     """     
rule figure_in_coding_regions:
    input:
     f=lambda wildcards: [f"../analysis/Graph/graph_construction/{sample}_graph_Alignment/{wildcards.var}.{wildcards.length}.variant_consequences_in_coding_areas.txt" for sample in config["samples"]]
    output:
     o= "../analysis/Graph/graph_construction/{var}.{length}.variant_consequences_in_coding_areas.pdf"
    shell:
     """
     python3 scripts/variant_barplot.py {output.o} {input.f} 
     """        