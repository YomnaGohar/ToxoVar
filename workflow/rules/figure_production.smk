ruleorder: statistics_from_vep > vep > zip_sorted_vcf > sort_separated_vcf > separate_vcf
ruleorder: getting_statistics_from_vep > overlapping_vep > overlapping_var
rule  generate_figures:
    input:
        expand("../analysis/Graph/graph_construction/variation_overlap/{var}/variant_consequences_in_coding_areas.txt", var=["SNP","INS","DEL"]),
        expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{var}.variant_consequences_in_coding_areas.txt", var=["SNP","INS","DEL"],sample=config["samples"])
rule separate_vcf:
    input:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter.vcf"
    output:
        var="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vcf"
    params: 
        vartype=lambda wildcards: wildcards.var.upper()  # Ensure the variable type is uppercase to match VCF conventions
    shell:
        """
        bcftools view -i 'INFO/VARTYPE="{params.vartype}" && (strlen(ALT)-strlen(REF))<=50' {input} -o {output.var}
        """
rule sort_separated_vcf:
    input:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vcf"
    output:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vcf"
    shell:
        """
        bcftools sort {input} > {output}
        """
rule zip_sorted_vcf:
    input:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vcf"
    output:
        "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vcf.gz"
    shell:
        """
        bgzip -c {input} > {output}
        tabix {output}
        """
rule vep:
    input:
      var="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vcf",
      gff=config["Files"]["gff"],
      ref=config["Files"]["ref"]
    output:
      vep="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vep",
      stat="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vep_summary.txt"
    params:
      vep_path=config["tools"]["Ensemble_vep"]    
    shell:
     """
     {params.vep_path} -i {input.var} -gff {input.gff} --fasta {input.ref} --stats_text -o {output.vep}
     """
rule statistics_from_vep:
    input:
     stat="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.sorted.vep_summary.txt"
    output:
     all_cons="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{var}.variant_consequences.all.txt",
     coding="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{var}.variant_consequences_in_coding_areas.txt"
    shell:
     """
     python3 scripts/vep.statistics.py {input.stat} {output.all_cons} {output.coding}
     """        
rule overlapping_var:
    input:
     var=lambda wildcards: [f"../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{wildcards.var}.sorted.vcf.gz" for sample in config["samples"]]
    output:
        var_ov=directory("../analysis/Graph/graph_construction/variation_overlap/{var}")
    shell:
        """
        bcftools isec -n3 -w1 -p {output.var_ov} {input.var}
        """        
rule overlapping_vep:
    input:
      var=lambda wildcards: f"../analysis/Graph/graph_construction/variation_overlap/{wildcards.var}/0000.vcf",
      gff=config["Files"]["gff"],
      ref=config["Files"]["ref"]
    output:
      vep="../analysis/Graph/graph_construction/variation_overlap/{var}/0000.vep",
      stat="../analysis/Graph/graph_construction/variation_overlap/{var}/0000.vep_summary.txt"
    params:
      vep_path=config["tools"]["Ensemble_vep"]    
    shell:
     """
     {params.vep_path} -i {input.var} -gff {input.gff} --fasta {input.ref} --stats_text -o {output.vep}
     """
rule getting_statistics_from_vep:
    input:
     stat=lambda wildcards: f"../analysis/Graph/graph_construction/variation_overlap/{wildcards.var}/0000.vep_summary.txt"
    output:
     all_cons="../analysis/Graph/graph_construction/variation_overlap/{var}/variant_consequences.all.txt",
     coding="../analysis/Graph/graph_construction/variation_overlap/{var}/variant_consequences_in_coding_areas.txt"
    shell:
     """
     python3 scripts/vep.statistics.py {input.stat} {output.all_cons} {output.coding}
     """