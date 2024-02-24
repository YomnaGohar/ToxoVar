rule generate_figures:
    input:
     expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/variant_consequences_in_coding_areas.{var}.txt", sample=config["samples"], var=["SNP","INS","DEL"])
rule separate_vcf:
    input:
      "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.vcf"
    output:
      "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vcf"
    params: 
     vartype=lambda wildcards: wildcards.var.upper()  # Ensure the variable type is uppercase to match VCF conventions
    shell:
     """
     bcftools view -i 'INFO/VARTYPE="{params.vartype}" && (strlen(ALT)-strlen(REF))<=50' {input} -o {output}
     """
rule separate_vep:
    input:
      var="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vcf",
      gff=config["Files"]["gff"],
      ref=config["Files"]["ref"]
    output:
      vep="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vep",
      stat="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vep_summary.txt"
    params:
      vep_path=config["tools"]["Ensemble_vep"]    
    shell:
     """
     {params.vep_path} -i {input.var} -gff {input.gff} --fasta {input.ref} --stats_text -o {output.vep}
     """
rule getting_statistics_from_vep:
    input:
     stat="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{var}.vep_summary.txt"
    output:
     all_cons="../analysis/Graph/graph_construction/{sample}_graph_Alignment/all_variant_consequences.{var}.txt",
     coding="../analysis/Graph/graph_construction/{sample}_graph_Alignment/variant_consequences_in_coding_areas.{var}.txt"
    params: 
     sample_name=lambda wildcards: wildcards.sample
    shell:
     """
     python3 scripts/vep.statistics.py {input.stat} {params.sample_name} {output.all_cons} {output.coding}
     """