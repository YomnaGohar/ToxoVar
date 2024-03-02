rule generate_figures_overlapping:
  input:
   expand("../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/variant_consequences_in_coding_areas.txt", var=["SNP","INS","DEL"],length=["less_than_50","greater_than_50"])
rule overlapping_var:
  input:
    var=lambda wildcards: [f"../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.{wildcards.var}.{wildcards.length}.sorted.vcf.gz" for sample in config["samples"]]
  output:
    var_ov=directory("../analysis/Graph/graph_construction/variation_overlap/{var}.{length}"),
    var="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/0000.vcf"
  shell:
    """
    bcftools isec -n3 -w1 -p {output.var_ov} {input.var}
    """        
rule overlapping_vep:
  input:
    var=lambda wildcards: f"../analysis/Graph/graph_construction/variation_overlap/{wildcards.var}.{wildcards.length}/0000.vcf",
    gff=config["Files"]["gff"],
    ref=config["Files"]["ref"]
  output:
    vep="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/0000.vep",
    stat="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/0000.vep_summary.txt"
  params:
    vep_path=config["tools"]["Ensemble_vep"]    
  shell:
    """
    {params.vep_path} -i {input.var} -gff {input.gff} --fasta {input.ref} --stats_text -o {output.vep}
    """
rule getting_statistics_from_vep:
  input:
    stat=lambda wildcards: f"../analysis/Graph/graph_construction/variation_overlap/{wildcards.var}.{wildcards.length}/0000.vep_summary.txt"
  output:
    all_cons="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/variant_consequences.all.txt",
    coding="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/variant_consequences_in_coding_areas.txt"
  shell:
    """
    python3 scripts/vep.statistics.py {input.stat} {output.all_cons} {output.coding}
    """
rule figure_overlap:
  input:
    f="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/variant_consequences.all.txt"
  output:
    o= "../analysis/Graph/graph_construction/{var}.{length}/variant_consequences.all.pdf"
  shell:
    """
    python3 scripts/ {output.o} {input.f} 
    """     
rule figure_in_coding_regions_ovelap:
  input:
    f="../analysis/Graph/graph_construction/variation_overlap/{var}.{length}/variant_consequences_in_coding_areas.txt"
  output:
    o= "../analysis/Graph/graph_construction/{var}.{length}./variant_consequences_in_coding_areas.pdf"
  shell:
    """
    python3 scripts/ {output.o} {input.f} 
    """            