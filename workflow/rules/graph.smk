
rule Graph:
    input:
      "../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref.csv"
     #expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_graph_calling_stats.txt", sample=config["samples"]) 

rule assign_id_SV:
    input:
        "../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.vcf"
    output:
        uncompressed="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.assignedID.vcf",
        compressed="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.assignedID.vcf.gz"
    shell:
        """
        cat {input} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
        bgzip -c {output.uncompressed} > {output.compressed}
        tabix -p vcf {output.compressed}
        """

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
         
rule merge_medaka_sniffles_sample:
    input: 
     medaka_file="../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_Qual_1_valid.newHead.sorted.assignedID.vcf.gz",
     sniffles_file="../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_valid_corrected.newHead.sorted.assignedID.vcf.gz"
    output: 
     "../analysis/Graph/{sample}/concat_medaka_sniffles.vcf"
    shell:
     "bcftools concat -a -o {output} {input.medaka_file} {input.sniffles_file}"

rule merge_concat_vcf:
    input:
      concat_input=expand("../analysis/Graph/{sample}/concat_medaka_sniffles.vcf", sample=config["samples"]),
      reference=config["Files"]["ref"]  
    output:
     uncompressedoutfile= "../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf",
     std_error="../analysis/Graph/merge_concat_medaka_sniffles/std_error.txt",
     compressedoutfile="../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz",
     tabix="../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz.tbi"
    shell:
     """
     python scripts/PanGenie_scripts/merge_vcfs_sniffles.py merge -r {input.reference} -vcf {input.concat_input} -ploidy 1 > {output.uncompressedoutfile} 2> {output.std_error} 
     bgzip -c {output.uncompressedoutfile} > {output.compressedoutfile}
     tabix -p vcf {output.compressedoutfile}
     """     
rule construct_graph_vg:
    input:
     ref=config["Files"]["ref"],
     infile="../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz",
     indexed_infile="../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz.tbi"
    output:
     outfile="../analysis/Graph/graph_construction/merged_graph.vg"
    params:
     vg_path=config["tools"]["vg"]
    shell: 
     """
     {params.vg_path} construct -r {input.ref} -v {input.infile} -a > {output.outfile}
     """     
rule vg_index:
    input:
     graph="../analysis/Graph/graph_construction/merged_graph.vg"
    output:
     xg_graph="../analysis/Graph/graph_construction/merged_graph.xg"
    params:
     vg_path=config["tools"]["vg"]
    shell:
     """
     {params.vg_path} index -t2 -x {output.xg_graph} -L {input.graph}
     """     
rule graphAligner:
    input:
     graph="../analysis/Graph/graph_construction/merged_graph.vg",
     reads=lambda wildcards: config["samples"][wildcards.sample]
    output:
     alignment="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam"
    conda:
     config["conda_env"]["GraphAligner"]
    shell:
     """
     GraphAligner -g {input.graph} -f {input.reads} -t 8 -a {output.alignment} -x vg --precise-clipping 0.75
     """

rule graphAlignment_stats:
  input:
    gam="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam"
  output:
    gam_stats="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_graphAlignment_stats.txt"
  params:
    vg=config["tools"]["vg"]
  shell:
    """
    {params.vg} stats -a {input.gam} > {output.gam_stats}
    """

rule vg_filter_gam_Q30:
  input: 
    alignment="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam"
  output:
    filtered_alignment="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_MQ30.gam"
  params: 
    vg=config["tools"]["vg"]
  shell:
    """
    {params.vg} filter -q 30 {input.alignment} > {output.filtered_alignment}
    """

rule vg_pack:
  input:
    graph="../analysis/Graph/graph_construction/merged_graph.xg",
    alignment="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_MQ30.gam"
  output:
    pack="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_MQ30_BQ20.pack"
  params: 
    vg=config["tools"]["vg"]
  shell:
    """
    {params.vg} pack -x {input.graph} -g {input.alignment} -Q 20 -o {output.pack}
    """

rule vg_call:
  input:
    graph="../analysis/Graph/graph_construction/merged_graph.xg",
    pack="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_MQ30_BQ20.pack",
    vcf= "../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz"
  output:
    calls="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20.vcf"
  params: 
    vg=config["tools"]["vg"]
  shell:
    """
    {params.vg} call {input.graph} -k {input.pack} -v {input.vcf} -d 1 > {output.calls}
    """

rule add_vartype_VG_calls:
  input: 
    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20.vcf"
  output:
    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.vcf"
  params:
    snpEff_path=config["tools"]["snpEff_path"]
  shell:
    "java -jar {params.snpEff_path} varType {input} > {output}" 

rule filter_vg_calls:
  input:
    vcf="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.vcf",
    agp=config["Files"]["agp"] 
  output:
    filter_vcf="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter.vcf",
    graph_calling_stats="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_graph_calling_stats.txt",
    mod_filter_vcf="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod.vcf",
    total_filtered="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_total_filtered.vcf"
  shell:
    """
    python scripts/snake_cl_filter_count_VG_variants_GTVCF.py {input.vcf} {input.agp} {output.filter_vcf} {output.mod_filter_vcf} {output.graph_calling_stats} {output.total_filtered}
    """

# rule sort_filtered_vcf:
#   input:
#     "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter.vcf"
#   output:
#    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_sorted.vcf"
#   shell:
#    """
#    bcftools sort {input} > {output}
#    """

# rule zip_sorted_vcf:
#   input:
#     "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_sorted.vcf"
#   output:
#     "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_sorted.vcf.gz"
#   shell:
#     """
#     bgzip -c {input} > {output}
#     tabix -p vcf {output}
#     """

rule make_result_table:
    input: 
        merged="../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf",
        vg_2015T_mod="../analysis/Graph/graph_construction/2015T_graph_Alignment/2015T_variants_MQ30_BQ20_vartype_total_filtered.vcf",
        vg_2020T_mod="../analysis/Graph/graph_construction/2020T_graph_Alignment/2020T_variants_MQ30_BQ20_vartype_total_filtered.vcf",
        vg_2000B_mod="../analysis/Graph/graph_construction/2000B_graph_Alignment/2000B_variants_MQ30_BQ20_vartype_total_filtered.vcf"
    output:
        out_table="../analysis/Graph/graph_construction/results/merged_vg_combined_table.csv"
    shell:
        "python scripts/snake_result_table_merged_vg.py {input.merged} {input.vg_2015T_mod} {input.vg_2020T_mod} {input.vg_2000B_mod} {output.out_table}"

rule reference_filtering_merged:
    input:
      combined_table="../analysis/Graph/graph_construction/results/merged_vg_combined_table.csv",
      medaka_2015T="../analysis/medaka/medaka_2015T/medaka.annotated_with_VarType.vcf",
      medaka_2020T="../analysis/medaka/medaka_2020T/medaka.annotated_with_VarType.vcf",
      medaka_2000B="../analysis/medaka/medaka_2000B/medaka.annotated_with_VarType.vcf",
      sniffles_2015T="../analysis/Sniffles/2015T/sniffles_2015T_with_reference.vcf",
      sniffles_2020T="../analysis/Sniffles/2020T/sniffles_2020T_with_reference.vcf",
      sniffles_2000B="../analysis/Sniffles/2000B/sniffles_2000B_with_reference.vcf",
      agp_file=config["Files"]["agp"]

    output:
      combined_table_out="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref.csv"
    shell:
      "python3 scripts/snake_medaka_sniffles_reference_filtering.py {input.combined_table} {input.medaka_2015T} {input.medaka_2020T} {input.medaka_2000B} {input.sniffles_2015T} {input.medaka_2020T} {input.medaka_2000B} {input.agp_file} {output.combined_table_out}"

#rule plot_venn_per_file:
#    input:
#        resulttable="../analysis/Graph/graph_construction/results/merged_vg_combined_table.csv"

#    output:
#        merged_png="../analysis/Graph/graph_construction/results/merged_venn_diagramm.png",
#        vg_png="../analysis/Graph/graph_construction/results/vg_venn_diagramm.png"
#    shell: 
#        "python scripts/snake_venn_diagramm_perFile.py {input.resulttable} {output.merged_png} {output.vg_png}"
