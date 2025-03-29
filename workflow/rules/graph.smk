num_samples = len(config["samples"])
def get_sniffles_input(sample):
    if config["manually_filtered_SV"].get(sample, False):
        return config["manual_vcf_paths"][sample]
    else:
        return f"{config['output_dir']}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID_sorted.vcf"
def get_sniffles_input2(sample):
    if config["manually_filtered_SV"].get(sample, False):
       return f"{config['output_dir']}/analysis/Sniffles/{sample}/sniffles_{sample}_manual.newHead_assignedID_sorted.vcf.gz"
    else:
       return f"{config['output_dir']}/analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID_sorted.vcf.gz"        

rule Graph:
  input: 
    expand("{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles.vcf", out=config["output_dir"],sample=config["samples"]),
    expand("{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf", out=config["output_dir"]),
    #"../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015.vcf",
    expand("{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_graph_calling_stats.txt",out=config["output_dir"], sample=config["graph_regenoyping_additional"]),
    #expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vep", sample=config["graph_regenoyping_additional"]),
    
rule assign_id_SV_manual:
  input:
      sniffles= lambda wildcards: get_sniffles_input(wildcards.sample),
  output:
      uncompressed="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_manual.newHead_assignedID.vcf",
      sorted = "{out}/analysis/Sniffles/{sample}/sniffles_{sample}_manual.newHead_assignedID_sorted.vcf",
      compressed="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_manual.newHead_assignedID_sorted.vcf.gz",
      idex="{out}/analysis/Sniffles/{sample}/sniffles_{sample}_manual.newHead_assignedID_sorted.vcf.gz.tbi"
  params:
     col_name="Toxo_ME49_{sample}"    
  shell:
      """
      cat {input.sniffles} | python3 scripts/PanGenie_scripts/assign-variant-ids.py > {output.uncompressed}
      bcftools sort {output.uncompressed} > {output.sorted}
      bgzip -c {output.sorted} > {output.compressed}
      tabix -p vcf {output.compressed}
      """                   
rule filter_overlaps_in_medaka:
  input:     
        medaka="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID.vcf", 
        sniffles= lambda wildcards: get_sniffles_input(wildcards.sample),
        bed=config["Files"]["positions_to_mask_small_variants_in"]
  output:
        vcf="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles.vcf",  
        compressed = "{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles.vcf.gz",
        idex="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles.vcf.gz.tbi",
        stat="{out}/analysis/medaka/medaka_{sample}/variants_excluded.txt",        
  shell:
      """          
      python3 scripts/remove_sniffles_overlap_from_medaka.py {input.medaka} {input.sniffles} {input.bed} {output.stat} {output.vcf}
      bgzip -c {output.vcf} > {output.compressed}
      tabix -p vcf {output.compressed}
      """
rule merge_medaka_sniffles_sample:
  input: 
    medaka_file="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles.vcf.gz",
    sniffles_file= lambda wildcards: get_sniffles_input2(wildcards.sample),
  output: 
    concat="{out}/analysis/Graph/{sample}/concat_medaka_sniffles.vcf",
  threads: 10
  shell:
    """
    bcftools concat -a --rm-dups all -o {output.concat} {input.medaka_file} {input.sniffles_file}
    """
rule merge_concat_vcf:
  input:
    concat_input=expand("{out}/analysis/Graph/{sample}/concat_medaka_sniffles.vcf",out=config["output_dir"], sample=config["samples"]),
    reference=config["Files"]["ref"]  
  threads: 10 
  output:
    uncompressedoutfile= "{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf",
    std_error="{out}/analysis/Graph/merge_concat_medaka_sniffles/std_error.txt",
    compressedoutfile="{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz",
    tabix="{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz.tbi"
  shell:
    """
    python scripts/PanGenie_scripts/merge_vcfs_sniffles.py merge -r {input.reference} -vcf {input.concat_input} -ploidy 1 > {output.uncompressedoutfile} 2> {output.std_error} 
    bgzip -c {output.uncompressedoutfile} > {output.compressedoutfile}
    tabix -p vcf {output.compressedoutfile}
    """ 
rule construct_graph_vg:
  input:
    ref=config["Files"]["ref"],
    infile="{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz",
    indexed_infile="{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz.tbi"
  output:
    outfile="{out}/analysis/Graph/graph_construction/merged_graph.vg"
  threads: 10
  shell: 
    """
    vg construct -r {input.ref} -v {input.infile} -a > {output.outfile}
    """     
rule vg_index:
  input:
    graph="{out}/analysis/Graph/graph_construction/merged_graph.vg"
  output:
    xg_graph="{out}/analysis/Graph/graph_construction/merged_graph.xg"
  threads: 10 
  shell:
    """
    vg index -t2 -x {output.xg_graph} -L {input.graph}
    """     
rule graphAligner:
  input:
    graph="{out}/analysis/Graph/graph_construction/merged_graph.vg",
    reads=lambda wildcards: config["graph_regenoyping_additional"][wildcards.sample]
  output:
    alignment="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam"
  threads: 10
  shell:
    """
    GraphAligner -g {input.graph} -f {input.reads} -t 8 -a {output.alignment} -x vg --precise-clipping 0.75
    """

rule graphAlignment_stats:
  input:
    gam="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam"
  output:
    gam_stats="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_graphAlignment_stats.txt"
  threads: 10  
  shell:
    """
    vg stats -a {input.gam} > {output.gam_stats}
    """

rule vg_filter_gam_Q30:
  input: 
    alignment="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam"
  output:
    filtered_alignment="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_MQ30.gam"
  shell:
    """
    vg filter -q 30 {input.alignment} > {output.filtered_alignment}
    """

rule vg_pack:
  input:
    graph="{out}/analysis/Graph/graph_construction/merged_graph.xg",
    alignment="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_MQ30.gam"
  output:
    pack="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_MQ30_BQ20.pack"
  shell:
    """
    vg pack -x {input.graph} -g {input.alignment} -Q 20 -o {output.pack}
    """

rule vg_call:
  input:
    graph="{out}/analysis/Graph/graph_construction/merged_graph.xg",
    pack="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_MQ30_BQ20.pack",
    vcf= "{out}/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf.gz"
  output:
    calls="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20.vcf"
  shell:
    """
    vg call {input.graph} -k {input.pack} -v {input.vcf} -d 1 > {output.calls}
    """

rule add_vartype_VG_calls:
  input: 
    "{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20.vcf"
  output:
    vcf="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.vcf",
    sort="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_sort.vcf",
    gz="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_sort.vcf.gz",
  params:
    snpEff_path=config["tools"]["snpEff_path"]
  shell:
    """
    java -jar {params.snpEff_path} varType {input} > {output.vcf}
    bcftools sort {output.vcf} > {output.sort}
    bgzip -c {output.sort} > {output.gz}
    tabix  {output.gz}
    """
               
rule filter_vg_calls:
  input:
    vcf="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype.vcf",
  output:
    filter_vcf="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter.vcf",
    graph_calling_stats="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_graph_calling_stats.txt",
    mod_filter_vcf="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod.vcf",
    total_filtered="{out}/analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_total_filtered.vcf"
  params:
     depth=config["Files"]["depth_graph"],
     allele_freq=config["Files"]["allele_freq_graph"],
     exclude_contig=config["Files"]["sequences_to_exclude"] 
  shell:
    """
    echo "Running: python scripts/snake_cl_filter_count_VG_variants_GTVCF.py {input.vcf} {output.filter_vcf} {output.mod_filter_vcf} {output.graph_calling_stats} {output.total_filtered} {params.exclude_contig} {params.depth} {params.allele_freq}"
    python scripts/snake_cl_filter_count_VG_variants_GTVCF.py {input.vcf} {output.filter_vcf} {output.mod_filter_vcf} {output.graph_calling_stats} {output.total_filtered} {params.exclude_contig} {params.depth} {params.allele_freq}
    """

rule sort_filtered_vcf:
  input:
    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod.vcf"
  output:
    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vcf"
  shell:
    """
    bcftools sort {input} > {output}
    """

rule zip_sorted_vcf:
  input:
    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vcf"
  output:
    "../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vcf.gz"
  shell:
    """
    bgzip -c {input} > {output}
    tabix -p vcf {output}
    """
rule make_result_table:
  input: 
      merged="../analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf",
      vg_mod=expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_total_filtered.vcf",sample=config["samples"]),
  output:
      out_table="../analysis/Graph/graph_construction/results/merged_vg_combined_table.vcf"
  shell:
      "python scripts/snake_result_table_merged_vg.py {input.merged} {input.vg_mod} {output.out_table}"


rule reference_filtering_merged:
  input:
    combined_table="../analysis/Graph/graph_construction/results/merged_vg_combined_table.vcf",
    medaka=expand("../analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_new_head.assignedID.vcf",sample=config["samples"]),
    sniffles=expand("../analysis/Sniffles/{sample}/sniffles_{sample}_with_reference_corrected.newHead_assignedID.vcf",sample=config["samples"]),

  output:
    combined_table_out="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref.vcf",
    for_igv="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf"
  params:
    num_samples = num_samples,
    samples = "\t".join(config["samples"].keys()) 
  shell:
    """ 
    python3 scripts/snake_medaka_sniffles_reference_filtering.py {input.combined_table} {params.num_samples} {params.samples} {input.medaka} {input.sniffles} {input.agp_file} {output.combined_table_out}
    awk 'BEGIN {{FS=OFS="\t"}} {{for (i=1; i<=NF; i++) if ($i == "-") $i="."; print}}' {output.combined_table_out} > {output.for_igv}    
    """
rule filter_by_2015:
     input:
          for_igv="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf",
          mod_filter_vcf="../analysis/Graph/graph_construction/2015T_graph_Alignment/2015T_variants_MQ30_BQ20_vartype_filter_mod.vcf", 
          ref="../analysis/Graph/graph_construction/2015T_graph_Alignment/trial_with_filtering_with_0.1/loci_to_be_filtered_from_ref.vcf"      
     output:
          for_igv_2015="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv_no_2015.vcf",
          gz="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015.vcf.gz"
     shell:
          """
          grep '^#' {input.for_igv} > {output.for_igv_2015}
          bedtools intersect -v -a {input.for_igv} -b  {input.mod_filter_vcf} {input.ref} >> {output.for_igv_2015}
          bgzip -c {output.for_igv_2015} > {output.gz}
          tabix {output.gz}
          """         
          
#rule plot_venn_per_file:
#    input:
#        resulttable="../analysis/Graph/graph_construction/results/merged_vg_combined_table.csv"

#    output:
#        merged_png="../analysis/Graph/graph_construction/results/merged_venn_diagramm.png",
#        vg_png="../analysis/Graph/graph_construction/results/vg_venn_diagramm.png"
#    shell: 
#        "python scripts/snake_venn_diagramm_perFile.py {input.resulttable} {output.merged_png} {output.vg_png}"


rule variant_annotation_graph:
  input:
    val="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vcf",
    gff=config["Files"]["gff"],
    ref=config["Files"]["ref"]
  output:
    vep="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vep" 
  threads: 4    
  shell:  
   """
   vep -i {input.val} -gff {input.gff} --fasta {input.ref} -o {output.vep}
   """
