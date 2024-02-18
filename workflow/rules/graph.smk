
rule Graph:
    input:
     expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_sequence_to_graph_alignment_newGA.gam", sample=config["samples"])
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
    shell:
     """
     GraphAligner -g {input.graph} -f {input.reads} -t 8 -a {output.alignment} -x vg --precise-clipping 0.75
     """

