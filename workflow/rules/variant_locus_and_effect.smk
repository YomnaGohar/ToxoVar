def get_graph_input(config):
    if config.get("manually_filtered_after_graph", False):
        return config["manual_vcf_paths_after_graph"]
    else:
        return f"{config['output_dir']}/analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf.gz"    
rule vcf_analysis:
     input:
        expand("{out}/analysis/Graph/graph_construction/results/normal.vcf", out=config["output_dir"]),
        expand("{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles_normal.vcf",sample=config["samples"], out=config["output_dir"]),
        expand("{out}/analysis/Graph/graph_construction/results/normalized.vep", out=config["output_dir"]),

rule get_intersection:
    input:
        vcf="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles.vcf.gz",
        tandem= config["Files"]["TandemRepeat"],
        homo= config["Files"]["homopolymer"],
        numt= config["Files"]["NUMT"]
    output:
        vcf_no_low_complexity="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles_normal.vcf",
        vcf_numt="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles_numt.vcf",
        vcf_homo="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_snifflesd_homo.vcf",
        vcf_tandem="{out}/analysis/medaka/medaka_{sample}/medaka.annotated_with_VarType_valid.newHead.sorted.assignedID_no_overlap_with_sniffles_tandem.vcf",
        temp=temp("{out}/analysis/medaka/medaka_{sample}/temp_intermediate.vcf"),
        header=temp("{out}/analysis/medaka/medaka_{sample}/header.txt")
    shell:
        """
        zgrep '^#' {input.vcf} > {output.header}
        cat {output.header} > {output.vcf_no_low_complexity}
        bedtools intersect -v -a {input.vcf} -b {input.tandem} {input.homo} {input.numt} >> {output.vcf_no_low_complexity}
        cat {output.header} > {output.vcf_numt}
        bedtools intersect -u -a {input.vcf} -b {input.numt} >> {output.vcf_numt}
        cat {output.header} > {output.vcf_homo}
        cat {output.header} > {output.temp}
        bedtools intersect -u -a {input.vcf} -b {input.homo} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} >> {output.vcf_homo}
        cat {output.header} > {output.temp}
        cat {output.header} > {output.vcf_tandem}
        bedtools intersect -u -a {input.vcf} -b {input.tandem} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} {input.homo} >> {output.vcf_tandem}
        """
rule get_intersection_merged_results:
    input:
        for_igv_gz=get_graph_input(config),
        tandem= config["Files"]["TandemRepeat"],
        homo= config["Files"]["homopolymer"],
        numt= config["Files"]["NUMT"]
    output:
        vcf_no_low_complexity="{out}/analysis/Graph/graph_construction/results/normal.vcf",
        vcf_numt="{out}/analysis/Graph/graph_construction/results/numt.vcf",
        vcf_homo="{out}/analysis/Graph/graph_construction/results/homo.vcf",
        vcf_tandem="{out}/analysis/Graph/graph_construction/results/tandem.vcf",
        temp=temp("{out}/analysis/Graph/graph_construction/results/temp_intermediate.vcf"),
        header=temp("{out}/analysis/Graph/graph_construction/results/header.txt")
    shell:
        """
        zgrep '^#' {input.for_igv_gz} > {output.header}
        cat {output.header} > {output.vcf_no_low_complexity}
        bedtools intersect -v -a {input.for_igv_gz} -b {input.tandem} {input.homo} {input.numt} >> {output.vcf_no_low_complexity}
        cat {output.header} > {output.vcf_numt}
        bedtools intersect -u -a {input.for_igv_gz} -b {input.numt} >> {output.vcf_numt}
        cat {output.header} > {output.vcf_homo}
        cat {output.header} > {output.temp}
        bedtools intersect -u -a {input.for_igv_gz} -b {input.homo} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} >> {output.vcf_homo}
        cat {output.header} > {output.temp}
        cat {output.header} > {output.vcf_tandem}
        bedtools intersect -u -a {input.for_igv_gz} -b {input.tandem} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} {input.homo} >> {output.vcf_tandem}
        """

rule normalize_filter_by_2015:
     input:
          for_igv_2015=get_graph_input(config),     
     output:
          for_igv_2015="{out}/analysis/Graph/graph_construction/results/normalized.vcf",
     shell:
          """
          bcftools norm -m- {input.for_igv_2015} | bcftools annotate -x ID -I +'%INFO/ID' -o {output.for_igv_2015}          
          """
        
rule normalize_filter_by_2015_vep:
     input:
          for_igv_2015="{out}/analysis/Graph/graph_construction/results/normalized.vcf",
          gff=config["Files"]["gff"],
          ref=config["Files"]["ref"]
     output:
          vep="{out}/analysis/Graph/graph_construction/results/normalized.vep" 
     threads: 4    
     shell:  
        """
         vep -i {input.for_igv_2015} -gff {input.gff} --pick --fasta {input.ref} -o {output.vep}
        """     
  
              	    

