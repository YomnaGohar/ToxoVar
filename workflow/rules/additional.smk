rule vcf_analysis:
     input:
          expand("{out}/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015.vcf", out=config["output_dir"]),
          "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed",


rule filter_by_2015:
     input:
          for_igv="{out}/analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf",
          mod_filter_vcf="/gpfs/project/yogah100/analysis/Graph/graph_construction/2015T_graph_Alignment/2015T_variants_MQ30_BQ20_vartype_filter_mod.vcf", 
          ref="/gpfs/project/yogah100/analysis/Graph/graph_construction/2015T_graph_Alignment/filter_for_ref/loci_to_be_filtered_from_ref.vcf"      
     output:
          for_igv_2015="{out}/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015.vcf",
          gz="{out}/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015.vcf.gz"
     shell:
          """
          grep '^#' {input.for_igv} > {output.for_igv_2015}
          bedtools intersect -v -a {input.for_igv} -b  {input.mod_filter_vcf} {input.ref} >> {output.for_igv_2015}
          bgzip -c {output.for_igv_2015} > {output.gz}
          tabix {output.gz}
          """         
rule extend_bed:
    input:
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta.bed",
        homopolymer= "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted.bed"
    output:
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_ex.bed",
        homopolymer_e= "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed"
    shell:
        """
        python scripts/extend_gamata_bed.py {input.tandem} {output.tandem}
        python scripts/extend_bed.py {input.homopolymer} {output.homopolymer_e}
        """


