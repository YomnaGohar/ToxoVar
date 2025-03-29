rule vcf_analysis:
     input:
          expand("../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted_tandem.vcf", sample=config["graph_regenoyping_additional"]),
          "../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_no_low_complexity.vcf",
          "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed",
          "../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem2.vcf",
          "../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized_vep_table.csv",
          #"../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized_vep_table.csv"

rule extend_bed:
    input:
        #tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/Gamata_less_than_4_no_JACEHA010000016.1_contig_1_contig_50_sorted.bed",
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta.bed",
        tandem2="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_2.bed",
        tandem3="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_3.bed",
        tandem4="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_4.bed",
        tandem5="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_greater_than_5.bed",
        homopolymer= "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted.bed"
    output:
        #tandem_e="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/Gamata_less_than_4_no_JACEHA010000016.1_contig_1_contig_50_sorted_ex.bed",
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_ex.bed",
        tandem2="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_2_ex.bed",
        tandem3="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_3_ex.bed",
        tandem4="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_4_ex.bed",
        tandem5="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_greater_than_5_ex.bed",
        homopolymer_e= "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed"
    shell:
        """
        python scripts/extend_gamata_bed.py {input.tandem} {output.tandem}
        python scripts/extend_gamata_bed.py {input.tandem2} {output.tandem2}
        python scripts/extend_gamata_bed.py {input.tandem3} {output.tandem3}
        python scripts/extend_gamata_bed.py {input.tandem4} {output.tandem4}
        python scripts/extend_gamata_bed.py {input.tandem5} {output.tandem5}
        python scripts/extend_bed.py {input.homopolymer} {output.homopolymer_e}
        """
rule get_intersection:
    input:
        vcf="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted.vcf.gz",
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_ex.bed",
        homo="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed",
        numt="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/numts/numts_only_no_contig_1_contig_50.bed"
    output:
        vcf_no_low_complexity="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted_no_low_complexity.vcf",
        vcf_numt="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted_numt.vcf",
        vcf_homo="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted_homo.vcf",
        vcf_tandem="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_variants_MQ30_BQ20_vartype_filter_mod_sorted_tandem.vcf",
        temp="../analysis/Graph/graph_construction/{sample}_graph_Alignment/{sample}_temp_intermediate.vcf",
        header="../analysis/Graph/graph_construction/{sample}_graph_Alignment/header.txt"
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
        for_igv_gz="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv_no_2015_manual.vcf",
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_ex.bed",
        homo="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed",
        numt="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/numts/numts_only_no_contig_1_contig_50.bed"
    output:
        for_igv_gz="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_manual.vcf.gz",
        vcf_no_low_complexity="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_no_low_complexity.vcf",
        vcf_numt="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_numt.vcf",
        vcf_homo="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_homo.vcf",
        vcf_tandem="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem.vcf",
        temp="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_temp_intermediate.vcf",
        header="../analysis/Graph/graph_construction/results/with_2015T_filtering/header.txt"
    shell:
        """
        bgzip -c {input.for_igv_gz} > {output.for_igv_gz}
        tabix {output.for_igv_gz}
        zgrep '^#' {output.for_igv_gz} > {output.header}
        cat {output.header} > {output.vcf_no_low_complexity}
        bedtools intersect -v -a {output.for_igv_gz} -b {input.tandem} {input.homo} {input.numt} >> {output.vcf_no_low_complexity}
        cat {output.header} > {output.vcf_numt}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.numt} >> {output.vcf_numt}
        cat {output.header} > {output.vcf_homo}
        cat {output.header} > {output.temp}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.homo} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} >> {output.vcf_homo}
        cat {output.header} > {output.temp}
        cat {output.header} > {output.vcf_tandem}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.tandem} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} {input.homo} >> {output.vcf_tandem}
        """


rule get_intersection_with_different_tandem_length:
    input:
        vcf_tandem="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem.vcf",
        tandem2="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_2_ex.bed",
        tandem3="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_3_ex.bed",
        tandem4="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_4_ex.bed",
        tandem5="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_greater_than_5_ex.bed"
    output:
        for_igv_gz="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem.vcf.gz",
        tandem2="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem2.vcf",
        tandem3="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem3.vcf",
        tandem4="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem4.vcf",
        tandem5="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem_more_than5.vcf",
        temp="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_tandem_intermediate.vcf",
        header="../analysis/Graph/graph_construction/results/with_2015T_filtering/header.txt"
    shell:
        """
        bgzip -c {input.vcf_tandem} > {output.for_igv_gz}
        tabix {output.for_igv_gz}
        zgrep '^#' {output.for_igv_gz} > {output.header}
        cat {output.header} > {output.tandem2}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.tandem2} >> {output.tandem2}
        cat {output.header} > {output.temp}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.tandem3} >> {output.temp}
        cat {output.header} > {output.tandem3}
        bedtools intersect -v -a {output.temp} -b {input.tandem2} >> {output.tandem3}
        cat {output.header} > {output.temp}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.tandem4} >> {output.temp}
        cat {output.header} > {output.tandem4}
        bedtools intersect -v -a {output.temp} -b {input.tandem2} {input.tandem3} >> {output.tandem4}
        cat {output.header} > {output.tandem5}
        bedtools intersect -v -a {output.for_igv_gz} -b {input.tandem2} {input.tandem3} {input.tandem4} >> {output.tandem5}
        """
rule get_intersection_merged_results_without_filtering:
    input:
        igv="../analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf",
        tandem="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/Gmata2/gmata_no_apicoblasta_ex.bed",
        homo="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/gmata_and_homolymers_for_the_paper/homopolymer_l4_o0_no_contig_1_contig_50_sorted_ex.bed",
        numt="/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/numts/numts_only_no_contig_1_contig_50.bed"
    output:
        for_igv_gz="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv.vcf.gz",
        vcf_no_low_complexity="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_low_complexity.vcf",
        vcf_numt="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_numt.vcf",
        vcf_homo="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_homo.vcf",
        vcf_tandem="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_tandem.vcf",
        temp="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_temp_intermediate.vcf",
        header="../analysis/Graph/graph_construction/results/without_2015T_filtering/header.txt"
    shell:
        """
        bgzip -c {input.igv} > {output.for_igv_gz}
        tabix {output.for_igv_gz}
        zgrep '^#' {output.for_igv_gz} > {output.header}
        cat {output.header} > {output.vcf_no_low_complexity}
        bedtools intersect -v -a {output.for_igv_gz} -b {input.tandem} {input.homo} {input.numt} >> {output.vcf_no_low_complexity}
        cat {output.header} > {output.vcf_numt}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.numt} >> {output.vcf_numt}
        cat {output.header} > {output.vcf_homo}
        cat {output.header} > {output.temp}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.homo} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} >> {output.vcf_homo}
        cat {output.header} > {output.temp}
        cat {output.header} > {output.vcf_tandem}
        bedtools intersect -u -a {output.for_igv_gz} -b {input.tandem} >> {output.temp}
        bedtools intersect -v -a {output.temp} -b {input.numt} {input.homo} >> {output.vcf_tandem}
        """    
rule normalize_filter_by_2015:
     input:
          for_igv_2015="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_manual.vcf.gz",     
     output:
          for_igv_2015="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vcf",
     shell:
          """
          bcftools norm -m- {input.for_igv_2015} | bcftools annotate -x ID -I +'%INFO/ID' -o {output.for_igv_2015}          
          """
rule normalize:
     input:
          for_igv_2015="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv.vcf.gz",     
     output:
          for_igv_2015="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized.vcf",
     shell:
          """
          bcftools norm -m- {input.for_igv_2015} | bcftools annotate -x ID -I +'%INFO/ID' -o {output.for_igv_2015}      
          """          
rule normalize_filter_by_2015_vep:
     input:
          for_igv_2015="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vcf",
          gff=config["Files"]["gff"],
          ref=config["Files"]["ref"]
     output:
          vep="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vep" 
     threads: 4    
     shell:  
        """
         vep -i {input.for_igv_2015} -gff {input.gff} --pick --fasta {input.ref} -o {output.vep}
        """     
rule normalize_vep:
     input:
          for_igv_gz="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized.vcf",          
          gff=config["Files"]["gff"],
          ref=config["Files"]["ref"]
     output:
          vep="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized.vep" 
     threads: 4    
     shell:  
        """
         vep -i {input.for_igv_gz} -gff {input.gff} --fasta {input.ref} -o {output.vep}
        """             
 
rule functional_annotaions_with_2015_filter:
     input:
          for_igv_2015="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vcf",          
          vep="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vep"
     output:
          table="../analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized_vep_table.csv" 
     threads: 4    
     shell:  
        """
         python scripts/Functional_annotation.py {input.vep} --pick {input.for_igv_2015} {output.table}
        """ 
rule functional_annotaions:
     input:
          for_igv_2015="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized.vcf",          
          vep="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized.vep"
     output:
          table="../analysis/Graph/graph_construction/results/without_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_normalized_vep_table.csv" 
     threads: 4    
     shell:  
        """
         python scripts/Functional_annotation.py {input.vep} {input.for_igv_2015} {output.table}
        """ 

                       
              	    

