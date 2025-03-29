import sys
input_bed= sys.argv[1] #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2000B_graph_Alignment/2000B_variants_MQ30_BQ20_vartype_filter_mod_homo_after_more_filtering_500_window.bed"
output_script = sys.argv[2] #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2000B_graph_Alignment/screenshot_homo_after_filtering.bat"
vcf= sys.argv[3] #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf"
genome_path = sys.argv[4] #"/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/scafs.fasta"  # Update this path
bam_file_1 = "/gpfs/project/yogah100/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/mapped_reads_to_scafs_sorted.bam"
# Update this path
numts= sys.argv[5] #"/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/numts/scafs.fasta.bed"
tand= sys.argv[6] #"/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/homopolymer/tandem.bed"
homo= sys.argv[7] #"/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/homopolymer/homopolymer.bed"
gtf= sys.argv[8] #"/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/scaffold.out.modified.with.gffread.gff3",
output_dir =  sys.argv[9] # "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2000B_graph_Alignment/screenshot_homo_after_filtering/"
medaka= sys.argv[10]
bam_file_= sys.argv[11:] #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/medaka/medaka_2000B/calls_to_ref.bam"


with open(input_bed, "r") as bed, open(output_script, "w") as script:
    # Initialize IGV session
    script.write(f"new\n")
    script.write(f"genome {genome_path}\n")
    script.write(f"load {vcf}\n")
    script.write(f"load {medaka}\n")
    script.write(f"load {bam_file_1}\n")
    for bam_file in bam_file_:
        script.write(f"load {bam_file}\n")
    script.write(f"load {tand}\n")
    script.write(f"load {homo}\n")
    script.write(f"load {numts}\n")
    script.write(f"load {gtf}\n")
    script.write(f"snapshotDirectory {output_dir}\n")

    for line in bed:
        chrom, start, end,label = line.strip().split("\t")
        start = max(0, int(start) )
        end = int(end) 
        script.write(f"goto {chrom}:{start}-{end}\n")
        script.write(f"snapshot {chrom}_{start}_{label}.png\n")

print(f"Batch script saved to: {output_script}")

