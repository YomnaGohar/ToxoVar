#!/usr/bin/env python3
from pysam import VariantFile
import sys

vcf_path = sys.argv[1]
fasta_path = sys.argv[2]  # Output FASTA file

with VariantFile(vcf_path) as vcf_in, open(fasta_path, 'w') as fasta_out:
    for rec in vcf_in.fetch():
        # Check if the variant is an insertion
        if len(rec.ref) < len(rec.alts[0]):  # Simple check for insertion
            seq = rec.alts[0]  # Assuming single ALT allele; adjust as necessary
            header = f">{rec.chrom}_{rec.pos}_insertion"
            fasta_out.write(f"{header}\n{seq}\n")

