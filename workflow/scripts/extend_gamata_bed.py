import sys

def extend_bed(input_file, output_file):
    """
    Extends BED file regions by one base pair to the left and right.

    Parameters:
        input_file (str): Path to the input BED file.
        output_file (str): Path to the output BED file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip empty lines or lines that start with a comment
            if not line.strip() or line.startswith("#"):
                continue

            # Split the line into columns
            cols = line.strip().split()
            chrom = cols[0]
            start = max(0, int(cols[1]) - 2)  # Ensure start is not negative
            end = int(cols[2]) + 2

            # Write the extended region to the output file
            outfile.write(f"{chrom}\t{start}\t{end}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extend_bed.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extend_bed(input_file, output_file)

