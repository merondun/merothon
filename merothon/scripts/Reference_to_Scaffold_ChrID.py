import argparse
import re
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Identify the most likely scaffold for each chromosome based on alignment data.")
    parser.add_argument("--paf", required=True, help="File containing alignment data.")
    parser.add_argument("--fai", required=True, help="FAI file containing scaffold lengths.")
    parser.add_argument("--out", required=True, help="Output file to write the scaffold to chromosome mapping.")
    parser.add_argument("--min_size", type=float, default=5.0, help="Minimum draft scaffold size to consider in Mb (default 5.0)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    alignment_file = args.paf
    fai_file = args.fai
    output_file = args.out
    min_size = args.min_size * 1e6  # Convert Mb to bp

    # Read scaffold lengths from the .fai file
    scaffold_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            scaffold = parts[0]
            length = int(parts[1])
            scaffold_lengths[scaffold] = length

    # Dictionary of total aligned length for each scaffold-chromosome pair
    alignment_lengths = defaultdict(lambda: defaultdict(int))
    strand_counts = defaultdict(lambda: defaultdict(lambda: {'+': 0, '-': 0}))

    # Read aln file and process each line
    with open(alignment_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            scaffold = fields[0]
            scaffold_start = int(fields[2])
            scaffold_end = int(fields[3])
            chromosome = fields[5]
            strand = fields[4]
            aligned_length = scaffold_end - scaffold_start

            # Update total aligned length for scaffold-chromosome pair
            alignment_lengths[scaffold][chromosome] += aligned_length
            strand_counts[scaffold][chromosome][strand] += aligned_length

    # Determine most related chromosome and predominant strand for each scaffold
    scaffold_to_chromosome = {}
    for scaffold, chrom_lengths in alignment_lengths.items():
        if scaffold_lengths[scaffold] > min_size:  # only consider scaffolds greater than min_size
            total_length = scaffold_lengths[scaffold]
            best_chromosome = max(chrom_lengths, key=chrom_lengths.get)
            best_percentage = (chrom_lengths[best_chromosome] / total_length) * 100
            predominant_strand = '+' if strand_counts[scaffold][best_chromosome]['+'] >= strand_counts[scaffold][best_chromosome]['-'] else '-'
            scaffold_to_chromosome[scaffold] = (best_chromosome, best_percentage, total_length, predominant_strand)

    # Write to the output file
    with open(output_file, 'w') as f:
        for scaffold, (chromosome, percentage, length, strand) in scaffold_to_chromosome.items():
            f.write(f"{scaffold}\t{chromosome}\t{percentage:.2f}%\t{length}\t{strand}\n")

    print(f"Scaffold to chromosome mapping has been written to {output_file}")

if __name__ == "__main__":
    main()