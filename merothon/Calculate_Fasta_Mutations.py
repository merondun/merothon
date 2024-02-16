import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate mutations in a FASTA file relative to a reference sequence.")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA file.")
    parser.add_argument("--reference", required=True, help="Reference sequence ID.")
    parser.add_argument("--out", required=True, help="Output file path.")
    return parser.parse_args()

def calculate_mutations(fasta_file, ref_seq_id, output_file_path):
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    ref_seq = next((seq for seq in sequences if seq.id == ref_seq_id), None)

    if ref_seq is None:
        print(f"Reference sequence with ID {ref_seq_id} not found.")
        return

    polymorphic_sites = set()

    # Identify polymorphic sites across all sequences
    for seq in sequences:
        for i, (ref_base, seq_base) in enumerate(zip(ref_seq.seq, seq.seq)):
            if ref_base != seq_base and seq_base not in {'N', '-'}:
                polymorphic_sites.add(i)

    with open(output_file_path, 'w') as output_file:
        output_file.write('sequence_id\tsequence_length\tnum_mutations\tnum_N_and_gaps\tnum_polymorphic\n')

        for seq in sequences:
            num_mutations = sum(1 for ref_base, seq_base in zip(ref_seq.seq, seq.seq) if ref_base != seq_base and seq_base not in {'N', '-'})
            num_N_and_gaps = sum(1 for base in seq.seq if base in {'N', '-'})
            sequence_length = len(seq.seq)

            output_file.write(f"{seq.id}\t{sequence_length}\t{num_mutations}\t{num_N_and_gaps}\t{len(polymorphic_sites)}\n")

def main():
    args = parse_args()
    calculate_mutations(args.fasta, args.reference, args.out)

if __name__ == "__main__":
    main()
