from Bio import SeqIO
from Bio import pairwise2
from multiprocessing import Pool
from tqdm import tqdm
import argparse

def create_parser():
    parser = argparse.ArgumentParser(
        description="Create pairwise alignments of all sequences in fasta1 vs all sequnences in fasta2. Calculates identity for best alignment of each pair and saves it to a csv file"
    )
    parser.add_argument(
        "fasta1",
        help="Path to fasta file with seqs to read",
    )
    parser.add_argument(
        "fasta2",
        help="Path to fasta file with seqs to read",
        )
    parser.add_argument(
        "out_csv",
        help="Path to resulting csv file",
        )

    return parser

def align_and_calculate_identity(pair):
    seq1, seq2 = pair
    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    best_alignment = max(alignments, key=lambda x: x.score)
    identity = best_alignment.score / max(len(seq1), len(seq2))
    return seq1.id, seq2.id, identity

def main(args):
    sequences_set1 = list(SeqIO.parse(args.fasta1, "fasta"))
    sequences_set2 = list(SeqIO.parse(args.fasta2, "fasta"))

    pairs = [(seq1, seq2) for seq1 in tqdm(sequences_set1) for seq2 in sequences_set2]

    with Pool() as pool:
        results = list(tqdm(pool.imap(align_and_calculate_identity, pairs), total=len(pairs)))

    with open(args.out_csv, "w") as out_file:
        for result in results:
            out_file.write(f"{result[0]},{result[1]},{result[2]*100}%\n")

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
