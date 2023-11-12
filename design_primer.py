from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import argparse
from Bio import SeqIO
import pandas as pd

def create_parser():
    parser = argparse.ArgumentParser(
        description="Design primers for chosen aa point mutants to exchange amino acids and calculate melting temperature with nearest-neighbor method"
    )
    parser.add_argument(
        "fasta_file",
        help="Nuclotide sequence to change",
    )
    parser.add_argument(
        "mutations",
        help="Plain txt file of mutations seperated by \n. Note the mutation should have the format WTidxAAnew e.g. if you want to change Q to R at position 3: Q3R",
    )
    parser.add_argument(
        "offset_5'",
        type = int,
        help="Length of primer in 5' direction",
    )
    parser.add_argument(
        "offset_3'",
        type = int,
        help="Length of primer in 5' direction",
    )
    parser.add_argument(
        "prefix",
        help="Prefix in primer to order table",
    )
    return parser

def substitute_aa(rbsc_seq, aa_pos, new_aa):
            
        
        codon_table = {'R':'CGT','H':'CAT','K':'AAA','D':'GAT','E':'GAA',
                    'S':'AGC','T':'ACC','N':'AAT','Q':'CAG','C':'TGC',
                    'G':'GGC','P':'CCG','A':'GCG','V':'GTG','I':'ATT',
                    'L':'CTG','M':'ATG','F':'TTT','Y':'TAT','W':'TGG'}
        dna_pos = (aa_pos - 1) * 3

        new_codon = codon_table[new_aa]

        new_rbsc_seq = rbsc_seq[:dna_pos] + new_codon + rbsc_seq[dna_pos + 3:]
        start = max(0, dna_pos - 18)
        end = min(len(rbsc_seq), dna_pos + 18)
        surrounding_seq = new_rbsc_seq[start:end]
        fwd_primer = Seq(surrounding_seq)
        rev_primer = fwd_primer.reverse_complement()
        
        return fwd_primer, rev_primer

def main(args):
    nclt_seq = str(list(SeqIO.parse(args.fasta_file, "fasta"))[0].seq)
    selected_mutants = [i.replace('\n','') for i in open(args.mutations).readlines()]
    primer_to_order = pd.DataFrame(columns  = ['NAME_FWD','FWD', 'NAME_REV','REV', 'TM_fwd', 'TM_rev'])
    for idx, mut in enumerate(selected_mutants):
        fwd_primer, rev_primer = substitute_aa(nclt_seq, int(mut[1:-1]), mut[-1])
        print(args.prefix+'_'+mut+str(idx),": Forward primer (5' -> 3'):", fwd_primer,f"Melting Temperature: {round(mt.Tm_NN(fwd_primer))} °C", "; Reverse primer(5' -> 3'):", rev_primer, f"Melting Temperature: {round(mt.Tm_NN(rev_primer))} °C\n")
        primer_to_order.loc[idx,:] = [args.prefix+'_'+mut+'_'+'fwd', str(fwd_primer), args.prefix+'_'+mut+'_'+'rev',str(rev_primer), round(mt.Tm_NN(fwd_primer)), round(mt.Tm_NN(rev_primer))]
    primer_to_order.to_csv(args.prefix+'_'+'primer_to_order.tsv', sep='\t')

if __name__ =="__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)