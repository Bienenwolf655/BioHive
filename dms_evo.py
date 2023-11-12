from Bio import SeqIO, AlignIO, Align
from Bio.Align.Applications import ClustalOmegaCommandline
from collections import Counter
import argparse
import numpy as np
import subprocess
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import math
from plotly.subplots import make_subplots
from tqdm import tqdm
import logomaker as lm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
from shutil import which
import os

def create_parser():
    parser = argparse.ArgumentParser(
        description="Pipeline that creates an MSA to a reference sequence (1. seq in provided fasta file)\
        and calculates the propensities for each position indexed by ref seq[off_by_n:off_by_c]. Uses the dms_enrichments \
        to provide radar/logo plots at each position that show the \
        propensity at each position and the difference enrichments - propensity therefore enrichemnts need to be \
        scaled between 0 (dead) and 1 (WT)"
    )
    parser.add_argument(
        "fasta_file",
        help="Path to fasta file with seqs to read and create the MSA ",
    )
    parser.add_argument(
        "dms_enrichment",
        help="Path to fasta file with seqs to read and create the MSA ",
    )
    parser.add_argument(
        "MSA_file",
        help="Path to which MSA file should be safed to",
    )
    parser.add_argument(
        "name_stem",
        help="To be used for the created files",
    )
    parser.add_argument(
        "off_by_n",
        type = int,
        help="How much indexing is off compared to the complete REF seq at the n-term",
    )
    parser.add_argument(
        "off_by_c",
        type = int,
        help="How much indexing is off compared to the complete REF seq at the c-term (minus indexing so if 1 shorter -1)",
    )

    return parser

def align_with_mafft(input_fasta, output_fasta):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    reference_seq = sequences[0]
    other_seqs = sequences[1:]
    
    ref_file = "reference.fasta"
    SeqIO.write([reference_seq], ref_file, "fasta")
    
    others_file = "others.fasta"
    SeqIO.write(other_seqs, others_file, "fasta")
    
    mafft_command = [
    'mafft',
    "--6merpair",
    "--thread", "-1",
    "--keeplength",
    "--addfragments", others_file,
    ref_file]

    print(' '.join(mafft_command))
    
    with open(output_fasta, "w") as output_handle:
        subprocess.run(mafft_command, stdout=output_handle)

def radarPlots(dat, file_name, ref_seq, off_by, order):
    labels = ['P', 'W', 'F', 'Y', 'I', 'L', 'M', 'V', 'A', 'G', 'S', 'T', 'C', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
    aa_color = ['mediumorchid'] + ['palegoldenrod']*3 + ['sandybrown']*4 + ['mediumseagreen']*7 + ['tomato']*2 + ['steelblue']*3
    rbsc_seq = ref_seq[off_by[0]:off_by[1]]
    num_plots = len(rbsc_seq)
    aa_label_dict = {'P':0, 'W':18, 'F':36, 'Y':54, 'I':72, 'L':90, 'M':108, 'V':126, 'A':144, 'G':162, 'S':180, 'T':198, 'C':216, 'N':234, 'Q':252, 'D':270, 'E':288, 'K':306, 'R':324, 'H':342}
    num_rows = math.ceil(num_plots / 5)
    fig = make_subplots(rows=num_rows, cols=5, subplot_titles=[f"{rbsc_seq[i]}{i+3}" for i in range(num_plots)], specs=[[{'type': 'polar'}]*5]*num_rows, horizontal_spacing=0.075)

    for i in tqdm(range(num_rows),desc=f'Creating the radar plot {file_name}'):
        for j in range(5):
            index = i * 5 + j
            if index < num_plots:
                ordered_data = [dat[index][order.index(amino)] for amino in labels]
                fig.add_trace(go.Barpolar(
                    r=ordered_data,
                    theta=labels,
                    marker_color = [aa_color[labels.index(amino)] for amino in order],
                    marker_line_color="black",
                    marker_line_width=1,
                    opacity=1,
                    name=f"Position {index+1}",
                    showlegend=False,
                ), row=i+1, col=j+1)

                fig.update_polars(radialaxis = dict(range=[-1, 1.25], tickvals=[-1, 0.5,0, .5, 1.25], ticks='', showticklabels=False), angularaxis = dict(showticklabels=True, ticks='', tickfont=dict(size=12)),radialaxis_angle = aa_label_dict[rbsc_seq[index]], row=i+1, col=j+1)   
                fig.update_layout(paper_bgcolor='rgba(0,0,0,0)')

    fig.update_layout(height=200*num_rows)

    pio.write_image(fig, file_name, format='pdf', scale=3, engine="kaleido")

def create_logoplot(matrix, file_name, off_by, amino_acids):
    matrix = np.where(matrix > 0, matrix, 20*matrix)
    segment_length = 50
    num_segments = matrix.shape[0] // segment_length

    df = pd.DataFrame(matrix, columns = amino_acids)

    with PdfPages(f'{file_name}.pdf') as pdf:
        fig, axes = plt.subplots(num_segments, 1, figsize=(20, 3*num_segments))
        
        for i in tqdm(range(num_segments), f'Creating the logo plot {file_name}'):
            start = i * segment_length
            end = start + segment_length
            segment_df = df[start:end]
            
            logo = lm.Logo(segment_df, ax=axes[i], color_scheme='NajafabadiEtAl2017', show_spines=False)
            logo.ax.set_xlabel('Position')
            logo.ax.set_ylim([-1, 20])
            logo.ax.set_yticks(np.concatenate([np.arange(-1.,0., .2)*20,np.arange(0.,25.,5.)]))
            logo.ax.set_xticks(np.arange(start,end,1))
            logo.ax.set_ylabel('Enrichment - Freq')
            logo.ax.set_yticklabels(np.concatenate([np.round(np.arange(-1.,0., .2),1),np.arange(0.,25.,5.)]))
            logo.ax.set_xticklabels(np.arange(start+off_by[0]+1,end+abs(off_by[1])+1,1))
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

def main(args):
    assert which('mafft'), "Please install mafft or move executable to name mafft: conda install -c bioconda mafft!"
    align_with_mafft(args.fasta_file, args.name_stem+'_'+args.MSA_file)
    alignment = AlignIO.read(args.name_stem+'_'+args.MSA_file, "fasta")
    num_positions = alignment.get_alignment_length()
    propensity_matrix = np.zeros((num_positions, 20))
    all_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_counts = {position: {aa: 0 for aa in all_amino_acids} for position in range(alignment.get_alignment_length())}

    for record in alignment:
        for position, aa in enumerate(record.seq):
            if aa in all_amino_acids:
                aa_counts[position][aa] += 1

    for position in aa_counts:
        total = sum(aa_counts[position].values())
        for aa in aa_counts[position]:
            aa_counts[position][aa] = aa_counts[position][aa] / total if total > 0 else 0
    out_df = pd.DataFrame(aa_counts)
    out_df.to_csv(args.name_stem+'_propensities.out')
    propensity_matrix = np.transpose(out_df.values)
    enrichments = pd.read_csv(args.dms_enrichment).values
    out = enrichments - propensity_matrix[args.off_by_n:args.off_by_c] 
    np.save(args.name_stem+'_diff_enrichment_propensities_raw', out)
    np.save(args.name_stem+'_propensities_raw', propensity_matrix)
    radarPlots(out, args.name_stem+"_diff_enrichment_propensities_radar.pdf", list(SeqIO.parse(args.fasta_file, "fasta"))[0], (args.off_by_n,args.off_by_c), all_amino_acids)
    radarPlots(propensity_matrix[args.off_by_n:args.off_by_c], args.name_stem+"_propensities_radar.pdf", list(SeqIO.parse(args.fasta_file, "fasta"))[0], (args.off_by_n,args.off_by_c), all_amino_acids)
    create_logoplot(out, args.name_stem+"_diff_enrichment_propensities_logo", (args.off_by_n,args.off_by_c), all_amino_acids)
    create_logoplot(propensity_matrix[args.off_by_n:args.off_by_c], args.name_stem+"_propensities_logo", (args.off_by_n,args.off_by_c), all_amino_acids)
    os.remove("reference.fasta"); os.remove("others.fasta")


if __name__ =="__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)



