from Bio import AlignIO
try:
    from itolapi import Itol, ItolExport
except:
    raise Exception("Please install the iTOL API Library: pip install itolapi")
import argparse
from pathlib import Path
import tempfile
import pandas as pd
from PyPDF2 import PdfReader, PdfWriter
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
from pdf2image import convert_from_path
import cv2
import numpy as np
from PIL import Image
import os



def create_parser():
    parser = argparse.ArgumentParser(
        description=""
    )
    parser.add_argument(
        "msa_file",
        help="Path to fasta file with seqs to read and create the MSA ",
    )
    parser.add_argument(
        "tree_file",
        type=Path,
        help="Path to fasta file with seqs to read and create the MSA ",
    )
    parser.add_argument(
        "position",
        type = int,
        help="Path to which MSA file should be safed to",
    )
    parser.add_argument(
        "output_image_path",
        help="To be used for the created files",
    )
    parser.add_argument(
        "ITOL_API_KEY",
        help="To be used for the created files",
    )

    return parser

def create_color_map_from_msa(msa_file, position):
    """
    Creates a color map based on amino acid at a specific position in the MSA.
    
    Parameters:
    msa_file (str): Path to the MSA file.
    position (int): Position in the alignment (0-indexed).
    
    Returns:
    dict: Mapping of species to colors.
    """
    global amino_acid_colors
    amino_acid_colors = {
        'A': '#ff9999', 'R': '#9999ff', 'N': '#ccccff', 'D': '#ffcc99',
        'C': '#ffff99', 'Q': '#99ff99', 'E': '#66ff66', 'G': '#ffcc66',
        'H': '#ff66ff', 'I': '#66ccff', 'L': '#99ccff', 'K': '#ff99ff',
        'M': '#99ffcc', 'F': '#ffccff', 'P': '#ff9999', 'S': '#cccccc',
        'T': '#ccffcc', 'W': '#cccc99', 'Y': '#cc99cc', 'V': '#ffcc00',
        '-': '#ffffff', 
    }

    alignment = AlignIO.read(msa_file, 'fasta')

    color_map = {}
    for record in alignment:
        amino_acid = str(record.seq[position]).upper()
        color = amino_acid_colors.get(amino_acid, '#ffffff')  # Default to white if amino acid not found
        color_map[record.id] = color

    return color_map

def upload_and_color_tree(tree_file, msa_file, position, output_image_path, ITOL_API_KEY):
    color_map = create_color_map_from_msa(msa_file, position)
    itol_uploader = Itol()
    itol_uploader.add_file(tree_file)
    print(f'API key: {ITOL_API_KEY}')
    itol_uploader.params['APIkey'] = ITOL_API_KEY
    itol_uploader.params['treeName'] = f'{position}_variation'
    itol_uploader.params['projectName'] = 'rbsc'
    status = itol_uploader.upload()
    if not status:
        print('Tree upload failed:', itol_uploader.comm.upload_output)
        return
    else:
        print(status)
    color_dataset = 'TREE_COLORS\nSEPARATOR TAB\nDATA\n'
    for species, color in color_map.items():
        color_dataset += f'{species}\trange\t{color}\t{species}\n'
    tree_id = itol_uploader.comm.tree_id
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as dataset_file:
        dataset_file.write(color_dataset)
        dataset_file_path = dataset_file.name
    itol_uploader.add_file(Path(dataset_file_path))
    status = itol_uploader.upload()
    if not status:
        print('Dataset upload failed:', itol_uploader.comm.upload_output)
        return
    else:
        print(status)
    itol_exporter = ItolExport()
    itol_exporter.set_export_param_value('tree', tree_id)
    itol_exporter.set_export_param_value('format', 'pdf')
    itol_exporter.set_export_param_value('datasets', 'colors')
    print('Exporting to pdf')
    itol_exporter = itol_uploader.get_itol_export()
    export_location = args.output_image_path
    itol_exporter.set_export_param_value('format', 'pdf')
    itol_exporter.set_export_param_value('datasets', 'dataset1')
    itol_exporter.export(export_location)
    _, ax = plt.subplots(figsize=(12, 1))
    ax.axis('off')
    for i, (aa, color) in enumerate(amino_acid_colors.items()):
        ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
        ax.text(i + 0.5, 0.5, aa, ha='center', va='center', fontsize=20)
    ax.set_xlim(0, len(amino_acid_colors))
    ax.set_ylim(0, 1)

    plt.savefig('colormap.png')
    pages = convert_from_path(args.output_image_path, 500)
    pages[0].save('full_page_image.jpg', 'JPEG')
    full_page_image = cv2.imread('full_page_image.jpg')
    image_to_be_added = cv2.imread('colormap.png')
    final_image = full_page_image.copy()
    final_image[300:400,:1200,:] = image_to_be_added[:,:,:]
    cv2.imwrite('final_image.jpg', final_image)
    final_image = Image.open('final_image.jpg')
    final_image = final_image.convert('RGB')
    final_image.save(args.output_image_path)
    os.remove('final_image.jpg')
    os.remove('full_page_image.jpg')
    os.remove('colormap.png')
    print('Exported tree to', export_location)

def main(args):
    upload_and_color_tree(args.tree_file, args.msa_file, args.position-1, args.output_image_path, args.ITOL_API_KEY)

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)

