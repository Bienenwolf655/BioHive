from tqdm import tqdm
import os
for i in tqdm(range(450)):
    os.system(f'python ITOL_tree_from_pos.py "edited_rub_msa-modified translation.fasta" rbsc_phydms_ExpCM_aa_preference_shortend_fr_tree.tree {i+1} phyloplots/rbsc_position_{i+1}_variation.png xLUDkNZKkrdXbZQMseD4zg')
