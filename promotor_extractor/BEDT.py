#!/usr/bin/env python

import os
import pandas as pd
import subprocess
import argparse

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bedfile",
        type=str,
     )

    parser.add_argument(
        "genelist",
        type=str,
    )

    parser.add_argument(
        "genomefasta",
    )
    
    parser.add_argument(
         "workingdir",
         type = str,
    )
    parser.add_argument(
         "shscript",
         type = str,
    )
    parser.add_argument(
        "--promotorregion",
        default=1000,
        type=int,
    )
    return parser


def main(args):
    os.makedirs(args.workingdir, exist_ok=True)
    with open(args.bedfile, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(os.path.join(args.workingdir,'interm.bed'), 'w') as fout:
        fout.writelines(data[1:])


    df = pd.read_csv(os.path.join(args.workingdir,'interm.bed'),sep='\t', names = ['chrom', 'chromStart' ,'chromEnd' ,'name' ,'score' ,'strand' ,'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])



    df = df.loc[:, 'chrom':'strand']



    with open(args.genelist,'r') as fn:
        GeneList = fn.read()
        GeneList = GeneList.split('\n')

    data = df['name'].apply(lambda x: any([k in x for k in GeneList]))
    entr = []
    for idx, i in enumerate(data.values):
        if i == True:
            entr.append(list(data.index)[idx])

    df2 = df.iloc[entr,:]

    df2[df2['strand']=='-']

    sel_gen = list(df2['name'])
    ingenes = [gene for gene in GeneList if sum(str(gene) in s for s in sel_gen) > 1]
    for idx,_ in enumerate(ingenes):
        sl = df2[df2['name'] in ingenes[idx]]
        sl.reset_index(inplace = True)
        df2.drop(labels= df2.index[df2['name']==ingenes[idx]].tolist(), axis=0, inplace=True)
        sub=pd.Series([sl.loc[len(sl)-1]['chrom'],sl.loc[len(sl)-1]['chromStart'], sl.loc[len(sl)-1]['chromEnd'],ingenes[idx],sl.loc[len(sl)-1]['score'],sl.loc[len(sl)-1]['strand']],index = df2.columns)
        df2=df2.append(sub,ignore_index=True)
        del sl


    pd.options.mode.chained_assignment = None
    for idx in df2['chromStart'].index:
        chromStart = df2['chromStart'][idx]
        chromEnde = df2['chromEnd'][idx]
        if chromStart >1000 and df2['strand'][idx] == '+':
            df2['chromStart'][idx] = chromStart-args.promotorregion
            df2['chromEnd'][idx] = chromStart
        if df2['strand'][idx] == '-':
            df2['chromStart'][idx] = chromEnde
            df2['chromEnd'][idx] = chromEnde+args.promotorregion

    df2['chromStart'] = df2['chromStart'].astype('int')
    df2['chromEnd'] = df2['chromEnd'].astype('int')

    df2.to_csv(os.path.join(args.workingdir,'finished.bed'), sep ='\t' , index=False, header=False)

    p = subprocess.Popen(["sh",f"{args.shscript}",f"{args.genomefasta}", f"{os.path.join(args.workingdir,'finished.bed')}",
                        f"{os.path.join(args.workingdir,'promotor.fasta')}"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


    stdout, stderr = p.communicate()
    print("stdout: '%s'" % stdout)
    print("stderr: '%s'" % stderr)


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)