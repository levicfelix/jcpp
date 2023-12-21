import os
import sys
import pandas as pd

def csv2dump(infile,outfile):
    df = pd.read_csv(infile)
    df = df.drop('neighbors', axis=1)
    df['id'] = df['id'].astype(int).astype(str)
    df['type'] = df['type'].astype(int).astype(str)
    df['coordination'] = df['coordination'].astype(int).astype(str)

    # Header
    Natoms = len(df)
    xlo = df['x'].min(); xhi = df['x'].max()
    ylo = df['y'].min(); yhi = df['y'].max()
    zlo = df['z'].min(); zhi = df['z'].max()
    box_bounds = ['ff','ff','ff']
    
    properties = df.columns.tolist(); properties = ' '.join(properties)
    
    lmptrj = open(outfile, 'w')
    lmptrj.write("ITEM: TIMESTEP\n")
    lmptrj.write("0\n")
    lmptrj.write("ITEM: NUMBER OF ATOMS\n")
    lmptrj.write(f"{Natoms}\n")
    lmptrj.write("ITEM: BOX BOUNDS"+' '.join(box_bounds)+"\n")
    lmptrj.write(f"{xlo} {xhi}\n")
    lmptrj.write(f"{ylo} {yhi}\n")
    lmptrj.write(f"{zlo} {zhi}\n")
    lmptrj.write("ITEM: ATOMS "+properties+"\n")
    
    for k, atom in df.iterrows():
        row = [str(value) for value in atom]
        lmptrj.write(' '.join(row)+'\n')
    
    lmptrj.close()