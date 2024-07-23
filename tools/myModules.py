import numpy as np

# Constants
Z=16.0*np.e/9.0

def average_bin(df,nbins,propt):
    xavg, yavg, savg = [], [], []
    count = 0
    x, y, s = 0, 0, 0
    x0 = df['x'].iloc[0]
    Lx = df['x'].iloc[-1]-x0
    dx = Lx/nbins
    for k, row in df.iterrows():
        if abs(row['x']-x0)<dx:
            x += row['x']
            y += row['y']
            s += row[propt]
            count += 1
        else:
            xavg.append(x/count)
            yavg.append(y/count)
            savg.append(s/count)
            x0 = row['x']
            x = x0
            y = row['y']
            s = row[propt]
            count = 1

    return np.array(xavg), np.array(yavg), np.array(savg)

def LinearFit(x,a,b):
    return a*x+b

def QuadFit(x,a,b,c):
    return a*x*x+b*x+c

def CubicFit(x,a,b,c,d):
    return a*x*x*x+b*x*x+c*x+d

def Needleman(x,a,b):
    return a*Z*(x/b)*np.exp(-Z*x/b+1)