import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad

# Input
infile = '../../crack/md/graphene/airebo/cohesion-zone/j_integral_dry/crack_tip_3010000.csv'
df = pd.read_csv(infile)

# Conversion units
Vatom = 145.155*146.659*3.35/8400
units = 'metal'

# Plot settings
plot_which=2
numbins = 50
ifit = 22 # Starting data point for fiting the traction-separation curve
###############################################################################

if units=='metal':
    strconv = 0.0001/Vatom # bar*Å³ to GPa
elif units=='real':
    strconv = 0.000101325/Vatom # atm*Å³ to GPa
else:
    strconv = 1.0
df['syy'] = df['syy']*strconv

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

def Needleman(x,a,b):
    z=16.0*np.e/9.0
    return a*z*(x/b)*np.exp(-z*x/b+1)

df1 = df.loc[( ( (df['pe']>-7.1) & (df['y']>74.5) ) | 
                 ( (df['x']< 71) & ( (df['y']>75-0.01*df['x']) & 
                                     (df['y']<77.8-0.01*df['x']) ) ) ) |
                ( (df['x']>71) & (df['y']>74.4) & 
                  (df['y']<76.1) )].sort_values(by="x")

df2 = df.loc[( ( (df['pe']>-7.1) & (df['y']<74.5) ) | 
               ( (df['x']< 71) & ( (df['y']>71.7+0.01*df['x']) & 
                                   (df['y']<73.8+0.01*df['x']) ) ) ) |
                ( (df['x']>71) & (df['y']>73.4) & 
                  (df['y']<74.5) )].sort_values(by="x")

xav1, yav1, sav1 = average_bin(df1, numbins, 'syy')
xav2, yav2, sav2 = average_bin(df2, numbins, 'syy')

if plot_which==0:
    df12 = pd.concat([df1, df2])
    df12.plot.scatter(x = 'x', y = 'y',c=None,cmap='bwr')

# Traction-position
elif plot_which==1:

    # Plot the results
    plt.scatter(df1['x'], df1['syy'], s = 50.0, alpha=0.2, label='Up (Raw Data)', color='blue')
    plt.scatter(df2['x'], df2['syy'], s = 50.0,  alpha=0.2, label='Down (Raw Data)', color='red')
    plt.plot(xav1,sav1, 'bs-', ms=5.0, label='Up (Average)')
    plt.plot(xav2,sav2, 'r^-', ms=5.0, label='Down (Average)')
    plt.plot(0.5*(xav1+xav2),0.5*(sav1+sav2), 'ko-', ms=5.0, label='Total Average')

    plt.xlabel(r'$x$ ($\AA$)')
    plt.ylabel(r'$\sigma_{yy}$ (GPa)')
    plt.legend()
    plt.show()

# Traction separation
elif plot_which==2:
    
    if len(xav1)!=len(xav2):
        sys.exit('Array size of upper and lower crack profiles are not compatible. Try changing the number of bins in the average or change how crack atoms are filtered...')
    
    u=yav1-yav2
    ty = 0.5*(sav1+sav2)
    
    sorted_pairs = sorted(zip(u, ty))
    ty = np.array([element[1] for element in sorted_pairs])
    u=np.array(sorted(u))
    
    u=u-u[0]
    ty=ty-ty[-1]
    
    # Perform curve fitting for the traction-separation law
    param, param_cov = curve_fit(Needleman, u[ifit:-1], ty[ifit:-1])
    # Calculate work of separation 
    # (integration from separation at maximum stress to ~inf )
    Ws, _ = quad(Needleman, 0.5, 8.0, args=(param[0], param[1]))
    Ws = Ws*1e-1 # GPa.Å to J/m²

    print("Fitting parameters:")
    print("σ_max = %.1f GPa" % param[0])
    print("δ_m = %.1f Å" % param[1])
    print("Work of separation = %.1f J/m²" % Ws)
    
    
    x = np.linspace(0, 8, 1000)
    plt.plot(x,Needleman(x, param[0], param[1]),'r-',label="Needleman's model")
    plt.plot(u,ty, 'bs-',label='MD atomistic')
    plt.xlabel(r'$\delta$ ($\AA$)')
    plt.ylabel(r'$T_{y}$ (GPa)')
    #plt.xlim(1.3,5)
    plt.legend()
    plt.show()

