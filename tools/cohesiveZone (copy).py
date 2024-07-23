import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad
from myModules import average_bin, Needleman, z

# Input
sys_label = 'cristo-dry'

if sys_label=='cristo-dry':
    infile = '../../j-integral/cristobalite/j_integral_dry/crack_tip_1390000.csv'
    Vatom = 142.176*151.589*5.15507/6912
    units = 'real'
    plot_which=0
    numbins = 60
    ifit = 15 # Starting data point for fiting the traction-separation curve
    x0, xf, xtip  = 20, 120.0, 78
    y1, y2 = 92.5, 58.05
    ylo, yup = 75, 77.9
    ym, dy = 70., 5
    mc = 0.156

if sys_label=='graphene':
    infile = '../../crack/md/graphene/airebo/cohesion-zone/j_integral_dry/crack_tip_3010000.csv'
    Vatom = 145.155*146.659*3.35/8400
    units = 'metal'
    plot_which = 2
    numbins = 60 # Adjust this to get compatible up and lower crack profiles
    ifit = 22 # Starting data point for fiting the traction-separation curve
    x0, xf, xtip  = 10, 300, 71
    y1, y2 = 80.44, 69.065 # Surface
    ylo, yup = 73.7, 76.1 # Bulk
    ym, dy = 74.5,  0.0
    mc = 0.050

if sys_label=='graphene-H':
    infile = '../../crack/md/graphene/airebo/cohesion-zone/j_integral_passivated/crack_tip_2720000.csv'
    Vatom = 145.155*146.659*3.35/8400
    units = 'metal'
    plot_which = 2
    numbins = 50 # Adjust this to get compatible up and lower crack profiles
    ifit = 22 # Starting data point for fiting the traction-separation curve
    x0, xf, xtip  = 25, 200, 71
    y1, y2 = 79.8, 68.7 # Surface
    ylo, yup = 72.9, 75.6# Bulk
    ym, dy = 74.2, 0.0
    mc = 0.05

if sys_label=='hematite-001_100':
    infile = '../../j-integral/hematite/j_integral_hematite-dry-001_100/crack_tip_1700000.csv'
    Vatom = 106.058*122.062*9.77174/10320
    units = 'real'
    plot_which=2
    numbins = 65
    ifit = 12 # Starting data point for fiting the traction-separation curve
    x0, xf, xtip  = 19, 75, 53.7
    y1, y2 = 77.2, 45.3 # Surface
    ylo, yup = 60.1, 63.4 # Bulk
    ym, dy = 61.65, 0.2
    mc = 0.20

if sys_label=='hematite-dry-small':
    infile = '../../j-integral/hematite/j_integral_hematite-dry-small/crack_tip_1390000.csv'
    Vatom = 51.7534*52.274*4.98011/1440
    units = 'real'
    plot_which=2
    numbins = 60
    ifit = 15 # Starting data point for fiting the traction-separation curve
    x0, xf, xtip  = 3, 44.7, 19.3
    y1, y2 = 35.8, 19.2
    ylo, yup = 25.7, 29
    ym, dy = 27.5, 0.5
    mc = 0.2

###############################################################################

if units=='metal':
    strconv = 0.0001/Vatom # bar*Å³ to GPa
elif units=='real':
    strconv = 0.000101325/Vatom # atm*Å³ to GPa
else:
    strconv = 1.0

df = pd.read_csv(infile, index_col=False)
df['syy'] = df['syy']*strconv

df1 = df.loc[( ( (df['x']>x0) & (df['x']<xtip) & (df['y']>ym) ) & (df['y']<y1-mc*df['x']) ) |
             ( (df['x']>xtip) & (df['x']<xf) & (df['y']>ym+dy) &  (df['y']<yup) )].sort_values(by="x")

df2 = df.loc[( ( (df['x']>x0) & (df['x']<xtip) & (df['y']<ym) ) & (df['y']>y2+mc*df['x']) ) |
             ( (df['x']>xtip) & (df['x']<xf) & (df['y']<ym-dy) &  (df['y']>ylo) )].sort_values(by="x")

if plot_which==0:
    df12 = pd.concat([df1, df2])
    plt.scatter(df12['x'], df12['y'], c=df12['type'], cmap='bwr', s=10)
    plt.xlabel(r'$x$ ($\AA$)')
    plt.ylabel(r'$y$ ($\AA$)')
    plt.show()
    sys.exit()

xav1, yav1, sav1 = average_bin(df1, numbins, 'syy')
xav2, yav2, sav2 = average_bin(df2, numbins, 'syy')

# Traction-position
if plot_which==1:
    
    #if len(xav1)!=len(xav2):
    #    sys.exit('Array size of upper and lower crack profiles are not compatible. Try changing the number of bins in the average or change how crack atoms are filtered...')

    # Plot the results
    plt.scatter(df1['x'], df1['syy'], s = 50.0, alpha=0.2, label='Up (Raw Data)', color='blue')
    plt.scatter(df2['x'], df2['syy'], s = 50.0,  alpha=0.2, label='Down (Raw Data)', color='red')
    #plt.plot(xav1,sav1, 'bs-', ms=5.0, label='Up (Average)')
    #plt.plot(xav2,sav2, 'r^-', ms=5.0, label='Down (Average)')
    #plt.plot(0.5*(xav1+xav2),0.5*(sav1+sav2), 'ko-', ms=5.0, label='Total Average')

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
    # (integration from separation at maximum stress d0 to ~inf )
    d0 = param[1]/z # get from dTy/dδ=0
    Ws, _ = quad(Needleman, d0,max(u), args=(param[0], param[1]))
    Ws = Ws*1e-1 # GPa.Å to J/m²

    print("Fitting parameters:")
    print("σ_max = %.1f GPa" % param[0])
    print("δ_m = %.1f Å" % param[1])
    print("")
    print("Work of separation = %.2f J/m²" % Ws)
    
    
    xfit = np.linspace(0, max(u), 1000)
    plt.plot(xfit,Needleman(xfit, param[0], param[1]),'r-',label="Needleman's model")
    plt.plot(u,ty, 'bs-',label='MD atomistic')
    plt.xlabel(r'$\delta$ ($\AA$)')
    plt.ylabel(r'$T_{y}$ (GPa)')
    #plt.xlim(1.3,5)
    plt.legend()
    plt.savefig(f'traction_separation_{sys_label}.png', dpi=900)
    plt.show()
