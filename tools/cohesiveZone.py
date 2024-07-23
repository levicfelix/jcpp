import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad
from myModules import average_bin, Needleman, Z

# Input
sys_label = 'cristo-dry'

if sys_label=='cristo-dry':
    infile = '/home/levi/Research/stress-corrosion/crack/md/cristobalite/dry/cz_1390000.csv'
    Vatom = 142.176*151.589*5.15507/6912
    units = 'real'
    plot_which=2
    numbins = 60
    ym = 78.0
    ifit = 10 # Starting data point for fiting the traction-separation curve

###############################################################################

if units=='metal':
    strconv = 0.0001/Vatom # bar*Å³ to GPa
elif units=='real':
    strconv = 0.000101325/Vatom # atm*Å³ to GPa
else:
    strconv = 1.0

df = pd.read_csv(infile)

if plot_which==0:
    plt.scatter(df['x'], df['y'], c=df['q'], cmap='bwr', s=10)
    plt.xlabel(r'$x$ ($\AA$)')
    plt.ylabel(r'$y$ ($\AA$)')
    plt.show()
    sys.exit()

# Extract data
dfup = df[(df['y']>ym) & (df['x']>20)].copy(); dfup.sort_values(by='x', inplace=True)
dflo = df[(df['y']<ym) & (df['x']>20)].copy(); dflo.sort_values(by='x', inplace=True)

# Average over upper and lower atoms
xup, yup, sup = average_bin(dfup,numbins,'syy')
xlo, ylo, slo = average_bin(dflo,numbins,'syy')

# Retain maximum stress
imax = np.where(sup==max(sup))[0][0]; sup[imax] = max(dfup['syy'])
imax = np.where(slo==max(slo))[0][0]; slo[imax] = max(dflo['syy'])


xavg = 0.5*(xup+xlo)
savg = 0.5*(sup+slo)

if plot_which==1:
    plt.scatter(dfup['x'], dfup['syy']*strconv, s = 100.0, alpha=0.3, label='Upper (Raw Data)')
    plt.scatter(dflo['x'], dflo['syy']*strconv, s = 100.0, alpha=0.3, label='Lower (Raw Data)')
    plt.plot(xup,sup*strconv, 's-', ms=5.0, label='Upper (Average)')
    plt.plot(xlo,slo*strconv, 's-', ms=5.0, label='Lower (Average)')
    plt.plot(xavg, savg*strconv, 'ko-', ms=5.0, label='Total Average')
    plt.xlabel(r'$x$ ($\AA$)')
    plt.ylabel(r'$\sigma_{yy}$ (GPa)')
    plt.legend()
    plt.show()

# Traction separation
elif plot_which==2:
    
    u=yup-ylo
    ty = savg*strconv
    
    sorted_pairs = sorted(zip(u, ty))
    ty = np.array([element[1] for element in sorted_pairs])
    u=np.array(sorted(u))
    
    u=u-u[0]
    ty=ty-ty[-1]
    
    # Perform curve fitting for the traction-separation law
    param, param_cov = curve_fit(Needleman, u[ifit:-1], ty[ifit:-1])
    # Calculate work of separation 
    # (integration from separation at maximum stress d0 to ~inf )
    d0 = param[1]/Z # get from dTy/dδ=0
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
