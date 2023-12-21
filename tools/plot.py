import numpy as np
import matplotlib.pyplot as plt

def plot_J(sys_label,outdir,q1,q2,shape,prefix,K0,dt,dKIt):
    
    outdir = outdir+f'/j_integral_{sys_label}/'
    
    if shape=='circular':
        fname = f'j_integral_R1_{q1}_R2_{q2}.txt'

    elif shape=='rectangular':
        fname = f'j_integral_Lx_{q1}_Ly_{q2}.txt'
    
    steps, Jen, Jst, J = np.genfromtxt(outdir+fname, usecols=(0,1,2,3), unpack=True)
    Karr = K0+steps*dt*dKIt
    Karr = Karr*1e-2 # GPa√Å to MPa√m

    J = Jen-Jst
    plt.plot(Karr, J, label=prefix)
    plt.legend()
    plt.xlabel(r'$K_I$ [MPa$\sqrt{m}$]')
    plt.ylabel(r'$J$-integral [J/m$^2$]')

    return


loop = 40.0
pltlabel = 'graphene-airebo-ei'

plot_J('dry','../examples/graphene-airebo/excluded-interactions', loop, loop, 
       'rectangular','dry',0.0,0.001,0.05)

plot_J('passivated','../examples/graphene-airebo/excluded-interactions', loop, loop,
       'rectangular','passivated',0.0,0.001,0.05)

plt.xlim(0.0,3.0)
plt.ylim(0.0,9.0)

plt.savefig(f'JKI_{pltlabel}.png',dpi=900)
plt.show()