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


#loop = 40.0
#pltlabel = 'graphene-airebo-ei'

#plot_J('dry','/home/levi/Research/stress-corrosion/crack/md/graphene/airebo/cohesion-zone', 10, 20, 
#       'circular','dry',0.0,0.001,0.05)

#plot_J('passivated','../examples/graphene-airebo/excluded-interactions', loop, loop,
#       'rectangular','passivated',0.0,0.001,0.05)

td, Jd = np.genfromtxt("/home/levi/Research/stress-corrosion/crack/md/graphene/airebo/cohesion-zone/j_integral_dry/j_integral_R1_10.0_R2_20.0.txt", 
                       usecols=(0,5), unpack=True)

tp, Jp = np.genfromtxt("/home/levi/Research/stress-corrosion/crack/md/graphene/airebo/cohesion-zone/j_integral_passivated/j_integral_R1_10.0_R2_20.0.txt", 
                       usecols=(0,5), unpack=True)
print(td)

plt.plot(td,Jd)
plt.plot(tp,Jp)


plt.xlim(0e6,4e6)
#plt.ylim(0.0,9.0)

#plt.savefig(f'JKI_{pltlabel}.png',dpi=900)
plt.show()