import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from myModules import LinearFit, QuadFit, CubicFit

# Implement later the the ramberg-osgood law for nonlinear materials!!!

# Input
sys_label = "corundum"

if sys_label=="hematite":
    infile = '../../basic-properties/hematite/young-modulus/001/thermo.txt'
    step, lx, ly, lz, pe, pzz = np.genfromtxt(infile,usecols=(0,1,2,3,4,7),
                                              unpack=True)
    # Conversion units
    units = 'real'
    # Crystal symmetry
    symm = 'cubic' # use cubic just as a test
    # Plot settings
    plot_which=3 # 1=energy-strain; 2=stress-strain; 3=strain-strain
    elasticity_model = "linear"
    ifit = 0
    ffit = 10
    xmin,xmax = 0.0, 6.0
    ymin,ymax = 0.0,16.0

elif sys_label=="corundum":
    infile = '../../basic-properties/corundum/young-modulus/001/thermo.txt'
    step, lx, ly, lz, pe, pzz = np.genfromtxt(infile,usecols=(0,1,2,3,4,7),
                                              unpack=True)
    # Conversion units
    units = 'real'
    # Crystal symmetry
    symm = 'cubic' # use cubic just as a test
    # Plot settings
    plot_which=1# 1=energy-strain; 2=stress-strain; 3=strain-strain
    elasticity_model = "quadratic"
    ifit = 0
    ffit = -25
    xmin,xmax = 0.0, 5.5
    ymin,ymax = 0.0, 50.0
    
elif sys_label=="cristobalite":
    infile = '../../basic-properties/cristobalite/young-modulus/thermo.txt'
    step, lx, ly, lz, pe, pzz = np.genfromtxt(infile,usecols=(0,1,3,2,4,6),
                                              unpack=True)
    # Conversion units
    units = 'real'
    # Crystal symmetry
    symm = 'cubic' # use cubic just as a test
    # Plot settings
    plot_which=3 # 1=energy-strain; 2=stress-strain; 3=strain-strain
    elasticity_model = "linear"
    ifit = 0
    ffit =15
    xmin,xmax = 0.0, 10.
    ymin,ymax = 0.0,20.0
else:
    sys.exit("ERROR: System label not recognized!")
    

###############################################################################

if units=='metal':
    eneconv = 1.0 # Use eV as standard energy unit
    strconv = 0.0001 # bar to GPa
elif units=='real':
    eneconv = 0.0433641153087705 # kcal/mol to eV
    strconv = 0.000101325 # atm to GPa
else:
    eneconv = 1.0
    strconv = 1.0
pe = (pe-pe[0])*eneconv
szz = -pzz*strconv

strain_x = 1e2*(lx-lx[0])/lx[0]
strain_y = 1e2*(ly-ly[0])/ly[0]
strain_z = 1e2*(lz-lz[0])/lz[0]

if elasticity_model is not None:
    elasticity_model = elasticity_model.lower()

### Calculate elastic moduli ###
if elasticity_model=="linear":
    # Young's Modulus (energy method)
    param2, param_cov2 = curve_fit(QuadFit, strain_z[ifit:ffit], pe[ifit:ffit])
    Yene = 2.0*param2[0]*1e2 # U = Y*eps**2/2
    print("Young's modulus (from energy) = %.2f GPa" % Yene)
    # Young's Modulus (stress method)
    param1, param_cov1 = curve_fit(LinearFit, strain_z[ifit:ffit], szz[ifit:ffit])
    Ystr = param1[0]*1e2 # U = Y*eps
    print("Young's modulus (from stress) = %.2f GPa" % (Ystr))
elif elasticity_model=="quadratic":
    param2, param_cov2 = curve_fit(CubicFit, strain_z[ifit:ffit], pe[ifit:ffit])
    Yene = 2.0*param2[1]*1e2 # U = D*eps**3/3 + Y*eps**2/2
    print("Young's modulus (from energy) = %.2f GPa" % Yene)
    # Young's Modulus (stress method)
    param1, param_cov1 = curve_fit(QuadFit, strain_z[ifit:ffit], szz[ifit:ffit])
    Ystr = param1[1]*1e2 # U = D*eps**2 + Y*eps
    print("Young's modulus (from stress) = %.2f GPa" % (Ystr))
elif elasticity_model==None:
    None
else:
    sys.exit("ERROR: Elasticity model not recognized!")

# Poisson's ratio
paramx, param_covx = curve_fit(LinearFit, strain_z[ifit:ffit],
                               strain_x[ifit:ffit])
paramy, param_covy = curve_fit(LinearFit, strain_z[ifit:ffit],
                               strain_y[ifit:ffit])
vzx = -paramx[0]
vzy = -paramy[0]
v = 0.5*(vzx+vzy)
print("Poisson's ratio (zx) = %.2f" % vzx)
print("Poisson's ratio (zy) = %.2f" % vzy)

if symm=='cubic':
    
    # For cubic symmetry
    B = Yene/(3.0*(1.0-2.0*v))
    G = Yene/(2.0*(1.0+v))
    print("")
    print("Cubic (averaged) Elastic moduli:")
    print("Bulk modulus = %.2f GPa" % B)
    print("Shear modulus = %.2f GPa" % G)


###############################################################################

# Energy-strain
if plot_which==1:
    xfit = np.linspace(0, max(strain_z)+1, 1000)
    plt.plot(strain_z, pe, 'bs', ms=3.0,label="MD atomistic")
    if elasticity_model=="linear":
        print("Fitting parameters (y=ax²+bx+c):")
        print("a = %.2f" % param2[0])
        print("b = %.2f" % param2[1])
        print("c = %.2f" % param2[2])
        plt.plot(xfit, QuadFit(xfit,param2[0],param2[1],param2[2]), 'r-', label="Quadratic fit")
    elif elasticity_model=="quadratic":
        print("Fitting parameters (y=ax³+bx²+cx+d):")
        print("a = %.2f" % param2[0])
        print("b = %.2f" % param2[1])
        print("c = %.2f" % param2[2])
        print("d = %.2f" % param2[3])
        plt.plot(xfit, CubicFit(xfit,param2[0],param2[1],param2[2],param2[3]), 'r-', label="Cubic fit")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xlabel("Normal strain (%)")
    plt.ylabel("Potential energy (eV)")
    plt.legend()
    plt.show()

# Stress-strain
elif plot_which==2:
    xfit = np.linspace(0, max(strain_z)+1, 1000)
    plt.plot(strain_z, szz, 'bs', ms=3.0,label="MD atomistic")
    if elasticity_model=="linear":
        print("Fitting parameters (y=ax+b):")
        print("a = %.2f" % param1[0])
        print("b = %.2f" % param1[1])
        plt.plot(xfit, LinearFit(xfit,param1[0],param1[1]), 'r-', label="Linear fit")
    elif elasticity_model=="quadratic":
        print("Fitting parameters (y=ax²+bx+c):")
        print("a = %.2f" % param2[0])
        print("b = %.2f" % param2[1])
        print("c = %.2f" % param2[2])
        plt.plot(xfit, QuadFit(xfit,param2[0],param2[1],param2[2]), 'r-', label="Quadratic fit")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xlabel("Normal strain (%)")
    plt.ylabel("Normal stress (GPa)")
    plt.legend()
    plt.show()

# Poisson's ratio
elif plot_which==3:    
    print("Fitting parameters (y=ax+b):")
    print("a(xx) = %.2f" % paramx[0])
    print("b(xx) = %.2f" % paramx[1])
    print("a(yy) = %.2f" % paramy[0])
    print("b(yy) = %.2f" % paramy[1])
    print("")

    xfit = np.linspace(0, max(strain_z)+1, 1000)
    plt.plot(strain_z, strain_x, 'bs', ms=4.0,label=r'$\epsilon_{xx}$ (MD)')
    plt.plot(strain_z, strain_y, 'g^', ms=4.0,label=r'$\epsilon_{yy}$ (MD)')
    plt.plot(xfit, LinearFit(xfit,paramx[0],paramx[1]), 'r-', 
             label=r'$\epsilon_{xx}$ (Fit)')
    plt.plot(xfit, LinearFit(xfit,paramy[0],paramy[1]), 'k--', 
             label=r'$\epsilon_{yy}$ (Fit)')
    #plt.xlim(xmin,xmax)
    #plt.ylim(ymin,ymax)
    plt.xlabel("Normal strain (%)")
    plt.ylabel("Transversal strain (%)")
    plt.legend()
    plt.show()
else:
    sys.exit("ERROR: Plot option not recognized! (plot_which=1,2,3)")
    