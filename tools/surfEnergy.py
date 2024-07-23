# Compute the crack surface energy Γ=2γ 

# Input
sys_label = 'graphene'

if sys_label=='cristobalite-dry':
    Eref = -1040690.5 # Bulk
    Esurf = -1037315.6
    Lx = 142.194 # Length along crack path of the surface
    Lz = 5.1531 # Thickness
    units = 'real'
    
if sys_label=='graphene':
    Eref = -62381.135 # Bulk
    Esurf = -62072.323
    Lx = 145.151 # Length along crack path of the surface
    Lz = 3.35 # Thickness
    units = 'metal'
    
if sys_label=='graphene-H':
    Eref = -62381.135-4.5062933*30 # Ebulk+EH2*NH2 (bulk + total energy of molecules)
    Esurf = -62304.921
    Lx = 147.57 # Length along crack path of the surface
    Lz = 3.35 # Thickness
    units = 'metal'
    
if sys_label=='hematite':
    Eref = -9393.9413
    Esurf = -9214.2448
    Lx = 8.70546 # Length along crack path of the surface
    Lz = 5.02636 # Thickness
    units = 'real'

###############################################################################

if units=='metal':
    eneconv = 16.02176565 # eV/Å² to J/m²
elif units=='real':
    eneconv = 0.6947693367887763 # kcal/mol/Å² to J/m²
else:
    eneconv = 1.0

A = Lx*Lz
Efrac = (Esurf-Eref)/A
Efrac = Efrac*eneconv

print("Crack surface energy per area = %.2f J/m²" % Efrac)