
// Calculate atomic von Mises stress
double vonMIses(const Atom& atom){
    double Mises2 = ( atom.sxx - atom.syy ) * ( atom.sxx - atom.syy ) + 
                    ( atom.syy - atom.szz ) * ( atom.syy - atom.szz ) +
                    ( atom.szz - atom.sxx ) * ( atom.szz - atom.sxx ) +
                    6.0*( atom.sxy*atom.sxy + atom.syz*atom.syz + atom.sxz*atom.sxz );
    double Mises = sqrt( 0.5 * Mises2 );
    return Mises;
}

int getAtomIdWithLargestMises(const std::unordered_map<int, Atom>& atoms) {
    int largestMisesAtomId = 0;
    double largestMisesValue = -1000.0; // begin with a large negative value
    
    for (const auto& entry : atoms) {
        const Atom& atom = entry.second;
        double atomMises = vonMIses(atom);
        if (atomMises > largestMisesValue) {
            
            largestMisesValue = atomMises;
            largestMisesAtomId = atom.id;
        }
    }

    return largestMisesAtomId;
}

int getAtomIdWithSecondLargestMises(const std::unordered_map<int, Atom>& atoms, int largestMisesAtomId) {
    int secondLargestMisesAtomId = 0;
    double secondLargestMisesValue = -1000.0; // begin with a large negative value

    for (const auto& entry : atoms) {
        const Atom& atom = entry.second;
        double atomMises = vonMIses(atom);

        if (atom.id == largestMisesAtomId) {
            // Skip the atom with the largest x
            continue;
        }

        if (atomMises > secondLargestMisesValue) {
            // Current x is larger than the second largest
            secondLargestMisesValue = atomMises;
            secondLargestMisesAtomId = atom.id;
        }
    }

    return secondLargestMisesAtomId;
}

// Function to calculate coordination number for each atom in atomsDomain
void coordinationAnalysis(std::unordered_map<int, Atom>& atomsDomain,
                               const std::unordered_map<int, Atom>& auxDomain,
                               const std::vector<int>& Nbulk,
                               std::vector<std::vector<double>> neighCutoff,
                               std::vector<std::vector<double>> neighScaling,
                               double lz, bool pbc) {
    //std::cout << neighCutoff[0][0] << std::endl;
    // Loop over atoms in atomsDomain
    for (auto& entryAtomsDomain : atomsDomain) {
        int atomId = entryAtomsDomain.first;
        Atom& atom = entryAtomsDomain.second;
        
        // Update the 'bulk_coordination' value for the current atom in atomsCrack
        int atype = atom.type;     
        atom.bulk_coordination = Nbulk[atype-1];
        
        // Count the number of neighbors within the specified cutoff
        int coordinationCount = 0;
        std::vector<int> neighIds;
        int type1 = atom.type;
        // Loop over atoms in auxDomain
        for (const auto& entryAuxDomain : auxDomain) {
            const Atom& auxAtom = entryAuxDomain.second;
            
            // Check if the neighbor distance is within the cutoff
            int type2 = auxAtom.type;
            if (neighborDistance(atom, auxAtom, pbc, lz) < neighScaling[type1-1][type2-1]*neighCutoff[type1-1][type2-1]
            && atom.id != auxAtom.id) {
                coordinationCount++;
                neighIds.push_back(auxAtom.id);
            }
        }

        // Update the coordination count for the current atom in atomsDomain
        atom.coordination = coordinationCount;
        atom.neighbor_list = vectorIntToString(neighIds);
        
    }
    //std::cout << "TEST!!!!" << std::endl;
}

// Compute the deformation gradient tensor
namespace deformationGradient {
    void minSquarError(const AtomMap& allAtoms, const AtomMap& bulkAtoms,
                        const AtomMap& surfAtoms, AtomMap& domainAtoms, 
                        double x0, double lz, bool pbc, bool checkfdim, double geomtol, int dimF){
        AtomMap refAtoms;
        // Loop over atoms in the contour domain
        for (auto& entry : domainAtoms) {
            int atomId = entry.first;
            Atom& atom1 = entry.second;
            std::vector<double> at1 = {atom1.x, atom1.y, atom1.z};
            
            if (atom1.coordination>0) { // Do not compute F for atoms with no neighbors
            // Choose reference configuration (bulk or surface) based on x position relative 
            // to the crack tip
            std::string refloc;
            if (atom1.x>x0){
                refloc="Bulk";
                refAtoms = bulkAtoms;
            } else {
                refloc="Surface";
                refAtoms = surfAtoms;
            }

            // Read data of atom in the reference configuration and attribute to vectors
            Atom atom10 = findAtomById(refAtoms,atomId);
            std::vector<double> at10 = {atom10.x, atom10.y, atom10.z};

            // Get neighbor list of ids
            std::vector<int> atomNeighs = stringToVectorInt(atom1.neighbor_list);

            // Introduce a condition here to skip this search if you are certain that it is not needed
            // Get dimensionality of F based on coplanarity (or colinearity) of atomic neighborhood
            std::vector<std::vector<double>> T(3, std::vector<double>(3, 0.0));
            if ( checkfdim==true ) {
                if (atom1.coordination==1) {
                    dimF=1;
                }
                else if (atom1.coordination==2) { 
                    if (atomsAreCollinear(atom1.id,atomNeighs,refAtoms,geomtol)) { 
                        dimF=1;
                    }
                    else { 
                        dimF=2;
                    }
                }
                else if (atom1.coordination>=3) {
                    if (atomsAreCoplanar(atom1.id,atomNeighs,refAtoms,geomtol)) {
                        dimF=2;
                    }
                    else {
                        dimF=3;
                    }
                }
            }
            
            // Get matrix transformation for atomic coordinates
            if (dimF==1) {
                T = transformMatrixLine(atom1.id,atomNeighs,refAtoms);
            } 
            else if (dimF==2) {
                T = transformMatrixPlane(atom1.id,atomNeighs,refAtoms);
            }
            else if (dimF==3) {
                T = identityMatrix;
            }
            else { 
                std::cout << "  Invalid F dimensionality..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // Test if T is singular
            double detT = determinantMatrix(T);
            if (detT==0) {
                std::cout << " ERROR: Singular coordinate transformation matrix found at atom " << atom1.id << "." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // Transform vectors
            at1 = transformVector(T,at1);
            at10 = transformVector(T,at10);
            
            // Define matrices used for the calculation of F 
            // (see Int. J. of Solids and Struct. 46 (2009) 238–253)
            std::vector<std::vector<double>> omega(dimF, std::vector<double>(dimF, 0.0));
            std::vector<std::vector<double>> eta(dimF, std::vector<double>(dimF, 0.0));
            std::vector<std::vector<double>> F = zeroesMatrix;

            // Iterate over neighbors of this atom
            for (const auto& neighborId : atomNeighs) {
                Atom atom2 = findAtomById(allAtoms,neighborId);
                Atom atom20 = findAtomById(refAtoms,neighborId);
                std::vector<double> at2 = {atom2.x, atom2.y, atom2.z};
                std::vector<double> at20 = {atom20.x, atom20.y, atom20.z};
            
                // Transform vectors in the reference state
                at2 = transformVector(T,at2);
                at20 = transformVector(T,at20);

                // Deformed (spatial) configuration
                std::vector<double>  xs = calculateRelativeVector(at1, at2, lz, pbc);
                // Reference (undeformed/material) configuration
                std::vector<double>  xm = calculateRelativeVector(at10, at20, lz, pbc);

                // Assign values to matrices omega and eta
                for (int i = 0; i < dimF; ++i) {
                    for (int j = 0; j < dimF; ++j) {
                        omega[i][j] += xs[i]*xm[j];
                        eta[i][j]   += xm[i]*xm[j];
                    }
                }
            }
            // End loop on atoms inside integration path

            // Test if eta is singular
            double deta = determinantMatrix(eta);
            if (deta==0) {
                std::cout << " ERROR: Singular matrix eta found at atom " << atom1.id << "." << std::endl;
                std::cout << " Try changing DefGradDimension or setting CheckFdim to True..." << std::endl;
                std::exit(EXIT_FAILURE);
            }

            // Calculate F based on neighborhood dimensionality
            if (dimF==3) { // F in 3D space
                std::vector<std::vector<double>> etaInv = inverseMatrix(eta);
                F = matrixMultiply(omega,etaInv);
            }
            else if (dimF==2) { // F in 2D space
                std::vector<std::vector<double>> etaInv = inverseMatrix(eta);
                std::vector<std::vector<double>> Fr = matrixMultiply(omega,etaInv);
                std::vector<std::vector<double>> Fp = matrixIn3D(Fr); // Project to 3D space           
                F = matrixMultiply(inverseMatrix(T),matrixMultiply(Fp,T)); // F = inv(T)F'T
            }
            else if (dimF==1) {
                std::vector<std::vector<double>> Fp = zeroesMatrix; 
                Fp[0][0] = omega[0][0]/eta[0][0]; // Assing only x'
                F = matrixMultiply(inverseMatrix(T),matrixMultiply(Fp,T)); // F = inv(T)F'T
            }

            // Assign deformation gradient values to atom
            atom1.Fxx = F[0][0]; atom1.Fxy = F[0][1]; atom1.Fxz = F[0][2];
            atom1.Fyx = F[1][0]; atom1.Fyy = F[1][1]; atom1.Fyz = F[1][2];
            atom1.Fzx = F[2][0]; atom1.Fzy = F[2][1]; atom1.Fzz = F[2][2];

            // Print atomic info
            //if (dimF!=2) {
            //std::cout << " ------ " << std::endl;
            //std::cout << "Atom ID: " << atom1.id << " Coord.: " << atom1.coordination << " dimF: " << dimF << " RefLoc: " << refloc << std::endl;
            //printMatrix(F);
            //}
            } // else { std::cout << "  WARNING: Atom with id: " << atom1.id << " has no neighbors!" << std::endl; }
        }
    } // end of minSquareError

    void finiteDifference(const AtomMap& allAtoms, const AtomMap& bulkAtoms,
                            const AtomMap& surfAtoms, AtomMap& domainAtoms, 
                            double x0, double lz, bool pbc){
        AtomMap refAtoms;
        // Loop over atoms in the contour domain
        for (auto& entry : domainAtoms) {
            int atomId = entry.first;
        
            // Read data of atom in the reference configuration and attribute to vectors
            Atom& atom1 = entry.second;
            
            std::string refloc;
            if (atom1.x>x0){
                refloc="Bulk";
                refAtoms = bulkAtoms;
            } else {
                refloc="Surface";
                refAtoms = surfAtoms;
            }
            Atom atom10 = findAtomById(refAtoms,atomId);
            
            std::vector<double> at1 = {atom1.x, atom1.y, atom1.z};
            std::vector<double> at10 = {atom10.x, atom10.y, atom10.z};
            
            // Choose reference configuration (bulk or surface) based on x position relative 
            // to the crack tip
            if (atom1.x>x0){
                refAtoms = bulkAtoms;
            } else {
                refAtoms = surfAtoms;
            }

            // Get neighbor list of ids
            std::vector<int> atomNeighs = stringToVectorInt(atom1.neighbor_list);
        
            // Define matrix that stores the deformation gradient tensor per atom
            std::vector<std::vector<double>> F = zeroesMatrix;
        
            // Iterate over neighbors of this atom
            for (const auto& neighborId : atomNeighs) {
                Atom atom2 = findAtomById(allAtoms,neighborId);
                Atom atom20 = findAtomById(refAtoms,neighborId);
                std::vector<double> at2 = {atom2.x, atom2.y, atom2.z};
                std::vector<double> at20 = {atom20.x, atom20.y, atom20.z};
            
                // Deformed (spatial) configuration
                std::vector<double>  xs = calculateRelativeVector(at1, at2, lz, pbc);
                // Reference (undeformed/material) configuration
                std::vector<double>  xm = calculateRelativeVector(at10, at20, lz, pbc);

                // Assign values to matrices omega and eta
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                            if (abs(xm[j])> 1.0e-5) { // To get rid of neighors along a given axis (x, y or z) to avoid nan F values
                                F[i][j] += (xs[i]-xm[i])/xm[j]/atom1.coordination+identityMatrix[i][j];
                            }                                
                        
                    }
                }
            }
            // End loop on atoms inside integration path
            
            // Assign deformation gradient values to atom
            atom1.Fxx = F[0][0]; atom1.Fxy = F[0][1]; atom1.Fxz = F[0][2];
            atom1.Fyx = F[1][0]; atom1.Fyy = F[1][1]; atom1.Fyz = F[1][2];
            atom1.Fzx = F[2][0]; atom1.Fzy = F[2][1]; atom1.Fzz = F[2][2];
            
            //if (dimF==2) {
            //std::cout << " ------ " << std::endl;
            //std::cout << "Atom ID: " << atom1.id << " Coord.: " << atom1.coordination << std::endl;
            //printMatrix(F);
            //}
        }
    } // end of finiteDifference

}