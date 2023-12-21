
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
                               double lz, int dim, bool pbc) {
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
            if (neighborDistance(atom, auxAtom, dim, pbc, lz) < neighScaling[type1-1][type2-1]*neighCutoff[type1-1][type2-1]
            && atom.id != auxAtom.id) {
                //std::cout << neighborDistance(atom, auxAtom, dim, pbc, lz) << " " << neighScaling[type1-1][type2-1]*neighCutoff[type1-1][type2-1] << std::endl;
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
void deformationGradient(const AtomMap& allAtoms, const AtomMap& bulkAtoms,
                         const AtomMap& surfAtoms, AtomMap& domainAtoms, 
                         double x0, double lz, int dim, bool pbc, std::string fmethod){
    AtomMap refAtoms;
    // Loop over atoms in the contour domain
    for (auto& entry : domainAtoms) {
        int atomId = entry.first;
        Atom& atom1 = entry.second;

        if (atom1.coordination>=dim){
    
            // Choose reference configuration (bulk or surface) based on x position relative to the crack tip
            if (atom1.x>x0){
                refAtoms = bulkAtoms;
            } else {
                refAtoms = surfAtoms;
            }
            
            Atom atom10 = findAtomById(refAtoms,atomId);

            // Define matrices used in the calculation of F (see Int. J. of Solids and Struct. 46 (2009) 238–253)
            std::vector<std::vector<double>> omega(dim, std::vector<double>(dim, 0.0));
            std::vector<std::vector<double>> eta(dim, std::vector<double>(dim, 0.0));
            //std::cout << atomId << std::endl;
            // Iterate over neighbors of this atom
            std::vector<int> atomNeighs = stringToVectorInt(atom1.neighbor_list);
            std::vector<std::vector<double>> F = zeroesMatrix;
            for (const auto& neighborId : atomNeighs) {
                Atom atom2 = findAtomById(allAtoms,neighborId);
                Atom atom20 = findAtomById(refAtoms,neighborId);
                // Spatial reference
                std::vector<double>  xs = calculateRelativeVector({atom1.x, atom1.y, atom1.z},
                                                            {atom2.x, atom2.y, atom2.z},lz, pbc, dim);
                // Material reference
                std::vector<double>  xm = calculateRelativeVector({atom10.x, atom10.y, atom10.z},
                                                            {atom20.x, atom20.y, atom20.z},lz, pbc, dim);

                // Assign values to matrices omega and eta
                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j) {
                        if (fmethod=="minsqerr") {
                            omega[i][j] += xs[i]*xm[j];
                            eta[i][j]   += xm[i]*xm[j];
                        }
                        else if (fmethod=="fdiff") {
                            if (abs(xm[j])> 1.0e-5) { // To get rid of neighors along a given axis (x, y or z) to avoid nan F values
                                F[i][j] += (xs[i]-xm[i])/xm[j]/atom1.coordination+identityMatrix[i][j];
                            }                                
                        }
                    }
                }
            }
            
            if (fmethod=="minsqerr") {
                // All neighbor atoms to a given atom in the integration domain 
                // have been colleted to arrays omega and eta
                std::vector<std::vector<double>> etaInv = inverseMatrix(eta);
                F = matrixMultiply(omega,etaInv);
            }
            // Assign deformation gradient values to atom
            atom1.Fxx = F[0][0]; atom1.Fxy = F[0][1]; 
            atom1.Fyx = F[1][0]; atom1.Fyy = F[1][1];
            if (dim ==3) {
                atom1.Fxz = F[0][2]; atom1.Fyz = F[1][2];
                atom1.Fzx = F[2][0]; atom1.Fzy = F[2][1]; 
                atom1.Fzz = F[2][2];
            }
            
            //std::cout << "TEST!!!" << std::endl;
            
        }
        else if (atom1.coordination>0 && atom1.coordination<dim) {
            std::cout << "WARNING: Some atoms with coordination number lower than system dimensionality are going to be excluded from the calculation of J. This may be avoided by increasing NeighborScaling for the specific atom type." << std:: endl;
        }
        //else if (atom1.coordination == 0 && atom1.bulk_coordination > 0) {
        //    std::cout << "Atomic coordination not updated! Exiting..." << std::endl;
        //    std::exit(EXIT_FAILURE);
        //}
    }
}