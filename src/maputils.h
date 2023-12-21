
// Filter atoms within a specified region
using AtomMap = std::unordered_map<int, Atom>;
void atomsCircularRegion(const AtomMap& allAtoms, AtomMap& filteredAtoms, 
                            double R1, double R2, double x0, double y0) {
    
    for (const auto& entry : allAtoms) {
        const Atom& atom = entry.second;
        
        if (calculateDistance(atom.x-x0, atom.y-y0) > R1 && calculateDistance(atom.x-x0, atom.y-y0) < R2) {
            //std::cout << atom.coordination << std::endl;
            filteredAtoms[atom.id] = atom;
        }
    }
}

void atomsRectangularRegion(const AtomMap& allAtoms, AtomMap& filteredAtoms,
                            double xlo, double xhi, double ylo, double yhi,
                            double x0, double y0) {
    
    for (const auto& entry : allAtoms) {
        const Atom& atom = entry.second;
        double xatom = atom.x-x0; double yatom = atom.y-y0; 
        if (xatom > xlo && xatom < xhi && yatom > ylo && yatom < yhi) {
            //std::cout << atom.coordination << std::endl;
            filteredAtoms[atom.id] = atom;
        }
    }
}

int getAtomIdWithLargestX(const std::unordered_map<int, Atom>& atoms) {
    int largestXAtomId = 0;
    double largestXValue = -1000.0; // begin with a large negative value
    
    for (const auto& entry : atoms) {
        const Atom& atom = entry.second;
        if (atom.x > largestXValue) {
            
            largestXValue = atom.x;
            largestXAtomId = atom.id;
        }
    }

    return largestXAtomId;
}

int getAtomIdWithSecondLargestX(const std::unordered_map<int, Atom>& atoms, int largestXAtomId) {
    int secondLargestXAtomId = 0;
    double secondLargestXValue = -1000.0; // begin with a large negative value

    for (const auto& entry : atoms) {
        const Atom& atom = entry.second;
        
        if (atom.id == largestXAtomId) {
            // Skip the atom with the largest x
            continue;
        }

        if (atom.x > secondLargestXValue) {
            // Current x is larger than the second largest
            secondLargestXValue = atom.x;
            secondLargestXAtomId = atom.id;
        }
    }

    return secondLargestXAtomId;
}

std::vector<double> getAtomXYById(const std::unordered_map<int, Atom>& atoms, int id) {
    auto it = atoms.find(id);
    //std::cout << atoms.x << std::endl;
    if (it != atoms.end()) {
        const Atom& atom = it->second;
        return std::vector<double>{atom.x, atom.y};
    } else {
        std::cerr << "ERROR: Atom with id " << id << " not found in function getAtomXYById. Try changing the crack search style or method. " << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void filterMapByIds(const AtomMap& completeMap, const AtomMap& inputMap, AtomMap& filteredMap) {
    // add only the ones with ids contained in inputMap
    for (const auto& entry : inputMap) {
        const Atom& atom = entry.second;
        int id = atom.id;
        auto it = completeMap.find(id);
        if (it != completeMap.end()) {
            filteredMap[atom.id] = it->second;
        }
    }
}

// Function to find an Atom in AtomMap by id
Atom findAtomById(const AtomMap& atomMap, int id) {
    auto it = atomMap.find(id);
    if (it != atomMap.end()) {
        return it->second;  // Return the Atom if found
    } else {
        // Handle the case when the id is not found
        std::cerr << "ERROR: Atom with id " << id << " not found in function findAtomById" << std::endl;
        std::exit(EXIT_FAILURE);
        // You might want to return a default Atom or throw an exception, depending on your needs
        return Atom();  // Return a default-constructed Atom
    }
}

void shiftPositions(AtomMap& atomMap, double shiftX, double shiftY, double shiftZ) {
    for (auto& entry : atomMap) {
        Atom& atom = entry.second;
        atom.x -= shiftX;
        atom.y -= shiftY;
        atom.z -= shiftZ;
    }
}

void shiftEnergy(AtomMap& targetMap, const AtomMap& bMap, const AtomMap& sMap, double x0) {
    
    for (auto& entry : targetMap) {
        int targetId = entry.first;
        Atom& targetAtom = entry.second;
        
        // Choose reference configuration (bulk or surface) based on x position relative to the crack tip
        //bool foundAtomBoth = false; // Whether atom is found in both maps
        //std::cout << targetId << "   " << targetAtom.type << "   " <<targetAtom.x << "   " << x0 << std::endl;
        if (targetAtom.x>x0 && targetAtom.bulk_coordination>0){
            auto it = bMap.find(targetId);
            //foundAtomBoth = it != bMap.end();
            if (it != bMap.end()) {
                targetAtom.ke -= it->second.ke;
                targetAtom.pe -= it->second.pe;
            } else {
                // Handle the case where the atom is not found in the reference map
                std::cerr << "ERROR: Atom with ID " << targetId << " not found in the bulk map on the shiftEnergy function!" << std::endl;
                std::exit(EXIT_FAILURE);
            }

        } else if (targetAtom.x<=x0 && targetAtom.bulk_coordination>0) {
            auto it = sMap.find(targetId);
            //foundAtomBoth = it != sMap.end();
            if (it != sMap.end()) {
                //std::cout << targetId << "   " << it->second.id << std::endl;
                targetAtom.ke -= it->second.ke;
                targetAtom.pe -= it->second.pe;
            } else {
                // Handle the case where the atom is not found in the reference map
                std::cerr << "ERROR: Atom with ID " << targetId << " not found in the surface map on the shiftEnergy function!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
}

void shiftStress(AtomMap& targetMap, const AtomMap& bMap, const AtomMap& sMap, double x0) {
    
    for (auto& entry : targetMap) {
        int targetId = entry.first;
        Atom& targetAtom = entry.second;
        
        // Choose reference configuration (bulk or surface) based on x position relative to the crack tip
        //bool foundAtomBoth = false; // Whether atom is found in both maps
        //std::cout << targetId << "   " << targetAtom.type << "   " <<targetAtom.x << "   " << x0 << std::endl;
        //if (targetAtom.x>x0){

        // Shift stress tensor relative to bulk only            
        auto it = bMap.find(targetId);
        if (it != bMap.end()) {
            targetAtom.sxx -= it->second.sxx;
            targetAtom.syy -= it->second.syy;
            targetAtom.szz -= it->second.szz;
            targetAtom.sxy -= it->second.sxy;
            targetAtom.sxz -= it->second.sxz;
            targetAtom.syz -= it->second.syz;
        } //else { // if crack is passivated, atoms at the surface are not present in the bulk (perfect) reference structure
            //auto it = sMap.find(targetId);
            //targetAtom.sxx -= it->second.sxx;
            //targetAtom.syy -= it->second.syy;
            //targetAtom.szz -= it->second.szz;
            //targetAtom.sxy -= it->second.sxy;
            //targetAtom.sxz -= it->second.sxz;
            //targetAtom.syz -= it->second.syz;
            // FOR LATER: Make user choose to print this warning (verbose)
            //std::cerr << "WARNING: Atom with ID " << targetId << " not found in the bulk map on the shiftStress function! If this atom is a passivating species, it should be ok..." << std::endl;
            //std::cerr << "ERROR: Atom with ID " << targetId << " not found in the bulk map on the shiftStress function!" << std::endl;
            //std::exit(EXIT_FAILURE);
        //}

        //} else {
        //    auto it = sMap.find(targetId);
        //    if (it != sMap.end()) {
                //std::cout << targetId << "   " << it->second.id << std::endl;
        //        targetAtom.sxx -= it->second.sxx;
        //        targetAtom.syy -= it->second.syy;
        //        targetAtom.szz -= it->second.szz;
        //        targetAtom.sxy -= it->second.sxy;
        //        targetAtom.sxz -= it->second.sxz;
        //        targetAtom.syz -= it->second.syz;
        //    } else {
        //        // Handle the case where the atom is not found in the reference map
        //        std::cerr << "ERROR: Atom with ID " << targetId << " not found in the surface map on the shiftStress function!" << std::endl;
        //        std::exit(EXIT_FAILURE);
        //    }
        //}
    }
}
