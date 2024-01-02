
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

        if (targetAtom.x>x0){

            // Shift stress tensor relative to bulk only            
            auto it = bMap.find(targetId);
            if (it != bMap.end()) {
                targetAtom.sxx -= it->second.sxx;
                targetAtom.syy -= it->second.syy;
                targetAtom.szz -= it->second.szz;
                targetAtom.sxy -= it->second.sxy;
                targetAtom.sxz -= it->second.sxz;
                targetAtom.syz -= it->second.syz;
            } 
        } 
        else {
            auto it = sMap.find(targetId);
            if (it != sMap.end()) {
                targetAtom.sxx -= it->second.sxx;
                targetAtom.syy -= it->second.syy;
                targetAtom.szz -= it->second.szz;
                targetAtom.sxy -= it->second.sxy;
                targetAtom.sxz -= it->second.sxz;
                targetAtom.syz -= it->second.syz;
            } else {
                // Handle the case where the atom is not found in the reference map
                std::cerr << "ERROR: Atom with ID " << targetId << " not found in the surface map on the shiftStress function!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
}

// Function to check if three points are collinear
bool atomsAreCollinear(int atomid, std::vector<int> Neighlist, const AtomMap& atomMap, double tol) {
    const Atom& atom = findAtomById(atomMap,atomid);
    size_t numPoints = Neighlist.size()+1;

    // Check if the number of points is at least 3
    if (numPoints < 3) {
        std::cerr << "ERROR: At least three points/atoms are required to test colinearity." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    Atom neighbor1 = findAtomById(atomMap,Neighlist[0]);
    Atom neighbor2 = findAtomById(atomMap,Neighlist[1]);

    // Check if vectors (p2 - p1) and (p3 - p1) are parallel
    std::vector<double> crossProduct = cross({neighbor1.x-atom.x, neighbor1.y-atom.y, neighbor1.z-atom.z},
                                             {neighbor2.x-atom.x, neighbor2.y-atom.y, neighbor2.z-atom.z});
    
    // Check if the cross product is close to the zero vector
    if (norm(crossProduct) < tol) {
        return true; // Points are collinear
    } else {
        return false; // Points are not collinear
    }
}

bool atomsAreCoplanar(int atomid, std::vector<int> Neighlist, const AtomMap& atomMap, double tol) {
    const Atom& atom = findAtomById(atomMap,atomid);
    size_t numPoints = Neighlist.size()+1;

    // Check if the number of points is at least 3
    if (numPoints < 3) {
        std::cerr << "ERROR: At least four points/atoms are required to test coplanarity." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::vector<double> xCoords(numPoints,0.0);
    std::vector<double> yCoords(numPoints,0.0);
    std::vector<double> zCoords(numPoints,0.0);

    // Feed coordinates of central atom to coord. vectors
    xCoords[0]=atom.x; yCoords[0]=atom.y; zCoords[0]=atom.z;

    for (int in = 1; in < numPoints; ++in) {
        Atom neighbor = findAtomById(atomMap,Neighlist[in-1]);
        xCoords[in] = neighbor.x;
        yCoords[in] = neighbor.y;
        zCoords[in] = neighbor.z;
    }
    // All points collected

    // Calculate the normal vector of the plane formed by the first three points
    std::vector<double> p01 = {xCoords[1]-xCoords[0],yCoords[1]-yCoords[0],zCoords[1]-zCoords[0]}; // p1-p0
    std::vector<double> p02 = {xCoords[2]-xCoords[0],yCoords[2]-yCoords[0],zCoords[2]-zCoords[0]}; // p2-p0
    std::vector<double> normalvec = cross(p01,p02);

    // Check if all the remaining points lie on the plane
    for (size_t i = 3; i < numPoints; ++i) {
        double pointCheck = (xCoords[i] - xCoords[0]) * normalvec[0] +
                            (yCoords[i] - yCoords[0]) * normalvec[1] +
                            (zCoords[i] - zCoords[0]) * normalvec[2];

        // Check if the point lies on the plane (within a small tolerance)
        if (std::abs(pointCheck) > tol) {
            return false; // Points are not coplanar
        }
    }
    return true; // All points lie on the same plane
}

// Function to find the direction vector of a line given two points
std::vector<double> directionVector(const Atom& p1, const Atom& p2) {
    return {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
}

// Function to create the transformation matrix
std::vector<std::vector<double>> transformMatrixLine(int atomid, std::vector<int> Neighlist, const AtomMap& atomMap) {
    
    // Get atom and its neighbors
    const Atom& atom = findAtomById(atomMap,atomid);
    Atom neighbor1 = findAtomById(atomMap,Neighlist[0]);
    Atom neighbor2 = findAtomById(atomMap,Neighlist[1]);

    // Find the direction vector of the line
    std::vector<double> lineDirection = directionVector(atom, neighbor1);

    //Normalize the direction vector and use it as the new X-axis
    std::vector<double> xAxis = normalize(lineDirection);

    // Choose two orthogonal vectors to complete the basis
    std::vector<double> yAxis = createOrthogonalVector(xAxis);
    std::vector<double> zAxis = normalize(cross(xAxis, yAxis));

    // Create the transformation matrix
    std::vector<std::vector<double>> transformationMatrix(3, std::vector<double>(3));

    // Set the columns of the matrix to the new basis vectors
    transformationMatrix[0] = {xAxis[0], xAxis[1], xAxis[2]};
    transformationMatrix[1] = {yAxis[0], yAxis[1], yAxis[2]};
    transformationMatrix[2] = {zAxis[0], zAxis[1], zAxis[2]};

    return transformationMatrix;
}

// Function to find the direction vector of a line given two points
std::vector<double> planeNormalVector(const Atom& p1, const Atom& p2, const Atom& p3) {
    
    // Calculate two vectors lying on the plane
    std::vector<double> v1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    std::vector<double> v2 = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};

    // Calculate the cross product to find the normal vector
    std::vector<double> v3 = cross(v1, v2);

    return normalize(v3);
}

// Function to create the transformation matrix for a planar case
std::vector<std::vector<double>> transformMatrixPlane(int atomid, std::vector<int> Neighlist, const AtomMap& atomMap) {

    // Get atom and two of its neighbors
    const Atom& atom = findAtomById(atomMap,atomid);
    Atom neighbor1 = findAtomById(atomMap,Neighlist[0]);
    Atom neighbor2 = findAtomById(atomMap,Neighlist[1]);

    // Find the normal vector of the plane
    std::vector<double> zAxis = planeNormalVector(atom, neighbor1, neighbor2);

    // Choose two orthogonal vectors to complete the basis
    // One possible orthogonal vector (X-axis) is the direction vector from atom to neighbor1
    std::vector<double> xAxis = {neighbor1.x - atom.x, neighbor1.y - atom.y, neighbor1.z - atom.z};
    xAxis = normalize(xAxis);

    // Choose the third basis vector as the cross product of the normal and X-axis vectors (Y-axis)
    std::vector<double> yAxis = normalize(cross(zAxis,xAxis));

    // Create the transformation matrix
    std::vector<std::vector<double>> transformationMatrix(3, std::vector<double>(3));

    // Set the columns of the matrix to the new basis vectors
    transformationMatrix[0] = {xAxis[0], xAxis[1], xAxis[2]};
    transformationMatrix[1] = {yAxis[0], yAxis[1], yAxis[2]};
    transformationMatrix[2] = {zAxis[0], zAxis[1], zAxis[2]};

    return transformationMatrix;
}