
// Find the crack tip
namespace findCrackTip {

    std::vector<double> underCoordination(
        const std::unordered_map<int, Atom>& allAtoms, double tz,
        std::vector<std::vector<double>> neighCutoff,
        double neighSkin, int ntypes, int dim,
        std::vector<double> crackRegion,
        const std::vector<int>& Nbulk, const std::string& outdir,
        int timestep, bool crackTipCSV, bool pbc, std::string searchstyle) {
    
        //std::cout << "Finding crack tip position... " << std::endl;
        double xlo = crackRegion[0]; double xhi = crackRegion[1];
        double ylo, yhi, ym, dy;

        if (searchstyle=="center") {
            ym = crackRegion[2]; dy = crackRegion[3];
            ylo = ym-dy; yhi = ym+dy;
        } else if (searchstyle=="range") {
            ylo = crackRegion[2]; yhi = crackRegion[3];
        } else {
            std::cerr << "ERROR: TipSearchStyle not recognized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        //std::cout << xlo << std::endl;

        std::unordered_map<int, Atom> atomsCrack;
        std::unordered_map<int, Atom> atomsAux;
        atomsRectangularRegion(allAtoms, atomsCrack, xlo, xhi, ylo, yhi, 0.0, 0.0);
        atomsRectangularRegion(allAtoms, atomsAux, xlo-neighSkin, xhi+neighSkin, ylo-neighSkin, yhi+neighSkin, 0.0, 0.0);

        // All neighScaling numbers should be equal to 1 since it is of interest
        // only deviation of actual coordination from bulk to identify a surface!
        coordinationAnalysis(atomsCrack, atomsAux, Nbulk, neighCutoff, onesMatrix, tz, dim, pbc);

        // Wheather to save csv for atoms near crack surfaces
        if (crackTipCSV) {
            std::string filename = outdir+"crack_tip_"+std::to_string(timestep)+".csv"; 
            writeCSV(atomsCrack, filename);
        }

        // Get crack tip position
        double crackTipX;
        double crackTipY;
    
        // Separate atomsCrack into two regions based on the y attribute
        if (searchstyle=="center") {
            std::unordered_map<int, Atom> upperRegion, lowerRegion;
        
            for (const auto& entry : atomsCrack) {
            
                const Atom& atom = entry.second;
            
                if (atom.y > ym && atom.coordination < atom.bulk_coordination) {
                    upperRegion[atom.id] = atom;
                } else if (atom.y < ym && atom.coordination < atom.bulk_coordination) {
                    lowerRegion[atom.id] = atom;
                }
            }
            int upperId = getAtomIdWithLargestX(upperRegion);
            int lowerId = getAtomIdWithLargestX(lowerRegion);
            //std::cout << upperId << " " << lowerId << std::endl;
            std::vector<double> upperXY = getAtomXYById(upperRegion, upperId);
            std::vector<double> lowerXY = getAtomXYById(lowerRegion, lowerId);
            crackTipX = std::max(upperXY[0],lowerXY[0]);
            crackTipY = 0.5*(upperXY[1]+lowerXY[1]);
        }

        // No need to separate between upper and lower crack faces
        else if (searchstyle=="range") {
            std::unordered_map<int, Atom> crackRegion;
            for (const auto& entry : atomsCrack) {
                const Atom& atom = entry.second;
        
                if (atom.coordination < atom.bulk_coordination) {
                    crackRegion[atom.id] = atom;
                }
            }
            int largestxId = getAtomIdWithLargestX(crackRegion);
            int secondlargestxId = getAtomIdWithSecondLargestX(crackRegion, largestxId);
        
            std::vector<double> firstXY = getAtomXYById(crackRegion, largestxId);
            std::vector<double> secondXY = getAtomXYById(crackRegion, secondlargestxId);
            crackTipX = std::max(firstXY[0],secondXY[0]);
            crackTipY = 0.5*(firstXY[1]+secondXY[1]);
        }
    
        //std::cout << crackTipX << std::endl;
        return {crackTipX, crackTipY};
    }

    std::vector<double> vonMisesMax(
        const std::unordered_map<int, Atom>& allAtoms,
        std::vector<double> crackRegion, const std::string& outdir,
        std::vector<std::vector<double>> neighCutoff, double neighSkin,
        const std::vector<int>& Nbulk, double tz, int dim, bool pbc,
        int timestep, bool crackTipCSV, std::string searchstyle) {
    
        //std::cout << "Finding crack tip position... " << std::endl;
        double xlo = crackRegion[0]; double xhi = crackRegion[1];
        double ylo = crackRegion[2]; double yhi = crackRegion[3];
    
        std::unordered_map<int, Atom> atomsCrack;
        atomsRectangularRegion(allAtoms, atomsCrack, xlo, xhi, ylo, yhi, 0.0, 0.0);
        
        // Wheather to save csv for atoms near crack surfaces
        if (crackTipCSV) {
            std::unordered_map<int, Atom> atomsAux;
            atomsRectangularRegion(allAtoms, atomsAux, xlo-neighSkin, xhi+neighSkin, ylo-neighSkin, yhi+neighSkin, 0.0, 0.0);
            coordinationAnalysis(atomsCrack, atomsAux, Nbulk, neighCutoff, onesMatrix, tz, dim, pbc);
            std::string filename = outdir+"crack_tip_"+std::to_string(timestep)+".csv"; 
            writeCSV(atomsCrack, filename);
        }

        // Get crack tip position
        double crackTipX;
        double crackTipY;
    
        int largestMisesId = getAtomIdWithLargestMises(atomsCrack);
        std::vector<double> firstXY = getAtomXYById(atomsCrack, largestMisesId);

        // Separate atomsCrack into two regions based on the y attribute
        if (searchstyle=="one") {    
            crackTipX = firstXY[0]; crackTipY = firstXY[1];
        }

        // No need to separate between upper and lower crack faces
        else if (searchstyle=="two") {
            int secondlargestMisesId = getAtomIdWithSecondLargestMises(atomsCrack, largestMisesId);
            std::vector<double> secondXY = getAtomXYById(atomsCrack, secondlargestMisesId);
            crackTipX = 0.5*(firstXY[0]+secondXY[0]);
            crackTipY = 0.5*(firstXY[1]+secondXY[1]);
        }  else {
            std::cerr << "ERROR: TipSearchStyle not recognized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    
        //std::cout << crackTipX << std::endl;
        return {crackTipX, crackTipY};
    }

}

// Compute q-function found in the integrand of J
void qFunction(const std::string& contourShape,
               const AtomMap& allAtoms, 
               AtomMap& domainAtoms, AtomMap& domainAux,
               std::vector<double> ctpos,
               double q1, double q2, double dq){
    
    if (contourShape=="circular"){
        atomsCircularRegion(allAtoms, domainAtoms, q1, q2, ctpos[0], ctpos[1]);
        //filterMapByIds(refAtoms, domainAtoms, domainRef);
        //filterMapByIds(allAtoms, domainAtoms, domainAux);
        //atomsCircularRegion(refAtoms, domainRef, q1, q2, 0.0, 0.0);
        atomsCircularRegion(allAtoms, domainAux, q1-dq, q2+dq, ctpos[0], ctpos[1]);

        // Loop over atoms in the contour domain
        for (auto& entry : domainAtoms) {
            int atomId = entry.first;
            Atom& atom = entry.second;
            double xatom = atom.x-ctpos[0]; double yatom = atom.y-ctpos[1]; 

            // Update q-function atomic info
            double r = calculateDistance(xatom,yatom);
            double q12 = q2-q1;
            atom.q = (r-q1)/q12;
            atom.dqx = xatom/r/q12;
            atom.dqy = yatom/r/q12;
            atom.dqz = atom.z/r/q12;
        }
    }
    else if (contourShape=="rectangular"){
        atomsRectangularRegion(allAtoms, domainAtoms, -q1, q1, -q2, q2, ctpos[0], ctpos[1]);
        //filterMapByIds(refAtoms, domainAtoms, domainRef);
        //filterMapByIds(allAtoms, domainAtoms, domainAux);
        //atomsRectangularRegion(refAtoms, domainRef, -q1, q1, -q2, q2, ctpos[0], ctpos[1]);
        atomsRectangularRegion(allAtoms, domainAux, -q1-dq, q1+dq, -q2-dq, q2+dq, ctpos[0], ctpos[1]);

        // Loop over atoms in the contour domain
        for (auto& entry : domainAtoms) {
            int atomId = entry.first;
            Atom& atom = entry.second;
            double xatom = atom.x-ctpos[0]; double yatom = atom.y-ctpos[1]; 

            // Update q-function atomic info
            double q1q2 = q1*q1*q2*q2;
            atom.q = 1.0 - (q1*q1-xatom*xatom)*(q2*q2-yatom*yatom)/q1q2;
            atom.dqx = 2.0*xatom*(q2*q2-yatom*yatom)/q1q2;
            atom.dqy = 2.0*yatom*(q1*q1-xatom*xatom)/q1q2;
            atom.dqz = 0.0;
        }
    }
    // Update atomic coordination
    //atomicCoordination(domainAtoms, domainAux, Nbulk, neighCutoff, neighScaling, lz, dim, pbc);
}

// Compute J-integral
std::vector<double> computeJintegral(std::unordered_map<int, Atom>& atomsDomain, 
                double tz, int dim, std::string units) {

    // Get LAMMPS units        
    std::vector<double> lmpunits = lammps::Units(units);
    double distanceUnits = lmpunits[0]; double energyUnits = lmpunits[1]; double stressUnits = lmpunits[2];

    tz = tz*distanceUnits;

    double J = 0.0; double JW = 0.0; double JT = 0.0;
    for (auto& entry : atomsDomain) {
        int atomId = entry.first;
        Atom& atom = entry.second;
        if (atom.coordination>dim) {
            
            // Strain energy contribution
            // q-function derivative
            std::vector<double> dq = {atom.dqx/distanceUnits, atom.dqy/distanceUnits, atom.dqz/distanceUnits};
            double w = (atom.ke+atom.pe)*energyUnits;
            double Jstrain = w*dq[0];

            // Traction contribution
            // Displacement gradient tensor Du=F-I, where I is the identity matrix
            std::vector<std::vector<double>> Du(dim, std::vector<double>(dim, 0.0));
            Du[0][0] = atom.Fxx-1.0; Du[0][1] = atom.Fxy;
            Du[1][0] = atom.Fyx; Du[1][1] = atom.Fyy-1.0;
            if (dim == 3) {
                Du[0][2] = atom.Fxz; Du[1][2] = atom.Fyz;
                Du[2][0] = atom.Fzx; Du[2][1] = atom.Fzy;
                Du[2][2] = atom.Fzz-1.0;
            }
            
            
            // Atomic stress tensor
            std::vector<std::vector<double>> Str(dim, std::vector<double>(dim, 0.0));
            Str[0][0] = atom.sxx*stressUnits; Str[0][1] = atom.sxy*stressUnits;
            Str[1][0] = atom.sxy*stressUnits; Str[1][1] = atom.syy*stressUnits;
            if (dim == 3){
                Str[0][2] = atom.sxz*stressUnits; Str[1][2] = atom.syz*stressUnits;
                Str[2][0] = atom.sxz*stressUnits; Str[2][1] = atom.syz*stressUnits;
                Str[2][2] = atom.szz*stressUnits;
            }
            
            double Jtraction = 0.0;
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    Jtraction += Str[i][j]*Du[i][0]*dq[j];
                }
            }
            JW += Jstrain/tz; JT += Jtraction/tz;
            J += (Jstrain - Jtraction)/tz;
        }
    }
    return {JW, JT, J};
}