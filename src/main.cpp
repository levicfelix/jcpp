#include "screen.h"
#include "dataformat.h"
#include "initparams.h"
#include "external.h"
#include "particleobjects.h"
#include "filemodules.h"
#include "math.h"
#include "maputils.h"
#include "mechanics.h"
#include "fracture.h"

int main(int argc, char* argv[]) {
    
    // Capture the start time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Print logo
    Logger logger;
    
    // Check if the correct number of command line arguments is provided
    if (argc != 2) {
        std::cerr << helpMessage << std::endl;
        return 1;
    }

    Params params; // Initialize with default values

    // Extract the file name from the command line arguments
    std::string paramsfName = argv[1];

    // Open the file
    std::ifstream paramsFile(paramsfName);

    // Check if the file is successfully opened
    if (!paramsFile.is_open()) {
        std::cerr << "ERROR: could not open file: " << paramsfName << std::endl;
        std::exit(EXIT_FAILURE);
        return 1;
    }

    // Read parameters from external file
    readParams(paramsfName, params);

    // Print for what system the calculation is being performed
    std::cout << "  Numerical discrete J-Integral computation for: " << params.systemlabel << std::endl;

    // Create output directory
    std::string outputdir = "j_integral_"+params.systemlabel+"/"; 
    createDirectory(outputdir);
    
    // Column ordering for atomic attribute labels
    std::vector<std::string> hLabels = {"ke", "pe", "sxx", "syy", "szz", "sxy", "sxz", "syz"};

    // Extract atomic info from LAMMPS dump files to csv for the reference configurations
    std::cout << "  Saving bulk reference state... " << std::endl;
    lammps::saveCSV(outputdir, params.bulktrajfile, params.bulkstep, 0, params.bulkstep+1, "bulk", hLabels);
    std::cout << "  Saving surface reference state... " << std::endl;
    lammps::saveCSV(outputdir, params.surfacetrajfile, params.surfacestep, 0, params.surfacestep+1, "surface", hLabels);

    // Delete temporary file with jsteps to avoid conflict with the current config. ones
    deleteFile(outputdir+"jsteps_tmp.txt");
    
    // Extract atomic info from LAMMPS dump files to csv for the current configuration(s)
    if (params.csvrestart==false) {
        std::cout << "  Saving crack configurations... " << std::endl;
        lammps::saveCSV(outputdir, params.cracktrajfile, params.crackstep, params.stepmin, params.stepmax, "crack", hLabels);
        std::filesystem::rename(outputdir+"jsteps_tmp.txt", outputdir+"jsteps.txt");
    }
    else if (params.csvrestart==true) {
        std::string directoryPath = outputdir;
        std::string filePattern = "*.csv";
        bool csvexist = false;
        for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
            if (entry.is_regular_file() && entry.path().extension() == ".csv") {
                csvexist = true;
            }
        }
        if (csvexist==true) {
            std::cout << " CSV files found!" << std::endl;
        }
        else {
            std::cerr << "ERROR: CSV files not found! Try setting CSVrestart to False." << std::endl;
            exit(EXIT_FAILURE);
        }
        bool jstepsexist = std::filesystem::exists(outputdir+"jsteps.txt");
        if (jstepsexist==false) {
            std::cerr << "ERROR: File jsteps.txt not found! Try setting CSVrestart to False." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    //-------------------- Start J-integral calculation --------------------//

    // Load the reference configurations data structure
    std::unordered_map<int, Atom> bulkMap; 
    std::unordered_map<int, Atom> surfMap; 
    readCSV(outputdir, params.bulkstep, "bulk", bulkMap);
    readCSV(outputdir, params.surfacestep, "surface", surfMap);

    // Load temporary file with the steps requested for computation
    std::ifstream jstepsFile(outputdir+"jsteps.txt");
    if (!jstepsFile.is_open()) {
        std::cerr << "ERROR: Could not open jsteps.txt file." << std::endl;
        return 1;
    }
    
    // Update domain label
    std::string domainLabel;
    if (params.contourshape=="circular"){
        domainLabel = "R1_"+formatFloatToString(params.contourparams[0],1)+"_R2_"+
        formatFloatToString(params.contourparams[1],1);
    }
    else if (params.contourshape=="rectangular"){
        domainLabel = "Lx_"+formatFloatToString(params.contourparams[0],1)+"_Ly_"+
        formatFloatToString(params.contourparams[1],1);
    }
    std::cout << "  Integration domain of "+params.contourshape+" shape: "+domainLabel << std::endl;

    // Save results to a file
    std::string jfname = outputdir + "j_integral_" + domainLabel + ".txt";
    std::ofstream jfile(jfname);
    // Check if the file is opened successfully
    if (!jfile.is_open()) {
    std::cerr << "ERROR: Could not open "+jfname+" file." << std::endl;
        return 1;  // Return an error code
    }
    // Write header of file
    jfile << "# Step     Jenergy     Jstress     Jtotal     CrackTipX     CrackTipY\n";
    // Print header of screen
    std::cout << "  Step   Tipx [Å]   Tipy [Å]   Jenergy [J/m²]   Jstress [J/m²]   Jtotal [J/m²]" << std::endl;
    //std::cout << params.neighborscaling[0][0] << std::endl;
    
    // -------------------- Loop over each timestep --------------------//
    std::string line;
    while (std::getline(jstepsFile, line)) {
        int jstep;
        if (std::istringstream(line) >> jstep) {

            // Condition on crack timestep
            // If there is already CSV for all steps and you want to run again but for a single step
            if (params.crackstep != -1) {
                if (jstep < params.crackstep){
                    continue;
                } else if (jstep > params.crackstep) {
                    break;
                } else if (params.crackstep<params.stepmin) {
                    std::cerr << "ERROR: CrackStep cannot be smaller than StepMin!" << std::endl;
                    std::exit(EXIT_FAILURE);
                }  else if (params.crackstep>params.stepmax) {
                    std::cerr << "ERROR: CrackStep cannot be greater than StepMax!" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else if (params.crackstep==-1){
                if (jstep < params.stepmin){
                    continue;
                } else if (jstep > params.stepmax){
                    break;
                }
            }
            std::unordered_map<int, Atom> curMap; // Particle object for the whole system
                                                  // at a given step

            // Update particle object with data from saved CSV files for the entire system in the:
            readCSV(outputdir, jstep, "crack", curMap); // Current configuration
            
            // Search crack tip position
            std::vector<double> crackTipPosition;
            if (params.findcracktip==true) {
                if (params.tipsearch=="Coordination"){
                    crackTipPosition = findCrackTip::underCoordination(curMap, params.thickness, 
                    params.neighborcutoff, params.neighborskin, params.numtypes, params.dimensions,
                    params.crackregion, params.bulkcoordination, outputdir, jstep, params.crackcsv, params.pbc,
                    params.tipstyle);
                }
                else if (params.tipsearch=="vonMises"){
                    crackTipPosition = findCrackTip::vonMisesMax(curMap, params.crackregion, outputdir,
                    params.neighborcutoff, params.neighborskin, params.bulkcoordination, params.thickness, params.dimensions,
                    params.pbc, jstep, params.crackcsv, params.tipstyle);
                }
                else {
                    std::cerr << "ERROR: TipSearchMethod not recognized!" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                
            } else {
                crackTipPosition = params.cracktippos;
            }
            //std::cout << "TEST!!!!" << std::endl;
            
            // Create particle maps for the integration domain with their q-funciton info in the:
            std::unordered_map<int, Atom> curDomain;                      // Current configuration
            std::unordered_map<int, Atom> refDomain;                      // Reference configuration
            std::unordered_map<int, Atom> auxDomain;                      // Neighbor skin
                        
            qFunction(params.contourshape, curMap, curDomain, auxDomain, crackTipPosition,
                      params.contourparams[0], params.contourparams[1], params.neighborskin);

            // Update atomic coordination
            coordinationAnalysis(curDomain, auxDomain, 
                               params.bulkcoordination, params.neighborcutoff, params.neighborscaling,
                               params.thickness, params.dimensions, params.pbc);

            // Compute deformation gradient tensor 
            // Based on Int. J. of Solids and Struct. 46 (2009) 238–253
            //std::cout << "Computing deformation gradient tensor, Fij ..." << std::endl;
            deformationGradient(curMap, bulkMap, surfMap, curDomain, crackTipPosition[0],
                                params.thickness, params.dimensions, params.pbc, params.defgradmethod);
            
            // Save CSV of domain info?                                
            if (params.domaincsv==true) {
                std::string filename = outputdir+"domain_"+domainLabel+"_timestep_"+std::to_string(jstep)+".csv"; 
                writeCSV(curDomain, filename);
            };
            
            // Shift energy and stress relative to their reference configurations
            shiftEnergy(curDomain,bulkMap,surfMap,crackTipPosition[0]);
            //shiftStress(curDomain,bulkMap,surfMap,crackTipPosition[0]);

            // Compute J-integral
            std::vector<double> J = computeJintegral(curDomain,params.thickness,params.dimensions, params.unitslammps);

            // Output to screen
            printScreenColumns(jstep, crackTipPosition, J);

            // Write to the file
            jfile << "  " << std::fixed << std::setprecision(4) << jstep << "        " << J[0] 
            << "           " << J[1] << "             " << J[2] << "          " 
            << crackTipPosition[0] << "             " << crackTipPosition[1] << "\n";
            
            std::cout.flush();
            
            } else {
        std::cerr << "ERROR: Could not convert line to integer: " << line << std::endl;
        std::exit(EXIT_FAILURE);
         }
    }
    // Close the file
    jfile.close();
    
    std::cout << std::endl;  // To make sure the last line keep appearing on terminal screen

    // Close the file
    jstepsFile.close();

    // Capture the end time
    auto endTime = std::chrono::high_resolution_clock::now();
    // Calculate the duration and print on screen
    printElapsedTime(startTime, endTime);

    std::cout << " All done!" << std::endl;
    // do not forget to add lines that clean all objects (maps, arrays, etc) to clean memory
    return 0;
}