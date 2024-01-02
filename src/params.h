
struct Params {
    std::string systemlabel = "crack";
    bool pbc = true;
    int numtypes = 2;
    std::vector<int> bulkcoordination = stringToVectorInt("4 4");
    std::vector<std::vector<double>> neighborcutoff = stringToMatrix("3.1 3.1 3.1", numtypes, numtypes);
    std::vector<std::vector<double>> neighborscaling = stringToMatrix("1.0 1.0 1.0", numtypes, numtypes);
    double neighborskin = 5.0;
    double thickness = 5.43065;
    bool csvrestart = false;
    int crackstep = 44791;
    int bulkstep = 0;
    int surfacestep = 0;
    int stepmin = 1;
    int stepmax = 9999999;
    std::vector<double> cracktippos = stringToVectorDouble("155.026 154.984");
    std::vector<double> crackregion = stringToVectorDouble("25.0 62.0 15.0 42.0");
    std::string unitslammps = "metal";
    std::string contourshape = "rectangular";
    std::vector<float> contourparams = stringToVectorFloat("10.0 10.0");
    std::string cracktrajfile = "./examples/silicon/critical_crack.lammpstrj";
    std::string bulktrajfile = "./examples/silicon/initial_crack.lammpstrj";
    std::string surfacetrajfile = "./examples/silicon/initial_crack.lammpstrj";
    int defgraddimension = 3;
    bool checkfdim = false;
    double geomtolerance = 0.5;
    std::string defgradmethod = "minsqerr";
    bool findcracktip = false;
    std::string tipsearch = "Rectangle";
    std::string tipstyle = "center";
    bool domaincsv = false;
    bool crackcsv = false;
    double czdelta = 1.0;
    double atomvol = 1.0;
};

// Function to read parameters from an external file
void readParams(const std::string& filename, Params& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::unordered_map<std::string, std::string> paramStrings;

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        // Split the line into key and value
        std::istringstream iss(line);
        std::string key, value;
        if (iss >> key) {
            std::getline(iss >> std::ws, value);
            paramStrings[key] = value;
            //std::cout << key << " " << value << std::endl;
        }
    }

    // Set parameters based on values from the file or defaults
    if (!paramStrings["SystemLabel"].empty()) {
        params.systemlabel = paramStrings["SystemLabel"];
    }
    if (!paramStrings["PBC"].empty()) {
        params.pbc = (paramStrings["PBC"] == "True");
    }
    if (!paramStrings["NumTypes"].empty()) {
        params.numtypes = std::stoi(paramStrings["NumTypes"]);
    }
    if (!paramStrings["BulkCoordination"].empty()) {
        params.bulkcoordination = stringToVectorInt(paramStrings["BulkCoordination"]);
    }
    if (!paramStrings["NeighborCutoff"].empty()) {
        params.neighborcutoff = vectorToMatrix(stringToVectorDouble(paramStrings["NeighborCutoff"]), params.numtypes);
    }
    if (!paramStrings["NeighborScaling"].empty()) {
        params.neighborscaling = vectorToMatrix(stringToVectorDouble(paramStrings["NeighborScaling"]), params.numtypes);
    }
    if (!paramStrings["NeighborSkin"].empty()) {
        params.neighborskin = std::stof(paramStrings["NeighborSkin"]);
    }
    if (!paramStrings["Thickness"].empty()) {
        params.thickness = std::stof(paramStrings["Thickness"]);
    }
    if (!paramStrings["CSVrestart"].empty()) {
        params.csvrestart = (paramStrings["CSVrestart"] == "True");
    }
    if (!paramStrings["CrackStep"].empty()) {
        params.crackstep = std::stoi(paramStrings["CrackStep"]);
    }
    if (!paramStrings["BulkStep"].empty()) {
        params.bulkstep = std::stoi(paramStrings["BulkStep"]);
    }
    if (!paramStrings["SurfaceStep"].empty()) {
        params.surfacestep = std::stoi(paramStrings["SurfaceStep"]);
    }
    if (!paramStrings["StepMin"].empty()) {
        params.stepmin = std::stoi(paramStrings["StepMin"]);
    }
    if (!paramStrings["StepMax"].empty()) {
        params.stepmax = std::stoi(paramStrings["StepMax"]);
    }
    if (!paramStrings["CrackTipPos"].empty()) {
        params.cracktippos = stringToVectorDouble(paramStrings["CrackTipPos"]);
    }
    if (!paramStrings["CrackRegion"].empty()) {
        params.crackregion = stringToVectorDouble(paramStrings["CrackRegion"]);
    }
    if (!paramStrings["UnitsLAMMPS"].empty()) {
        params.unitslammps = paramStrings["UnitsLAMMPS"];
    }
    if (!paramStrings["ContourShape"].empty()) {
        params.contourshape = paramStrings["ContourShape"];
    }
    if (!paramStrings["ContourParams"].empty()) {
        params.contourparams = stringToVectorFloat(paramStrings["ContourParams"]);
    }
    if (!paramStrings["CrackTrajFile"].empty()) {
        params.cracktrajfile = paramStrings["CrackTrajFile"];
    }
    if (!paramStrings["BulkTrajFile"].empty()) {
        params.bulktrajfile = paramStrings["BulkTrajFile"];
    }
    if (!paramStrings["SurfaceTrajFile"].empty()) {
        params.surfacetrajfile = paramStrings["SurfaceTrajFile"];
    }
    if (!paramStrings["DefGradDimension"].empty()) {
        params.defgraddimension = std::stoi(paramStrings["DefGradDimension"]);
    }
    if (!paramStrings["CheckFdim"].empty()) {
        params.checkfdim = (paramStrings["CheckFdim"] == "True");
    }
    if (!paramStrings["GeomTol"].empty()) {
        params.geomtolerance = std::stof(paramStrings["GeomTol"]);
    }
    if (!paramStrings["DefGradMethod"].empty()) {
        params.defgradmethod = paramStrings["DefGradMethod"];
    }
    if (!paramStrings["FindCrackTip"].empty()) {
        params.findcracktip = (paramStrings["FindCrackTip"] == "True");
    }
    if (!paramStrings["CrackCSV"].empty()) {
        params.crackcsv = (paramStrings["CrackCSV"] == "True");
    }
    if (!paramStrings["TipSearchMethod"].empty()) {
        params.tipsearch = paramStrings["TipSearchMethod"];
    }
    if (!paramStrings["TipSearchStyle"].empty()) {
        params.tipstyle = paramStrings["TipSearchStyle"];
    }
    if (!paramStrings["DomainCSV"].empty()) {
        params.domaincsv = (paramStrings["DomainCSV"] == "True");
    }
    if (!paramStrings["CZoneSpacing"].empty()) {
        params.czdelta = std::stof(paramStrings["CZoneSpacing"]);
    }
    if (!paramStrings["AtomVolume"].empty()) {
        params.atomvol = std::stof(paramStrings["AtomVolume"]);
    }
}