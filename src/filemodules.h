namespace fs = std::filesystem;

bool directoryExists(const std::string& path) {
    return fs::exists(path) && fs::is_directory(path);
}

bool createDirectory(const std::string& path) {
    if (!directoryExists(path)) {
        try {
            fs::create_directory(path);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
            return false;
        }
    } else {
        std::cout << "  Folder '"+path+"' already exists." << std::endl;
        return true;
    }
}

// Function to delete a file
void deleteFile(const std::string& fileName) {
    std::remove(fileName.c_str());
}

// Function to read CSV and update the map with Atom entries
void readCSV(const std::string& outdir, int timestep, const std::string& prefix, std::unordered_map<int, Atom>& atomMap) {
    // Read data from CSV file
    std::string filepath = outdir + prefix + "_" + std::to_string(timestep) + ".csv";
    
std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error opening CSV file." << std::endl;
        return;
    }

    // Read the header line
    std::string line;
    std::getline(file, line);
    std::istringstream headerStream(line);
    std::vector<std::string> headers;

    // Split the header line into column names
    std::string header;
    while (std::getline(headerStream, header, ',')) {
        headers.push_back(header);
    }

    // Read data lines
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        Atom atom;

        // Read each column and update the corresponding attribute in the Atom struct
        int colIndex = 0;
        std::string value;
        while (std::getline(lineStream, value, ',')) {
            if (headers[colIndex] == "id") {
                atom.id = std::stoi(value);
            } else if (headers[colIndex] == "type") {
                atom.type = std::stoi(value);
            } else if (headers[colIndex] == "x") {
                atom.x = std::stod(value);
            } else if (headers[colIndex] == "y") {
                atom.y = std::stod(value);
            } else if (headers[colIndex] == "z") {
                atom.z = std::stod(value);
            } else if (headers[colIndex] == "ke") {
                atom.ke = std::stod(value);
            } else if (headers[colIndex] == "pe") {
                atom.pe = std::stod(value);
            } else if (headers[colIndex] == "sxx") {
                atom.sxx = std::stod(value);
            } else if (headers[colIndex] == "syy") {
                atom.syy = std::stod(value);
            } else if (headers[colIndex] == "szz") {
                atom.szz = std::stod(value);
            } else if (headers[colIndex] == "sxy") {
                atom.sxy = std::stod(value);
            } else if (headers[colIndex] == "sxz") {
                atom.sxz = std::stod(value);
            } else if (headers[colIndex] == "syz") {
                atom.syz = std::stod(value);
            }

            // Move to the next column
            colIndex++;
        }

        // Update the atom in the map
        atomMap[atom.id] = atom;
    }

    //std::cout << "CSV file read successfully." << std::endl;
}

template <typename KeyType, typename ValueType>
void writeCSV(const std::unordered_map<KeyType, ValueType>& myMap, const std::string& filename) {
    std::ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Print header
    outputFile << "id,type,x,y,z,ke,pe,sxx,syy,szz,sxy,sxz,syz,bulk,coordination,neighbors,q,dqx,dqy,dqz,Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz\n";

    // Print data
    for (const auto& entry : myMap) {
        outputFile << entry.first << ',' << entry.second.type << ',' 
        << entry.second.x << ',' << entry.second.y << ',' << entry.second.z << ',' 
        << entry.second.ke << ',' << entry.second.pe << ','
        << entry.second.sxx << ',' << entry.second.syy << ',' << entry.second.szz << ',' 
        << entry.second.sxy << ',' << entry.second.sxz << ',' << entry.second.syz << ',' 
        << entry.second.bulk_coordination << ',' << entry.second.coordination << ','
        << entry.second.neighbor_list << ","
        << entry.second.q << ',' << entry.second.dqx << ',' << entry.second.dqy << ',' << entry.second.dqz << ',' 
        << entry.second.Fxx << ',' << entry.second.Fxy << ',' << entry.second.Fxz << ',' 
        << entry.second.Fyx << ',' << entry.second.Fyy << ',' << entry.second.Fyz << ','
        << entry.second.Fzx << ',' << entry.second.Fzy << ',' << entry.second.Fzz << ','  << '\n';
    }

    outputFile.close();
}