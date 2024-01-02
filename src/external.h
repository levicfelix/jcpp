#include <tuple>
#include <algorithm>  // for std::find

// Function to split a string into tokens
std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

namespace lammps {

    // Function to convert LAMMPS units to SI units
    std::vector<double> Units(const std::string& units) {
        double distance, energy, stress;

        if (units == "metal") {
            distance = 1.0e-10;  // Å to m
            energy = 1.6e-19;    // eV to J
            stress = 1e-25;      // bar*Å³ to J
        } else if (units == "real") {
            distance = 1.0e-10;                     // Å to m
            energy = 0.0433641153087705 * 1.6e-19;  // kcal/mol to J
            stress = 9.869232667160128e-26;         // atm*Å³ to J
        } else {
            std::cerr << "ERROR: LAMMPS units system not implemented yet or does not exist!" << std::endl;
            exit(EXIT_FAILURE);
        }

        return {distance, energy, stress};
    }

    // Function to save atom data for a given timestep to a CSV file
    void saveCSVTimestep(const std::string& outdir, const std::vector<std::vector<std::string>>& atomData, int timestep, const std::vector<std::string>& header, const std::string& prefix, const std::vector<std::string>& customLabels) {
        std::ofstream csvFile(outdir + prefix + "_" + std::to_string(timestep) + ".csv");
        if (!csvFile.is_open()) {
            std::cerr << "ERROR: Could not open CSV file." << std::endl;
            exit(EXIT_FAILURE);
            return;
        }

        // Write header with "index" as the first column
        csvFile << "index";
        for (size_t i = 0; i < header.size(); ++i) {
            csvFile << ",";

            // If the current header label is one of the specified ones, use the corresponding custom label
            if (i >= 5 && i < 14) {
                size_t customIndex = i - 5;  // Adjust the index for customLabels
                if (customIndex < customLabels.size()) {
                    csvFile << customLabels[customIndex];
                    continue;  // Move to the next iteration
                }
            }

            // Otherwise, use the original header label
            csvFile << header[i];
        }
        csvFile << "\n";

        // Write atom data with index values
        for (size_t index = 0; index < atomData.size(); ++index) {
            csvFile << index + 1;  // Adding 1 to start the index from 1
            for (size_t i = 0; i < atomData[index].size(); ++i) {
                csvFile << "," << atomData[index][i];
            }
            csvFile << "\n";
        }
        std::cout << "\r    CSV file saved for timestep: " << timestep << std::flush;
    }

    // Function to convert LAMMPS dump to a series of CSV files
    void saveCSV(const std::string& outdir, const std::string& lmptrj, int instep, int stpmin, int stpmax, const std::string& prefix, const std::vector<std::string>& customLabels) {

        std::ifstream inputFile(lmptrj);

        if (!inputFile.is_open()) {
            std::cerr << "  File not found! ("+lmptrj+")" << std::endl;
            exit(EXIT_FAILURE);
            return;
        }

        int targetTimestep = instep;

        std::string line;
        std::vector<std::vector<std::string>> atomData;
        int currentTimestep = -1;

        // Define header manually
        std::vector<std::string> header;
        std::ofstream jstepsfile(outdir+"jsteps_tmp.txt"); // Stores all timesteps use to calculate J
        while (std::getline(inputFile, line)) {
            if (line.find("ITEM: TIMESTEP") != std::string::npos) {
                // Update current timestep
                std::getline(inputFile, line);
                currentTimestep = std::stoi(line);

                // Read and skip the header line
                std::getline(inputFile, line);

                // Read the number of atoms
                int numAtoms;
                std::getline(inputFile, line);
                numAtoms = std::stoi(line);

                // Read the simulation box dimensions
                std::getline(inputFile, line);
                std::getline(inputFile, line);
                std::getline(inputFile, line);
                std::getline(inputFile, line);

                // Read the manually defined header line
                std::getline(inputFile, line);
                header = split(line.substr(line.find("ITEM: ATOMS") + 12));

                // Read atom data for the current timestep
                for (int i = 0; i < numAtoms; ++i) {
                    std::getline(inputFile, line);
                    atomData.push_back(split(line));
                }

                // Save the CSV file for the current timestep
                if ((instep == -1 && currentTimestep >= stpmin && currentTimestep <= stpmax)) {
                    jstepsfile << currentTimestep << "\n";
                    saveCSVTimestep(outdir, atomData, currentTimestep, header, prefix, customLabels);
                }
                else if (currentTimestep == targetTimestep) {
                    jstepsfile << currentTimestep << "\n";
                    saveCSVTimestep(outdir, atomData, currentTimestep, header, prefix, customLabels);
                    break;
                }
                else if (instep == -1 && currentTimestep >= stpmax) {
                    break;
                }

                // Clear the atomData vector for the next timestep
                atomData.clear();
            }
        }
        std::cout << std::endl;  // To make sure the last line keep appearing on terminal screen

        jstepsfile.close();
        inputFile.close();
        return;
    }
}  // end namespace lammps