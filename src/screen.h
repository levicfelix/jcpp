#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <chrono>

#ifndef LOGGER_H
#define LOGGER_H
class Logger {
private:
    std::string message;

public:
    Logger() {
        (*this)();
        (*this)("              _____   ____                  ");
        (*this)("             |_   _| | ___| _      _        ");
        (*this)("               | |  | |   _| |_  _| |_      ");
        (*this)("            _  | |  | |  |_   _||_   _|     ");
        (*this)("           | |_| |  | |___ |_|    |_|       ");
        (*this)("            |___|    |____|                 ");
        (*this)();
        (*this)(" Version: 10Dec2023 ");
        (*this)(" Author: Levi Felix, PhD ");
        (*this)(" Git address: https://github.com/levicfelix ");
        (*this)(" Contact info: levi.felix@rice.edu ");
        (*this)("               levifelix1@gmail.com ");
        (*this)();
    }

    void operator()(const std::string& msg = "") {
        message = msg;
        std::cout << std::setw(20) << std::left << message << std::endl;
    }

    const std::string& getMessage() const {
        return message;
    }
};
#endif // LOGGER_H

// Define a multi-line string using raw string literals
const char* helpMessage = R"(
Program usage:

    ./lattmech <filename>.jparams

where a control file (.jparams) should be in this folder.
)";

// Print in column format
void printScreenColumns(int step, std::vector<double> ctpos, std::vector<double> J){

    std::cout << "  " << step << "   " 
                      << ctpos[0] << "   " 
                      << ctpos[1] << "   " 
                      << J[0] << "   " 
                      << J[1] << "   " 
                      << J[2] << std::endl; //" \n";
}

// Print elapsed time at the end of the run
void printElapsedTime(std::chrono::_V2::system_clock::time_point t0, std::chrono::_V2::system_clock::time_point t) {

    auto HRS = std::chrono::duration_cast<std::chrono::hours>(t - t0).count();
    auto TOTMINS = std::chrono::duration_cast<std::chrono::minutes>(t - t0).count();
    auto TOTSECS = std::chrono::duration_cast<std::chrono::seconds>(t - t0).count();
    
    auto MINS = int(TOTMINS - HRS*60);
    auto SECS = int(TOTSECS - HRS*3600 - TOTMINS*60);

    // Print total walltime to screen
    std::cout << " Total time (hh:mm:ss): " << HRS  << ":" 
                                            << MINS << ":" 
                                            << SECS << std::endl;
}