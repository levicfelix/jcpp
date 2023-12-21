#include <vector>
#include <string>
#include <unordered_map>

// Function to format a double with fixed-point notation and a specified precision
std::string formatFloatToString(float value, int precision = 1) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

// Convert string to integer vector
std::vector<int> stringToVectorInt(const std::string& input) {
    std::istringstream iss(input);
    std::vector<int> result;

    int value;
    while (iss >> value) {
        result.push_back(value);
    }

    return result;
}

// Convert string to float vector
std::vector<float> stringToVectorFloat(const std::string& input) {
    std::istringstream iss(input);
    std::vector<float> result;

    float value;
    while (iss >> value) {
        result.push_back(value);
    }

    return result;
}

// Convert string to double precision vector
std::vector<double> stringToVectorDouble(const std::string& input) {
    std::istringstream iss(input);
    std::vector<double> result;

    double value;
    while (iss >> value) {
        result.push_back(value);
    }

    return result;
}

// Convert integer vector to string
std::string vectorIntToString(const std::vector<int>& myVector) {
    std::ostringstream oss;
    for (size_t i = 0; i < myVector.size(); ++i) {
        oss << myVector[i];
        if (i < myVector.size() - 1) {
            oss << " ";  // Add a space if it's not the last element
        }
    }
    return oss.str();
}

// Convert vector containing neighbor cuttofs to matrix of Ntype x Ntypes
// used to define both coordination number and crack tip position
 std::vector<std::vector<double>> stringToMatrix(const std::string& input, int rows, int cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols, 0.0));

    std::istringstream iss(input);
    double value;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols && iss >> value; ++j) {
            matrix[i][j] = value;
        }
    }

    return matrix;
}