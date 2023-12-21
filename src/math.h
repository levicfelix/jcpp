// Calculate distance between two atoms
double calculateDistance(double x, double y) {
    return std::sqrt(x * x + y * y);
}

// Calculate the relative vector between two atoms
std::vector<double> calculateRelativeVector(std::vector<double> r1, std::vector<double> r2, double lz, bool pbc, int dim) {
    double dx = r1[0]-r2[0];
    double dy = r1[1]-r2[1];
    double dz = r1[2]-r2[2];

    //if (dim == 2) {
        // Calculate 2D distance between atoms
    //    return {dx, dy};
    //} else 
    //if (dim == 3) {
        if (pbc) {
            // Apply minimum image convention for periodic boundary conditions
            dz = dz - std::round(dz / lz) * lz;
        }
        // Calculate 3D distance
        return {dx, dy, dz};
    //} else {
        // Invalid dimensionality, print error message and exit
        //std::cerr << "ERROR: Invalid system dimensionality!" << std::endl;
        //std::exit(EXIT_FAILURE);
    //}
}

// Calculate neighbor distance based on dimensionality and PBC
double neighborDistance(const Atom& atom1, const Atom& atom2, int dim, bool pbc, double lz) {
    //if (dim == 2) {
        // Calculate 2D distance between atoms
    //    return std::sqrt(std::pow(atom1.x - atom2.x, 2) + std::pow(atom1.y - atom2.y, 2));
    //} else if (dim == 3) {
        double dx = atom1.x - atom2.x;
        double dy = atom1.y - atom2.y;
        double dz = atom1.z - atom2.z;

        if (pbc) {
            // Apply minimum image convention for periodic boundary conditions
            dz = dz - std::round(dz / lz) * lz;
        }

        // Calculate 3D distance
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    //} else {
        // Invalid dimensionality, print error message and exit
    //    std::cerr << "ERROR: Invalid system dimensionality!" << std::endl;
    //    std::exit(EXIT_FAILURE);
    //}
}

// Define identity matrix
std::vector<std::vector<double>> identityMatrix = {{1.0, 0.0, 0.0},
                                                   {0.0, 1.0, 0.0},
                                                   {0.0, 0.0, 1.0}};

// Define matrix of zeroes
std::vector<std::vector<double>> zeroesMatrix = {{0.0, 0.0, 0.0},
                                               {0.0, 0.0, 0.0},
                                               {0.0, 0.0, 0.0}};

// Define matrix of ones
std::vector<std::vector<double>> onesMatrix = {{1.0, 1.0, 1.0},
                                               {1.0, 1.0, 1.0},
                                               {1.0, 1.0, 1.0}};


// Function to calculate the inverse of a matrix
std::vector<std::vector<double>> inverseMatrix(const std::vector<std::vector<double>>& matrix) {
    int matrixSize = matrix.size();

    // Function to calculate the inverse of a 2x2 matrix
    auto inverseMatrix2x2 = [](double a11, double a12, double a21, double a22,
                               double& invA11, double& invA12, double& invA21, double& invA22) {
        double det = a11 * a22 - a12 * a21;

        if (det != 0.0) {
            double invDet = 1.0 / det;
            invA11 = a22 * invDet;
            invA12 = -a12 * invDet;
            invA21 = -a21 * invDet;
            invA22 = a11 * invDet;
        }
    };

    // Function to calculate the inverse of a 3x3 matrix
    auto inverseMatrix3x3 = [](double a11, double a12, double a13,
                               double a21, double a22, double a23,
                               double a31, double a32, double a33,
                               double& invA11, double& invA12, double& invA13,
                               double& invA21, double& invA22, double& invA23,
                               double& invA31, double& invA32, double& invA33) {
        double det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);

        if (det != 0.0) {
            double invDet = 1.0 / det;
            invA11 = (a22 * a33 - a23 * a32) * invDet;
            invA12 = (a13 * a32 - a12 * a33) * invDet;
            invA13 = (a12 * a23 - a13 * a22) * invDet;
            invA21 = (a23 * a31 - a21 * a33) * invDet;
            invA22 = (a11 * a33 - a13 * a31) * invDet;
            invA23 = (a13 * a21 - a11 * a23) * invDet;
            invA31 = (a21 * a32 - a22 * a31) * invDet;
            invA32 = (a12 * a31 - a11 * a32) * invDet;
            invA33 = (a11 * a22 - a12 * a21) * invDet;
        }
    };

    // Check if the matrix size is valid
    if (matrixSize != 2 && matrixSize != 3) {
        std::cerr << "Invalid matrix size. Please enter 2 or 3." << std::endl;
        return {};
    }

    // Calculate the inverse
    std::vector<std::vector<double>> invMatrix(matrixSize, std::vector<double>(matrixSize, 0.0));
    if (matrixSize == 2) {
        inverseMatrix2x2(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1],
                          invMatrix[0][0], invMatrix[0][1], invMatrix[1][0], invMatrix[1][1]);
    } else {
        inverseMatrix3x3(matrix[0][0], matrix[0][1], matrix[0][2],
                          matrix[1][0], matrix[1][1], matrix[1][2],
                          matrix[2][0], matrix[2][1], matrix[2][2],
                          invMatrix[0][0], invMatrix[0][1], invMatrix[0][2],
                          invMatrix[1][0], invMatrix[1][1], invMatrix[1][2],
                          invMatrix[2][0], invMatrix[2][1], invMatrix[2][2]);
    }

    return invMatrix;
}

// Function to multiply two matrices
std::vector<std::vector<double>> matrixMultiply(const std::vector<std::vector<double>>& A,
                                                const std::vector<std::vector<double>>& B) {
    int m = A.size();
    int n = B[0].size();
    int p = B.size();

    std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < p; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

