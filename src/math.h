// Calculate absolute value
double absv(double x) {
    return std::sqrt(x*x);
}

// Calculate norm of vector
double norm (const std::vector<double> vec) {
    return std::sqrt( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

// Function to normalize a vector
std::vector<double> normalize(const std::vector<double> vec) {
    return {vec[0] / norm(vec), vec[1] / norm(vec), vec[2] / norm(vec)};
}

// Function to project a vector B into the direction of another vector A
std::vector<double> project(const std::vector<double>& A, const std::vector<double>& B) {
    double dotProduct = 0.0;
    for (size_t i = 0; i < A.size(); ++i) {
        dotProduct += A[i] * B[i];
    }

    std::vector<double> projection;
    for (double component : A) {
        projection.push_back(dotProduct / (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]) * component);
    }

    return projection;
}

// Funtion to 
std::vector<double> orthogonalize(const std::vector<double>& A, const std::vector<double>& B) {
    std::vector<double> projection = project(A, B);

    std::vector<double> result;
    for (size_t i = 0; i < A.size(); ++i) {
        result.push_back(B[i] - projection[i]);
    }

    return result;
}

// Uses the Gran-Schmidt process to construct an orthogonal vector to an inputVector
std::vector<double> createOrthogonalVector(const std::vector<double>& A) {
    if (A.size() != 3) {
        std::cerr << "Input vector must have exactly three components.\n";
        return {};
    }

    std::vector<double> orthogonalVector = orthogonalize(A, {1.0,1.0,1.0}); // Choose any arbitrary vector ({1,1,1} in this case)
    return normalize(orthogonalVector);
}

// Function to calculate the cross product of two vectors
std::vector<double> cross(const std::vector<double> vec1, const std::vector<double> vec2) {
    std::vector<double> vec3 = {vec1[1] * vec2[2] - vec1[2] * vec2[1],
                                vec1[2] * vec2[0] - vec1[0] * vec2[2],
                                vec1[0] * vec2[1] - vec1[1] * vec2[0]};
    return vec3;
}

std::vector<double> transformVector(const std::vector<std::vector<double>>& matrix, const std::vector<double> vector) {
    // Check if the matrix and vector are compatible for multiplication
    size_t matrixSize = matrix.size();
    if (matrixSize == 0 || matrix[0].size() != matrixSize || vector.size() != matrixSize) {
        std::cerr << "Error: Incompatible matrix and vector sizes for multiplication." << std::endl;
        return {};
    }

    // Perform matrix-vector multiplication
    std::vector<double> result(matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; ++i) {
        for (size_t j = 0; j < matrixSize; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

// Calculate distance between two atoms
double calculateDistance(double x, double y) {
    return std::sqrt(x * x + y * y);
}

// Calculate the relative vector between two atoms
std::vector<double> calculateRelativeVector(std::vector<double> r1, std::vector<double> r2, double lz, bool pbc) {
    double dx = r1[0]-r2[0];
    double dy = r1[1]-r2[1];
    double dz = r1[2]-r2[2];

    if (pbc) {
        // Apply minimum image convention for periodic boundary conditions
        dz = dz - std::round(dz / lz) * lz;
    }
    // Calculate 3D distance
    return {dx, dy, dz};
}

// Calculate neighbor distance based on dimensionality and PBC
double neighborDistance(const Atom& atom1, const Atom& atom2, bool pbc, double lz) {
    double dx = atom1.x - atom2.x;
    double dy = atom1.y - atom2.y;
    double dz = atom1.z - atom2.z;

    if (pbc) {
        // Apply minimum image convention for periodic boundary conditions
        dz = dz - std::round(dz / lz) * lz;
    }

    // Return neighbors distance
    return std::sqrt(dx * dx + dy * dy + dz * dz);
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

// Function to calculate the determinant of a matrix
double determinantMatrix(const std::vector<std::vector<double>>& matrix) {
    int matrixSize = matrix.size();
    double det;
    // Function to calculate the determinant of a 2x2 matrix
    if (matrixSize==2) {
        det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    } else if (matrixSize==3) {
        double block1 = matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];
        double block2 = matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0];
        double block3 = matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0];
        det = matrix[0][0] * block1 - matrix[0][1] * block2 + matrix[0][2] * block3;
    }
    return det;
}

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
        } else {
            std::cout << "ERROR: Singular matrix found!" << std::endl;
            std::exit(EXIT_FAILURE);
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
        } else {
            std::cout << "Singular matrix found at the computation of F!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    };

    // Check if the matrix size is valid
    if (matrixSize != 2 && matrixSize != 3) {
        std::cerr << "  ERROR: Invalid matrix size. Exiting..." << std::endl;
        std::exit(EXIT_FAILURE);
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

// Function to ptoject 1D or 2D matrices to 3D ones
std::vector<std::vector<double>> matrixIn3D(std::vector<std::vector<double>>& m) {
    size_t N = m.size();

    std::vector<std::vector<double>> C(3, std::vector<double>(3, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i][j] = m[i][j];
        }
    }
    return C;
}

// Function to create a reduced matrix by excluding a specified row and column
std::vector<std::vector<double>> reduceMatrix(const std::vector<std::vector<double>>& matrix, int excludeRow) {
    // Check if the matrix is a 3x3 matrix
    if (matrix.size() != 3 || matrix[0].size() != 3 || matrix[1].size() != 3 || matrix[2].size() != 3) {
        std::cerr << "Error: Input matrix must be a 3x3 matrix." << std::endl;
        return {};
    }

    // Create a 2x2 matrix to store the reduced matrix
    std::vector<std::vector<double>> reducedMatrix(2, std::vector<double>(2, 0));

    // Copy elements to the reduced matrix while excluding the specified row and column
    int newRow = 0;
    for (int i = 0; i < 3; ++i) {
        // Skip the excluded row
        if (i == excludeRow) {
            continue;
        }

        int newCol = 0;
        for (int j = 0; j < 3; ++j) {
            // Skip the excluded column
            if (j == excludeRow) {
                continue;
            }

            // Copy the element to the reduced matrix
            reducedMatrix[newRow][newCol] = matrix[i][j];
            ++newCol;
        }

        ++newRow;
    }

    return reducedMatrix;
}

// Function to print a matrix
void printMatrix(const std::vector<std::vector<double>>& matrix) {
    // Get the number of rows and columns in the matrix
    size_t rows = matrix.size();
    size_t cols = (rows > 0) ? matrix[0].size() : 0;
    
    // Iterate through the matrix and print its elements
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}
