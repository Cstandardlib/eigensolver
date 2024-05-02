#include "utils.h"
#include "lobpcg.h"
#include "ortho.h"
#include <Eigen/Dense>

// intermediate tool functions
// check_init_guess

// check_guess - 检查给定的矩阵evec是否为零矩阵或者是否正交。
// 如果evec为零矩阵，将用随机数填充并正交化。
// 如果evec不正交，将进行正交化处理。
/**
 * @brief Checks the validity of the given guess matrix and performs orthogonalization if necessary.
 * 
 * If evec is zero, provide a random guess, filled by uniform [0,1).
 * 
 * @param n The number of rows in the guess matrix.
 * @param m The number of columns in the guess matrix.
 * @param evec The guess matrix to be checked and orthogonalized if needed.
 */
void check_init_guess(int n, int m, Eigen::MatrixXd& evec) {
    double enorm = evec.norm(); // 计算矩阵的范数
// #ifdef DEBUG_LOBPCG
//     std::cout << "-- start check_init_guess --" << std::endl;
// #endif
    if (std::abs(enorm) < LOBPCG_CONSTANTS::tol_zero) {
        // if evec is all zero(which means no initial guess is provided by user), fill it with random values in U[0,1)
        std::default_random_engine engine;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n*m; ++i) evec(i) = dist(engine);
        ortho(n, m, evec);// orthogonalize evec
    } else {
        // compute overlao matrix M = X^T X (size m*m), and do ortho
        // first check for orthogonality: whether X^TX == Identity
        Eigen::MatrixXd overlap = evec.transpose() * evec;
        double diag_norm = overlap.diagonal().array().square().sum();   // square sum of overlap diagonal
        double out_norm = (overlap.array().square()).sum() - diag_norm; // square sum of non-diagonal elements
        
        // 检查对角线元素是否为1，非对角线元素是否为0
        if (std::abs(diag_norm - m) > 1e-10 || std::abs(out_norm) > 1e-10) {
            // 如果不正交，进行正交化处理
            ortho(n, m, evec);
        }
    }
// #ifdef DEBUG_LOBPCG
//     std::cout << "-- end check_init_guess --" << std::endl;
// #endif
}


// timing tools
/*
 * wrap for std::chrono::high_resolution_clock::now();
 * auto start = get_current_time();
 * auto end = get_current_time();
 * std::chrono::duration<double> duration = end - start;
 */
std::chrono::time_point<std::chrono::high_resolution_clock> get_current_time() {
    return std::chrono::high_resolution_clock::now();
}


// dense eigen solver
/* overwrite A and eig
 * solves eigenvalue and eigenvectors of a selfadjoint(i.e. inverse equal to its adjoint) matrix A.
 * solve A_nn(n, n) topleft corner of size n input
*/


int selfadjoint_eigensolver(Eigen::MatrixXd &A_to_be_eigenvecs, Eigen::VectorXd &eig, int n)
{
#ifdef DEBUG_LOBPCG
    std::cout << "--- starts eigensolver ---" << std::endl;
    std::cout << "A size: " << A_to_be_eigenvecs.size() <<", ("<<A_to_be_eigenvecs.rows()<<", "<<A_to_be_eigenvecs.cols()<<")" <<std::endl;
#endif
    if(A_to_be_eigenvecs.rows() != n || A_to_be_eigenvecs.cols()!=n){
        std::cerr << "input matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"<<std::endl; return LOBPCG_CONSTANTS::eig_fail;
    }
    if(eig.size() < n){
        std::cerr << "eigenvalues vector must not be smaller than n = "<<n<<std::endl; return LOBPCG_CONSTANTS::eig_fail;
    }
#ifdef DEBUG_EIGENSOLVER
    std::cout << "A size: " << A_to_be_eigenvecs.size() <<", ("<<A_to_be_eigenvecs.rows()<<", "<<A_to_be_eigenvecs.cols()<<")" <<std::endl;
    // std::cout << "original A = \n" << A_to_be_eigenvecs << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "solving A(" << n << "," << n << ") matrix" << std::endl;
    // Eigen::MatrixXd A_nn = A_to_be_eigenvecs.topLeftCorner(n,n);// topleft corner of A_to_be_eigenvecs
    std::cout << "A to be solved = \n" << A_to_be_eigenvecs << std::endl;
#endif
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A_to_be_eigenvecs);
    if (eigensolver.info() == Eigen::Success)
    {
        eig.head(n) = eigensolver.eigenvalues();
        // A_to_be_eigenvecs.topLeftCorner(n,n) = eigensolver.eigenvectors();
        A_to_be_eigenvecs = eigensolver.eigenvectors();
#ifdef DEBUG_EIGENSOLVER
        std::cout << "The eigenvalues are:\n" << eig << std::endl;
        std::cout << "The eigenvectors are:\n" << A_to_be_eigenvecs << std::endl;
#endif
    }
    else
    {
        std::cerr << "Eigen decomposition failed." << std::endl;
        return LOBPCG_CONSTANTS::eig_fail;
    }
    return LOBPCG_CONSTANTS::eig_success;
#ifdef DEBUG_LOBPCG
    std::cout << "--- end eigensolver ---" << std::endl;
#endif
}