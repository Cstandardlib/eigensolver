#ifndef LOBPCG_H
#define LOBPCG_H

/* includes:
    * LOBPCG_CONSTANTS          -- defines constants(tolerance, success/fail flags)
    * logpcg_solve              -- main driver and entry point of LOBPCG
*/

// Include necessary headers
#include <Eigen/Dense>
#include <chrono>
#include <random>

#include "test.h"

// #define LOCKING_BY_TRACE // locking by subspace residual norm, else by individual eigenpairs

namespace LOBPCG_CONSTANTS {
    // 使用static constexpr以确保类型安全和编译时评估
    static constexpr double tol_zero = 1.0e-10;
    static constexpr int success = 0;       // success flag for lobpcg_solve
    static constexpr int fail = 1;          // fail flag for lobpcg_solve
    static constexpr int eig_success = 0;   // success flag for eigensolver
    static constexpr int eig_fail = 1;      // fail flag for eigensolver
    static constexpr double tol_subspace_res = 1.0e-2;
}

// Declare function prototypes
int lobpcg_solve(
    void (*avec)(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs),             // AX, external operator for X nxm
    void (*precnd)(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift),   // TX, external operator for X nxm
    void (*bvec)(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs),             // BX, external operator for X nxm
    Eigen::VectorXd& eig,   // lambdas, should be allocated size n_max_subspace
    Eigen::MatrixXd& evec,  // X, should be allocated size (n,n_max_subspace)
    int n,                  // size of A, B
    int n_eigenpairs,       // X(n, n_eigenpairs)
    int n_max_subspace,     // maximum size of the search subspace
    bool solving_generalized, // true for generalized, B=I if false
    int max_iter,           // maximum number of iterations
    double tol,             // tolerance for convergence test
    double shift,           // shift for Cholesky & precnd
    bool verbose            // print intermediate results if true
);
// template <typename Derived>
// void lobpcg_solve(
//     void (*avec)(int n, int m, Eigen::MatrixBase<Derived>& vecs, Eigen::MatrixBase<Derived> avecs),             // AX, external operator for X nxm
//     void (*precnd)(int n, int m, Eigen::MatrixBase<Derived>& vecs, Eigen::MatrixBase<Derived>& tvecs),   // TX, external operator for X nxm
//     void (*bvec)(int n, int m, Eigen::MatrixBase<Derived>& vecs, Eigen::MatrixBase<Derived>& bvecs),             // BX, external operator for X nxm
//     Eigen::VectorXd& eig,   // lambdas, should be allocated size n_max_subspace
//     Eigen::MatrixXd& evec,  // X, should be allocated size (n,n_max_subspace)
//     int n,                  // size of A, B
//     int n_eigenpairs,       // X(n, n_eigenpairs)
//     int n_max_subspace,     // maximum size of the search subspace
//     bool solving_generalized, // true for generalized, B=I if false
//     int max_iter,           // maximum number of iterations
//     double tol,             // tolerance for convergence test
//     double shift,           // shift for Cholesky & precnd
//     bool& ok,               // true if converged
//     bool verbose            // print intermediate results if true
// )


// moved to ortho.h
// basic tool functions -- linear algebra routines

// ortho
// void ortho(int n, int m, Eigen::MatrixXd& u);
// void ortho_qr(int n, int m, Eigen::MatrixXd &u);

// void b_ortho(int n, int m, Eigen::MatrixXd& u, Eigen::MatrixXd& bu);
// orthogonalize x against given orthonormal set y and orthonormalize x
// void ortho_against_y(int n, int m, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& y);

// void b_ortho_against_y(int n, int m, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& y, Eigen::MatrixXd& by);


// moved to utils.h
// intermediate tool functions
// check_init_guess
// void check_init_guess(int n, int m, Eigen::MatrixXd& evec);
// // timing tools
// std::chrono::time_point<std::chrono::high_resolution_clock> get_current_time();

// // dense eigen solver
// int selfadjoint_eigensolver(Eigen::MatrixXd& A_to_be_eigenvecs,
//  Eigen::VectorXd &eig,
//  int n // dim of A
// );


#endif // LOBPCG_H