#include "lobpcg.h" // lobpcg_solve
#include "matvec.h" // avec, bvec, precnd
#include "utils.h"  //  selfadjoint_eigensolver

#include <iostream>
#include <fstream>
#include <fast_matrix_market/app/Eigen.hpp> // construct sparse from mtx

void _no_matvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& still_vecs){
    // do nothing but copy, still_vecs = I*vesc
    still_vecs = vecs;
}

void _diag_2_matvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& two_vecs){
    // diagonal 2 diagonal matrix-vector product
    two_vecs = 3 * vecs;
}

// sparseA.mtx 20*20 diag
void test_sparse_diag_A() {
    int n = 20;//100;
    int n_eigenpairs = 2;//5;
    int n_max_subspace = 4;//10; // n_eig = min(2*n_want, n_want + 5)
    bool solving_generalized = false;
    int max_iter = 10;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;

    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(sparseAvec, precnd, bvec, eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);
    // lobpcg_solve(avec, precnd, bvec, eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, ok, verbose);

    // Perform assertions to check if the results are correct
    // ...
    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
}



// a1.mtx 9*9 diag
void test_a1_9() {
    int n = 9;//100;
    int n_eigenpairs = 2;//5;
    int n_max_subspace = 3;//10; // n_eig = min(2*n_want, n_want + 5)
    bool solving_generalized = false;
    int max_iter = 100;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;

    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        _diag_2_matvec,/*_diag_2_matvec*//*a1vec*/
        _no_matvec,
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig << std::endl;
    // std::cout << "eigenvectors = \n"<< evec << std::endl;

    std::ifstream f("../../../a1.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    Eigen::SparseMatrix<double> mat;
    fast_matrix_market::read_matrix_market_eigen(f, mat); // std::cout << mat << std::endl;
    mat.diagonal() <<1,2,3,4,5,6,7,8,9;
    //2,2,2,2,2,2,2,2,2;
    Eigen::MatrixXd A_dense = mat.toDense();
    // std::cout << "A_dense = \n"<< A_dense << std::endl;
    Eigen::VectorXd eig_real(9);
    selfadjoint_eigensolver(A_dense, eig_real, 9);
    std::cout << "eig_real = \n"<< eig_real.transpose() << std::endl;
}



int main(){
    // test_sparse_diag_A();
    test_a1_9();
    return 0;
}