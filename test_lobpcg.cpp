#include "lobpcg.h" // lobpcg_solve
#include "matvec.h" // avec, bvec, precnd
#include "sparseMv.h"
// #include "utils.h"  //  selfadjoint_eigensolver

#include <iostream>
#include <fstream>
#include <fast_matrix_market/app/Eigen.hpp> // construct sparse from mtx

// identity matrix-vector product
void _no_matvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& still_vecs){
    // do nothing but copy, still_vecs = I*vesc
    still_vecs.topLeftCorner(n,m) = vecs.topLeftCorner(n,m);
}

void _no_precnd(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& still_vecs, double shift){
    // do nothing but copy, still_vecs = I*vesc
    still_vecs.topLeftCorner(n,m) = vecs.topLeftCorner(n,m);
}

void _diag_matvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& after_vecs){
    // diagonal matrix-vector product
    Eigen::SparseMatrix<double> mat(9,9);
    mat.diagonal() <<5,5,5,5,5,5,5,5,5;
    after_vecs = mat * vecs;
}

// sparseA.mtx 20*20 diag
void test_sparse_diag_A() {
    int n = 20;//100;
    int n_eigenpairs = 3;//5;
    int n_max_subspace = 6;//10; // n_eig = min(2*n_want, n_want + 5)
    bool solving_generalized = false;
    int max_iter = 10;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;

    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        sparseAvec,/*_diag_matvec*//*a1vec*/
        _no_precnd,
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);
    // lobpcg_solve(avec, precnd, bvec, eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, ok, verbose);

    // Perform assertions to check if the results are correct
    // ...
    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig << std::endl;
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
        a1vec,/*_diag_matvec*//*a1vec*/
        _no_precnd,/*_no_matvec*/ /*mprec*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
    std::cout << "LOBPCG eigenvectors = \n"<< evec.leftCols(n_eigenpairs) << std::endl;
}

void test_large_1000(){
    int n = 1000;
    int n_eigenpairs = 20;//5;
    int n_max_subspace = std::min(2*n_eigenpairs, n_eigenpairs + 5);
    bool solving_generalized = false;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        avec,/*_diag_matvec*//*a1vec*/
        _no_precnd,/*_no_matvec*/ /*mprec*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}
void test_large_5000(){
    int n = 5000;
    int n_eigenpairs = 200;//5;
    int n_max_subspace = std::min(2*n_eigenpairs, n_eigenpairs + 25);
    bool solving_generalized = true;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        avec,/*_diag_matvec*//*a1vec*/
        mprec,/*_no_matvec*/ /*mprec*/
        bvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void test_gen_9() {
    int n = 9;//100;
    int n_eigenpairs = 2;//5;
    int n_max_subspace = 3;//10; // n_eig = min(2*n_want, n_want + 5)
    bool solving_generalized = true; //false;
    int max_iter = 100;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;

    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        a1vec,/*_diag_matvec*//*a1vec*/
        _no_precnd,/*_no_matvec*/ /*mprec*/
        b1vec,/*bvec*//*b1vec*/
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void test_gen_1000(){
    int n = 1000;
    int n_eigenpairs = 20;//5;
    int n_max_subspace = std::min(2*n_eigenpairs, n_eigenpairs + 5);
    bool solving_generalized = true;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        avec,/*_diag_matvec*//*a1vec*/
        mprec,/*_no_matvec*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void run_dense_Si2(){
    int n = 769;
    int n_eigenpairs = 8;//5;
    int n_max_subspace = std::min(2*n_eigenpairs, n_eigenpairs + 5);
    bool solving_generalized = false;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        avec_Si2,/*_diag_matvec*//*a1vec*/
        precnd_Si2,/*_no_matvec*/ /*mprec*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void run_sparse_Si2(){
    int n = 769;
    int n_eigenpairs = 2;//5;
    int n_max_subspace = std::min(2*n_eigenpairs, n_eigenpairs + 5);
    bool solving_generalized = false;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        sparse_avec_Si2,/*_diag_matvec*//*a1vec*/
        _no_precnd,/*_no_matvec*/ /*mprec*//*sparse_precnd_Si2*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void run_sparse_Na5(){
    int n = 5832;
    int n_eigenpairs = 1;//5;
    int n_max_subspace = 2;//std::min(2*n_eigenpairs, n_eigenpairs + 5);
    bool solving_generalized = false;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        sparse_avec_Na5,/*_diag_matvec*//*a1vec*/
        _no_precnd,/*_no_matvec*/ /*mprec*/ /*sparse_precnd_Na5*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void run_sparse_Si5H12(){
    int n = 19896;
    int n_eigenpairs = 199;//200;//5;
    int n_max_subspace = std::min(2*n_eigenpairs, n_eigenpairs + 5); //225;
    bool solving_generalized = false;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        sparse_avec_Si5H12,/*_diag_matvec*//*a1vec*/
        _no_precnd,/*_no_precnd*/ /*mprec*/ /*sparse_precnd_Si5H12*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

void run_sparse_Ga3As3H12(){
    int n = 61349;
    int n_eigenpairs = 5;//5;
    int n_max_subspace = 10;//std::min(2*n_eigenpairs, n_eigenpairs + 50);
    bool solving_generalized = false;
    int max_iter = 1000;
    double tol = 1e-6;
    double shift = 0.0;
    bool verbose = true;
    Eigen::VectorXd eig(n_max_subspace);
    Eigen::MatrixXd evec(n, n_max_subspace);
    eig.setZero(); evec.setZero();
    int ok = lobpcg_solve(
        sparse_avec_Ga3As3H12,/*_diag_matvec*//*a1vec*/
        sparse_precnd_Ga3As3H12,/*_no_matvec*/ /*mprec*//*sparse_precnd_Ga3As3H12*/
        _no_matvec/*bvec*/,
        eig, evec, n, n_eigenpairs, n_max_subspace, solving_generalized, max_iter, tol, shift, verbose);

    if(ok != LOBPCG_CONSTANTS::success) std::cerr<< "not ok! "<< std::endl;
    else std::cout << "ok! converged "<< std::endl;
    std::cout << "------- final -------" << std::endl;
    std::cout << "LOBPCG eigenvalues = \n"<< eig.head(n_eigenpairs) << std::endl;
}

int main(){
    // test_sparse_diag_A();
    // test_a1_9();
    // test_large_1000();
    // test_large_5000();
    // test_gen_9();
    // run_dense_Si2();
    // run_sparse_Si2();
    // run_sparse_Na5();
    // run_sparse_Si5H12();
    run_sparse_Ga3As3H12();
    return 0;
}