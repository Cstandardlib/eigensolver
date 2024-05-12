#include "sparseMv.h"
#include <iostream>
#include <fstream>      // ifstream
#include <fast_matrix_market/app/Eigen.hpp>

// -----------------------------------------------------------------------
// Na5

Eigen::SparseMatrix<double> sparseA_Na5;
void init_SparseA_Na5(){
    std::ifstream f("../../../Na5.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building sparseA_Na5, ";
    fast_matrix_market::read_matrix_market_eigen(f, sparseA_Na5);
    std::cout <<"size = " << sparseA_Na5.rows() << ", " << sparseA_Na5.cols() << std::endl;
    std::cout << "sparseA_Na5 size = " << sparseA_Na5.rows() << ", " << sparseA_Na5.cols() << std::endl;
}

void sparse_avec_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    if(sparseA_Na5.rows() == 0) init_SparseA_Na5();

    if(sparseA_Na5.rows() != n || sparseA_Na5.cols() != n){std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;}
    if(vecs.rows() != n || vecs.cols() != m){std::cerr << "vecs must be of size (n,m)"; return;}
    if(avecs.rows() != n || avecs.cols() != m){std::cerr << "avecs must be of size (n,m)"; return;}

    avecs = sparseA_Na5 * vecs;
}

void sparse_precnd_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(sparseA_Na5.rows() == 0) init_SparseA_Na5();

    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            // 检查分母是否为0，避免除以0的错误
            // if (abs(a(i, i) + fac) > 1.0e-5) {
            if (abs(sparseA_Na5.coeff(i, i)+shift) > 1.0e-5) {
                tvecs(i, icol) = vecs(i, icol) / (sparseA_Na5.coeff(i, i)+shift);
            }
        }
    }
}

Eigen::SparseMatrix<double> tridiagA_Na5;
void init_tridiagA_Na5(double shift){
    if(sparseA_Na5.rows() == 0) init_SparseA_Na5();

    int n = sparseA_Na5.rows();
    tridiagA_Na5.resize(n, n);
    tridiagA_Na5.reserve(3*n); // reserve space for tridiagonal elements

    for(int i = 0; i < n; i++){
        tridiagA_Na5.insert(i, i) = sparseA_Na5.coeff(i, i)+shift;
        if(i > 0)
            tridiagA_Na5.insert(i, i-1) = sparseA_Na5.coeff(i, i-1);
        if(i < n-1)
            tridiagA_Na5.insert(i, i+1) = sparseA_Na5.coeff(i, i+1);
    }

    tridiagA_Na5.makeCompressed();
}

void tridiagA_precnd_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(tridiagA_Na5.rows() == 0) init_tridiagA_Na5(shift);
    for(int i = 0; i < n; i++){
        tridiagA_Na5.coeffRef(i, i) = sparseA_Na5.coeff(i, i)+shift;
    }
    // solve tridiagA_Na5 * tvecs = vecs // 计算LU分解
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(tridiagA_Na5);
    solver.factorize(tridiagA_Na5);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Decomposition failed when doing tridiagonal precnd" << std::endl;}
    // 使用LU分解来解线性系统
    tvecs = solver.solve(vecs);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Solving failed when doing tridiagonal precnd" << std::endl;}
    // 输出解矩阵
    // std::cout << "The solution matrix tvecs is:\n" << tvecs << std::endl;
}




// -----------------------------------------------------------------------
// Si2




Eigen::SparseMatrix<double> sparseA_Si2;
void init_SparseA_Si2(){
    std::ifstream f("../../../Si2.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building sparseA_Si2, ";
    fast_matrix_market::read_matrix_market_eigen(f, sparseA_Si2);
    std::cout <<"size = " << sparseA_Si2.rows() << ", " << sparseA_Si2.cols() << std::endl;
    std::cout << "sparseA_Si2 size = " << sparseA_Si2.rows() << ", " << sparseA_Si2.cols() << std::endl;
}

void sparse_avec_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    if(sparseA_Si2.rows() == 0) init_SparseA_Si2();

    if(sparseA_Si2.rows() != n || sparseA_Si2.cols() != n){std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;}
    if(vecs.rows() != n || vecs.cols() != m){std::cerr << "vecs must be of size (n,m)"; return;}
    if(avecs.rows() != n || avecs.cols() != m){std::cerr << "avecs must be of size (n,m)"; return;}

    avecs = sparseA_Si2 * vecs;
}

void sparse_precnd_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(sparseA_Si2.rows() == 0) init_SparseA_Si2();

    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            // 检查分母是否为0，避免除以0的错误
            // if (abs(a(i, i) + fac) > 1.0e-5) {
            if (abs(sparseA_Si2.coeff(i, i)+shift) > 1.0e-5) {
                tvecs(i, icol) = vecs(i, icol) / (sparseA_Si2.coeff(i, i)+shift);
            }
        }
    }
}

Eigen::SparseMatrix<double> tridiagA_Si2;
void init_tridiagA_Si2(double shift){
    if(sparseA_Si2.rows() == 0) init_SparseA_Si2();

    int n = sparseA_Si2.rows();
    tridiagA_Si2.resize(n, n);
    tridiagA_Si2.reserve(3*n); // reserve space for tridiagonal elements

    for(int i = 0; i < n; i++){
        tridiagA_Si2.insert(i, i) = sparseA_Si2.coeff(i, i)+shift;
        if(i > 0)
            tridiagA_Si2.insert(i, i-1) = sparseA_Si2.coeff(i, i-1);
        if(i < n-1)
            tridiagA_Si2.insert(i, i+1) = sparseA_Si2.coeff(i, i+1);
    }

    tridiagA_Si2.makeCompressed();
}

void tridiagA_precnd_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(tridiagA_Si2.rows() == 0) init_tridiagA_Si2(shift);
    for(int i = 0; i < n; i++){
        tridiagA_Si2.coeffRef(i, i) = tridiagA_Si2.coeff(i, i)+shift;
    }
    // solve tridiagA_Si2 * tvecs = vecs // 计算LU分解
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(tridiagA_Si2);
    solver.factorize(tridiagA_Si2);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Decomposition failed when doing tridiagonal precnd" << std::endl;}
    // 使用LU分解来解线性系统
    tvecs = solver.solve(vecs);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Solving failed when doing tridiagonal precnd" << std::endl;}
    // 输出解矩阵
    // std::cout << "The solution matrix tvecs is:\n" << tvecs << std::endl;
}


// Eigen::SparseMatrix<double> sparseA_Si5H12;
// void init_SparseA_Si2(){
//     std::ifstream f("../../../Si2.mtx");
//     if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
//     std::cout <<"building sparseA_Si2, ";
//     fast_matrix_market::read_matrix_market_eigen(f, sparseA_Si2);
//     std::cout <<"size = " << sparseA_Si2.rows() << ", " << sparseA_Si2.cols() << std::endl;
//     std::cout << "sparseA_Si2 size = " << sparseA_Si2.rows() << ", " << sparseA_Si2.cols() << std::endl;
// }

// void sparse_avec_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);

// void sparse_precnd_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift);





Eigen::SparseMatrix<double> sparseA_Si5H12;
void init_SparseA_Si5H12(){
    std::ifstream f("../../../Si5H12.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building sparseA_Si5H12, ";
    fast_matrix_market::read_matrix_market_eigen(f, sparseA_Si5H12);
    std::cout <<"size = " << sparseA_Si5H12.rows() << ", " << sparseA_Si5H12.cols() << std::endl;
    std::cout << "sparseA_Si5H12 size = " << sparseA_Si5H12.rows() << ", " << sparseA_Si5H12.cols() << std::endl;
}

void sparse_avec_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    if(sparseA_Si5H12.rows() == 0) init_SparseA_Si5H12();

    if(sparseA_Si5H12.rows() != n || sparseA_Si5H12.cols() != n){std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;}
    if(vecs.rows() != n || vecs.cols() != m){std::cerr << "vecs must be of size (n,m)"; return;}
    if(avecs.rows() != n || avecs.cols() != m){std::cerr << "avecs must be of size (n,m)"; return;}

    avecs = sparseA_Si5H12 * vecs;
}

void sparse_precnd_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(sparseA_Si5H12.rows() == 0) init_SparseA_Si5H12();

    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            if (abs(sparseA_Si5H12.coeff(i, i)+shift) > 1.0e-5) {
                tvecs(i, icol) = vecs(i, icol) / (sparseA_Si5H12.coeff(i, i)+shift);
            }
        }
    }
}


Eigen::SparseMatrix<double> tridiagA_Si5H12;
void init_tridiagA_Si5H12(double shift){
    if(sparseA_Si5H12.rows() == 0) init_SparseA_Si5H12();

    int n = sparseA_Si5H12.rows();
    tridiagA_Si5H12.resize(n, n);
    tridiagA_Si5H12.reserve(3*n); // reserve space for tridiagonal elements

    for(int i = 0; i < n; i++){
        tridiagA_Si5H12.insert(i, i) = sparseA_Si5H12.coeff(i, i)+shift;
        if(i > 0)
            tridiagA_Si5H12.insert(i, i-1) = sparseA_Si5H12.coeff(i, i-1);
        if(i < n-1)
            tridiagA_Si5H12.insert(i, i+1) = sparseA_Si5H12.coeff(i, i+1);
    }

    tridiagA_Si5H12.makeCompressed();
}

void tridiagA_precnd_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(tridiagA_Si5H12.rows() == 0) init_tridiagA_Si5H12(shift);
    for(int i = 0; i < n; i++){
        tridiagA_Si5H12.coeffRef(i, i) = tridiagA_Si5H12.coeff(i, i)+shift;
    }
    // solve tridiagA_Si5H12 * tvecs = vecs // 计算LU分解
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(tridiagA_Si5H12);
    solver.factorize(tridiagA_Si5H12);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Decomposition failed when doing tridiagonal precnd" << std::endl;}
    // 使用LU分解来解线性系统
    tvecs = solver.solve(vecs);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Solving failed when doing tridiagonal precnd" << std::endl;}
    // 输出解矩阵
    // std::cout << "The solution matrix tvecs is:\n" << tvecs << std::endl;
}


















// Ga3As3H12
Eigen::SparseMatrix<double> sparseA_Ga3As3H12;
void init_SparseA_Ga3As3H12(){
    std::ifstream f("../../../Ga3As3H12.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building sparseA_Ga3As3H12, ";
    fast_matrix_market::read_matrix_market_eigen(f, sparseA_Ga3As3H12);
    std::cout <<"size = " << sparseA_Ga3As3H12.rows() << ", " << sparseA_Ga3As3H12.cols() << std::endl;
    std::cout << "sparseA_Ga3As3H12 size = " << sparseA_Ga3As3H12.rows() << ", " << sparseA_Ga3As3H12.cols() << std::endl;
}

void sparse_avec_Ga3As3H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    if(sparseA_Ga3As3H12.rows() == 0) init_SparseA_Ga3As3H12();

    if(sparseA_Ga3As3H12.rows() != n || sparseA_Ga3As3H12.cols() != n){std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;}
    if(vecs.rows() != n || vecs.cols() != m){std::cerr << "vecs must be of size (n,m)"; return;}
    if(avecs.rows() != n || avecs.cols() != m){std::cerr << "avecs must be of size (n,m)"; return;}

    avecs = sparseA_Ga3As3H12 * vecs;
}

void sparse_precnd_Ga3As3H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(sparseA_Ga3As3H12.rows() == 0) init_SparseA_Ga3As3H12();

    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            if (abs(sparseA_Ga3As3H12.coeff(i, i)+shift) > 1.0e-5) {
                tvecs(i, icol) = vecs(i, icol) / (sparseA_Ga3As3H12.coeff(i, i)+shift);
            }
        }
    }
}

Eigen::SparseMatrix<double> tridiagA_Ga3As3H12;
void init_tridiagA_Ga3As3H12(double shift){
    if(sparseA_Ga3As3H12.rows() == 0) init_SparseA_Ga3As3H12();

    int n = sparseA_Ga3As3H12.rows();
    tridiagA_Ga3As3H12.resize(n, n);
    tridiagA_Ga3As3H12.reserve(3*n); // reserve space for tridiagonal elements

    for(int i = 0; i < n; i++){
        tridiagA_Ga3As3H12.insert(i, i) = sparseA_Ga3As3H12.coeff(i, i)+shift;
        if(i > 0)
            tridiagA_Ga3As3H12.insert(i, i-1) = sparseA_Ga3As3H12.coeff(i, i-1);
        if(i < n-1)
            tridiagA_Ga3As3H12.insert(i, i+1) = sparseA_Ga3As3H12.coeff(i, i+1);
    }

    tridiagA_Ga3As3H12.makeCompressed();
}

void tridiagA_precnd_Ga3As3H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(tridiagA_Ga3As3H12.rows() == 0) init_tridiagA_Ga3As3H12(shift);
    for(int i = 0; i < n; i++){
        tridiagA_Ga3As3H12.coeffRef(i, i) = tridiagA_Ga3As3H12.coeff(i, i)+shift;
    }
    // solve tridiagA_Ga3As3H12 * tvecs = vecs // 计算LU分解
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(tridiagA_Ga3As3H12);
    solver.factorize(tridiagA_Ga3As3H12);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Decomposition failed when doing tridiagonal precnd" << std::endl;}
    // 使用LU分解来解线性系统
    tvecs = solver.solve(vecs);
    if (solver.info() != Eigen::Success) {std::cerr << "SparseLU Solving failed when doing tridiagonal precnd" << std::endl;}
    // 输出解矩阵
    // std::cout << "The solution matrix tvecs is:\n" << tvecs << std::endl;
}