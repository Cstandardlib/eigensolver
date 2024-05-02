#include "matvec.h"
#include "test.h"
#include <iostream>
#include <numeric>  // iota
#include <fstream>

#include <Eigen/Sparse>
#include <fast_matrix_market/app/Eigen.hpp>

// test avec function, A*vecs(n,m) -> avecs(n,m)
// incremental diagonal elements from 1 to n
void avec_diag_n(int n, int m, Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    Eigen::VectorXd vec(n);
    std::iota(vec.data(), vec.data() + n, 1); // 从1开始填充

    // 使用向量vec创建对角矩阵
    Eigen::MatrixXd a = vec.asDiagonal();
    for (int icol = 0; icol < m; ++icol) {
        avecs.col(icol) = a * vecs.col(icol);
    }
}

// a1.mtx 9*9 diag
void a1vec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
#ifdef DEBUG_LOBPCG
    std::cout << "-- starts a1vec --" << std::endl;
    std::cout << "size of vecs(n, m) = " << n << ", " << m << std::endl;
#endif
    std::ifstream f("../../../a1.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    Eigen::SparseMatrix<double> mat;
    fast_matrix_market::read_matrix_market_eigen(f, mat); // std::cout << mat << std::endl;
    if(mat.rows() != n || mat.cols() != n){
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;
    }
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "avecs must be of size (n,m)"; return;
    }
    // std::cout << "A-vec a1 = \n" << mat << std::endl;
    // std::cout << "A-vec vecs = \n"<< vecs << std::endl;
    avecs = mat * vecs;
    // std::cout << "A-vec avecs = \n"<< avecs << std::endl;
#ifdef DEBUG_LOBPCG
    std::cout << "-- end a1vec --" << std::endl;
#endif
}

void sparseAvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    // std::ifstream f("../../../../matrix/sparseA.mtx");
    std::ifstream f("../../../sparseA.mtx");
    if(!f.is_open()){
        std::cerr << "failed to open file" << std::endl;
        return;
    }

    Eigen::SparseMatrix<double> mat;

    fast_matrix_market::read_matrix_market_eigen(f, mat);
    // std::cout << mat << std::endl;
    if(mat.rows() != n || mat.cols() != n){
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;
    }
    // check vecs and avecs are of size(n,m)
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "bvecs must be of size (n,m)"; return;
    }

    avecs = mat * vecs;
}



void avec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    // check vecs and avecs are of size(n,m)
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "bvecs must be of size (n,m)"; return;
    }
    // Assume 'a' is a global variable or class member matrix that is already defined and initialized.
    static Eigen::MatrixXd a(n, n);
    for (int i = 0; i < n; ++i) {
        a(i, i) = static_cast<double>(i + 1) + 1.0;
        for (int j = 0; j < i; ++j) {
            a(j,i) = 1.0 / static_cast<double>(i + j);
            a(i,j) = a(j,i);
        }
    }

    // Perform matrix-vector multiplication for each column of x
    for (int icol = 0; icol < m; ++icol) {
        avecs.col(icol) = a * vecs.col(icol);
    }
}


void bvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs){
    // check vecs and bvecs are of size(n,m)
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(bvecs.rows() != n || bvecs.cols() != m){
        std::cerr << "bvecs must be of size (n,m)"; return;
    }

    bvecs = vecs;
}

// identity preconditioner
void precnd(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs){
    // check vecs and tvecs are of size(n,m)
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(tvecs.rows() != n || tvecs.cols() != m){
        std::cerr << "tvecs must be of size (n,m)"; return;
    }
    tvecs = vecs;
}