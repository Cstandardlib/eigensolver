#include "matvec.h"
#include "test.h"
#include <iostream>
#include <numeric>  // iota
#include <fstream>

#include <Eigen/Sparse>
#include <fast_matrix_market/app/Eigen.hpp>

// namespace dense_a{}
static Eigen::MatrixXd a;
void initializeMatrix(int n) {
    std::cout << "initializing a..." << std::endl;
    a.resize(n, n);
    for (int i = 0; i < n; ++i) {
        a(i, i) = static_cast<double>(i + 1) + 1.0;
        for (int j = 0; j < i; ++j) {
            a(j,i) = 1.0 / static_cast<double>(i+1 + j+1);
            a(i,j) = a(j,i);
        }
    }
    std::cout << "a size = " << a.rows() << ", " << a.cols() << std::endl;
#ifdef DEBUG_MATVEC
    std::cout << "a initialized = \n" << a << std::endl;
#endif
}


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
    avecs = mat * vecs;
}
// b1.mtx 9*9 diag
void b1vec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs){
#ifdef DEBUG_LOBPCG
    std::cout << "-- starts b1vec --" << std::endl;
    std::cout << "size of vecs(n, m) = " << n << ", " << m << std::endl;
#endif
    std::ifstream f("../../../b1.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    Eigen::SparseMatrix<double> mat;
    fast_matrix_market::read_matrix_market_eigen(f, mat); // std::cout << mat << std::endl;
    if(mat.rows() != n || mat.cols() != n){
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;
    }
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(bvecs.rows() != n || bvecs.cols() != m){
        std::cerr << "bvecs must be of size (n,m)"; return;
    }
    bvecs = mat * vecs;
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
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<"), but "
                    << "size ("<<mat.rows()<<", "<<mat.cols()<<") is provided" << std::endl;
        return;
    }
    // check vecs and avecs are of size(n,m)
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m) = "<<"("<<n<<", "<<m<<"), but "
                    << "size ("<<vecs.rows()<<", "<<vecs.cols()<<") is provided" << std::endl;
        return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "avecs must be of size (n,m) = "<<"("<<n<<", "<<m<<"), but "
                    << "size ("<<avecs.rows()<<", "<<avecs.cols()<<") is provided" << std::endl;
        return;
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

    if (a.rows() == 0) {
        initializeMatrix(n);
    }

    // Perform matrix-vector multiplication for each column of x
    for (int icol = 0; icol < m; ++icol) {
        avecs.col(icol) = a * vecs.col(icol);
    }
}

// b = diag(2,2, ...)
void bvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs){
    // check vecs and bvecs are of size(n,m)
    std::cout << "bvec: n, m = " << n << ", " << m << std::endl;
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(bvecs.rows() != n || bvecs.cols() != m){
        std::cerr << "bvecs must be of size (n,m)"; return;
    }
    bvecs = 2*vecs;
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



void mprec(int n, int m, const Eigen::MatrixXd& x, Eigen::MatrixXd& px) {
    // double fac, 
    // 检查输入矩阵x的维度是否正确
    assert(x.rows() == n && x.cols() == m);
    
    // 检查输出矩阵px是否已经正确初始化
    assert(px.rows() == n && px.cols() == m);

    if (a.rows() == 0) {
        initializeMatrix(n);
    }
    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            // 检查分母是否为0，避免除以0的错误
            // if (abs(a(i, i) + fac) > 1.0e-5) {
            if (abs(a(i, i)) > 1.0e-5) {
                px(i, icol) = x(i, icol) / (a(i, i));
            }
            // else {
                // 如果分母接近0，可以设置一个错误值或者抛出异常
                // 这里简单地设置为0，实际情况可能需要更复杂的错误处理
                // px(i, icol) = 0.0;
                // do nothing
            // }
        }
    }
}


void avec_Si5H12(int n, int m, const Eigen::MatrixXd &vecs, Eigen::MatrixXd &avecs)
{
}