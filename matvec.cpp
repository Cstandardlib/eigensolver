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
        std::cerr << "avec failed! vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){std::cerr << "avec failed! avecs must be of size (n,m)"; return;}
    if (a.rows() == 0) initializeMatrix(n);
    // Perform matrix-vector multiplication for each column of x
    for (int icol = 0; icol < m; ++icol) avecs.col(icol) = a * vecs.col(icol);
}

// b = diag(2,2, ...)
void bvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs){
    // check vecs and bvecs are of size(n,m)
    // std::cout << "bvec: n, m = " << n << ", " << m << std::endl;
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "bvecs failed! vecs must be of size (n,m)"; return;
    }
    if(bvecs.rows() != n || bvecs.cols() != m){
        std::cerr << "bvecs failed! bvecs must be of size (n,m)"; return;
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



void mprec(int n, int m, const Eigen::MatrixXd& x, Eigen::MatrixXd& px, double shift) {
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
            if (abs(a(i, i)+shift) > 1.0e-5) {
                px(i, icol) = x(i, icol) / (a(i, i)+shift);
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

static Eigen::MatrixXd denseA_Si2;
void initializeMatrixSi2(){
    std::ifstream f("../../../Si2.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building denseA_Si2, ";
    fast_matrix_market::read_matrix_market_eigen_dense(f, denseA_Si2);
    std::cout <<"size = " << denseA_Si2.rows() << ", " << denseA_Si2.cols() << std::endl;
    std::cout << "denseA_Si2 size = " << denseA_Si2.rows() << ", " << denseA_Si2.cols() << std::endl;
}

void avec_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    if(denseA_Si2.rows() == 0){
        initializeMatrixSi2();
    }
    if(denseA_Si2.rows() != n || denseA_Si2.cols() != n){
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;
    }
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "avecs must be of size (n,m)"; return;
    }
    avecs = denseA_Si2 * vecs;
}

void precnd_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(denseA_Si2.rows() == 0){
        initializeMatrixSi2();
    }
    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            // 检查分母是否为0，避免除以0的错误
            // if (abs(a(i, i) + fac) > 1.0e-5) {
            if (abs(denseA_Si2(i, i)+shift) > 1.0e-5) {
                tvecs(i, icol) = vecs(i, icol) / (denseA_Si2(i, i)+shift);
            }
        }
    }
}

// build Na5 avec and precnd
static Eigen::MatrixXd denseA_Na5;
void initializeMatrixNa5(){
    std::ifstream f("../../../Na5.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building denseA_Na5, ";
    fast_matrix_market::read_matrix_market_eigen_dense(f, denseA_Na5);
    std::cout <<"size = " << denseA_Na5.rows() << ", " << denseA_Na5.cols() << std::endl;
    std::cout << "denseA_Si2 size = " << denseA_Na5.rows() << ", " << denseA_Na5.cols() << std::endl;
}

void avec_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    if(denseA_Na5.rows() == 0){
        initializeMatrixNa5();
    }
    if(denseA_Na5.rows() != n || denseA_Na5.cols() != n){
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;
    }
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "avecs must be of size (n,m)"; return;
    }
    avecs = denseA_Na5 * vecs;
}

void precnd_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift){
    if(denseA_Na5.rows() == 0){
        initializeMatrixNa5();
    }
    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            // 检查分母是否为0，避免除以0的错误
            // if (abs(a(i, i) + fac) > 1.0e-5) {
            if (abs(denseA_Na5(i, i)+shift) > 1.0e-5) {
                tvecs(i, icol) = vecs(i, icol) / (denseA_Na5(i, i)+shift);
            }
        }
    }
}






// --------------------------------------------
static Eigen::SparseMatrix<double> sparseA;
void initializeSparse(){
    std::ifstream f("../../../Si5H12.mtx");
    if(!f.is_open()) {std::cerr << "failed to open file" << std::endl; return;}
    std::cout <<"building spareseA, ";
    fast_matrix_market::read_matrix_market_eigen(f, sparseA);
    std::cout <<"size = " << sparseA.rows() << ", " << sparseA.cols() << std::endl;
    std::cout << "sparseA size = " << sparseA.rows() << ", " << sparseA.cols() << std::endl;
}

void avec_Si5H12(int n, int m, const Eigen::MatrixXd &vecs, Eigen::MatrixXd &avecs)
{
    if(sparseA.rows() == 0){
        initializeSparse();
    }
    if(sparseA.rows() != n || sparseA.cols() != n){
        std::cerr << "input sparse matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"; return;
    }
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "avecs must be of size (n,m)"; return;
    }
    avecs = sparseA * vecs;
}

// void precnd_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs){
//     if(sparseA.rows() == 0){
//         initializeSparse();
//     }
//     for (int icol = 0; icol < m; ++icol) {
//         for (int i = 0; i < n; ++i) {
//             // 检查分母是否为0，避免除以0的错误
//             // if (abs(a(i, i) + fac) > 1.0e-5) {
//             if (abs(sparseA.coeffDiagonal(i)) > 1.0e-5) {
//                 tvecs(i, icol) = vecs(i, icol) / (sparseA(i, i));
//             }
//         }
//     }
// }