/*
test for sparse matrix A read from .mtx
to perform matrix-vector multiplication
on dense block vectors v
that is A*v
where A is sparse
and v is a dense block vector
*/

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>

#include <fast_matrix_market/app/Eigen.hpp>
#include <fast_matrix_market/fast_matrix_market.hpp>

void test_read_mtx(){
    std::ifstream f("../../../../../matrix/Si5H12/Si5H12.mtx");
    // std::ifstream f("../../../../input.mtx");
    std::cout << "read." << std::endl;
    

    Eigen::SparseMatrix<double> mat;
    std::cout << "here." << std::endl;
    fast_matrix_market::read_matrix_market_eigen(f, mat);
    std::cout << "here." << std::endl;
    std::cout << mat << std::endl;

    // std::ifstream input_stream;//
    // input_stream.open("input.txt", std::ios::in);

    // if (!input_stream.is_open()) {
    //     std::cerr << "Error opening file." << std::endl;
    //     return;
    // }

    // // 创建Eigen的稀疏矩阵对象
    // Eigen::SparseMatrix<double> mat;

    // // 使用fast_matrix_market读取文件到Eigen的稀疏矩阵
    // fast_matrix_market::read_matrix_market_eigen(input_stream, mat);

    // // 此时mat包含了Matrix Market文件中的数据

    // // 打印矩阵信息
    // std::cout << "Matrix is " << mat.rows() << " by " << mat.cols()
    //           << " with " << mat.nonZeros() << " non-zero entries" << std::endl;
}

template <typename Derived>
void avec(int n, int m, Eigen::MatrixBase<Derived>& vecs, Eigen::MatrixBase<Derived>& avecs){
// void avec(int n, int m, Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
    std::ifstream f("../../../../input.mtx");    

    Eigen::SparseMatrix<double> mat;

    fast_matrix_market::read_matrix_market_eigen(f, mat);

    std::cout << mat << std::endl;

    if(mat.rows() != n || mat.cols() != n){
        std::cerr << "input sparse mat must be of size (n,n)"; return;
    }

    // check vecs and avecs are of size(n,m)
    if(vecs.rows() != n || vecs.cols() != m){
        std::cerr << "vecs must be of size (n,m)"; return;
    }
    if(avecs.rows() != n || avecs.cols() != m){
        std::cerr << "bvecs must be of size (n,m)"; return;
    }

    avecs = mat * vecs;
    // Assume 'a' is a global variable or class member matrix that is already defined and initialized.
    // static Eigen::MatrixXd a(n, n);
    // for (int i = 0; i < n; ++i) {
    //     a(i, i) = static_cast<double>(i + 1) + 1.0;
    //     for (int j = 0; j < i; ++j) {
    //         a(j,i) = 1.0 / static_cast<double>(i + j);
    //         a(i,j) = a(j,i);
    //     }
    // }

    // // Perform matrix-vector multiplication for each column of x
    // for (int icol = 0; icol < m; ++icol) {
    //     avecs.col(icol) = a * vecs.col(icol);
    // }
}

void sparse_demo()
{
    int rows = 3;
    int cols = 3;
    // 创建并初始化一个稀疏矩阵
    Eigen::SparseMatrix<double> sparseMat(rows, cols);
    // 插入非零元素
    sparseMat.insert(0, 1) = 1.0; // 位于(0,1)位置的元素
    sparseMat.insert(1, 2) = 2.0; // 位于(1,2)位置的元素
    sparseMat.insert(2, 0) = 3.0; // 位于(2,0)位置的元素
    sparseMat.makeCompressed();    // 应用插入，优化存储
    // 输出稀疏矩阵
    std::cout << "Sparse Matrix:\n" << sparseMat << std::endl;
    Eigen::MatrixXd denseMat(cols, cols);
    denseMat.setIdentity();
    denseMat = sparseMat * denseMat;
    // std::cout << "sparse to Dense Matrix:\n" << denseMat << std::endl;
    // denseMat(3, 3);
    denseMat << 1, 0, 4,
                0, 5, 0,
                7, 0, 6;

    // 输出密集矩阵
    std::cout << "Dense Matrix:\n" << denseMat << std::endl;
    // 进行矩阵乘法
    Eigen::MatrixXd result = sparseMat * denseMat;
    // 输出结果矩阵
    std::cout << "Result of multiplication:\n" << result << std::endl;
}

void test_sparsemv(){
    int n = 3;
    int m = 2;
    Eigen::MatrixXd vecs = Eigen::MatrixXd::Ones(n, m);
    std::cout << "vecs = \n" << vecs << std::endl;
    Eigen::MatrixXd avecs = vecs;
    std::cout << "avecs = \n" << avecs << std::endl;
    avec(n, m, vecs, avecs);
    std::cout <<"doing avecs" << std::endl;
    std::cout << "avecs = \n" << avecs << std::endl;
}

int main(){
    Eigen::SparseMatrix<double> sparseMat(3, 3);
    sparseMat.insert(0, 0) = 1.0;
    sparseMat.insert(0, 1) = 1.0; // 位于(0,1)位置的元素
    sparseMat.insert(1, 2) = 2.0; // 位于(1,2)位置的元素
    sparseMat.insert(2, 0) = 3.0; // 位于(2,0)位置的元素
    sparseMat.insert(2, 2) = -1.0;
    std::cout << "sparseMat = \n"<< sparseMat << std::endl;
    Eigen::VectorXd diagvec= sparseMat.diagonal(); 
    std::cout << "at (2,1) lies " << sparseMat.coeff(2, 1) << std::endl;
    std::cout << "at (1,2) lies " << sparseMat.coeff(1, 2) << std::endl;
    for(int i=0; i<3; ++i) std::cout << "diag("<<i<<") is " << sparseMat.coeff(i,i) << std::endl;
    std::cout << diagvec << std::endl;
}