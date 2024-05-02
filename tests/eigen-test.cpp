#include <Eigen/Dense>
#include <iostream>
#include "utils.h"

void test()
{
    // 创建一个4x4的矩阵
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(4, 4);
    mat.diagonal() << 8,7,6,5;

    // 使用EigenSolver来计算特征值和特征向量
    // Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(mat);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mat);
    
    // 检查是否成功
    if (eigensolver.info() == Eigen::Success)
    {
        // 获取特征值
        Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
        // 获取特征向量
        Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
        
        // 输出特征值
        std::cout << "The eigenvalues are:\n" << eigenvalues << std::endl;
        // std::cout << "The eigenvalues are:\n" << eigensolver.eigenvalues() << std::endl;
        
        // 输出特征向量
        std::cout << "The eigenvectors are:\n" << eigenvectors << std::endl;
        // std::cout << "The eigenvectors are:\n" << eigensolver.eigenvectors() << std::endl;
        std::cout << "orthonormal check: \n" << eigenvectors.transpose() * eigenvectors << std::endl;
    }
    else
    {
        std::cerr << "Eigen decomposition failed." << std::endl;
    }
}
// Test for selfadjoint_eigensolver
void test_selfadjoint_eigensolver() {
    int n = 4; // number of rows
    Eigen::MatrixXd A(n, n); // initialize A with zeros
    A << 1, 2, 3, 4,
         2, 5, 6, 7,
         3, 6, 8, 9,
         4, 7, 9, 10;

    auto A_copy = A;

    Eigen::VectorXd eig(n);
    eig << 1, 2, 3, 4;
    std::cout << eig << std::endl;
    Eigen::MatrixXd eigd = eig.asDiagonal();
    std::cout << "eig_as_diagonal: \n" << eigd << std::endl;
    // Eigen::Vector4i M; M<<1, 1, 3, 4;
    // Eigen::MatrixXi V(M.asDiagonal());
    // std::cout<<V<<std::endl;

    std::cout << "A: \n" << A << std::endl;

    selfadjoint_eigensolver(A, eig, n);

    std::cout << "A: \n" << A << std::endl;
    std::cout << "eig: \n" << eig << std::endl;
    // Check if eigenvalues and eigenvectors are computed correctly
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A_copy);
    // std::cout << "A * eigenvectors = \n" << A_copy*A << std::endl;
    // std::cout << "eig * eigenvectors = \n" << eig.asDiagonal() * A << std::endl;
    assert(eigensolver.info() == Eigen::Success);
    assert(eig.isApprox(eigensolver.eigenvalues()));
    assert(A.isApprox(eigensolver.eigenvectors()));

    A = A_copy;
    std::cout << "A: \n" << A << std::endl;

    selfadjoint_eigensolver(A, eig, 2);

    std::cout << "A: \n" << A << std::endl;
    std::cout << "eig: \n" << eig << std::endl;
    eigensolver.compute(A_copy.topLeftCorner(2,2));

    assert(eigensolver.info() == Eigen::Success);
    assert(eig.head(2).isApprox(eigensolver.eigenvalues()));
    std::cout << "std eig = \n" << eigensolver.eigenvalues() << std::endl;
    assert(A.topLeftCorner(2,2).isApprox(eigensolver.eigenvectors()));

}
void test_vav(){
    int n = 10;
    int nmax=2;
    int size_space = 3*nmax;
    
    Eigen::MatrixXd evec(n,nmax);

    Eigen::MatrixXd v(n,size_space);
    Eigen::MatrixXd A(n,n);
    A.diagonal() << 1,2,3,4,5,6,7,8,9,10;
    Eigen::MatrixXd vav(size_space, size_space); // 3*nmax, 3*nmax
    int cnt=1;
    for(int i=0; i<n; ++i){
        for(int j=0; j<nmax; ++j){
            evec(i,j) = ++cnt;
        }
    }
    std::cout << evec << std::endl;
    std::cout << "v.size = " <<v.rows() << "," << v.cols() << std::endl;
    v.topLeftCorner(n,nmax) = evec;
    std::cout << "v.size = " << v.rows() << "," << v.cols() << std::endl;
    std::cout << v << std::endl;
    Eigen::MatrixXd av = A * v;
    vav = v.transpose() *av;
    std::cout << "vav = \n" << vav << std::endl;
    Eigen::VectorXd eig(size_space);
    selfadjoint_eigensolver(vav, eig, nmax);
    /* eigenvalues and eigenvectors should be
    eigenvalue
Out[7]: array([2.81480609e+00, 2.57921852e+04])

In [8]: eigenvectors
Out[8]: 
array([[-0.72864504, -0.68489153],
       [ 0.68489153, -0.72864504]])
    according to numpy
    */
   std::cout << "eigenvalues are:\n" << eig << std::endl;
   std::cout << "eigenvectors are: \n" << vav << std::endl;
}

#include <chrono>
void test_v_Areduced(){
    // test_selfadjoint_eigensolver();
    // test_vav();
    // test for different size of matrices
    int n = 10;
    int n_max_subspace = 2;
    // Eigen::MatrixXd evec = Eigen::MatrixXd::Random(n,n_max_subspace);
    int size_space = 3*n_max_subspace;          // Matrix<typename _Scalar, int _Rows, int _Cols, ...>
    Eigen::MatrixXd v(n, size_space);
    v.setOnes();
    Eigen::MatrixXd A_reduced(size_space, size_space); // 6 * 6
    A_reduced.setConstant(3.0);
    Eigen::MatrixXd evec = Eigen::MatrixXd::Ones(n,n_max_subspace);
    // selfadjoint_eigensolver(A_reduced, evec, size_space);
    // std::cout << "evec = \n" << evec << std::endl;
    std::cout << "v = \n" << v << std::endl;
    std::cout << "A_reduced = \n" << A_reduced << std::endl;

    // evec = (v * A_reduced).leftCols(n_max_subspace);
    v.leftCols(n_max_subspace) = v.leftCols(n_max_subspace) * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
    // evec = (v * A_reduced);
    // std::cout << "evec = \n" << evec << std::endl;
    std::cout << "v = \n" << v << std::endl;
    // std::cout << "A_reduced = \n" << A_reduced << std::endl;
    
    // std::cout << evec << std::endl;
    // std::chrono::duration<double> t_diag, tt, ts;
    // std::cout << t_diag.count() << std::endl;
    // std::cout << tt.count() << std::endl;
    // std::cout << ts.count() << std::endl;
}

int main(){
    test();
    // test_selfadjoint_eigensolver();
    // Eigen::MatrixXd A = Eigen::MatrixXd::Identity(3,3);
    // std::cout << A.topLeftCorner(3,3) << std::endl;
    // std::cout << A.size() << std::endl;
}