#include <Eigen/Dense>
#include <iostream>
#include "matvec.h"

// 函数原型，修改矩阵 A 中从 (iRow, iCol) 开始的 size x size 大小的块
void modifyBlock(Eigen::MatrixXd& A, int iRow, int iCol, int size, const Eigen::MatrixXd& block) {
    // 使用 Eigen::Block 来获取 A 中相应的块并修改它
    A.block(iRow, iCol, size, size) = block;
}

void set_block(Eigen::Block<Eigen::MatrixXd> A, int num){
    A.setConstant(num);
}

void mprec(int n, int m,  const Eigen::MatrixXd &x, Eigen::MatrixXd &px) {
    double fac = 2;
    // 假设 'a' 是一个已经定义并初始化的矩阵，其类型与Eigen::MatrixXd兼容
    Eigen::MatrixXd a(n, n);
    for (int i = 0; i < n; ++i) {
        a(i, i) = static_cast<double>(i + 1) + 1.0;
        for (int j = 0; j < i; ++j) {
            a(j,i) = 1.0 / static_cast<double>(i + j);
            a(i,j) = a(j,i);
        }
    }

    // 输出生成的对称矩阵a
    std::cout << "Symmetric matrix a:\n" << a << std::endl;
    // Eigen::MatrixXd xx(n, m);
    // xx.setConstant(3);
    // std::cout << "x = \n"<< xx << std::endl;
    // std::cout << "a * x = \n"<< a * xx << std::endl;

    // 检查 'a' 矩阵是否已经正确定义和初始化
    if (a.rows() != n || a.cols() != n) {
        std::cerr << "a size != n" << std::endl;
        return;
    }

    // Eigen库中的矩阵是列主序存储的，这与Fortran的行主序不同
    // 因此，我们需要确保按照正确的顺序访问元素
/*
 do icol = 1, m
      do i = 1, n
        if (abs(a(i,i)+fac).gt.tol) px(i,icol) = x(i,icol)/(a(i,i) + fac)
      end do
    end do
    return
*/

    double tol = 1.0e-5;
    for (int icol = 0; icol < m; ++icol) {
        for (int i = 0; i < n; ++i) {
            // 注意Eigen中使用()运算符来访问元素，而不是Fortran中所用的逗号
            if (std::abs(a(i, i) + fac) > tol) {
                px(i, icol) = x(i, icol) / (a(i, i) + fac);
            } else {
                // 如果分母接近零，可以设置px(i, icol)为零或者进行其他处理
                px(i, icol) = 0.0;
            }
        }
    }
}

void v_w_test(){
    int n = 5;
    int n_max_subspace = 3;
    int size_space = 3*n_max_subspace;          // Matrix<typename _Scalar, int _Rows, int _Cols, ...>
    Eigen::MatrixXd v(n, size_space);
    Eigen::MatrixXd r(n, n_max_subspace);   // residuals
    Eigen::MatrixXd tr(n, n_max_subspace);  // preconditioned residuals
    v.setConstant(3);
    r.setOnes();
    tr.setZero();
    std::cout << "v = \n" << v << std::endl;

    int index_w = n_max_subspace;
    /* v(n, n_max_subspace) = t * r(n, n_max_subspace)*/
    mprec(n, n_max_subspace, r, tr);   // tr = t * r, r = ones
    std::cout << "tr = \n" << tr << std::endl;
    v.middleCols(index_w, n_max_subspace) = tr;
    std::cout << "v = \n" << v << std::endl;

    // Eigen::MatrixXd a(n, n);
    // for (int i = 0; i < n; ++i) {
    //     a(i, i) = static_cast<double>(i + 1) + 1.0;
    //     for (int j = 0; j < i; ++j) {
    //         a(j,i) = 1.0 / static_cast<double>(i + j);
    //         a(i,j) = a(j,i);
    //     }
    // }

    double tol = 1.0e-5;
    Eigen::MatrixXd px(n,n_max_subspace);
    Eigen::MatrixXd x = Eigen::MatrixXd::Ones(n,n_max_subspace);
    std::cout << "x = \n" << x << std::endl;
    mprec(n, n_max_subspace, x, px);
    std::cout << "preconditioned px = \n" << px << std::endl;

}

void pass_block_test() {
    int rows = 4;
    int cols = 4;
    int block_size = 2;

    // 创建一个 4x4 的矩阵 A，并用随机数填充
    Eigen::MatrixXd A(rows, cols);
    A = Eigen::MatrixXd::Identity(rows, cols);//Random(rows, cols);

    // 创建一个 2x2 的矩阵 block，并用随机数填充
    Eigen::MatrixXd block(block_size, block_size); // 2x2 block
    block = Eigen::MatrixXd::Random(block_size, block_size);

    // 将 block 赋值到 A 的左上角部分
    modifyBlock(A, 0, 0, block_size, block);

    // 输出修改后的矩阵 A
    std::cout << "Matrix A after modifying a block:\n" << A << std::endl;

    Eigen::Block<Eigen::MatrixXd> uBlock = A.block(2, 0, 2, 2);

    set_block(uBlock, 3);

    std::cout << "Matrix A after modifying a block(2, 0, 2, 2):\n" << A << std::endl;

    // Eigen::Block<Eigen::MatrixXd> rightCols = A.middleCols(2,2);
    // std::cout << rightCols << std::endl;
    auto rightCols = A.block(0, 2, 4, 2);

    set_block(rightCols, 5);

    std::cout << "Matrix A after modifying a block(0, 2, 4, 2):\n" << A << std::endl;

}

void test_make_reduce(){
    int n=5;
    int n_max_subspace = 2;
    int size_space = 3*n_max_subspace;
    Eigen::MatrixXd v(n, size_space); v.leftCols(n_max_subspace).setConstant(3);
    
    std::cout << "v = \n" << v << std::endl;
    Eigen::MatrixXd A(n,n);
    A.diagonal() << 1,2,3,4,5;

    Eigen::MatrixXd av(n, size_space);
    av = A*v;
    std::cout << "av = \n" << av << std::endl;

    Eigen::MatrixXd A_reduced(size_space, size_space); A_reduced.setZero();
    A_reduced = v.transpose() * av;
    std::cout << "A_reduced =\n" << A_reduced << std::endl;
    std::cout << "v'Av = \n" << v.transpose()*A*v<<std::endl;
    Eigen::MatrixXd v2=v.leftCols(n_max_subspace);
    std::cout << "v2 = \n" << v2 << std::endl;
    Eigen::MatrixXd v3(n,n_max_subspace);
    v3 << 3,3,
          6,6,
          9,9,
          12,12,
          15,15;
    std::cout << "v3 = \n" << v3 << std::endl;
    std::cout << "ans = \n" << v2.transpose()*v3 << std::endl;
}

void corner_test(){
    // v_w_test();
    int n=5;
    int m=3;
    int size_space = 3*m;
    Eigen::MatrixXd v(n,size_space);
    v.setConstant(3);
    Eigen::MatrixXd A(size_space,size_space);
    A.setZero();
    A.diagonal()<< 1,2,3,4,5,6,7,8,9;
    std::cout << "A = \n" << A << std::endl;

    v.leftCols(m) = v.leftCols(m) * A.topLeftCorner(m,m);
    std::cout << "v = \n" << v << std::endl;
}

void change_dynamic_test(){
    // MatrixXd can be changed at runtime
    Eigen::MatrixXd mat(3, 4); mat.setOnes(); std::cout << mat << std::endl;
    std::cout << "size = " <<  mat.size() << std::endl;
    Eigen::MatrixXd anotherMat(2, 2); anotherMat.setConstant(2); std::cout << anotherMat << std::endl;
    mat = anotherMat; // 这将产生编译错误，因为尺寸不匹配
    std::cout << mat << std::endl;
    std::cout << "size = " <<  mat.size() << std::endl;
}

template <typename Derived>
void print_block(const Eigen::DenseBase<Derived>& b, int x, int y, int r, int c)
{
  std::cout << "block: \n" << b.block(x,y,r,c) << std::endl;
}

void base_arguments_test(){
    int n = 6;
    Eigen::MatrixXd mat(n,n);
    for(int i=0; i<mat.size();++i){
        mat(i) = i;
    }
    std::cout <<"mat = \n" << mat << std::endl;
    std::cout <<"mat(1, 1, 2, 2) = \n" << std::endl;
    print_block(mat, 1, 1, 2, 2);
    std::cout <<"mat(2, 2, 2, 2) = \n" << std::endl;
    print_block(mat.block(1,1,3,3), 1, 1, 2, 2);

}

template <typename Derived>
void local_avec(int n, int m, Eigen::DenseBase<Derived>& vecs, Eigen::DenseBase<Derived>& avecs){
// void avec(int n, int m, Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
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

void base_argument_avec_test(){
    int n = 3;
    int m = 2;
    Eigen::MatrixXd vecs = Eigen::MatrixXd::Ones(n, m);
    std::cout << "vecs = \n" << vecs << std::endl;
    Eigen::MatrixXd avecs = vecs;
    std::cout << "avecs = \n" << avecs << std::endl;
    local_avec(n, m, vecs, avecs);
    std::cout <<"doing avecs" << std::endl;
    std::cout << "avecs = \n" << avecs << std::endl;
}

void xwp_test(){
    int n=10;
    int n_max_subspace = 3;
    int n_active = 2;
    int size_space = 3*n_max_subspace;
    Eigen::MatrixXd v,x,w,p;
    v = Eigen::MatrixXd::Constant(n, size_space, 0);
    x = Eigen::MatrixXd::Constant(n, n_max_subspace, 1);
    w = Eigen::MatrixXd::Constant(n, n_active, 2);
    p = Eigen::MatrixXd::Constant(n, n_active, 3);
    std::cout << "v = \n" << v << std::endl;
    std::cout << "x = \n" << x << std::endl;
    std::cout << "w = \n" << w << std::endl;
    std::cout << "p = \n" << p << std::endl;
    v.leftCols(n_max_subspace) = x;
    v.middleCols(n_max_subspace, n_active) = w;
    v.middleCols(n_max_subspace + n_active, n_active) = p;
    std::cout << "v = \n" << v << std::endl;
}
// test for substract a block by I
void coeff_test(){
    int n_max_subspace = 4;
    int n_active = 2;
    int size_space = 3*n_max_subspace;
    int n_working_space = n_max_subspace + 2*n_active;
    Eigen::MatrixXd A_reduced = Eigen::MatrixXd::Constant(size_space, size_space, 2);
    Eigen::MatrixXd coeff = A_reduced.topLeftCorner(size_space, n_max_subspace);
    std::cout << "coeff = \n" << coeff << std::endl;
    auto c_p = coeff.block(n_max_subspace-n_active, n_max_subspace-n_active,n_active,n_active);
    c_p -= Eigen::MatrixXd::Identity(n_active, n_active);
    std::cout << "coeff = \n" << coeff << std::endl;
}
// 拼接 x and p
void xp_from_x_and_p_test(){
     // 假设 x 是一个 m x n 的矩阵，p 是一个 m x k 的矩阵
    int m = 3; // 行数
    int n = 2; // x 的列数
    int k = 1; // p 的列数

    // 创建矩阵 x 和 p
    Eigen::MatrixXd x(m, n);
    Eigen::MatrixXd p(m, k);

    x << 1, 2,
         3, 4,
         5, 6;
    p << 7,
         8,
         9;

    assert(x.rows() == p.rows());
    std::cout << "x = \n" << x << std::endl;
    std::cout << "p = \n" << p << std::endl;

    // 创建 xp 矩阵，其左半部分是 x，右半部分是 p
    Eigen::MatrixXd xp(x.rows(), x.cols() + p.cols());
    xp << x, p;

    std::cout << "Matrix xp:\n" << xp << std::endl;
}

int main(){
    // test_make_reduce();
    // base_arguments_test();
    // base_argument_avec_test();
    // xwp_test();
    // coeff_test();
    xp_from_x_and_p_test();
    // Eigen::VectorXd ddd;
    // std::cout << ddd.size() << std::endl; // will be zero
}