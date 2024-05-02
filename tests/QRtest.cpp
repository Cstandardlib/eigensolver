#include "Eigen/Dense"
#include <iostream>
#include <Eigen/QR>
#include "lobpcg.h"

/*
official demo
https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html#a89bc9bbb69d140898c5f5e8de11f8066
thin Q

MatrixXf A(MatrixXf::Random(5,3)), thinQ(MatrixXf::Identity(5,3)), Q;
A.setRandom();
HouseholderQR<MatrixXf> qr(A);
Q = qr.householderQ();
thinQ = qr.householderQ() * thinQ;
std::cout << "The complete unitary matrix Q is:\n" << Q << "\n\n";
std::cout << "The thin matrix Q is:\n" << thinQ << "\n\n";
*/


void test_householderQR(){
    /*
    official demo
    https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html#a89bc9bbb69d140898c5f5e8de11f8066
    thin Q
    */
    Eigen::MatrixXf A(Eigen::MatrixXf::Random(5,3)), thinQ(Eigen::MatrixXf::Identity(5,3)), Q;
    A.setRandom();
    Eigen::HouseholderQR<Eigen::MatrixXf> qr(A);
    Q = qr.householderQ();
    thinQ = qr.householderQ() * thinQ;
    std::cout << "The complete unitary matrix Q is:\n" << Q << "\n\n";
    std::cout << "The thin matrix Q is:\n" << thinQ << "\n\n";
    std::cout << "thinQ' thinQ = \n" << thinQ.transpose() * thinQ << std::endl;
}

void test_AQR_HouseholderQR(){
    // to check whether the QR decomposition is in correct order
    // i.e., whether there is a permutation matrix P
    std::cout << "--- HouseholderQR A = QR ---" << std::endl;
    int n=5;
    int m=3;

    Eigen::MatrixXd A(Eigen::MatrixXd::Random(n,m)), thinQ(Eigen::MatrixXd::Identity(n,m)), Q;
    A.setRandom();
    // A << 11,5,4,3,2,9,7,8,6,10,11,12,13,14,15;
    std::cout << "A = \n" << A << std::endl;
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();
    std::cout << "The matrixQR upper is:\n" << R << "\n\n";
    for(int i=0; i<R.cols(); ++i){
        std::cout << "R column "<< i <<" =\n" << R.col(i) << std::endl;
    }
    Eigen::MatrixXd upperR = R.topLeftCorner(m,m);
    // Eigen::MatrixXd upperR = R.block(0,0,m,m);
    std::cout << "The upper matrix R is:\n" << upperR << "\n\n";
    thinQ = qr.householderQ() * thinQ;
    std::cout << "The complete unitary matrix Q is:\n" << Q << "\n\n";
    std::cout << "The thin matrix Q is:\n" << thinQ << "\n\n";
    // std::cout << "thinQ' thinQ = \n" << thinQ.transpose() * thinQ << std::endl;
    std::cout << "The upper matrix R is:\n" << upperR << "\n\n";
    std::cout << " Q1 * R = \n" << thinQ * upperR << std::endl;
    std::cout << "A = \n" << A << std::endl;
    std::cout << "A - Q1*R =  " << (A - thinQ*upperR) <<std::endl;
    // solve XR=U
    Eigen::MatrixXd Uortho = A * upperR.reverse();
    std::cout << "Uortho = \n" << Uortho << std::endl;
    std::cout << "Uortho * R = \n" << Uortho * upperR << std::endl;
    std::cout << "R-1 * R = " << upperR.reverse() * upperR << std::endl;
}
void test_APQR_ColPivHouseholderQR(){
    // to check whether the QR decomposition is in correct order
    // i.e., whether there is a permutation matrix P
    std::cout << "--- ColPivHouseholderQR AP = QR ---" << std::endl;
    int n=5;
    int m=3;

    Eigen::MatrixXd A(Eigen::MatrixXd::Random(n,m)), thinQ(Eigen::MatrixXd::Identity(n,m)), Q;
    A.setRandom();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    Q = qr.matrixQ();
    Q = Q.topLeftCorner(n,m);
    Eigen::MatrixXd R = qr.matrixR().template triangularView<Eigen::Upper>();
    thinQ = qr.householderQ() * thinQ;
    std::cout << "The complete unitary matrix Q is:\n" << Q << "\n\n";
    std::cout << "The thin matrix Q is:\n" << thinQ << "\n\n";
    std::cout << "thinQ' thinQ = \n" << thinQ.transpose() * thinQ << std::endl;
    std::cout << " Q1 * R = \n" << thinQ * R << std::endl;

}

void test_2qr(){
    // test_householderQR();
    // return 0;
    int n = 4;
    int m = 3;
    Eigen::MatrixXd u = Eigen::MatrixXd::Random(n, m);
    u << 1,2,3,4,5,6,7,8,9,10,11,12;
    auto v = u;
    auto ucopy = u;
    std::cout << "u = \n" << u << std::endl;
    std::cout << "v = \n" << v << std::endl;
    std::cout << std::endl;
    ortho_qr(n,m,u); std::cout << "--- ortho HouseholderQR ---" << std::endl;
    std::cout << "u = \n" << u << std::endl;
    std::cout << "u'u = \n" << u.transpose() * u << std::endl;
    ortho(n, m, v); std::cout << "--- ortho ColPivHouseholderQR ---" << std::endl;
    std::cout << "v = \n" << v << std::endl;
    std::cout << "v'v = \n" << v.transpose() * v << std::endl;
    // solve Qx = U
    std::cout << "----- Solving Qx = U-----" << std::endl;
    std::cout << "U = \n" << ucopy << std::endl;
    std::cout << "----- HouseholderQR -----" << std::endl;
    Eigen::MatrixXd r1 = u.colPivHouseholderQr().solve(ucopy);
    std::cout << "q1 = \n" << u << std::endl;
    std::cout << "r1 = \n" << r1 << std::endl;
    std::cout << "q1 * r1 = \n" << u * r1 << std::endl;
    std::cout << "----- ColPivHouseholderQR -----" << std::endl;
    Eigen::MatrixXd r2 = v.colPivHouseholderQr().solve(ucopy);
    std::cout << "q2 = \n" << v << std::endl;
    std::cout << "r2 = \n" << r2 << std::endl;
    std::cout << "q2 * r2 = \n" << v * r2 << std::endl;
}

int main(){
    test_AQR_HouseholderQR();
    // test_APQR_ColPivHouseholderQR();
}

/* scipy
u = np.array([[1,2,3,4],
[5,6,7,8]])
q, r = qr(u, mode='complete')
normalized_u = q / np.linalg.norm(q, axis=0)
*/

/*
official demo
https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html#a89bc9bbb69d140898c5f5e8de11f8066
thin Q

MatrixXf A(MatrixXf::Random(5,3)), thinQ(MatrixXf::Identity(5,3)), Q;
A.setRandom();
HouseholderQR<MatrixXf> qr(A);
Q = qr.householderQ();
thinQ = qr.householderQ() * thinQ;
std::cout << "The complete unitary matrix Q is:\n" << Q << "\n\n";
std::cout << "The thin matrix Q is:\n" << thinQ << "\n\n";
*/

// /*

// */

int test_basis() {
    // // 创建一个3x3的矩阵，这里仅为示例
    // Eigen::MatrixXd B(3,3);
    // for(int i=0; i<9; ++i){
    //     B(i) = i-4;
    //     // B[i] = i; // bugs, not the proper index
    // } 
    // std::cout << B << std::endl;
    // // 可以看到，Eigen的矩阵存储是列优先的！！

    Eigen::MatrixXd A(3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    // // 进行Householder QR分解
    // Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);

    // // 获取Q矩阵
    // // Eigen::MatrixXd Q = qr.matrixQR();

    // // 获取R矩阵
    // Eigen::MatrixXd R = qr.matrixQR();

    // // 输出Q和R矩阵
    // // std::cout << "Q matrix:\n" << Q << std::endl;
    // std::cout << "R matrix:\n" << R << std::endl;

    // 创建一个动态大小的矩阵A
    // Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 3); // 4x3随机矩阵

    // 对矩阵A进行列主元Householder QR分解
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);

    // 获取Q矩阵
    Eigen::MatrixXd Q = qr.matrixQ();

    // 获取R矩阵
    Eigen::MatrixXd R = qr.matrixR().template triangularView<Eigen::Upper>();

    // 输出Q和R矩阵
    std::cout << "The matrix Q is:\n" << Q << std::endl;
    std::cout << "The matrix R is:\n" << R << std::endl;

    // 验证QR分解
    std::cout << "The product QR is:\n" << Q * R << std::endl;
    std::cout << "The original matrix A is:\n" << A << std::endl;

    Eigen::MatrixXd u = A;

    Eigen::HouseholderQR<Eigen::MatrixXd> qr_solver = u.householderQr();
    Eigen::MatrixXd Q2 = qr_solver.householderQ();
    Eigen::MatrixXd R2 = qr_solver.matrixQR();

    std::cout << "The matrix U is:\n" << u << std::endl;
    std::cout << "The matrix Q2 is:\n" << Q2 << std::endl;
    std::cout << "The matrix R2 is:\n" << R2 << std::endl;
    std::cout << "The matrix produce Q*R is:\n" << Q2*R2 << std::endl;
    
    // Apply the orthogonalization to u
    u = Q2 * u;

    std::cout << "The matrix u is:\n" << u << std::endl;

    std::cout << " ------ " << std::endl;
    u << 2, 3, 1,
         5, 6, 4,
         8, 9, 7;

    std::cout << " A = \n" << A << std::endl;
    std::cout << " U = \n" << u << std::endl;
    // ortho(true, 3, 3, A, u);
    std::cout << " A = \n" << A << std::endl;
    std::cout << " U = \n" << u << std::endl;


    return 0;
}