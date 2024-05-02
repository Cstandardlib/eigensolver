#include <iostream>
#include <Eigen/Eigen>

int main() {
    // 创建一个3x3的矩阵
    Eigen::MatrixXd matrix(3, 3);
    matrix << 4, 2, 0,
              2, 5, 1,
              0, 1, 6;

    // 计算矩阵的特征值
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matrix);
    Eigen::VectorXcd eigenvalues = eigensolver.eigenvalues();
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

    std::cout << "The eigenvalues are:" << std::endl << eigenvalues << std::endl;
    std::cout << "The eigenvectors are:" << std::endl << eigenvectors << std::endl;

    return 0;
}
