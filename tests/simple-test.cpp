#include "lobpcg.h"
#include <iostream>

int main(){
    // Eigen::MatrixXd vecs = Eigen::MatrixXd::Constant(4,2,3);
    // Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(4,4);
    // mat.diagonal() << 2,2,2,2;
    // std::cout << mat << std::endl;
    // std::cout << "*" << std::endl;
    // std::cout << vecs << std::endl;
    // vecs = mat * vecs;
    // std::cout << "=" << std::endl;
    // std::cout << vecs << std::endl;
    // vecs.col(0) -= Eigen::VectorXd::Constant(4,1);
    // std::cout << vecs << std::endl;

    Eigen::VectorXd activeMask = Eigen::VectorXd::Zero(4);
    activeMask(2) = 1;
    std::cout << activeMask.transpose() << std::endl;
    int n_eigenpairs = 3;
    std::cout << (activeMask.head(n_eigenpairs).array() == 0).all() << std::endl;
}