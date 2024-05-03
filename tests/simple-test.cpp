#include "lobpcg.h"
#include <iostream>

void _strange_avec_change_size(int n, int m, const Eigen::MatrixXd vec, Eigen::MatrixXd& after_vec){
    after_vec = vec;
}

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

    // Eigen::VectorXd activeMask = Eigen::VectorXd::Zero(4);
    // activeMask(2) = 1;
    // std::cout << activeMask.transpose() << std::endl;
    // int n_eigenpairs = 3;
    // std::cout << (activeMask.head(n_eigenpairs).array() == 0).all() << std::endl;

    // test for MatrixXd arguments
    Eigen::MatrixXd vecs = Eigen::MatrixXd::Constant(4,2,3); // vec(4,2)
    Eigen::MatrixXd after_vecs = Eigen::MatrixXd::Constant(4,3,5); // vec(4,3)
    std::cout << "after_vecs: \n" << after_vecs << std::endl;
    _strange_avec_change_size(4,2,vecs,after_vecs);
    std::cout << "after_vecs: \n" << after_vecs << std::endl;
}