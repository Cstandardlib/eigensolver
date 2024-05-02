#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <chrono>

// intermediate tool functions

void check_init_guess(int n, int m, Eigen::MatrixXd& evec);

std::chrono::time_point<std::chrono::high_resolution_clock> get_current_time();

int selfadjoint_eigensolver(Eigen::MatrixXd &A_to_be_eigenvecs, Eigen::VectorXd &eig, int n);

#endif // UTILS_H