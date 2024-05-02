#ifndef UTILS_H
#define UTILS_H
/* includes:
    * get_current_time          -- get current time by <chrono>
    * check_init_guess          -- make init guess and ortho
    * selfadjoint_eigensolver   -- wrapper for selfadjoint_eigensolver
*/

#include <Eigen/Dense>
#include <chrono>

// intermediate tool functions

void check_init_guess(int n, int m, Eigen::MatrixXd& evec);

std::chrono::time_point<std::chrono::high_resolution_clock> get_current_time();

int selfadjoint_eigensolver(Eigen::MatrixXd &A_to_be_eigenvecs, Eigen::VectorXd &eig, int n);

#endif // UTILS_H