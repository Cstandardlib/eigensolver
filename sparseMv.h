#ifndef SPARSEMV_H
#define SPARSEMV_H

// Include any necessary headers here
#include <Eigen/Dense>  // MatrixXd
#include <Eigen/Sparse> // SparseMatrix

// Define your class or function declarations here
void sparse_avec_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);
void sparse_precnd_Si2(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift);

void sparse_avec_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);
void sparse_precnd_Na5(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift);

void sparse_avec_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);
void sparse_precnd_Si5H12(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs, double shift);

#endif // SPARSEMV_H