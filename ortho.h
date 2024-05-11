#ifndef ORTHO_H
#define ORTHO_H

// Add any necessary includes here
#include <Eigen/Dense>
#include <iostream>

// # define DEBUG_ORTHO

// basic tool functions -- linear algebra routines
// ortho
/* contains:
 * ortho(QR)
 * ortho(Cholesky)
 * b_ortho
 * ortho_against_y
 * b_ortho_against_y
 */

// deprecated, QR
// #define USE_QR

#define USE_THIN_QR

// else use cholesky

void ortho(int n, int m, Eigen::MatrixXd& u);
void ortho_qr(int n, int m, Eigen::MatrixXd &u);
// void ortho_cho(int n, int m, Eigen::MatrixXd &x);

void b_ortho(int n, int m, Eigen::MatrixXd& u, Eigen::MatrixXd& bu);

void ortho_against_y(int n, int m, int k, Eigen::MatrixXd& x, const Eigen::MatrixXd& y);

void b_ortho_against_y(int n, int m, int k, Eigen::MatrixXd& x, const Eigen::MatrixXd& y, const Eigen::MatrixXd& by);

#endif // ORTHO_H