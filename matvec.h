#ifndef MATVEC_H
#define MATVEC_H

#include <Eigen/Dense>


void sparseAvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);
// a1.mtx 9*9 diag
void a1vec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);

// template <typename Derived>
void avec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs);

// template <typename Derived>
void bvec(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs);

// template <typename Derived>
void precnd(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs);

// preconditioner for avec
void mprec(int n, int m, const Eigen::MatrixXd& x, Eigen::MatrixXd& px);

// template <typename Derived>
// void avec(int n, int m, Eigen::DenseBase<Derived>& vecs, Eigen::DenseBase<Derived>& avecs){
// // void avec(int n, int m, Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs){
//     // check vecs and avecs are of size(n,m)
//     if(vecs.rows() != n || vecs.cols() != m){
//         std::cerr << "vecs must be of size (n,m)"; return;
//     }
//     if(avecs.rows() != n || avecs.cols() != m){
//         std::cerr << "bvecs must be of size (n,m)"; return;
//     }
//     // Assume 'a' is a global variable or class member matrix that is already defined and initialized.
//     static Eigen::MatrixXd a(n, n);
//     for (int i = 0; i < n; ++i) {
//         a(i, i) = static_cast<double>(i + 1) + 1.0;
//         for (int j = 0; j < i; ++j) {
//             a(j,i) = 1.0 / static_cast<double>(i + j);
//             a(i,j) = a(j,i);
//         }
//     }

//     // Perform matrix-vector multiplication for each column of x
//     for (int icol = 0; icol < m; ++icol) {
//         avecs.col(icol) = a * vecs.col(icol);
//     }
// }

// template <typename Derived>
// void bvec(int n, int m, Eigen::DenseBase<Derived>& vecs, Eigen::DenseBase<Derived>& bvecs){
// // void bvec(int n, int m, Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs){
//     // check vecs and bvecs are of size(n,m)
//     if(vecs.rows() != n || vecs.cols() != m){
//         std::cerr << "vecs must be of size (n,m)"; return;
//     }
//     if(bvecs.rows() != n || bvecs.cols() != m){
//         std::cerr << "bvecs must be of size (n,m)"; return;
//     }

//     bvecs = vecs;
// }

// // identity preconditioner
// template <typename Derived>
// void precnd(int n, int m, Eigen::DenseBase<Derived>& vecs, Eigen::DenseBase<Derived>& tvecs){
// // void precnd(int n, int m, Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs){
//     // check vecs and tvecs are of size(n,m)
//     if(vecs.rows() != n || vecs.cols() != m){
//         std::cerr << "vecs must be of size (n,m)"; return;
//     }
//     if(tvecs.rows() != n || tvecs.cols() != m){
//         std::cerr << "tvecs must be of size (n,m)"; return;
//     }
//     tvecs = vecs;
// }



#endif // MATVEC_H