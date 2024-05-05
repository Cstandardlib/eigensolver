
#include "ortho.h"

#include <Eigen/Dense>
#include <iostream>

// basic tool functions -- linear algebra routines
// ortho
/* contains:
 * ortho(QR)
 * ortho(Cholesky)
 * b_ortho
 * ortho_against_y
 * b_ortho_against_y
 */


/**
 * @brief Orthonormalize a set of m vectors (of size n) using QR decomposition.
 *
 * First compute the QR decomposition of the input matrix
 * U = Q * (R)
 *         (0)
 * to get R,
 * and then solve the upper triangular system U_{ortho}R = U.
 * We get U_{ortho} = UR^{-1} = Q, and Q is an orthogonal matrix.
 *
 * @param n Size of the vectors (number of rows in u and w).
 * @param m Number of vectors to orthogonalize (number of columns in u and w).
 * @param u Input matrix of vectors, size of n x m ,to be orthogonalized (will be overwritten with Q).
 */
void ortho(int n, int m, Eigen::MatrixXd &u)
{
    // direct approach, using QR
    // Thin QR Factorization!
    /*Eigen::ColPivHouseholderQR performs
      a rank-revealing QR decomposition of a matrix A into matrices P, Q and R
      such that AP=QR*/
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(u);
    if (qr.info() == Eigen::Success){
        Eigen::MatrixXd Qfull = qr.matrixQ();
        u.topLeftCorner(n,m) = Qfull.topLeftCorner(n,m);
    }
    else
        std::cerr << "QR decomposition was not successful." << std::endl;
}


void ortho_qr(int n, int m, Eigen::MatrixXd &u)
{
    // direct approach, using QR
    // Thin QR Factorization! U=Q1R1, see Matrix Computations Page 237
    /* Eigen::HouseholderQR performs
        a QR decomposition of a matrix A into matrices Q and R
        such thatA=QR
    */
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(u);
    // Eigen::MatrixXd thinQ(u); // same shape as u
    u.setIdentity();
    u = qr.householderQ() * u;  // now u=Q
}


void ortho_cho(int n, int m, Eigen::MatrixXd &x)
{
    // direct approach, using Cholesky
    // overlap = x'x
    Eigen::MatrixXd overlap = x.transpose() * x;
    // overlap = LL^T
    Eigen::LLT<Eigen::MatrixXd> llt(overlap);
    if (llt.info() == Eigen::Success)
    {
        Eigen::MatrixXd Lt = llt.matrixL().transpose();
        Lt = Lt.inverse();
        x = x * Lt;
    }
    else
        std::cerr << "Cholesky decomposition failed." << std::endl;
}

void ortho_qr_two_phase(int n, int m, Eigen::MatrixXd &x){
    // ortho using QR, two-phase
    // 1. QR decomposition, get R
    /*
     * First compute the QR decomposition of the input matrix
     * U = Q * (R)
     *         (0)
    */
    Eigen::MatrixXd overlap = x.transpose() * x;
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(overlap);
    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();

    // 2. Solve X_ortho * R = X to get X_ortho = XR^{-1} = Q
    x = x * R.inverse();
}

/**
 * @brief Performs b-orthogonalization of the given matrices `u`, given `bu`.
 * 
 * b-orthogonalize m vectors of lenght n using the Cholesky decomposition
 * of the overlap matrix u^tbu
 * 
 * @param n The number of rows in the matrices `u` and `bu`.
 * @param m The number of columns in the matrices `u` and `bu`.
 * @param u The input matrix `u` to be b-orthonormalized.
 * @param bu The input matrix `bu` to be b-orthonormalized.
 * 
 * @note This function modifies the input matrices `u` and `bu` in-place.
 * 
 * works only if u is already orthonormal.
 * overlap matrix u^tbu matrix can be very ill-conditioned
 * 
 */
void b_ortho(int n, int m, Eigen::MatrixXd& u, Eigen::MatrixXd& bu) {
    // 局部变量
    Eigen::MatrixXd ubu(m, m);
    Eigen::VectorXd sigma(m);
    Eigen::MatrixXd u_svd(m, m), vt_svd(m, m), temp(n, m);
    double tol_svd = 1.0e-5;
    bool use_svd = false; //true;

    // 计算重叠矩阵
    ubu = u.transpose() * bu;
    // ubu = (u.transpose() * bu).eval();
    std::cout << "overlap ubu = \n" << ubu << std::endl;

    if (use_svd) {
        // 使用SVD进行b-正交化
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(ubu, Eigen::ComputeThinU | Eigen::ComputeThinV);
        sigma = svd.singularValues();
        u_svd = svd.matrixU();
        vt_svd = svd.matrixV().transpose();

        // 计算sigma的逆平方根
        for (int i = 0; i < sigma.size(); ++i) {
            if (sigma(i) > tol_svd)  sigma(i) = 1 / sqrt(sigma(i));
            else sigma(i) = 0;
        }

        ubu = Eigen::MatrixXd::Zero(m, m);
        for (int i = 0; i < m; ++i) {
            ubu += sigma(i) * vt_svd.col(i) * u_svd.row(i).transpose();
        }

        // 应用变换
        temp = u * ubu;
        u = temp;
        temp = bu * ubu;
        bu = temp;
    } else {
        // 使用Cholesky分解
        // ubu = LL^T
        Eigen::LLT<Eigen::MatrixXd> llt(ubu);
        if (llt.info() == Eigen::Success) {
            // 计算u * L^-T 和 bu * L^-T
            Eigen::MatrixXd Lt = llt.matrixL().transpose();
            // std::cout << "L^T = \n" << Lt << std::endl;
            Lt = Lt.inverse();
            // std::cout << "L^-T = \n" << Lt << std::endl;
            // std::cout << "u = \n" << u << std::endl;
            // std::cout << "bu = \n" << bu << std::endl;
            u = u * Lt;
            bu = bu * Lt;
            // std::cout << "u = \n" << u << std::endl;
            // std::cout << "bu = \n" << bu << std::endl;
        } else {
            std::cerr << "Cholesky decomposition failed." << std::endl;
        }
    }
}


/**
 * Orthogonalizes block vector `x` against a given \b orthonormal set `y`
 * and orthonormalizes `x`. If `y` is not orthonormal, extra computation
 * is needed inside this function.
 * 
 * `y(n, m)` is a given block vector, assumed to be orthonormal
 * `x(n, k)` is the block vector to be orthogonalized
 * first against x and then internally.
 * 
 * if y is not orthonormal, we will do by solving (Y'Y)y_coeff = Y'X
 * X = X - Y(Y'Y)^{-1}Y'X = X - Y * y_coeff
 *
 * @param n The number of rows in `x` and `y`.
 * @param m The number of columns in `y`.
 * @param k The number of columns in `x`.
 * @param x The matrix to be orthogonalized. size of n x k.
 * @param y The matrix to orthogonalize against. size of n x m. assumed to be orthonormal
 */
void ortho_against_y(int n, int m, int k, Eigen::MatrixXd& x, const Eigen::MatrixXd& y){
#ifdef DEBUG_LOBPCG
    // std::cout << "----- start ortho against y -----" << std::endl;
#endif
    // orthogonalize x(n, k) against y(n, m)
    double tol_ortho = 1.0e-10;//1.0e-13;

    // first check input y is orthonormal
    bool is_y_orthonormal = true;
    Eigen::MatrixXd yby = y.transpose() * y;
    Eigen::MatrixXd ybx(m,m);
    double diag_norm = yby.diagonal().array().square().sum();
    double out_norm = (yby.array().square()).sum() - diag_norm;
    if (!( (std::abs(diag_norm - m) <= tol_ortho) && (std::abs(out_norm) <= tol_ortho) )){
        std::cout << "Input y is not orthonormal! continue by Solving Y^TY" << std::endl;
#ifdef DEBUG_ORTHO
        std::cout << "y = \n" << y << std::endl << "yby = \n" << yby << std::endl;
        std::cerr << "diag norm:" << diag_norm << " | and out norm:" <<out_norm << std::endl;
#endif
        // return;
        is_y_orthonormal = false;
    }

    // start with initial orthogonalization
    ortho(n, k, x);

    // X = X - Y * y_coeff
    // Y'BY y_coeff = Y'BX
    // X = X - Y(Y'Y)^{-1}Y'X
    
    Eigen::MatrixXd y_coeff = Eigen::MatrixXd::Identity(m, k); // to store Y'X or (Y'Y)^{-1}Y'X
    // since y will not change, we will compute Cholesky factor of yby only once
    Eigen::LLT<Eigen::MatrixXd> factYBY(yby);
    if(!is_y_orthonormal) factYBY.compute(yby); //only solve Y'Y y_coeff = Y'X if y'y is not I

    // do ortho while overlap norm > tol_ortho
    double norm_overlap = 10.0; // big init value
    const int ITER_MAX = 10;//10;
    int iter_cnt = ITER_MAX; // max 10 iters
    while(norm_overlap >= tol_ortho){
#ifdef DEBUG_LOBPCG
        // std::cout <<"--- "<< ITER_MAX - iter_cnt << "th iter" <<" ---" << std::endl;
#endif
        // orthogonalize x against y, X = X - Y(Y'Y)^{-1}Y'X
        // X = X - Y * y_coeff
        // Y'Y y_coeff = Y'X
        ybx = y.transpose() * x;

        if(!is_y_orthonormal){
            y_coeff = factYBY.solve(ybx);
            x = x - y * y_coeff; // (Y'Y)^{-1}Y'X
        }else{
            x = x - y * ybx; // Y'X
        }
        // compute overlap = Y^TX, size(m, k)
        Eigen::MatrixXd overlap;

        // ortho_qr(n, k, x);
        ortho(n, k, x);

        // compute overlap norm ||y^T x|| after x orthonormalization
        overlap = y.transpose() * x;
        norm_overlap = overlap.norm();
        // 判定条件：norm_overlap >= tol_ortho
#ifdef DEBUG_LOBPCG
        // std::cout << " norm = " << norm_overlap << std::endl;
#endif
        --iter_cnt;
        if(iter_cnt < 0){ // to many ortho iterations, exit 
            std::cerr << "Too many iterations in ortho_against_x. Failed to reach tolerance." << std::endl;
            break;
        }
    }
#ifdef DEBUG_LOBPCG
    // std::cout << "----- end ortho against y -----" << std::endl;
#endif
}



/**
 * B-Orthogonalizes the matrix `x` against the matrix `y` given y and by.
 * 
 * b-orthogonalize x(n, k) against y(n, m)
 * 
 * @param n The number of rows in `x` and `y`.
 * @param m The number of columns in `y`.
 * @param k The number of columns in `x`.
 * @param x The matrix `x` to be b-orthogonalized.
 * @param y The matrix `y` used for b-orthogonalization.
 * @param by The matrix `by` used for b-orthonormalization.
 */
void b_ortho_against_y(int n, int m, int k, Eigen::MatrixXd& x, const Eigen::MatrixXd& y, const Eigen::MatrixXd& by){
    // b-orthogonalize x(n, k) against y(n, m)
#ifdef DEBUG_LOBPCG
    // std::cout << "----- start b-ortho against y -----" << std::endl;
#endif
    double tol_ortho = 1.0e-10;//1.0e-13;

    // first check input y is b-orthonormal
    Eigen::MatrixXd yby = by.transpose() * y;
    Eigen::MatrixXd ybx(m,m);
    bool is_y_orthonormal = true;
    double diag_norm = yby.diagonal().array().square().sum();
    double out_norm = (yby.array().square()).sum() - diag_norm;
    if (!( (std::abs(diag_norm - m) <= tol_ortho) && (std::abs(out_norm) <= tol_ortho) )){
        std::cerr << "Input y is not b-orthonormal!" << std::endl;
        std::cerr << "diag norm:" << diag_norm << " | and out norm:" <<out_norm << std::endl;
        // return;
        std::cout << "continue by Solving Y^TBY" << std::endl;
        is_y_orthonormal = false;
    }
    // X = X - Y * y_coeff
    // Y'BY y_coeff = Y'BX
    // X = X - Y(Y'BY)^{-1} Y'BX
    Eigen::LLT<Eigen::MatrixXd> factYBY(yby);
    Eigen::MatrixXd y_coeff = Eigen::MatrixXd::Identity(m, k);
    // since y will not change, we will compute Cholesky factor of yby only once
    if(!is_y_orthonormal) factYBY.compute(yby);

    // do ortho while overlap norm > tol_ortho
    double norm_overlap = 10.0; // big init value
    const int ITER_MAX = 10;//10;
    int iter_cnt = ITER_MAX; // max 10 iters
    while(norm_overlap >= tol_ortho){
#ifdef DEBUG_LOBPCG
        // std::cout <<"--- "<< ITER_MAX - iter_cnt << "th iter" <<" ---" << std::endl;
#endif
        // orthogonalize x against y, X = X - Y(Y'BY)^{-1}(YB)'X
        // X = X - Y * y_coeff
        // Y'BY y_coeff = Y'BX
        // X = X - Y(Y'BY)^{-1}(YB)'X
        ybx = by.transpose() * x;

        if(!is_y_orthonormal){
            y_coeff = factYBY.solve(ybx);
            x = x - y * y_coeff; // (Y'BY)^{-1}Y'X
        }else{
            x = x - y * ybx; // Y'BX
        }

        // compute overlap = (BY)^TX, size(m, k)
        Eigen::MatrixXd overlap;

        ortho(n, k, x);

        overlap = by.transpose() * x;
        norm_overlap = overlap.norm();

        --iter_cnt;
        if(iter_cnt < 0){ // to many ortho iterations, exit
            std::cerr << "Too many iterations in ortho_against_x. Failed to reach tolerance." << std::endl;
            break;
        }
    }
#ifdef DEBUG_LOBPCG
    // std::cout << "----- end b-ortho against y -----" << std::endl;
#endif
}
