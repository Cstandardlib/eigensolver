#include "lobpcg.h"
#include "ortho.h"

#include <Eigen/Dense>

#include <chrono>

#include <random> // make a random init guess
#include <iostream>



/**
 * @brief Main driver for the Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG) method.
 *
 * This function is the primary driver for the LOBPCG algorithm, which is an iterative method for solving
 * large sparse eigenvalue problems. The function aims to compute the specified number of eigenpairs
 * of a given matrix or matrix pencil.
 * 
 * solve Ax = \lambda Bx, block version AX = \lambda BX, where B is optional
 *
 * @param[in] avec Ax - Function pointer for the matrix-vector multiplication with the main operator A.
 * @param[in] precnd Function pointer for applying the preconditioner. The function should accept the size of the matrix,
 *               the number of vectors, the shift value, the input vectors, and the vectors to be preconditioned.
 * @param[in] bvec Bx - Function pointer for the matrix-vector multiplication with the B operator (only used in the generalized case).
 * @param[in,out] eig Reference to a vector where the computed eigenvalues will be stored, in ascending order if converged.
 * @param[in,out] evec Reference to a blockvector with an initial guess for the eigenvectors; upon successful convergence, it will contain
 *             the computed eigenvectors.
 * @param[in] n Size of the matrix to be diagonalized.
 * @param[in] n_eigenpairs Number of required eigenpairs to be computed.
 * @param[in] n_max_subspace Maximum size of the search subspace. should be greater than n_eigenpairs.
 *                           for better convergence, a value larger than n_eigenpairs (eg., n_eigenpairs + 10) is recommended.           
 * @param[in] solving_generalized Boolean flag indicating whether the problem is a generalized eigenvalue problem.
 * @param[in] max_iter Maximum allowed number of iterations for convergence.
 * @param[in] tol Convergence threshold for the eigenvalues.
 * @param[in] shift Diagonalshift value for the eigenvalue problem.
 * @param[in] verbose Logical flag to control the verbosity, whether to print intermediate results
 *                  at each iteration (eigenvalues, residuals...)
 */
int lobpcg_solve(
    void (*avec)(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& avecs),             // AX, external operator for X nxm
    void (*precnd)(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& tvecs),   // TX, external operator for X nxm
    void (*bvec)(int n, int m, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& bvecs),             // BX, external operator for X nxm
    Eigen::VectorXd& eig,   // lambdas, should be allocated size n_max_subspace
    Eigen::MatrixXd& evec,  // X, should be allocated size (n,n_max_subspace)
    int n,                  // size of A, B
    int n_eigenpairs,       // X(n, n_eigenpairs)
    int n_max_subspace,     // maximum size of the search subspace
    bool solving_generalized, // true for generalized, B=I if false
    int max_iter,           // maximum number of iterations
    double tol,             // tolerance for convergence test
    double shift,           // shift for Cholesky & precnd
    bool verbose            // print intermediate results if true
){
    /* ! note that eig and evec should be allocated n_max_subspace and (n,n_max_subspace),
        not n_eigenpairs and (n,n_eigenpairs) */
    
    // --- 0. allocate and initialize ---
#ifdef DEBUG_LOBPCG
    std::cout << "--- 0. allocate and initialize ---" << std::endl;
#endif

    /* 0.0 allocate memory for lapack */
    /* no need at present, use other wrapped API instead*/
    // lwork
    // allocate work, used in dsyev

    /* 0.1 allocate memory for expansion space, and corresponding Matrix-Vector results and residuals
       i.e. X(n, n_max_subspace), W(n, n_active), P(n, n_active)
       if using one unified block of memory, it will need V(n, 3*n_max_subspace)
    2 ways of doing x, w, p
       1. we can construct v, av, A_reduced = V'AV only when needed by Rayleigh-Ritz instead
       of stroring x, w, p together in v
       2. we can use v as [x p w] and access different parts of v as x, p, w

       https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
       As it is often confusing to use Eigen's templated base classes,
       we will use <1> (i.e. to construct A_reduced as needed) instead of <2>.
    */
    int size_space = 3*n_max_subspace;  // Matrix<typename _Scalar, int _Rows, int _Cols, ...>
    Eigen::MatrixXd v(n, size_space);   // dynamic searching space V = [x p w], n x (n_max_subspace+2*n_active)
    Eigen::MatrixXd av(n, size_space);
    /* since we only use bx, bp for b-orthogonalization of w(in need of [x p], [wx wp])
        we will not store bv but to get b[x p] when needed(in 2.7) instead
    */
    // Eigen::MatrixXd bv(n, size_space);
    // due to Lazy Evaluation and Aliasing of Eigen, we will not need to
    // set values of v,av,bv until we need them to construct A_reduced=v'av
    v.setZero(); av.setZero(); // bv.setZero();

    Eigen::MatrixXd r(n, n_max_subspace);       // residuals
    Eigen::MatrixXd w(n, n_max_subspace);       // preconditioned residuals
    Eigen::MatrixXd p(n, n_max_subspace); 
    r.setZero(); w.setZero(); p.setZero();
 
    /* 0.2 allocate memory for temporary copies of X, AX, BX */
    Eigen::MatrixXd x(n, n_max_subspace);  
    Eigen::MatrixXd ax(n, n_max_subspace); 
    Eigen::MatrixXd bx(n, n_max_subspace);
    x.setZero(); ax.setZero(); bx.setZero();
    Eigen::MatrixXd aw(n, n_max_subspace);
    Eigen::MatrixXd ap(n, n_max_subspace); 
    aw.setZero(); ap.setZero();
    Eigen::MatrixXd bw(n, n_max_subspace); // only used for b-ortho w
    Eigen::MatrixXd bp(n, n_max_subspace); 
    bw.setZero(); bp.setZero();

    /* 0.3 allocate memory for the reduced matrix and its eigenvalues */
    // fixed size A_reduced(size_space, size_space), eig_reduced(size_space)
    // the subspace expands from [X] to [x p w]
    // where X has the width of n_max_subspace, while W and P have the width of n_active
    Eigen::MatrixXd A_reduced(size_space, size_space); A_reduced.setZero();
    Eigen::VectorXd eig_reduced(size_space); eig_reduced.setZero();

    /* 0.4 allocate memory for convergence check */
    // std::vector<bool> done(n_max_subspace);
    Eigen::VectorXd r_norm_2(n_max_subspace); // 2-norm

    /* 0.5 initialize time variables */
    auto tp_start = get_current_time(); // tp for time point, different from duration
    auto tp_end = get_current_time();
    auto tp_1 = get_current_time(); // temp time point 1 for inner durations
    auto tp_2 = get_current_time(); // temp time point 2 for inner durations
    // duration default 0
    std::chrono::duration<double> t_solveRR;
    std::chrono::duration<double> t_avec;
    std::chrono::duration<double> t_ortho;
    std::chrono::duration<double> t_total;

    int eig_flag = LOBPCG_CONSTANTS::eig_success; // flag for eigensolver

    // --- 1. first iter ---  explicitly do the first Rayleigh-Ritz of X'AX
    tp_start = get_current_time();      // start timing for the whole iteration algorithm
    /* 1.0 check for guess */
    /* check for initial guess. If not set(all zeros), generate a random
        guess from values in Uniform[0,1), and then orthogonalize */
    check_init_guess(n, n_max_subspace, evec);
    // now evec contains orthogonal initial guess

    /* if solving generalized problem, compute bvec and do b-ortho */
    if(solving_generalized) {
        bvec(n, n_max_subspace, evec, bx);      // bx <- b*evec
        b_ortho(n, n_max_subspace, evec, bx);   // b-ortho bx against evec
    }

    /* --- Rayleigh-Ritz --- */
    /* 1.1 construct the reduced matrix and diagonalization to get first-round eigenpairs */
    /* x <- evec */
    x = evec;
    // v.leftCols(n_max_subspace) = evec;
#ifdef DEBUG_LOBPCG
    // std::cout << "init guess x = \n" << x << std::endl << "and x'x = \n" << x.transpose()*x << std::endl;
#endif
    /* ax <- a*evec */
    tp_1 = get_current_time();
    avec(n, n_max_subspace, evec, ax);
    tp_2 = get_current_time();
    t_avec += tp_2 - tp_1;
    // av.leftCols(n_max_subspace) = ax;
    /* BV <- bx, bx generated before in 1.0 */
    // if(solving_generalized){bv.leftCols(n_max_subspace) = bx;}

    /* first round, A_reduced <- X' * AX */
    A_reduced = x.transpose() * ax;
    /* do eigen */
    /* first round, solves only A_reduced = x'ax (n_max_subspace, n_max_subspace)
       get to use bigger A_reduced later
       as the subspace expands from [X] to [x p w]
    */
#ifdef DEBUG_LOBPCG
    std::cout << "--- first round, do eigen & solves only A_reduced = x'ax (n_max_subspace, n_max_subspace) ---" << std::endl;
#endif
    tp_1 = get_current_time();
    /* A_reduced u = \lambda u */
    eig_flag = selfadjoint_eigensolver(A_reduced, eig_reduced, n_max_subspace); // V=[X], use(n, n_max_subspace) of (n, 3*n_max_subspace)
    if(eig_flag == LOBPCG_CONSTANTS::eig_fail){
        std::cerr << "eigensolver failed in first round" << std::endl;
        return LOBPCG_CONSTANTS::fail;
    }
    /* now A_reduced(n_max_subspace, n_max_subspace) contains eigenvectors and eig_reduced contains eigenvalues */
    tp_2 = get_current_time();
    t_solveRR += tp_2 - tp_1;
    eig = eig_reduced.head(n_max_subspace);

    /* 1.2 compute the Ritz vectors, X, AX and if required, BX */
    /* update new guess X_1 = X_0 u, update corresponding V[X] left n_max_subspace columns */
    /* x <- x * A_reduced*/
    x = x * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
    // v.leftCols(n_max_subspace) = v.leftCols(n_max_subspace) * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
    /* ax <- ax * A_reduced*/
    ax = ax * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
    // av.leftCols(n_max_subspace) = av.leftCols(n_max_subspace) * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
    /* bx <- bx *A_reduced*/
    if(solving_generalized){
        bx = bx * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
        // bv.leftCols(n_max_subspace) = bv.leftCols(n_max_subspace) * A_reduced.topLeftCorner(n_max_subspace,n_max_subspace);
    }

    /* 1.3 compute the residuals and norms */
    /* r <- Ax - eig Bx */
    r = ax; // r = av.leftCols(n_max_subspace);
    if(solving_generalized){
        for(int i=0; i<n_max_subspace; ++i) r.col(i) -= eig(i) * bx.col(i);
    }else{ // b = I
        for(int i=0; i<n_max_subspace; ++i) r.col(i) -= eig(i) * x.col(i);
    }

    /* 1.4 compute the preconditioned residuals W = TR */
    precnd(n, n_max_subspace, r, w);   // w = t * r
    /* [X W]  0~n_max_subspace-1 | n_max_subspace~2*n_max_subspace-1 */
    

    /* 1.5 orthogonalize W; and then orthonormalize it */
    tp_1 = get_current_time(); // t_ortho

    if(solving_generalized) { /* compute and b-orthogonalize bw */
        b_ortho_against_y(n, n_max_subspace, n_max_subspace, x, bx, w);
        bvec(n, n_max_subspace, w, bw);     // bw = b*w
        b_ortho(n, n_max_subspace, w, bw);  // b-orthogonalize w
    } else {
        // ortho(n, n_max_subspace, w);
        ortho_against_y(n, n_max_subspace, n_max_subspace, w, x); // orthogonalize w to x
    }
    tp_2 = get_current_time();
    t_ortho += tp_2 - tp_1;

    /* 1.6 build first round v and av */
    /* x */
    v.leftCols(n_max_subspace) = x;
    av.leftCols(n_max_subspace) = ax;
    int index_w = n_max_subspace;
    v.middleCols(index_w, n_max_subspace) = w; // now v = [X w], w of width n_max_subspace
#ifdef DEBUG_LOBPCG
    std::cout << "before main loop: v = \n" << v << std::endl << "before main loop: av = \n" << av << std::endl;
#endif
    // --- 2. main loop ---
    /* prepare for the main loop, initialize parameters */
    // X size(n, n_max_subspace), W size(n, n_active), P size(n, n_active) in loop
    // active searching size of w and p, changed when checking residuals and locking

    int n_active = n_max_subspace;
    const int ACTIVE = 1;
    const int INACTIVE = 0;
    Eigen::VectorXi activeMask(n_max_subspace); /* active mask can be 0 or 1, indicating whether x_i, w_i and p_i are active or not*/
    activeMask.setConstant(ACTIVE); // at first, all vecs are active

    if(verbose){
        std::cout << "\nLOBPCG iter starts with max_iter = " << max_iter
            << " and tolerance = " << tol << std::endl;
    }
/* ----------------------------------------------------------------------*/
    /* start main loop */
for(int iter = 0; iter < max_iter; ++iter){


    /* --- Rayleigh-Ritz --- */
#ifdef DEBUG_LOBPCG
    std::cout << "----- starts Loop #" << iter << " -----" << std::endl;
#endif
    /* 2.1 matrix-blockvector multiplication, calculate aw = a*w */
    tp_1 = get_current_time(); // avec
    avec(n, n_active, w, aw);   // aw = a*w
    tp_2 = get_current_time();
    t_avec += tp_2 - tp_1;
    av.middleCols(index_w, n_max_subspace) = aw;

    /* 2.2 construct the reduced matrix and diagonalization */
    /* v = [x p w]*/
    /* size[
        x [0..n_max_subspace), n_max_subspace in total
        w [n_max_subspace..n_max_subspace+n_active) n_active in total
        p [n_max_subspace+n_active..n_max_subspace+2*n_active) n_active in total
       ]*/
    /* notice: here w and p should be of size(n, n_active),
       or the assignment of v will fail
    */
    int n_working_space = n_max_subspace + 2*n_active; // current valid v size v(n, n_v)
    // v.resize(n_working_space)

    // v.leftCols(n_max_subspace) = x;
    // v.middleCols(index_w, n_active) = w;
    // if(iter > 0) v.middleCols(n_max_subspace + n_active, n_active) = p; // first iter there is no p

    /* av = a*v, v(n, n_working_space) */
    // avec(n, n_working_space, v, av); av is made the same time as v by merging [x p w]
    /* notice that X full(active and locked) still engage in Rayleigh-Ritz */
    A_reduced = v.transpose() * av; // full v'av, (n_working_space x n_working_space)
#ifdef DEBUG_LOBPCG
    std::cout << "A_reduced size = " << A_reduced.rows() << " x " << A_reduced.cols() << std::endl;
#endif
    tp_1 = get_current_time(); // t_solveRR
    eig_flag = selfadjoint_eigensolver(A_reduced, eig_reduced , n_working_space);
    if(eig_flag == LOBPCG_CONSTANTS::eig_fail){
        std::cerr << "eigensolver failed in round " << iter << std::endl;
        return LOBPCG_CONSTANTS::fail;
    }
    tp_2 = get_current_time();
    t_solveRR += tp_2 - tp_1;
    eig = eig_reduced.head(n_max_subspace);

    /* check the eigenpairs */

    /* 2.3 update X, AX and, if required BX */
    // from now on x, ax, bx store new x^{k+1}, ax^{k+1} and bx^{k+1}
    x = v.leftCols(n_working_space) * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    ax = av.leftCols(n_working_space) * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    if(solving_generalized){
        // bx = v.leftCols(n_working_space) * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
// ???
// how to deal with bx = bv(n, n_working_space) * A_reduced(n_working_space, n_max_subspace)
        bx = bx * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    }
#ifdef DEBUG_LOBPCG
    std::cout << "--- starts compute residuals & norms ---" << std::endl;
#endif
    /* 2.4 compute residuals & norms */
    r = ax; // r(n, n_max_subspace)
    /* loop for eigenpairs */
    for(int i = 0; i < n_max_subspace; ++i){
         // converged vecs should not engage in computing residuals
        if(activeMask(i) != ACTIVE) continue; // locked, skip

        if(solving_generalized){
            r.col(i) -= eig(i) * bx.col(i); // r = ax - eig*bx
        } else {
            r.col(i) -= eig(i) * x.col(i); // r = ax - eig*x
        }
        r_norm_2(i) = r.col(i).norm(); // r_norm_2(i) = ||r_i||_2
    }
#ifdef DEBUG_LOBPCG
    std::cout << "--- starts checking convergence and locking ---" << std::endl;
#endif
    /* 2.5 check convergence and locking */
    // i range from 0 to n_max_subspace-1,   total n_max_subspace vectors
    for(int i = 0; i < n_max_subspace; ++i){
        // alreadhy locked
        if(activeMask(i) != ACTIVE) continue;

        if(r_norm_2(i) < tol*std::sqrt(n) && iter >0){
            activeMask(i) = INACTIVE; // lock the vector
        }
        if(activeMask(i) == ACTIVE){
            // if not locked, update the active size
            // every vec after this one should be active
            activeMask.segment(i+1, n_max_subspace-1-i) = Eigen::VectorXi::Constant(n_max_subspace-1-i, ACTIVE);
        }
    }
#ifdef DEBUG_LOBPCG
    std::cout << "--- checking overall convergence ---" << std::endl;
#endif
    // checking overall convergence for n_eigenpairs needed
    if((activeMask.head(n_eigenpairs).array() == INACTIVE).all()){
        // all vectors are locked, i.e. converged
        evec = x; // x(n, n_max_subspace)
        return LOBPCG_CONSTANTS::success;
    }

    /* 2.6 check active eigenvalues and update blockvectors X, P, W */
    /* 2.6.1 compute the number of active eigenvalues */
    n_active = activeMask.count(); // count number of active eigenvalues
    /* x can be partitioned into [X_locked | X_active]
       and we care the latter part: X_active*/

    /* 2.6.2 compute the new P and AP blockvector,
        by computing the coefficients u_p and then orthonalizing to u_x */
    /* [x p w]*/
    /* coefficients of p can be computed as follows:
        !!! notice that
        !!! we only deal with active part of X, and thus W and P
        !!! of size(n, n_active)
        0. we get coefficients from A_reduced(n_working_space, n_working_space)
        A_reduced [ c_x
                    c_w
                    c_p ]
        1. x^{k+1} = x^k * c_x + w^k * c_w + p^k * c_p
        2. p^k = x^{k+1} - x^k = (x^k - I) * c_x + w^k * c_w + p^k * c_p
        3. size of p: (n, n_active)
        4. we get active x^k here from x(:, n_max_subspace-n_active to n_max_subspace-1)
            corresponding coeff block(n_max_subspace-n_active, n_max_subspace-n_active,n_active,n_active)
    */
#ifdef DEBUG_LOBPCG
    std::cout << "--- start computing the coefficients of p ---" << std::endl;
#endif
    // -- compute the coefficients of p
    /* coeff consists of c_x, c_w, c_p */
    Eigen::MatrixXd coeff = A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    auto c_p = coeff.block(n_max_subspace-n_active, n_max_subspace-n_active,n_active,n_active);
    c_p -= Eigen::MatrixXd::Identity(n_active, n_active);
    // ortho coeff p against coeff x
    // size(n_working_space, n_active)
    Eigen::MatrixXd coeff_p = coeff.middleCols(n_max_subspace-n_active,n_active);
    ortho_against_y(n_working_space, n_max_subspace, n_active, coeff_p, coeff);
    // -- end of computing the coefficients of p
    /* p(n, n_active) = v(n, n_working_space) * coeff_p(n_working_space, n_active)
       where v(n, n_working_space) contains old [x p w]
    */
    p = v * coeff_p;
    ap = v * coeff_p;
    if(solving_generalized) bp = v * coeff_p;
    // -- now p, ap, bp contains new values of step k+1
    
    /* 2.6.4 compute the preconditioned residuals W */
    precnd(n, n_active, r, w); // w = tr

    /* 2.7.1 orthogonalize W to X, P; and then orthonormalize it */
    const Eigen::MatrixXd xp = v.leftCols(n_max_subspace + n_active);
    tp_1 = get_current_time();
    if(solving_generalized){ // get orthogonalized w and bw
        // const Eigen::MatrixXd b_xp = bv.leftCols(n_max_subspace + n_active);
        Eigen::MatrixXd b_xp(n, n_max_subspace + n_active);
        // assert(bx.rows() == bp.rows());
        b_xp << bx, bp;
        b_ortho_against_y(n, n_max_subspace + n_active, n_active, w, xp, b_xp); // b-ortho w to [x p]
        bvec(n, n_active, w, bw); // bw = b*w
        b_ortho(n, n_active, w, bw); // b-ortho
    } else {
        ortho_against_y(n, n_max_subspace+n_active, n_active, w, xp);// ortho w to [x p]
    }
    
    tp_2 = get_current_time(); // t_ortho
    t_ortho += tp_2 - tp_1;

    /* 2.7 update corresponding blockvectors V from X, P, W*/
    /* -- these will be used in the next iteration for constructing
        reduced matrix A_reduced = v'av for Rayleigh-Ritz */
    /* now that we have new p, ap and bp, we shall put them in v */
    /* v(n, n_working_space) [x(n, n_max_subspace) p(n, n_active) w(n, n_active)]*/
    /* x */
    v.leftCols(n_max_subspace) = x;
    av.leftCols(n_max_subspace) = ax;
    /* p */
    v.middleCols(n_max_subspace, n_active) = p;
    av.middleCols(n_max_subspace, n_active) = ap;
    // if(solving_generalized) bv.middleCols(n_max_subspace, n_active) = bp;
    /* w */
    v.middleCols(n_max_subspace+n_active, n_active) = w;
    av.middleCols(n_max_subspace+n_active, n_active) = aw;

} // end -  main loop for

    /* verbose output */
    tp_end = get_current_time();
    t_total = tp_end - tp_start;
    std::cout << "LOBPCG total time: " << t_total.count() << std::endl;

    // --- 3. clean up ---
    // deallocate and check

}

// intermediate tool functions
// check_init_guess

// check_guess - 检查给定的矩阵evec是否为零矩阵或者是否正交。
// 如果evec为零矩阵，将用随机数填充并正交化。
// 如果evec不正交，将进行正交化处理。
/**
 * @brief Checks the validity of the given guess matrix and performs orthogonalization if necessary.
 * 
 * If evec is zero, provide a random guess, filled by uniform [0,1).
 * 
 * @param n The number of rows in the guess matrix.
 * @param m The number of columns in the guess matrix.
 * @param evec The guess matrix to be checked and orthogonalized if needed.
 */
void check_init_guess(int n, int m, Eigen::MatrixXd& evec) {
    double enorm = evec.norm(); // 计算矩阵的范数
// #ifdef DEBUG_LOBPCG
//     std::cout << "-- start check_init_guess --" << std::endl;
// #endif
    if (std::abs(enorm) < LOBPCG_CONSTANTS::tol_zero) {
        // if evec is all zero(which means no initial guess is provided by user), fill it with random values in U[0,1)
        std::default_random_engine engine;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n*m; ++i) evec(i) = dist(engine);
        ortho(n, m, evec);// orthogonalize evec
    } else {
        // compute overlao matrix M = X^T X (size m*m), and do ortho
        // first check for orthogonality: whether X^TX == Identity
        Eigen::MatrixXd overlap = evec.transpose() * evec;
        double diag_norm = overlap.diagonal().array().square().sum();   // square sum of overlap diagonal
        double out_norm = (overlap.array().square()).sum() - diag_norm; // square sum of non-diagonal elements
        
        // 检查对角线元素是否为1，非对角线元素是否为0
        if (std::abs(diag_norm - m) > 1e-10 || std::abs(out_norm) > 1e-10) {
            // 如果不正交，进行正交化处理
            ortho(n, m, evec);
        }
    }
// #ifdef DEBUG_LOBPCG
//     std::cout << "-- end check_init_guess --" << std::endl;
// #endif
}


// timing tools
/*
 * wrap for std::chrono::high_resolution_clock::now();
 * auto start = get_current_time();
 * auto end = get_current_time();
 * std::chrono::duration<double> duration = end - start;
 */
std::chrono::time_point<std::chrono::high_resolution_clock> get_current_time() {
    return std::chrono::high_resolution_clock::now();
}


// dense eigen solver
/* overwrite A and eig
 * solves eigenvalue and eigenvectors of a selfadjoint(i.e. inverse equal to its adjoint) matrix A.
 * solve A_nn(n, n) topleft corner of size n input
*/


int selfadjoint_eigensolver(Eigen::MatrixXd &A_to_be_eigenvecs, Eigen::VectorXd &eig, int n)
{
#ifdef DEBUG_LOBPCG
    std::cout << "--- starts eigensolver ---" << std::endl;
    std::cout << "A size: " << A_to_be_eigenvecs.size() <<", ("<<A_to_be_eigenvecs.rows()<<", "<<A_to_be_eigenvecs.cols()<<")" <<std::endl;
#endif
    if(A_to_be_eigenvecs.rows() != n || A_to_be_eigenvecs.cols()!=n){
        std::cerr << "input matrix must be of size (n,n) = "<<"("<<n<<", "<<n<<")"<<std::endl; return LOBPCG_CONSTANTS::eig_fail;
    }
    if(eig.size() < n){
        std::cerr << "eigenvalues vector must not be smaller than n = "<<n<<std::endl; return LOBPCG_CONSTANTS::eig_fail;
    }
#ifdef DEBUG_EIGENSOLVER
    std::cout << "A size: " << A_to_be_eigenvecs.size() <<", ("<<A_to_be_eigenvecs.rows()<<", "<<A_to_be_eigenvecs.cols()<<")" <<std::endl;
    // std::cout << "original A = \n" << A_to_be_eigenvecs << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "solving A(" << n << "," << n << ") matrix" << std::endl;
    // Eigen::MatrixXd A_nn = A_to_be_eigenvecs.topLeftCorner(n,n);// topleft corner of A_to_be_eigenvecs
    std::cout << "A to be solved = \n" << A_to_be_eigenvecs << std::endl;
#endif
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A_to_be_eigenvecs);
    if (eigensolver.info() == Eigen::Success)
    {
        eig.head(n) = eigensolver.eigenvalues();
        // A_to_be_eigenvecs.topLeftCorner(n,n) = eigensolver.eigenvectors();
        A_to_be_eigenvecs = eigensolver.eigenvectors();
#ifdef DEBUG_EIGENSOLVER
        std::cout << "The eigenvalues are:\n" << eig << std::endl;
        std::cout << "The eigenvectors are:\n" << A_to_be_eigenvecs << std::endl;
#endif
    }
    else
    {
        std::cerr << "Eigen decomposition failed." << std::endl;
        return LOBPCG_CONSTANTS::eig_fail;
    }
    return LOBPCG_CONSTANTS::eig_success;
#ifdef DEBUG_LOBPCG
    std::cout << "--- end eigensolver ---" << std::endl;
#endif
}