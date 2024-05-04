#include "lobpcg.h"
#include "ortho.h"
#include "utils.h"

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cassert>


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

    // 设置输出格式，这里设置浮点数的精度为5位小数
    Eigen::IOFormat fmt(5);
    // Eigen::IOFormat fmt(5, Eigen::DontAlignCols, ", ", ";\n", "", "", "[", "]");
    
// --- 0. allocate and initialize ---
    /* 0.1 allocate memory for expansion space, and corresponding Matrix-Vector results and residuals
       i.e. X(n, n_max_subspace), P(n, n_active), W(n, n_active)
       if using one unified block of memory, it will need V(n, 3*n_max_subspace)
    2 ways of doing x, w, p
       1. we can construct v, av, A_reduced = V'AV only when needed by Rayleigh-Ritz instead of stroring x, w, p together in v
       2. we can use v as [x p w] and access different parts of v as x, p, w

       https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html  block params need to be passed as MatrixBase<>
       As it is often confusing to use Eigen's templated base classes, we will use <1> (i.e. to construct A_reduced as needed).
    */
    int size_space = 3*n_max_subspace;  // total max space size of V = [x p w]
    // dynamic searching space V(n, n_working_space) = [x p w]
    // n_working_space=n_max_subspace+2*n_active will change in the main loop
    Eigen::MatrixXd v(n, size_space);
    Eigen::MatrixXd av(n, size_space);
    /* since we only use bx, bp for b-orthogonalization of w(in need of [x p], [wx wp])
        we will not store bv but to get b[x p] when needed(in 2.7) instead */
    // Eigen::MatrixXd bv(n, size_space);
    v.setZero(); av.setZero(); // bv.setZero();
 
    /* 0.2 allocate memory for temporary copies of X, AX, BX */
    Eigen::MatrixXd x(n, n_max_subspace);       // X(n, n_max_subspace)
    Eigen::MatrixXd ax(n, n_max_subspace); 
    Eigen::MatrixXd bx(n, n_max_subspace);
    x.setZero(); ax.setZero(); bx.setZero();

    Eigen::MatrixXd r(n, n_max_subspace);       // residuals, r(n, n_active)
    Eigen::MatrixXd w(n, n_max_subspace);       // preconditioned residuals w(n, n_active)
    Eigen::MatrixXd aw(n, n_max_subspace);
    Eigen::MatrixXd bw(n, n_max_subspace);      // only used for b-ortho w
    r.setZero();  w.setZero(); aw.setZero(); bw.setZero();

    Eigen::MatrixXd p(n, n_max_subspace);       // implicit difference between current and previous eigvec approx, P(n, n_active)
    Eigen::MatrixXd ap(n, n_max_subspace);
    Eigen::MatrixXd bp(n, n_max_subspace); 
    p.setZero(); ap.setZero(); bp.setZero();

    /* 0.3 allocate memory for the reduced matrix and its eigenvalues */
    /* soft locking by Knyazev, all X, active p, active w will engage in Rayleigh-Ritz */
    // dynamic size A_reduced(n_working_space, n_working_space), eig_reduced(n_working_space)
    // the subspace expands from V=[X] to v=[x p w], A_reduced = v'av
    // where X has the width of n_max_subspace, while W and P have the width of n_active
    Eigen::MatrixXd A_reduced(size_space, size_space); A_reduced.setZero();
    Eigen::VectorXd eig_reduced(size_space); eig_reduced.setZero();

    /* 0.4 allocate memory for convergence check */
    // Eigen::VectorXi activeMask(n_max_subspace) will be defined before we start main loop
    Eigen::VectorXd r_norm_2(n_max_subspace); // 2-norm of preconditioned-residuals

    /* 0.5 initialize time variables */
    auto tp_start   = get_current_time(); // tp means 'time point', different type from duration
    auto tp_end     = get_current_time();
    auto tp_1       = get_current_time(); // temp time point 1 for inner durations
    auto tp_2       = get_current_time(); // temp time point 2 for inner durations
    // duration = tp_2 - tp_1, duration.count() returns second elapsed
    std::chrono::duration<double> t_solveRR, t_avec, t_ortho, t_total;

    int eig_flag = LOBPCG_CONSTANTS::eig_success; // flag for eigensolver

// --- 1. first iter ---  explicitly do the first Rayleigh-Ritz of X'AX
    tp_start = get_current_time();      // start timing for the whole iteration algorithm

    /* 1.0 check for guess */
    /* check for initial guess. If not set(all zeros), generate a random
        guess from values in Uniform[0,1), and then orthogonalize */
    evec << -0.045008846777,-0.012487056830, 0.451774729742,
            -0.595457985511,-0.302271780449, 0.262108946955,
            -0.683658739782,-0.026606629965,-0.547438361389,
            -0.214624268522,-0.097081891179, 0.451904306374,
            -0.326352098345, 0.603950319551, 0.108625784707,
            -0.115209128264, 0.356526801138, 0.228854930509,
            0.075833974684, 0.399936355871,-0.300634194996,
            0.045339444223,-0.459789809603,-0.236122012805,
            -0.048640002418,-0.187403125974, 0.113945458005;
    check_init_guess(n, n_max_subspace, evec); // now evec contains orthogonal initial guess
#ifdef DEBUG_LOBPCG
std::cout << "initial guess = \n" << evec << std::endl;
#endif
    /* if solving generalized problem, compute bvec and do b-ortho */
    if(solving_generalized) {
        bvec(n, n_max_subspace, evec, bx);      // bx <- b*evec
        b_ortho(n, n_max_subspace, evec, bx);   // b-ortho bx against evec
    }
#ifdef DEBUG_LOBPCG
if(solving_generalized){
    std::cout << "init bx = \n" << bx << std::endl;
    std::cout << "guess (bx)'b(bx) = \n" << bx.transpose() * bx << std::endl;
}
#endif   
    /* --- Rayleigh-Ritz --- */
    /* 1.1 construct the reduced matrix and diagonalization to get first-round eigenpairs */
    x = evec; /* x <- evec */

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
       as the subspace expands from [X] to [x p w] */

    tp_1 = get_current_time();
    /* A_reduced u = \lambda u */
#ifdef DEBUG_LOBPCG
std::cout << "first round A_reduced = \n" << A_reduced << std::endl;
#endif
    eig_flag = selfadjoint_eigensolver(A_reduced, eig_reduced, n_max_subspace); // V=[X], use(n, n_max_subspace) of (n, 3*n_max_subspace)
    if(eig_flag == LOBPCG_CONSTANTS::eig_fail){
        std::cerr << "eigensolver failed in first round" << std::endl;
        return LOBPCG_CONSTANTS::fail;
    }
#ifdef DEBUG_LOBPCG
std::cout << "first round A_reduced diag = \n" << A_reduced << std::endl;
#endif
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
        b_ortho_against_y(n, n_max_subspace, n_max_subspace, w, x, bx);
        bvec(n, n_max_subspace, w, bw);     // bw = b*w
        b_ortho(n, n_max_subspace, w, bw);  // b-orthogonalize w and bw
    } else {
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

// --- 2. main loop ---
#ifdef DEBUG_LOBPCG
    std::cout << "before main loop: v = \n" << v << std::endl << "before main loop: av = \n" << av << std::endl;
#endif
    /* prepare for the main loop, initialize parameters */
    // X size(n, n_max_subspace), W size(n, n_active), P size(n, n_active) in loop

    // active searching size of w and p n_active, changed when checking residuals and locking
    int n_active = n_max_subspace;
    int n_working_space = 2*n_max_subspace; // init [x w], later will be [x p w] size
    const int ACTIVE = 1;
    const int INACTIVE = 0;
    Eigen::VectorXi activeMask(n_max_subspace); /* active mask can be 0 or 1, indicating whether x_i, w_i and p_i are active or not*/
    activeMask.setConstant(ACTIVE); // at first, all vecs are active

    if(verbose){
        std::cout << "\nLOBPCG iter starts with max_iter = " << max_iter
            << " and tolerance = " << tol << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        // std::cout << "iter#  root         eigenvalues         residuals         converged" << std::endl;
        std::cout << std::left;
        std::cout << std::setw(12) << "iter#"
                << std::setw(13) << "root"
                << std::setw(20) << "eigenvalues"
                << std::setw(20) << "residuals"
                << std::setw(11) << "converged"
                << std::endl;
    }


/* ----------------------------------------------------------------------*/
    /* start main loop */
for(int iter = 0; iter < max_iter; ++iter){



    /* --- Rayleigh-Ritz --- */
#ifdef DEBUG_LOBPCG
    std::cout << "\n\n----- starts Loop #" << iter << " -----" << std::endl;
#endif
    /* 2.1 matrix-blockvector multiplication, calculate aw = a*w */
    tp_1 = get_current_time(); // avec
    avec(n, n_active, w, aw);   // aw = a*w(n, n_active)
    tp_2 = get_current_time();
    t_avec += tp_2 - tp_1;
    av.middleCols(index_w, n_active) = aw;
    // av =  [ax         ap          aw]
    // [n_max_subspace | n_active | n_active]
    //                              ^ index_w = n_max_subspace + n_active

#ifdef DEBUG_LOBPCG
    // std::cout << "v and av for constructing reduced matrix: \n";
    // std::cout << "v = \n" << v << std::endl;
    // std::cout << "av = \n" << av << std::endl;
#endif
    /* 2.2 construct the reduced matrix and diagonalization */
    /* v = [x p w]*/
    /* size[
        x [0..n_max_subspace), n_max_subspace in total
        p [n_max_subspace..n_max_subspace+n_active) n_active in total
        w [n_max_subspace+n_active..n_max_subspace+2*n_active) n_active in total
       ]*/
    /* notice: here w and p should be of size(n, n_active),
       or the assignment of v will fail
    */
    
    
    // v.leftCols(n_max_subspace) = x;
    // v.middleCols(index_w, n_active) = w;
    // if(iter > 0) v.middleCols(n_max_subspace + n_active, n_active) = p; // first iter there is no p

    /* av = a*v, v(n, n_working_space) */
    // avec(n, n_working_space, v, av); av is made the same time as v by merging [x p w]
    /* notice that X full(active and locked) still engage in Rayleigh-Ritz */

    //!!!!! when iter #0, v = [x w], no p here, we would not have blank columns in A_reduced, 
    // which leads to incorrect size of A_reduced(n_max_subspace blank columns, 0 eigenvalues and unit eigenvectors)
    n_working_space = n_max_subspace + 2*n_active; // current valid v size v(n, n_working_space)
    if(0 == iter){
        // std::cout << "0th A_reduced"<<std::endl;
        n_working_space = 2*n_max_subspace;
    }
    A_reduced = v.leftCols(n_working_space).transpose() * av.leftCols(n_working_space);

#ifdef DEBUG_LOBPCG
    std::cout << "A_reduced before = \n" << A_reduced << std::endl;
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
#ifdef DEBUG_LOBPCG
    std::cout << "A_reduced after = \n" << A_reduced << std::endl;
#endif
    /* check the eigenpairs */

    /* 2.3 update X, AX and, if required BX */
    // from now on x, ax, bx store new x^{k+1}, ax^{k+1} and bx^{k+1}
    x = v.leftCols(n_working_space) * A_reduced.leftCols(n_max_subspace); // (n, n_max_subspace)
    ax = av.leftCols(n_working_space) * A_reduced.leftCols(n_max_subspace);
    // x = v.leftCols(n_working_space) * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    // ax = av.leftCols(n_working_space) * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    if(solving_generalized){
        // bx = v.leftCols(n_working_space) * A_reduced.topLeftCorner(n_working_space, n_max_subspace);
// ???
// how to deal with bx = bv(n, n_working_space) * A_reduced(n_working_space, n_max_subspace)
        bx = bx * A_reduced.leftCols(n_max_subspace);
    }
#ifdef DEBUG_LOBPCG
// std::cout << "updated A_reduced = \n" << A_reduced << std::endl;
// std::cout << "updated x = \n" << x << std::endl;
// std::cout << "updated ax = \n" << ax << std::endl;
    // std::cout << "--!!-- updated eig = \n" << eig << std::endl;
    // std::cout << "--- starts compute residuals & norms ---" << std::endl;
#endif
    /* 2.4 compute residuals & norms */
/*???
here R[k] = AX[k] - X[k]eig[K]
but ax is AX[k+1]
now v, av is old [k]
x, ax is new[k+1]
*/
    r = ax; // r(n, n_max_subspace), ax is new ax[k]*A_reduced
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
    // std::cout << "--- starts checking convergence and locking ---" << std::endl;
#endif
    /* 2.5 check convergence and locking */
    // i range from 0 to n_max_subspace-1,   total n_max_subspace vectors
    for(int i = 0; i < n_max_subspace; ++i){
        // alreadhy locked
        if(activeMask(i) != ACTIVE) continue;

        if(r_norm_2(i) < tol*std::sqrt(n)
/*???*/
            && iter >0
        ){
            activeMask(i) = INACTIVE; // lock the vector
#ifdef DEBUG_LOBPCG
            // std::cout << "eigvec " << (i) << " is locked" << std::endl;
#endif
        }
        if(activeMask(i) == ACTIVE){
            // if not locked, update the active size
            // every vec after this one should be active
#ifdef DEBUG_LOBPCG
            // std::cout << "activeMask.segment("<< i+1<<", "<< n_max_subspace-1-i<<") =\n" << activeMask.segment(i+1, n_max_subspace-1-i)
            //     <<"\nis now Eigen::VectorXi::Constant("<<n_max_subspace-1-i<<", "<< ACTIVE <<") = \n"
            //     << Eigen::VectorXi::Constant(n_max_subspace-1-i, ACTIVE) << std::endl;
#endif 
            activeMask.segment(i+1, n_max_subspace-1-i) = Eigen::VectorXi::Constant(n_max_subspace-1-i, ACTIVE);
#ifdef DEBUG_LOBPCG
            // std::cout << "eigvecs from " << (i+1) << " to " << n_max_subspace-1 << " are re-activated" << std::endl;
#endif            
        }
    }
    if(verbose){
        for(int i=0; i<n_eigenpairs; ++i){
            // std::cout << "iter#  root         eigenvalues         residuals         converged" << std::endl;
            std::cout  << std::setw(12) << iter
              << std::setw(13) << i
              << std::setw(20) << eig(i)
              << std::setw(20) << r_norm_2(i)
              << std::setw(11) << !activeMask(i) << std::endl;
        }
    }
#ifdef DEBUG_LOBPCG
    // std::cout << "--- checking overall convergence ---" << std::endl;
#endif
    // checking overall convergence for n_eigenpairs needed
    if((activeMask.head(n_eigenpairs).array() == INACTIVE).all()){
        // all vectors are locked, i.e. converged
        std::cout << "> all of required " << n_eigenpairs << " eigenpairs converged! <" << std::endl;
        evec = x; // x(n, n_max_subspace)
        break; // exit main loop
    }

    /* 2.6 check active eigenvalues and update blockvectors X, P, W */
    /* 2.6.1 compute the number of active eigenvalues */
    n_active = activeMask.count(); // count number of active eigenvalues
    int index_x_active = n_max_subspace-n_active;
    // index_p = n_max_subspace
    index_w = n_max_subspace+n_active;
#ifdef DEBUG_LOBPCG
    std::cout << "number of active eigenvalues: " << n_active << std::endl;
#endif
    /* x can be partitioned into [X_locked | X_active] and we care the latter part: X_active*/

    /* 2.6.2 compute the new P and AP blockvector,
        by computing the coefficients u_p and then orthonalizing to u_x */
    /* [x p w]*/
    /* coefficients of p can be computed as follows:
        !!! notice that we only deal with active part of X, and thus W and P of size(n, n_active)
        0. we get coefficients from A_reduced(n_working_space, n_working_space)
        A_reduced [ c_x
                    c_p
                    c_w ]
        1. x^{k+1} = x^k * c_x + p^k * c_p + w^k * c_w
        2. p^k = x^{k+1} - x^k = (x^k - I) * c_x + w^k * c_w + p^k * c_p
        3. size of p: (n, n_active)
        4. we get active x^k here from x(:, n_max_subspace-n_active to n_max_subspace-1)
            corresponding coeff block(n_max_subspace-n_active, n_max_subspace-n_active,n_active,n_active)
    */
#ifdef DEBUG_LOBPCG
    std::cout << "--- start computing the coefficients of p ---" << std::endl;
#endif
    // -- compute the coefficients of p
    /* coeff consists of c_x, c_p, c_w */
    Eigen::MatrixXd coeff = A_reduced.topLeftCorner(n_working_space, n_max_subspace);
    // in c_p: c_x - I where x is active, i.e. n_max_subspace-n_active to n_max_subspace-1
    // auto c_p = coeff.block(n_max_subspace-n_active, n_max_subspace-n_active,n_active,n_active);
    // c_p -= Eigen::MatrixXd::Identity(n_active, n_active);
    /* do not change coeff of x, but change coeff of p itself
        because coeff_p needs to be ortho against original c_x */
    
    
    Eigen::MatrixXd coeff_p = coeff.middleCols(index_x_active,n_active);
    // for(int i=0; i<n_active; ++i){coeff_p(n_max_subspace-n_active+i, i) -= 1;}
    auto c_active = coeff_p.middleRows(index_x_active, n_active);
#ifdef DEBUG_LOBPCG
    std::cout << "coeff = \n" << coeff << std::endl;
    std::cout << "coeff_p = \n" << coeff_p << std::endl;
    // std::cout << "c_active = \n" << c_active << std::endl;
#endif
    c_active -= Eigen::MatrixXd::Identity(n_active, n_active);
#ifdef DEBUG_LOBPCG
    // std::cout << "c_active = \n" << c_active << std::endl;
    std::cout << "coeff_p-I = \n" << coeff_p << std::endl;
#endif
    // ortho coeff p(n_working_space, n_active) against coeff x(n_working_space, n_max_subspace)
    ortho_against_y(n_working_space, n_max_subspace, n_active, coeff_p, coeff);
#ifdef DEBUG_LOBPCG
    // std::cout << "coeff = \n" << coeff << std::endl;
    std::cout << "coeff = \n" << coeff << std::endl;
    std::cout << "coeff_p-ortho = \n" << coeff_p << std::endl;
    std::cout << "coeff_p-ortho overlap u_p'u_p = \n" << coeff_p.transpose() * coeff_p << std::endl;
    std::cout << "coeff_p-ortho overlap u_x'u_p = \n" << coeff.transpose() * coeff_p << std::endl;
#endif
    // -- end of computing the coefficients of p
    /* p(n, n_active) = v(n, n_working_space) * coeff_p(n_working_space, n_active)
       where v(n, n_working_space) contains old [x p w]
    */
#ifdef DEBUG_LOBPCG
    // std::cout << "v size: " << v.rows() << " x " << v.cols() << std::endl;
    // std::cout << "coeff_p size: " << coeff_p.rows() << " x " << coeff_p.cols() << std::endl;
#endif
    // iter#0, v = [x w], n_working_space = 2*n_max_subspace
    // else  v = [x p w], n_working_space = n_max_subspace + 2*n_active
    // p(n, n_active) = x(n, n_working_space) * coeff_p(n_working_space, n_active)
    p = v.leftCols(n_working_space) * coeff_p;
    ap = v.leftCols(n_working_space) * coeff_p;
/*???
how to maintain bv
*/    
    if(solving_generalized) bp = v.leftCols(n_working_space) * coeff_p;
    /* now that we have new p, ap and bp, we shall put them in v */
    /* v(n, n_working_space) [x(n, n_max_subspace) p(n, n_active) w(n, n_active)]*/
    /* update p in v[x p w] */
    v.middleCols(n_max_subspace, n_active) = p;
    av.middleCols(n_max_subspace, n_active) = ap;
    /* update x in v[x p w] */
    v.leftCols(n_max_subspace) = x;
    av.leftCols(n_max_subspace) = ax;
#ifdef DEBUG_LOBPCG
    std::cout << "x, p already updated, w to be done" << std::endl;
    std::cout << "updated v[x p w] = \n" << v << std::endl;
    std::cout << "updated av[x p w] = \n" << av << std::endl;
#endif

    // -- now p, ap, bp contains new values of step k+1
    
    /* 2.6.4 compute the preconditioned residuals W */
    /* only residuals of active x will be taken into W */
    Eigen::MatrixXd r_active = r.rightCols(n_active);
    w = r_active; // resize w to (n, n_active)
    precnd(n, n_active, r_active, w); // w(n, n_active) = t * r(n, n_active)
        // eigen will change w size if needed during precnd

    /* 2.7.1 orthogonalize W to X, P; and then orthonormalize it */
    /* w size(n, n_active), x size(n, n_max_subspace), p size(n, n_active)*/
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
    
    
    // if(solving_generalized) bv.middleCols(n_max_subspace, n_active) = bp;
    /* w */
    
    v.middleCols(index_w, n_active) = w;
    // av.middleCols(index_w, n_active) = aw;

    std::cout << "updated v[x p w_new] = \n" << v << std::endl;

} // end -  main loop for

// --- 3. clean up ---
    /* verbose output */
    tp_end = get_current_time();
    t_total = tp_end - tp_start;
    std::cout << std::setw(35) << "LOBPCG total time: " << t_total.count() << std::endl;
    std::cout << std::setw(35) << "LOBPCG time for avec: " << t_avec.count() << std::endl;
    std::cout << std::setw(35) << "LOBPCG time for ortho: " << t_ortho.count() << std::endl;
    std::cout << std::setw(35) << "LOBPCG time for Rayleigh-Ritz: " << t_solveRR.count() << std::endl;
    return LOBPCG_CONSTANTS::success;
} // end - lobpcg_solve
