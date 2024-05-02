"""
Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).

References
----------
.. [1] A. V. Knyazev (2001),
       Toward the Optimal Preconditioned Eigensolver: Locally Optimal
       Block Preconditioned Conjugate Gradient Method.
       SIAM Journal on Scientific Computing 23, no. 2,
       pp. 517-541. :doi:`10.1137/S1064827500366124`

.. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov (2007),
       Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX)
       in hypre and PETSc.  :arxiv:`0705.2626`

.. [3] A. V. Knyazev's C and MATLAB implementations:
       https://github.com/lobpcg/blopex
"""

import warnings
import numpy as np
from scipy.linalg import (inv, eigh, cho_factor, cho_solve,
                          cholesky, LinAlgError)
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import issparse

__all__ = ["lobpcg"]


def _report_nonhermitian(M, name):
    """如果M的共轭转置不等于M(使用二者之差的范数判断)，发出warning
    Report if `M` is not a Hermitian matrix given its type.
    """
    from scipy.linalg import norm

    md = M - M.T.conj()
    nmd = norm(md, 1)
    tol = 10 * np.finfo(M.dtype).eps
    tol = max(tol, tol * norm(M, 1))
    if nmd > tol:
        warnings.warn(
              f"Matrix {name} of the type {M.dtype} is not Hermitian: "
              f"condition: {nmd} < {tol} fails.",
              UserWarning, stacklevel=4
         )

def _as2d(ar):
    """  2d: n*m = (100,3) 这样的叫2D
    If the input array is 2D return it, if it is 1D, append a dimension,
    making it a column vector.
    """
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = np.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux

# 将m转换成可调用的函数 矩阵m -> 函数，返回 m @ v
def _makeMatMat(m):
    if m is None:
        return None
    elif callable(m):
        return lambda v: m(v)
    else:
        return lambda v: m @ v


def _matmul_inplace(x, y, verbosityLevel=0):
    """Perform 'np.matmul' in-place if possible.
x = x @ y  结果存在x中
    If some sufficient conditions for inplace matmul are met, do so.
    Otherwise try inplace update and fall back to overwrite if that fails.
    """
    if x.flags["CARRAY"] and x.shape[1] == y.shape[1] and x.dtype == y.dtype:
        # conditions where we can guarantee that inplace updates will work;
        # i.e. x is not a view/slice, x & y have compatible dtypes, and the
        # shape of the result of x @ y matches the shape of x.
        np.matmul(x, y, out=x)
    else:
        # ideally, we'd have an exhaustive list of conditions above when
        # inplace updates are possible; since we don't, we opportunistically
        # try if it works, and fall back to overwriting if necessary
        try:
            np.matmul(x, y, out=x)
        except Exception:
            if verbosityLevel:
                warnings.warn(
                    "Inplace update of x = x @ y failed, "
                    "x needs to be overwritten.",
                    UserWarning, stacklevel=3
                )
            x = x @ y
    return x


def _applyConstraints(blockVectorV, factYBY, blockVectorBY, blockVectorY):
    """Changes blockVectorV in-place."""
    YBV = blockVectorBY.T.conj() @ blockVectorV # YBV = Y*BV
    tmp = cho_solve(factYBY, YBV) # YBY x = YBV  tmp = x  给定A的Cholesky分解，求解方程
    blockVectorV -= blockVectorY @ tmp # V = V - Y @ tmp


def _b_orthonormalize(B, blockVectorV, blockVectorBV=None,
                      verbosityLevel=0):
    """in-place B-orthonormalize the given block vector using Cholesky."""
    if blockVectorBV is None:
        if B is None:
            blockVectorBV = blockVectorV
        else:
            try:
                blockVectorBV = B(blockVectorV)
            except Exception as e:
                if verbosityLevel:
                    warnings.warn(
                        f"Secondary MatMul call failed with error\n"
                        f"{e}\n",
                        UserWarning, stacklevel=3
                    )
                    return None, None, None
            if blockVectorBV.shape != blockVectorV.shape:
                raise ValueError(
                    f"The shape {blockVectorV.shape} "
                    f"of the orthogonalized matrix not preserved\n"
                    f"and changed to {blockVectorBV.shape} "
                    f"after multiplying by the secondary matrix.\n"
                )

    VBV = blockVectorV.T.conj() @ blockVectorBV  # VBV = V^T B V
    try:
        # VBV is a Cholesky factor from now on...
        VBV = cholesky(VBV, overwrite_a=True) # 此刻起VBV变量存放V^T B V的Cholesky分解上三角
        VBV = inv(VBV, overwrite_a=True)  # 求逆 K = (L^T) ^-1
        blockVectorV = _matmul_inplace(  # V = VK 更新后的 V 满足各列 B-正交化
            blockVectorV, VBV,
            verbosityLevel=verbosityLevel
        )
        if B is not None:
            blockVectorBV = _matmul_inplace(
                blockVectorBV, VBV,
                verbosityLevel=verbosityLevel
            )
        return blockVectorV, blockVectorBV, VBV
    except LinAlgError:
        if verbosityLevel:
            warnings.warn(
                "Cholesky has failed.",
                UserWarning, stacklevel=3
            )
        return None, None, None


def _get_indx(_lambda, num, largest):
    """Get `num` indices into `_lambda` depending on `largest` option."""
    ii = np.argsort(_lambda)
    if largest:
        ii = ii[:-num - 1:-1]
    else:
        ii = ii[:num]

    return ii


def _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel):
    if verbosityLevel:
        _report_nonhermitian(gramA, "gramA")
        _report_nonhermitian(gramB, "gramB")


def lobpcg(
    A,  # 对称矩阵，常是稀疏矩阵
    X,  # X 前k个特征向量初猜
    B=None, # B : 可选，广义特征值问题 默认B=I，标准特征值问题
    M=None, # M 约等于A^-1 预条件
    Y=None, # n*sizeY 约束，迭代在Y的B-正交补中进行，要求Y满秩
    tol=None, # tolerance， 停止迭代标准
    maxiter=None, # 最大迭代次数，默认20
    largest=True, # 解最大or最小特征值
    verbosityLevel=0, # 控制输出信息量
    retLambdaHistory=False, # Whether to return eigenvalue history.
    retResidualNormsHistory=False, # Whether to return history of residual norms.
    restartControl=20,
    # 如果与“retresidualnormshhistory”中记录的最小值相比，
    # 残差跳跃“2**restartControl”倍，则迭代重新开始。
    # 默认值为restartControl=20，使得重启很少见。
):
    """Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).

    LOBPCG is a preconditioned eigensolver for large real symmetric and complex
    Hermitian definite generalized eigenproblems.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator, callable object}
        The Hermitian linear operator of the problem, usually given by a
        sparse matrix.  Often called the "stiffness matrix".
    X : ndarray, float32 or float64
        Initial approximation to the ``k`` eigenvectors (non-sparse).
        If `A` has ``shape=(n,n)`` then `X` must have ``shape=(n,k)``.
    B : {sparse matrix, ndarray, LinearOperator, callable object}
        Optional. By default ``B = None``, which is equivalent to identity.
        The right hand side operator in a generalized eigenproblem if present.
        Often called the "mass matrix". Must be Hermitian positive definite.
    M : {sparse matrix, ndarray, LinearOperator, callable object}
        Optional. By default ``M = None``, which is equivalent to identity.
        Preconditioner aiming to accelerate convergence.
    Y : ndarray, float32 or float64, default: None
        An ``n-by-sizeY`` ndarray of constraints with ``sizeY < n``.
        The iterations will be performed in the ``B``-orthogonal complement
        of the column-space of `Y`. `Y` must be full rank if present.
    tol : scalar, optional
        The default is ``tol=n*sqrt(eps)``.
        Solver tolerance for the stopping criterion.
    maxiter : int, default: 20
        Maximum number of iterations.
    largest : bool, default: True
        When True, solve for the largest eigenvalues, otherwise the smallest.
    verbosityLevel : int, optional
        By default ``verbosityLevel=0`` no output.
        Controls the solver standard/screen output.
    retLambdaHistory : bool, default: False
        Whether to return iterative eigenvalue history.
    retResidualNormsHistory : bool, default: False
        Whether to return iterative history of residual norms.
    restartControl : int, optional.
        Iterations restart if the residuals jump ``2**restartControl`` times
        compared to the smallest recorded in ``retResidualNormsHistory``.
        The default is ``restartControl=20``, making the restarts rare for
        backward compatibility.

    Returns
    -------
    lambda : ndarray of the shape ``(k, )``.
        Array of ``k`` approximate eigenvalues.
    v : ndarray of the same shape as ``X.shape``.
        An array of ``k`` approximate eigenvectors.
    lambdaHistory : ndarray, optional.
        The eigenvalue history, if `retLambdaHistory` is ``True``.
    ResidualNormsHistory : ndarray, optional.
        The history of residual norms, if `retResidualNormsHistory`
        is ``True``.

    Notes
    -----
    The iterative loop runs ``maxit=maxiter`` (20 if ``maxit=None``)
    iterations at most and finishes earler if the tolerance is met.
    Breaking backward compatibility with the previous version, LOBPCG
    now returns the block of iterative vectors with the best accuracy rather
    than the last one iterated, as a cure for possible divergence.
    循环最多跑maxit=maxiter（默认20）次，如果提前达到误差限则停止
    此版本返回最高精度而非最后一次迭代的迭代向量组

    If ``X.dtype == np.float32`` and user-provided operations/multiplications
    by `A`, `B`, and `M` all preserve the ``np.float32`` data type,
    all the calculations and the output are in ``np.float32``.
    A B M 和用户提供的操作均为float32类型，则计算和输出也全部采用此精度

    The size of the iteration history output equals to the number of the best
    (limited by `maxit`) iterations plus 3: initial, final, and postprocessing.
    迭代历史记录数量 迭代次数+3(初始，最终，后处理)

    If both `retLambdaHistory` and `retResidualNormsHistory` are ``True``,
    the return tuple has the following format
    ``(lambda, V, lambda history, residual norms history)``.
    返回值

    In the following ``n`` denotes the matrix size and ``k`` the number
    of required eigenvalues (smallest or largest).
    n是矩阵大小（列向量维数），k是需要的特征值数目

    The LOBPCG code internally solves eigenproblems of the size ``3k`` on every
    iteration by calling the dense eigensolver `eigh`, so if ``k`` is not
    small enough compared to ``n``, it makes no sense to call the LOBPCG code.
    Moreover, if one calls the LOBPCG algorithm for ``5k > n``, it would likely
    break internally, so the code calls the standard function `eigh` instead.
    It is not that ``n`` should be large for the LOBPCG to work, but rather the
    ratio ``n / k`` should be large. It you call LOBPCG with ``k=1``
    and ``n=10``, it works though ``n`` is small. The method is intended
    for extremely large ``n / k``.
    每个循环内部调用稠密矩阵特征值求解器eigh求解大小 3k 的特征问题，
    如果 k 和 n 相比不够小，lobpcg方法没有意义
    如果 5k > n 很可能内部中断，直接调用标准的eigh
    要求 n 和 k 的比值 n/k 足够大，本方法用于解决 n/k 非常大的问题

    The convergence speed depends basically on three factors:
收敛速度取决于
    1. Quality of the initial approximations `X` to the seeking eigenvectors.
       Randomly distributed around the origin vectors work well if no better
       choice is known.
       对特征向量初猜 X 的质量。
       如果没有已知的高质量初始值，选择零向量周围的随机分布。

    2. Relative separation of the desired eigenvalues from the rest
       of the eigenvalues. One can vary ``k`` to improve the separation.
       
       待求特征值与其他特征值的相对分离程度。可以试着改变“k”来提高分离程度

    3. Proper preconditioning to shrink the spectral spread.
       For example, a rod vibration test problem (under tests
       directory) is ill-conditioned for large ``n``, so convergence will be
       slow, unless efficient preconditioning is used. For this specific
       problem, a good simple preconditioner function would be a linear solve
       for `A`, which is easy to code since `A` is tridiagonal.
       
       恰当的预条件，减小谱分布范围（？） spectral spread
       对于这个特定的问题，一个好的简单的预条件函数应该是a的线性解（？）， linear solve
       因为a是三对角线的，所以很容易代码实现。

    References
    ----------
    .. [1] A. V. Knyazev (2001),
           Toward the Optimal Preconditioned Eigensolver: Locally Optimal
           Block Preconditioned Conjugate Gradient Method.
           SIAM Journal on Scientific Computing 23, no. 2,
           pp. 517-541. :doi:`10.1137/S1064827500366124`

    .. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov
           (2007), Block Locally Optimal Preconditioned Eigenvalue Xolvers
           (BLOPEX) in hypre and PETSc. :arxiv:`0705.2626`

    .. [3] A. V. Knyazev's C and MATLAB implementations:
           https://github.com/lobpcg/blopex

    Examples
    --------
    Our first example is minimalistic - find the largest eigenvalue of
    a diagonal matrix by solving the non-generalized eigenvalue problem
    ``A x = lambda x`` without constraints or preconditioning.
    第一个例子 标准特征值问题 没有预条件和正交化约束 Ax = lambda x

    >>> import numpy as np
    >>> from scipy.sparse import spdiags
    >>> from scipy.sparse.linalg import LinearOperator, aslinearoperator
    >>> from scipy.sparse.linalg import lobpcg

    The square matrix size is

    >>> n = 100

    and its diagonal entries are 1, ..., 100 defined by
    对角矩阵

    >>> vals = np.arange(1, n + 1).astype(np.int16)

    The first mandatory input parameter in this test is
    the sparse diagonal matrix `A`
    of the eigenvalue problem ``A x = lambda x`` to solve.
    第一个强制参数A

    >>> A = spdiags(vals, 0, n, n)
    >>> A = A.astype(np.int16)
    >>> A.toarray()
    array([[  1,   0,   0, ...,   0,   0,   0],
           [  0,   2,   0, ...,   0,   0,   0],
           [  0,   0,   3, ...,   0,   0,   0],
           ...,
           [  0,   0,   0, ...,  98,   0,   0],
           [  0,   0,   0, ...,   0,  99,   0],
           [  0,   0,   0, ...,   0,   0, 100]], dtype=int16)

    The second mandatory input parameter `X` is a 2D array with the
    row dimension determining the number of requested eigenvalues.
    `X` is an initial guess for targeted eigenvectors.
    `X` must have linearly independent columns.
    If no initial approximations available, randomly oriented vectors
    commonly work best, e.g., with components normally distributed
    around zero or uniformly distributed on the interval [-1 1].
    Setting the initial approximations to dtype ``np.float32``
    forces all iterative values to dtype ``np.float32`` speeding up
    the run while still allowing accurate eigenvalue computations.
    第二个强制参数X，对待求特征向量数目的初猜，列数决定待求解特征向量数，要求X各列线性无关
    如果没有可用的初始值，一般情况下选取随机正态分布的向量即可
    例如选择以0为中心两侧[-1 1]正太或均匀分布的向量
    强制所有元素为float32类型加速计算，并且仍然能获得精确结果

    >>> k = 1
    >>> rng = np.random.default_rng()
    >>> X = rng.normal(size=(n, k))
    >>> X = X.astype(np.float32)

    >>> eigenvalues, _ = lobpcg(A, X, maxiter=60)
    >>> eigenvalues
    array([100.])
    >>> eigenvalues.dtype
    dtype('float32')

    LOBPCG needs only access the matrix product with `A` rather
    then the matrix itself. Since the matrix `A` is diagonal in
    this example, one can write a function of the product
    ``A @ X`` using the diagonal values ``vals`` only, e.g., by
    element-wise multiplication with broadcasting
    对于矩阵A，只要求能够计算 A @ X 乘积，可以通过其他方式给出函数，不一定需要严格矩阵相乘

    >>> A_f = lambda X: vals[:, np.newaxis] * X

    and use the handle ``A_f`` to this callable function as an input

    >>> eigenvalues, _ = lobpcg(A_f, X, maxiter=60)
    >>> eigenvalues
    array([100.])

    The next example illustrates computing 3 smallest eigenvalues of
    the same matrix given by the function handle ``A_f`` with
    constraints and preconditioning.
    复杂的算例，计算3个最小特征值，有约束和预条件

    >>> k = 3
    >>> X = rng.normal(size=(n, k))

    Constraints - an optional input parameter is a 2D array comprising
    of column vectors that the eigenvectors must be orthogonal to
    约束 可选，以矩阵形式给出列向量组，所有特征向量必须与向量组正交

    >>> Y = np.eye(n, 3)

    The preconditioner acts as the inverse of `A` in this example, but
    in the reduced precision ``np.float32`` even though the initial `X`
    and thus all iterates and the output are in full ``np.float64``.
    预条件作为'A'的逆出现，采用低精度（虽然X和全程计算均为64位高精度）

    >>> inv_vals = 1./vals
    >>> inv_vals = inv_vals.astype(np.float32)
    >>> M = lambda X: inv_vals[:, np.newaxis] * X

    Let us now solve the eigenvalue problem for the matrix `A` first
    without preconditioning requesting 80 iterations
    求解，无预条件，80次迭代

    >>> eigenvalues, _ = lobpcg(A_f, X, Y=Y, largest=False, maxiter=80)
    >>> eigenvalues
    array([4., 5., 6.])
    >>> eigenvalues.dtype
    dtype('float64')

    With preconditioning we need only 20 iterations from the same `X`
    有预条件，同样的 X 只需要20次迭代

    >>> eigenvalues, _ = lobpcg(A_f, X, Y=Y, M=M, largest=False, maxiter=20)
    >>> eigenvalues
    array([4., 5., 6.])

    Note that the vectors passed in `Y` are the eigenvectors of the 3
    smallest eigenvalues. The results returned above are orthogonal to those.
    Y给出了3个最小特征值的对应特征向量。结果和它们正交。

    Finally, the primary matrix `A` may be indefinite, e.g., after shifting
    ``vals`` by 50 from 1, ..., 100 to -49, ..., 50, we still can compute
    the 3 smallest or largest eigenvalues.
    A不要求是正定的。

    >>> vals = vals - 50
    >>> X = rng.normal(size=(n, k))
    >>> eigenvalues, _ = lobpcg(A_f, X, largest=False, maxiter=99)
    >>> eigenvalues
    array([-49., -48., -47.])
    >>> eigenvalues, _ = lobpcg(A_f, X, largest=True, maxiter=99)
    >>> eigenvalues
    array([50., 49., 48.])

    """
### 步骤 算法2.2对应的编号
### 初始化 1~3
### 1 给变量分配空间
    blockVectorX = X
    bestblockVectorX = blockVectorX
    blockVectorY = Y
    residualTolerance = tol
    if maxiter is None:
        maxiter = 20

    bestIterationNumber = maxiter

    sizeY = 0
    if blockVectorY is not None:
        if len(blockVectorY.shape) != 2:
            warnings.warn(
                f"Expected rank-2 array for argument Y, instead got "
                f"{len(blockVectorY.shape)}, "
                f"so ignore it and use no constraints.",
                UserWarning, stacklevel=2
            )
            blockVectorY = None
        else:
            sizeY = blockVectorY.shape[1] # Y 是 (n,sizeY)矩阵

    # Block size.
    if blockVectorX is None:
        raise ValueError("The mandatory initial matrix X cannot be None")
    if len(blockVectorX.shape) != 2:
        raise ValueError("expected rank-2 array for argument X")

    n, sizeX = blockVectorX.shape   # X (n,sizeX)

    # Data type of iterates, determined by X, must be inexact
    # inexact 不精确表示的数值类型 如浮点数 https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.inexact
    if not np.issubdtype(blockVectorX.dtype, np.inexact): 
        warnings.warn(
            f"Data type for argument X is {blockVectorX.dtype}, "
            f"which is not inexact, so casted to np.float32.",
            UserWarning, stacklevel=2
        )
        blockVectorX = np.asarray(blockVectorX, dtype=np.float32)
# 初始化history 分配空间
    if retLambdaHistory:
        lambdaHistory = np.zeros((maxiter + 3, sizeX),
                                 dtype=blockVectorX.dtype)
    if retResidualNormsHistory:
        residualNormsHistory = np.zeros((maxiter + 3, sizeX),
                                        dtype=blockVectorX.dtype)
# aux 输出用字符串
    '''
    Solving standard|generalized eigenvalue problem with|without preconditioning
    matrix size %d block size %d  {No constraints}|%d constraints
    '''
    if verbosityLevel:
        aux = "Solving "
        if B is None:
            aux += "standard"
        else:
            aux += "generalized"
        aux += " eigenvalue problem with"
        if M is None:
            aux += "out"
        aux += " preconditioning\n\n"
        aux += "matrix size %d\n" % n
        aux += "block size %d\n\n" % sizeX
        if blockVectorY is None:
            aux += "No constraints\n\n"
        else:
            if sizeY > 1:
                aux += "%d constraints\n\n" % sizeY
            else:
                aux += "%d constraint\n\n" % sizeY
        print(aux)
# 待求解特征向量个数相对向量空间维数来说太大了  5k>n 调用标准的eigh求解
    if (n - sizeY) < (5 * sizeX):
        warnings.warn(
            f"The problem size {n} minus the constraints size {sizeY} "
            f"is too small relative to the block size {sizeX}. "
            f"Using a dense eigensolver instead of LOBPCG iterations."
            f"No output of the history of the iterations.",
            UserWarning, stacklevel=2
        )

        sizeX = min(sizeX, n)

        if blockVectorY is not None:
            raise NotImplementedError(
                "The dense eigensolver does not support constraints."
            )

        # Define the closed range of indices of eigenvalues to return.
        if largest:
            eigvals = (n - sizeX, n - 1)
        else:
            eigvals = (0, sizeX - 1)
# 将算子A和B转换成矩阵传给eigh求解
        try:
            if isinstance(A, LinearOperator):
                A = A(np.eye(n, dtype=int))
            elif callable(A):
                A = A(np.eye(n, dtype=int))
                if A.shape != (n, n):
                    raise ValueError(
                        f"The shape {A.shape} of the primary matrix\n"
                        f"defined by a callable object is wrong.\n"
                    )
            elif issparse(A):
                A = A.toarray() # 返回稀疏矩阵的ndarray表示
            else:
                A = np.asarray(A)
        except Exception as e:
            raise Exception(
                f"Primary MatMul call failed with error\n"
                f"{e}\n")

        if B is not None:
            try:
                if isinstance(B, LinearOperator):
                    B = B(np.eye(n, dtype=int))
                elif callable(B):
                    B = B(np.eye(n, dtype=int))
                    if B.shape != (n, n):
                        raise ValueError(
                            f"The shape {B.shape} of the secondary matrix\n"
                            f"defined by a callable object is wrong.\n"
                        )
                elif issparse(B):
                    B = B.toarray()
                else:
                    B = np.asarray(B)
            except Exception as e:
                raise Exception(
                    f"Secondary MatMul call failed with error\n"
                    f"{e}\n")
# 调用标准的eigh解特征问题
        try:
            vals, vecs = eigh(A,
                              B,
                              subset_by_index=eigvals, # 定义需要求特征值的索引范围，下标从0开始
                              check_finite=False)
            if largest:
                # Reverse order to be compatible with eigs() in 'LM' mode.
                vals = vals[::-1]
                vecs = vecs[:, ::-1]

            return vals, vecs
        except Exception as e:
            raise Exception(
                f"Dense eigensolver failed with error\n"
                f"{e}\n"
            )
# 初始化tolerance
    if (residualTolerance is None) or (residualTolerance <= 0.0):
        residualTolerance = np.sqrt(np.finfo(blockVectorX.dtype).eps) * n
# 转换成函数
    A = _makeMatMat(A)
    B = _makeMatMat(B)
    M = _makeMatMat(M)
# 2. Apply the constraints to X: 关于Y B-正交化
    # Apply constraints to X.
    if blockVectorY is not None:

        if B is not None:
            blockVectorBY = B(blockVectorY)
            if blockVectorBY.shape != blockVectorY.shape:
                raise ValueError(
                    f"The shape {blockVectorY.shape} "
                    f"of the constraint not preserved\n"
                    f"and changed to {blockVectorBY.shape} "
                    f"after multiplying by the secondary matrix.\n"
                )
        else:
            blockVectorBY = blockVectorY
    # 关于Y B-正交化
        # gramYBY is a dense array.
        gramYBY = blockVectorY.T.conj() @ blockVectorBY
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = cho_factor(gramYBY, overwrite_a=True)
        except LinAlgError as e:
            raise ValueError("Linearly dependent constraints") from e

        _applyConstraints(blockVectorX, gramYBY, blockVectorBY, blockVectorY)
# 3. B-orthonormalize X
    ##X列向量组 B-正交化
    # B-orthonormalize X.
    blockVectorX, blockVectorBX, _ = _b_orthonormalize(
        B, blockVectorX, verbosityLevel=verbosityLevel)
    if blockVectorX is None:
        raise ValueError("Linearly dependent initial approximations")
# 4. Compute the initial Ritz vectors: solve the eigenproblem
# 计算初始的Ritz向量组  
    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = A(blockVectorX)
    if blockVectorAX.shape != blockVectorX.shape:
        raise ValueError(
            f"The shape {blockVectorX.shape} "
            f"of the initial approximations not preserved\n"
            f"and changed to {blockVectorAX.shape} "
            f"after multiplying by the primary matrix.\n"
        )

    gramXAX = blockVectorX.T.conj() @ blockVectorAX
    
    # 求解 X^T A X v = lambda v
    # X = X * v
    # AX = AX * v, BX = BX * v

    _lambda, eigBlockVector = eigh(gramXAX, check_finite=False)
    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    if retLambdaHistory:
        lambdaHistory[0, :] = _lambda

    eigBlockVector = np.asarray(eigBlockVector[:, ii])
    blockVectorX = _matmul_inplace(
        blockVectorX, eigBlockVector,
        verbosityLevel=verbosityLevel
    )
    blockVectorAX = _matmul_inplace(
        blockVectorAX, eigBlockVector,
        verbosityLevel=verbosityLevel
    )
    if B is not None:
        blockVectorBX = _matmul_inplace(
            blockVectorBX, eigBlockVector,
            verbosityLevel=verbosityLevel
        )
# 5. Define the index set J of active iterates to be {1, . . . , m}.
# 定义接下来需要求解的列向量集合
# 初始为全部(1~m)
# 当部分向量的残差模长小于tol时，从index set J中排除，此后不需要
    ##
    # Active index set.
    activeMask = np.ones((sizeX,), dtype=bool)
    
# ----------------------------------------------------------------------------
# 主循环
# 6. for k = 0, . . . , MaxIterations:
    ##
    # Main iteration loop.

    blockVectorP = None  # set during iteration
    blockVectorAP = None
    blockVectorBP = None

    smallestResidualNorm = np.abs(np.finfo(blockVectorX.dtype).max)

    iterationNumber = -1
    restart = True
    forcedRestart = False
    explicitGramFlag = False
    while iterationNumber < maxiter:
        iterationNumber += 1
# 7. Compute the residuals: Wj = AXj − BXj ∗ Λj.
        # W = AX - BX*_lambda
        if B is not None:
            aux = blockVectorBX * _lambda[np.newaxis, :]
        else:
            aux = blockVectorX * _lambda[np.newaxis, :]

        blockVectorR = blockVectorAX - aux # 残差 R = AX - BX*_lambda

        aux = np.sum(blockVectorR.conj() * blockVectorR, 0) # 平方和，列和  elementwise product
        residualNorms = np.sqrt(np.abs(aux)) 
        if retResidualNormsHistory:
            residualNormsHistory[iterationNumber, :] = residualNorms
        residualNorm = np.sum(np.abs(residualNorms)) / sizeX # 残差范数
# 残差太大 restart
        if residualNorm < smallestResidualNorm:
            smallestResidualNorm = residualNorm
            bestIterationNumber = iterationNumber
            bestblockVectorX = blockVectorX
        elif residualNorm > 2**restartControl * smallestResidualNorm:
# forcedRestart = True
            forcedRestart = True
            blockVectorAX = A(blockVectorX)
            if blockVectorAX.shape != blockVectorX.shape:
                raise ValueError(
                    f"The shape {blockVectorX.shape} "
                    f"of the restarted iterate not preserved\n"
                    f"and changed to {blockVectorAX.shape} "
                    f"after multiplying by the primary matrix.\n"
                )
            if B is not None:
                blockVectorBX = B(blockVectorX)
                if blockVectorBX.shape != blockVectorX.shape:
                    raise ValueError(
                        f"The shape {blockVectorX.shape} "
                        f"of the restarted iterate not preserved\n"
                        f"and changed to {blockVectorBX.shape} "
                        f"after multiplying by the secondary matrix.\n"
                    )
        '''
8. Exclude from the index set J the indices that correspond to residual
vectors for which the norm has become smaller than the tolerance.
If J then becomes empty, exit loop
        '''
        ii = np.where(residualNorms > residualTolerance, True, False)
        activeMask = activeMask & ii
        currentBlockSize = activeMask.sum()

        if verbosityLevel:
            print(f"iteration {iterationNumber}")
            print(f"current block size: {currentBlockSize}")
            print(f"eigenvalue(s):\n{_lambda}")
            print(f"residual norm(s):\n{residualNorms}")
        # If J then becomes empty, exit loop
        if currentBlockSize == 0:
            break

        activeBlockVectorR = _as2d(blockVectorR[:, activeMask]) # 取当前需要迭代的列

        if iterationNumber > 0:
            activeBlockVectorP = _as2d(blockVectorP[:, activeMask])
            activeBlockVectorAP = _as2d(blockVectorAP[:, activeMask])
            if B is not None:
                activeBlockVectorBP = _as2d(blockVectorBP[:, activeMask])
# 9. Apply the preconditioner T to the residuals: WJ = T ∗ WJ
        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M(activeBlockVectorR) # R = T @ R
#10 Apply the constraints to the preconditioned residuals WJ: activeBlockVectorR
# W对Y B-正交化
        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            _applyConstraints(activeBlockVectorR,
                              gramYBY,
                              blockVectorBY,
                              blockVectorY)
# 11. Compute BWJ and B-orthonormalize WJ: 
# BWJ = B ∗ WJ; R = chol(WJT ∗ BWJ); WJ = WJ ∗ R−1; BWJ = BWJ ∗ R−1.

# ? 似乎是加的
# 预条件残差W 关于X B-正交化

        ##
        # B-orthogonalize the preconditioned residuals to X.
        if B is not None:
            activeBlockVectorR = activeBlockVectorR - (
                blockVectorX @
                (blockVectorBX.T.conj() @ activeBlockVectorR)
            )
        else:
            activeBlockVectorR = activeBlockVectorR - (
                blockVectorX @
                (blockVectorX.T.conj() @ activeBlockVectorR)
            )
# 预条件残差W (矩阵activeBlockVectorR) B-正交化
        ##
        # B-orthonormalize the preconditioned residuals.
        aux = _b_orthonormalize(
            B, activeBlockVectorR, verbosityLevel=verbosityLevel)
        activeBlockVectorR, activeBlockVectorBR, _ = aux
# 12. Compute AWJ: AWJ = A ∗ WJ. 计算AW
        if activeBlockVectorR is None:
            warnings.warn(
                f"Failed at iteration {iterationNumber} with accuracies "
                f"{residualNorms}\n not reaching the requested "
                f"tolerance {residualTolerance}.",
                UserWarning, stacklevel=2
            )
            break
        activeBlockVectorAR = A(activeBlockVectorR)
        '''
13. if k > 0
14. B-orthonormalize PJ: R = chol(PJT ∗ BPJ); PJ = PJ ∗ R−1;
15. Update APJ = APJ ∗ R−1; BPJ = BPJ ∗ R−1.
16. end if
        '''
        if iterationNumber > 0:
            if B is not None:
                aux = _b_orthonormalize(
                    B, activeBlockVectorP, activeBlockVectorBP,
                    verbosityLevel=verbosityLevel
                )
                activeBlockVectorP, activeBlockVectorBP, invR = aux
            else:
                aux = _b_orthonormalize(B, activeBlockVectorP,
                                        verbosityLevel=verbosityLevel)
                activeBlockVectorP, _, invR = aux
            # Function _b_orthonormalize returns None if Cholesky fails
            if activeBlockVectorP is not None:
                activeBlockVectorAP = _matmul_inplace(
                    activeBlockVectorAP, invR,
                    verbosityLevel=verbosityLevel
                )
                restart = forcedRestart
            else:
                restart = True
# 16. end if

#--------------------------------------------------------------------#
# Perform the Rayleigh Ritz Procedure:

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:
# 是否需要显示Gram？ 设置flag
        if activeBlockVectorAR.dtype == "float32":
            myeps = 1
        else:
            myeps = np.sqrt(np.finfo(activeBlockVectorR.dtype).eps)

        if residualNorms.max() > myeps and not explicitGramFlag:
            explicitGramFlag = False
        else:
            # Once explicitGramFlag, forever explicitGramFlag.
            explicitGramFlag = True
#
        # Shared memory assingments to simplify the code 共享内存，如果B=I，BX=X，BR=R，BP=P
        if B is None:
            blockVectorBX = blockVectorX
            activeBlockVectorBR = activeBlockVectorR
            if not restart:
                activeBlockVectorBP = activeBlockVectorP

        # Common submatrices:
        gramXAR = np.dot(blockVectorX.T.conj(), activeBlockVectorAR)
        gramRAR = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAR)

        gramDtype = activeBlockVectorAR.dtype
        if explicitGramFlag:
            gramRAR = (gramRAR + gramRAR.T.conj()) / 2
            gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)
            gramXAX = (gramXAX + gramXAX.T.conj()) / 2
            gramXBX = np.dot(blockVectorX.T.conj(), blockVectorBX)
            gramRBR = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorBR)
            gramXBR = np.dot(blockVectorX.T.conj(), activeBlockVectorBR)
        else:
            gramXAX = np.diag(_lambda).astype(gramDtype)
            gramXBX = np.eye(sizeX, dtype=gramDtype)
            gramRBR = np.eye(currentBlockSize, dtype=gramDtype)
            gramXBR = np.zeros((sizeX, currentBlockSize), dtype=gramDtype)
# 17. if k>0
        if not restart:
            gramXAP = np.dot(blockVectorX.T.conj(), activeBlockVectorAP)
            gramRAP = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAP)
            gramPAP = np.dot(activeBlockVectorP.T.conj(), activeBlockVectorAP)
            gramXBP = np.dot(blockVectorX.T.conj(), activeBlockVectorBP)
            gramRBP = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorBP)
            if explicitGramFlag:
                gramPAP = (gramPAP + gramPAP.T.conj()) / 2
                gramPBP = np.dot(activeBlockVectorP.T.conj(),
                                 activeBlockVectorBP)
            else:
                gramPBP = np.eye(currentBlockSize, dtype=gramDtype)
# 18. gramA
            gramA = np.block(
                [
                    [gramXAX, gramXAR, gramXAP],
                    [gramXAR.T.conj(), gramRAR, gramRAP],
                    [gramXAP.T.conj(), gramRAP.T.conj(), gramPAP],
                ]
            )
# 19. gramB
            gramB = np.block(
                [
                    [gramXBX, gramXBR, gramXBP],
                    [gramXBR.T.conj(), gramRBR, gramRBP],
                    [gramXBP.T.conj(), gramRBP.T.conj(), gramPBP],
                ]
            )

            _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel)
            '''
24. Solve the generalized eigenvalue problem:
gramA ∗ C = gramB ∗ C ∗ Λ, where the first m eigenvalues in
increasing order are in the diagonal matrix Λ and the
gramB-orthonormalized eigenvectors are the columns of C.
            '''
            try:
                _lambda, eigBlockVector = eigh(gramA,
                                               gramB,
                                               check_finite=False)
            except LinAlgError as e:
                # raise ValueError("eigh failed in lobpcg iterations") from e
                if verbosityLevel:
                    warnings.warn(
                        f"eigh failed at iteration {iterationNumber} \n"
                        f"with error {e} causing a restart.\n",
                        UserWarning, stacklevel=2
                    )
                # try again after dropping the direction vectors P from RR
                restart = True
        '''
20. else (if k>0 else)
21. gramA
22. gramB
        '''
        if restart:
            gramA = np.block([[gramXAX, gramXAR], [gramXAR.T.conj(), gramRAR]])
            gramB = np.block([[gramXBX, gramXBR], [gramXBR.T.conj(), gramRBR]])

            _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel)
            '''
24. Solve the generalized eigenvalue problem:
gramA ∗ C = gramB ∗ C ∗ Λ, where the first m eigenvalues in
increasing order are in the diagonal matrix Λ and the
gramB-orthonormalized eigenvectors are the columns of C.
            '''
            try:
                _lambda, eigBlockVector = eigh(gramA,
                                               gramB,
                                               check_finite=False)
            except LinAlgError as e:
                # raise ValueError("eigh failed in lobpcg iterations") from e
                warnings.warn(
                    f"eigh failed at iteration {iterationNumber} with error\n"
                    f"{e}\n",
                    UserWarning, stacklevel=2
                )
                break

        ii = _get_indx(_lambda, sizeX, largest)
        _lambda = _lambda[ii]
        eigBlockVector = eigBlockVector[:, ii]
        if retLambdaHistory:
            lambdaHistory[iterationNumber + 1, :] = _lambda
# Compute Ritz vectors:

        # Compute Ritz vectors.
        if B is not None:
# 25. if k > 0
            if not restart:
# 26. Partition C = [Cx Cw Cp]  分解成X R（残差） P
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:
                                                 sizeX + currentBlockSize]
                eigBlockVectorP = eigBlockVector[sizeX + currentBlockSize:]
# 27. Compute P
                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                pp += np.dot(activeBlockVectorP, eigBlockVectorP)
# 27. Compute AP
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                app += np.dot(activeBlockVectorAP, eigBlockVectorP)
# 27. Compute BP
                bpp = np.dot(activeBlockVectorBR, eigBlockVectorR)
                bpp += np.dot(activeBlockVectorBP, eigBlockVectorP)
# 29. else (if k>0 else)
            else:
# 30. Partition C  = [Cx Cw]
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:]
# 31. Compute P, AP, BP
                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                bpp = np.dot(activeBlockVectorBR, eigBlockVectorR)
# 28. 32. 更新 X, AX, BX
            blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
            blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app
            blockVectorBX = np.dot(blockVectorBX, eigBlockVectorX) + bpp

            blockVectorP, blockVectorAP, blockVectorBP = pp, app, bpp
            
# 33. end if

# B=I
# 25. if k > 0
# ...
# 33. end if
        else:
            if not restart:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:
                                                 sizeX + currentBlockSize]
                eigBlockVectorP = eigBlockVector[sizeX + currentBlockSize:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                pp += np.dot(activeBlockVectorP, eigBlockVectorP)

                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                app += np.dot(activeBlockVectorAP, eigBlockVectorP)
            else:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)

            blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
            blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app

            blockVectorP, blockVectorAP = pp, app
# 33. end if
    if B is not None:
        aux = blockVectorBX * _lambda[np.newaxis, :]
    else:
        aux = blockVectorX * _lambda[np.newaxis, :]


#-----------------------------------------------------------------------#
# while iterationNumber < maxiter
# ends here
#-----------------------------------------------------------------------#

#-----------------------------------------------------------------------#
# lobpcg postprocessing starts
#-----------------------------------------------------------------------#

# 7. Compute the residuals: Wj = AXj − BXj ∗ Λj.

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))
    # Use old lambda in case of early loop exit.
    if retLambdaHistory:
        lambdaHistory[iterationNumber + 1, :] = _lambda
    if retResidualNormsHistory:
        residualNormsHistory[iterationNumber + 1, :] = residualNorms
    residualNorm = np.sum(np.abs(residualNorms)) / sizeX
    if residualNorm < smallestResidualNorm:
        smallestResidualNorm = residualNorm
        bestIterationNumber = iterationNumber + 1
        bestblockVectorX = blockVectorX
# reach iterationNumber not reaching the requested tolerance
# 迭代次数达到后未满足给定精度，选取最佳结果
    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited at iteration {iterationNumber} with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.\n"
            f"Use iteration {bestIterationNumber} instead with accuracy \n"
            f"{smallestResidualNorm}.\n",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final iterative eigenvalue(s):\n{_lambda}")
        print(f"Final iterative residual norm(s):\n{residualNorms}")
# 2. Apply the constraints to X:
    blockVectorX = bestblockVectorX
    # Making eigenvectors "exactly" satisfy the blockVectorY constrains
    if blockVectorY is not None:
        _applyConstraints(blockVectorX,
                          gramYBY,
                          blockVectorBY,
                          blockVectorY)
# 最后返回值包括：
# _lambda, blockVectorX
# 最后一轮Rayleigh-Ritz
    # Making eigenvectors "exactly" othonormalized by final "exact" RR
    blockVectorAX = A(blockVectorX)
    if blockVectorAX.shape != blockVectorX.shape:
        raise ValueError(
            f"The shape {blockVectorX.shape} "
            f"of the postprocessing iterate not preserved\n"
            f"and changed to {blockVectorAX.shape} "
            f"after multiplying by the primary matrix.\n"
        )
    gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)

    blockVectorBX = blockVectorX
    if B is not None:
        blockVectorBX = B(blockVectorX)
        if blockVectorBX.shape != blockVectorX.shape:
            raise ValueError(
                f"The shape {blockVectorX.shape} "
                f"of the postprocessing iterate not preserved\n"
                f"and changed to {blockVectorBX.shape} "
                f"after multiplying by the secondary matrix.\n"
            )

    gramXBX = np.dot(blockVectorX.T.conj(), blockVectorBX)
    _handle_gramA_gramB_verbosity(gramXAX, gramXBX, verbosityLevel)
    gramXAX = (gramXAX + gramXAX.T.conj()) / 2
    gramXBX = (gramXBX + gramXBX.T.conj()) / 2
    try:
        _lambda, eigBlockVector = eigh(gramXAX,
                                       gramXBX,
                                       check_finite=False)
    except LinAlgError as e:
        raise ValueError("eigh has failed in lobpcg postprocessing") from e

    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    eigBlockVector = np.asarray(eigBlockVector[:, ii])

    blockVectorX = np.dot(blockVectorX, eigBlockVector)
    blockVectorAX = np.dot(blockVectorAX, eigBlockVector)

    if B is not None:
        blockVectorBX = np.dot(blockVectorBX, eigBlockVector)
        aux = blockVectorBX * _lambda[np.newaxis, :]
    else:
        aux = blockVectorX * _lambda[np.newaxis, :]

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))
# log
    if retLambdaHistory:
        lambdaHistory[bestIterationNumber + 1, :] = _lambda
    if retResidualNormsHistory:
        residualNormsHistory[bestIterationNumber + 1, :] = residualNorms

    if retLambdaHistory:
        lambdaHistory = lambdaHistory[
            : bestIterationNumber + 2, :]
    if retResidualNormsHistory:
        residualNormsHistory = residualNormsHistory[
            : bestIterationNumber + 2, :]
# 后处理过程中eigh失败报错
    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited postprocessing with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final postprocessing eigenvalue(s):\n{_lambda}")
        print(f"Final residual norm(s):\n{residualNorms}")
# 返回值，输出log处理
    if retLambdaHistory:
        lambdaHistory = np.vsplit(lambdaHistory, np.shape(lambdaHistory)[0])
        lambdaHistory = [np.squeeze(i) for i in lambdaHistory]
    if retResidualNormsHistory:
        residualNormsHistory = np.vsplit(residualNormsHistory,
                                         np.shape(residualNormsHistory)[0])
        residualNormsHistory = [np.squeeze(i) for i in residualNormsHistory]

    if retLambdaHistory:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, lambdaHistory, residualNormsHistory
        else:
            return _lambda, blockVectorX, lambdaHistory
    else:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, residualNormsHistory
        else:
            return _lambda, blockVectorX
