import numpy as np
from scipy.linalg import lapack,blas, qr, solve_triangular

def ortho(n, m, u):
    """
    Orthogonalize m vectors of length n using the QR decomposition with dgeqrf and dtrsm.
    
    Parameters:
    - do_other: bool, if True also orthogonalize the matrix w
    - n: int, number of rows (length of vectors)
    - m: int, number of columns (number of vectors)
    - u: numpy.ndarray with shape (n, m), the matrix to be orthogonalized in-place
    - w: numpy.ndarray with shape (n, m), a second matrix to be orthogonalized if do_other is True
    """
    
    v = u @ np.eye(m)
    
    # Compute the QR decomposition of u using dgeqrf
    u = lapack.dgeqrf(u) #, overwrite_a=True)

    # Now u_copy contains the R matrix in the upper triangle and the Householder
    # vectors in the lower triangle. tau contains the scalar factors for the Householder vectors.

    # Solve the triangular system U^T * R = V^T for R using dtrsm
    # We are interested in U * R = V, so we solve for R
    # Since dtrsm solves AX = B or XA = B, we use the right side (R) and
    # transpose 'T' to apply the transformation from the right.
    q = blas.dtrsm(1.0, u, v)

    # Update u with the result of the triangular solve
    u[:] = v


# Example usage:
n, m = 4, 3
u = np.random.rand(n, m)


ortho(n, m, u)

print("Orthonormalized matrix u:\n", u)