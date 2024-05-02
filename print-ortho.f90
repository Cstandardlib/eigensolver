! orthogonalization routines:
! 正交化例程：
! ==========================
!
  subroutine ortho(do_other,n,m,u,w)
    implicit none
!
!   orthogonalize m vectors of length n using the QR decomposition.
!   使用QR分解对长度为n的m个向量进行正交化。
!   this is done by computing U = QR and then by solving the upper
!   triangular system U(ortho)R = U.
!   这是通过计算U = QR，然后通过解上三角系统U(ortho)R = U来完成的。
!
!   using this strategy allows to apply the same linear transformation
!   that orthogonalizes U to a second set of vectors that usually contain
!   the product AU, where A is some matrix that has already been applied
!   to U.
!   使用这种策略允许将正交化U的相同线性变换应用于通常包含AU乘积的第二组向量，其中A是已经应用于U的某个矩阵。
!   this is useful when U and AU are built together without explicitly
!   performing the matrix vector multiplication.
!   当U和AU一起构建而没有显式执行矩阵向量乘法时，这是有用的。
!
!   arguments:
!   参数：
!
    logical,                     intent(in)    :: do_other
    integer,                     intent(in)    :: n, m
    real(dp),  dimension(n,m),   intent(inout) :: u, w
!
!   local scratch
!   局部暂存：
!
    real(dp), allocatable :: v(:,:)
!
!   external functions:
!   外部函数：
!
    external dgeqrf, dtrsm
!
    allocate (v(n,m))
    v = u
    call dgeqrf(n,m,u,n,tau,work,lwork,info) ! now u: upper contains R, lower contains Householders
    !
    ! U * R(upper) = v, result stores on v
    !
    call dtrsm('r','u','n','n',n,m,one,u,n,v,n)
!
    if (do_other) call dtrsm('r','u','n','n',n,m,one,u,n,w,n)
!
    u = v
!
    deallocate (v)
    return
  end subroutine ortho
!
  subroutine b_ortho(n,m,u,bu)
    implicit none
!
!   b-orthogonalize m vectors of length n using the Cholesky decomposition
!   of the overlap matrix.
!   使用重叠矩阵的Cholesky分解对长度为n的m个向量进行b-正交化。
!   this is in principle not a good idea, as the u'bu matrix can be very
!   ill-conditioned, independent of how bas is b, and only works if x is
!   already orthonormal.
!   原则上这不是一个好主意，因为u'bu矩阵可能病态，与b的好坏无关，只有在x已经正交时才有效。
!
!   arguments:
!   参数：
!
    integer,                     intent(in)    :: n, m
    real(dp),  dimension(n,m),   intent(inout) :: u, bu
!
!   local variables
!   局部变量
!
    integer               :: info, i, j
    real(dp), allocatable :: metric(:,:), sigma(:), u_svd(:,:), vt_svd(:,:), &
                             temp(:,:)
    real(dp), parameter   :: tol_svd = 1.0e-5_dp
    logical,  parameter   :: use_svd = .true.
!
!   external functions:
!   外部函数：
!
    external dpotrf, dtrsm, dgemm
!
    allocate (metric(m,m))
    !
    ! b**(-1/2)
    !
    call dgemm('t','n',m,m,n,one,u,n,bu,n,zero,metric,m)
    if (use_svd) then
      allocate (sigma(m), u_svd(m,m), vt_svd(m,m), temp(n,m))
      call dgesvd('a','a',m,m,metric,m,sigma,u_svd,m,vt_svd,m,work,lwork,info)
      !
      ! sigma**(-1/2)
      !
      do i = 1, m
        if (sigma(i) .gt. tol_svd) then
          sigma(i) = 1/sqrt(sigma(i))
        else
          sigma(i) = zero
        end if
      end do
      !
      ! compute metric ** (-1/2). first, compute sigma ** (-1/2) vt
      !
      metric = zero
      do i = 1, m
        do j = 1, m
          metric(j,i) = metric(j,i) + sigma(j)*vt_svd(j,i)
        end do
      end do
      !
      ! now, multiply for u:
      !
      vt_svd = metric
      call dgemm('n','n',m,m,m,one,u_svd,m,vt_svd,m,zero,metric,m)
      !
      ! metric contains s ** (-1/2), and projects out directions corresponding
      ! to pathological singular values.
      ! orthogonalize u and bu:
      !
      call dgemm('n','n',n,m,m,one,u,n,metric,m,zero,temp,n)
      u = temp
      call dgemm('n','n',n,m,m,one,bu,n,metric,m,zero,temp,n)
      bu = temp
      deallocate (sigma, u_svd, vt_svd, temp)
    else
      !
      ! compute the Cholesky factorization of the metric.
      !
      call dpotrf('l',m,metric,m,info)
      !
      ! get u * l^-T and bu * l^-T
      !
      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,bu,n)
    end if
    deallocate (metric)
    return
  end subroutine b_ortho
!
  subroutine ortho_cd(do_other,n,m,u,w,ok)
    implicit none
!
!   orthogonalize m vectors of length n using the Cholesky factorization
!   of their overlap.
!   使用它们重叠的Cholesky分解对长度为n的m个向量进行正交化。
!   this is done by metric = U^t U and then by computing its Cholesky
!   decompositoin metric = L L^t.
!   这是通过metric = U^t U完成的，然后计算其Cholesky分解metric = L L^t。
!   The orthogonal vectors are obtained then by solving the triangular linear system
!   通过解三角线性系统获得正交向量
!     U(ortho)L^T = U
!
!   using this strategy allows to apply the same linear transformation
!   that orthogonalizes U to a second set of vectors that usually contain
!   the product AU, where A is some matrix that has already been applied
!   to U.
!   使用这种策略允许将正交化U的相同线性变换应用于通常包含AU乘积的第二组向量，其中A是已经应用于U的某个矩阵。
!   this is useful when U and AU are built together without explicitly
!   performing the matrix vector multiplication.
!   当U和AU一起构建而没有显式执行矩阵向量乘法时，这是有用的。
!
!   as Cholesky decomposition is not the most stable way of orthogonalizing
!   a set of vectors, the orthogonalization is refined iteratively.
!   由于Cholesky分解不是正交化一组向量最稳定的方法，因此正交化是迭代细化的。
!   still, this routine can fail.
!   尽管如此，这个例程仍然可能失败。
!   a logical flag is then set to false, so that the calling program can
!   call a more robust orthogonalization routine without aborting.
!   然后，一个逻辑标志被设置为false，这样调用程序就可以调用一个更健壮的正交化例程而不中断。
!
!   arguments:
!   参数：
!
    logical,                   intent(in)    :: do_other
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(inout) :: u, w
    logical,                   intent(inout) :: ok
!
!   local variables
!   局部变量
!
    integer               :: it, it_micro
    real(dp)              :: metric_norm, dnrm2, alpha, unorm, shift
    logical               :: macro_done, micro_done
    integer, parameter    :: maxit = 10
!
!   local scratch
!   局部暂存：
!
    real(dp), allocatable :: metric(:,:), msave(:,:)
!
!   external functions:
!   外部函数：
!
    external              :: dgemm, dgeqrf, dtrsm, dnrm2
!
!   get memory for the metric.
!   为metric分配内存。
!
    allocate (metric(m,m), msave(m,m))
    macro_done = .false.
!
!   assemble the metric
!   组装metric
!
    it = 0
    do while(.not. macro_done)
      it = it + 1
      if (it .gt. maxit) then
!
!       ortho_cd failed. return with an error message
!       ortho_cd失败。返回带有错误消息
!
        ok = .false.
        write(6,100) ' maximum number of iterations reached.'
        return
      end if
      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
      msave = metric
!
!   compute the Cholesky factorization of the metric.
!   计算metric的Cholesky分解。
!
      call dpotrf('l',m,metric,m,info)
!
!     if dpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!     如果dpotrf失败，尝试在对metric的对角线进行水平移动后再次尝试。
!
      if (info.ne.0) then
!
        alpha      = 100.0_dp
        unorm      = dnrm2(n*m,u,1)
        it_micro   = 0
        micro_done = .false.
!
!       add larger and larger shifts to the diagonal until dpotrf manages to factorize it.
!       在dpotrf能够分解它之前，向对角线添加越来越大的偏移。
!
        do while (.not. micro_done)
          it_micro = it_micro + 1
          if (it_micro.gt.maxit) then
!
!           something went very wrong. return with an error status, the orthogonalization
!           will be carried out using a different algorithm.
!           出了一些问题。返回错误状态，正交化将使用不同的算法执行。
!
            ok = .false.
            write(6,100) ' maximum number of iterations for factorization reached.'
            return
          end if
!
          shift = max(epsilon(one)*alpha*unorm,tol_ortho)
          metric = msave
          call diag_shift(m,shift,metric)
          call dpotrf('l',m,metric,m,info)
          alpha = alpha * 10.0_dp
          micro_done = info.eq.0
        end do
!
      end if
!
!     orthogonalize u by computing a solution to u(ortho) l^t = u
!     如果需要，将相同的变换应用于w。
!
      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
      if (do_other) call dtrsm('r','l','t','n',n,m,one,metric,m,w,n)
!
!     check that the vectors in v are really orthogonal
!     检查v中的向量是否真的正交
!
      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
      metric_norm = abs(dnrm2(m*m,metric,1) - sqrt(dble(m)))
      macro_done = metric_norm .lt. tol_ortho
    end do
!
    100 format(t3,'ortho_cd failed with the following error:',a)
!
    ok = .true.
!
    deallocate (metric)
    return
  end subroutine ortho_cd
!
  subroutine ortho_vs_x(do_other,n,m,k,x,u,ax,au)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x
!   is assumed to be orthogonal, orthogonalize u against x.
!   如果需要，正交化au到ax，其中ax和au是矩阵a对x和u的应用结果。
!   furthermore, orthonormalize u and, if required, apply the
!   same linear transformation to au.
!   此外，对u进行正交归一化，如果需要，将相同的线性变换应用于au。
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold.
!   本例程执行u与x的正交化及其随后的u正交归一化，迭代进行，直到x与正交化u之间的重叠小于一个（紧密的）阈值。
!
!   arguments:
!   参数：
!
    logical,                   intent(in)    :: do_other
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x, ax
    real(dp),  dimension(n,k), intent(inout) :: u, au
!
!   local variables:
!   局部变量：
!
    logical                :: done, ok
    integer                :: it
    real(dp)               :: xu_norm
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   外部函数：
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
    !
    integer, parameter     :: maxit = 10
    logical, parameter     :: useqr = .false.
!
!   allocate space for the overlap between x and u.
!   为x和u之间的重叠分配空间。
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it   = 0
!
!   start with an initial orthogonalization to improve conditioning.
!   以初始正交化开始以改善条件。
!
    if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
    if (.not. ok .or. useqr) call ortho(do_other,n,k,u,au)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!   迭代地对u进行正交化，然后对u进行正交归一化。
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^t u)
!     u = u - (x^t u) x
!
      call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
      if (do_other) call dgemm('n','n',n,k,m,-one,ax,n,xu,m,one,au,n)
!
!     now, orthonormalize u.
!     现在，对u进行正交归一化。
!
      if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
      if (.not. ok .or. useqr) call ortho(do_other,n,k,u,au)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!     计算正交归一化的u和x之间的重叠，并决定正交化过程是否收敛。
!
      call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
      xu_norm = dnrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!     如果事情真的出了问题，中止。
!
      if (it.gt.maxit) stop ' catastrophic failure of ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine ortho_vs_x
!
  subroutine b_ortho_vs_x(n,m,k,x,bx,u)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x
!   is assumed to be orthogonal, b-orthogonalize u against x.
!   如果需要，对u进行正交归一化。
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold.
!   本例程执行u与x的正交化及其随后的u正交归一化，迭代进行，直到x与正交化u之间的重叠小于一个（紧密的）阈值。
!
!   arguments:
!   参数：
!
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x, bx
    real(dp),  dimension(n,k), intent(inout) :: u
!
!   local variables:
!   局部变量：
!
    logical                :: done, ok
    integer                :: it
    real(dp)               :: xu_norm, xx(1)
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   外部函数：
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
    !
    integer, parameter     :: maxit = 10
    logical, parameter     :: useqr = .true.
!
!   allocate space for the overlap between x and u.
!   为x和u之间的重叠分配空间。
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it   = 0
!
!   start with an initial orthogonalization to improve conditioning.
!   以初始正交化开始以改善条件。
!
    if (.not. useqr) call ortho_cd(.false.,n,k,u,xx,ok)
    if (.not. ok .or. useqr) call ortho(.false.,n,k,u,xx)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!   迭代地对u进行正交化，然后对u进行正交归一化。
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (bx^t u)
!     u = u - (bx^t u) x
!
      call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
!     now, orthonormalize u.
!     现在，对u进行正交归一化。
!
      if (.not. useqr) call ortho_cd(.false.,n,k,u,xx,ok)
      if (.not. ok .or. useqr) call ortho(.false.,n,k,u,xx)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!     计算正交归一化的u和x之间的重叠，并决定正交化过程是否收敛。
!
      call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
      xu_norm = dnrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!     如果事情真的出了问题，中止。
!
      if (it.gt.maxit) stop ' catastrophic failure of b_ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine b_ortho_vs_x
!
! utilities
! 实用工具
! =========
!
  subroutine diag_shift(n,shift,a)
    implicit none
!  
!   add shift to the diagonal elements of the matrix a
!   向矩阵a的对角线元素添加偏移
!  
    integer,                  intent(in)    :: n
    real(dp),                 intent(in)    :: shift
    real(dp), dimension(n,n), intent(inout) :: a
!  
    integer :: i
!  
    do i = 1, n
      a(i,i) = a(i,i) + shift
    end do
!  
    return
  end subroutine diag_shift
!
  subroutine get_coeffs(len_a,len_u,n_max,n_act,a_red,u_x,u_p)
    implicit none
!
!   given the eigenvectors of the reduced matrix in a_red, extract
!   the expansion coefficients for x_new (u_x) and assemble the
!   ones for p_new in u_p.
!   给定简化矩阵a_red的特征向量，提取x_new（u_x）的展开系数，并为p_new在u_p中组装它们。
!   the coefficients u_p are computed as the difference between the
!   coefficients for x_new and x_old, and only the columns associated
!   with active eigenvectors are considered.
!   u_p的系数计算为x_new和x_old的系数之间的差，并且只考虑与活动特征向量相关的列。
!   u_p is then orthogonalized to u_x: this not only guarantees that
!   the p_new vectors will be orthogonal to x_new, but also allows one
!   to reuse the ax, aw, and ap vectors to compute ap_new, without
!   loosing numerical precision.
!   然后对u_p进行正交化以对u_x：这不仅保证了p_new向量将与x_new正交，而且还允许重新使用ax、aw和ap向量来计算ap_new，而不会失去数值精度。
!
    integer,                          intent(in)    :: len_a, len_u, n_max, n_act
    real(dp), dimension(len_a,len_a), intent(in)    :: a_red
    real(dp), dimension(len_u,n_max), intent(inout) :: u_x
    real(dp), dimension(len_u,n_act), intent(inout) :: u_p
!
    integer               :: ind_x, off_x, i_eig
    real(dp)              :: xx(1)
!
    off_x = n_max - n_act
    ind_x = off_x + 1
!
    u_x(1:len_u,1:n_max) = a_red(1:len_u,1:n_max)
!
!   u_p = u_x for the active vectors only
!   u_p = 仅针对活动向量的u_x
!
    u_p = u_x(:,ind_x:n_max)
!
!   remove the coefficients for x from u_p
!   从u_p中移除x的系数
!
    do i_eig = 1, n_act
      u_p(off_x + i_eig,i_eig) = u_p(off_x + i_eig,i_eig) - one
    end do
!
!   orthogonalize:
!   正交化：
!
    call ortho_vs_x(.false.,len_u,n_max,n_act,u_x,u_p,xx,xx)
!
!   all done.
!   全部完成。
!
    return
!
  end subroutine get_coeffs