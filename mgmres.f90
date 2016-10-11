!  MGMRES is a FORTRAN90 library which applies the restarted Generalized Minimum
!  Residual (GMRES) algorithm to solve a sparse linear system.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    * A(1:NZ_NUM),  the value of the entry;
!    * IA(1:N+1),    row I values occur in entries IA(I) to IA(I+1)-1;
!    * JA(1:NZ_NUM), the column of the entry;
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!    The routine REARRANGE_CR guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!  Last modified:
!
!    28 August 2012
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.

subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )
!*****************************************************************************80
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end

subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
!*****************************************************************************80
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end

subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )
!*****************************************************************************80
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) tl
  integer ( kind = 4 ) ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end

subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )
!*****************************************************************************80
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end

subroutine mult_givens ( c, s, k, g )
!*****************************************************************************80
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer ( kind = 4 ) K, indicates the location of the first
!    vector entry.
!
!    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = 8 ) c
  real ( kind = 8 ) g(1:k+1)
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end

subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel )
!*****************************************************************************80
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(mr+1)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(mr+1)
  real ( kind = 8 ) h(mr+1,mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) l(ia(n+1)+1)
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) s(mr+1)
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) v(n,mr+1);
  logical, parameter :: verbose = .true.
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(mr+1)

  itr_used = 0

  call rearrange_cr ( n, nz_num, ia, ja, a )

  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end

subroutine rearrange_cr ( n, nz_num, ia, ja, a )
!*****************************************************************************80
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
!
!    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) i4temp
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end
