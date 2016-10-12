! -*- mode: fortran; -*-
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
           module mgmres
           implicit none
           private
           public pmgmres_ilu_cr
           integer, parameter :: si = selected_int_kind (9)      ! short int
           integer, parameter :: dp = selected_real_kind (15, 307)  ! double
           contains

           subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )
!*****************************************************************************80
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Parameters:
!
!    N, the order of the system.
!
!    NZ_NUM, the number of nonzeros.
!
!    IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    A(NZ_NUM), the matrix values.
!
!    X(N), the vector to be multiplied by A.
!
!    W(N), the value of A*X.
!
      implicit none

      integer (kind=si), intent (in)  :: n
      integer (kind=si), intent (in)  :: nz_num
      integer (kind=si), intent (in)  :: ia(n+1)
      integer (kind=si), intent (in)  :: ja(nz_num)
      real    (kind=dp), intent (in)  :: a(nz_num)
      real    (kind=dp), intent (in)  :: x(n)
      real    (kind=dp), intent (out) :: w(n)

      integer (kind=si) :: i
      integer (kind=si) :: k
      integer (kind=si) :: k1
      integer (kind=si) :: k2

      w(1:n) = 0.0D+00

      do i = 1, n
        k1 = ia(i)
        k2 = ia(i+1) - 1
        w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
      end do

      return
      end subroutine ax_cr

      subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
!*****************************************************************************80
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Parameters:
!
!    Input N, the order of the system.
!
!    Input NZ_NUM, the number of nonzeros.
!
!    IA(N+1), JA(NZ_NUM), the row and column indices of the matrix values.
!    The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    UA(N), the index of the diagonal element of each row.
!
      implicit none

      integer (kind=si), intent (in)  :: n
      integer (kind=si), intent (in)  :: nz_num
      integer (kind=si), intent (in)  :: ia(n+1)
      integer (kind=si), intent (in)  :: ja(nz_num)
      integer (kind=si), intent (out) :: ua(n)

      integer (kind=si) :: i
      integer (kind=si) :: k

      ua(1:n) = -1

      do i = 1, n
        do k = ia(i), ia(i+1) - 1
          if ( ja(k) == i ) then
            ua(i) = k
          end if
        end do
      end do

      return
      end subroutine diagonal_pointer_cr

      subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )
!*****************************************************************************80
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Parameters:
!
!    N, the order of the system.
!
!    NZ_NUM, the number of nonzeros.
!
!    IA(N+1), JA(NZ_NUM), the row and column indices of the matrix values.
!    The row vector has been compressed.
!
!    A(NZ_NUM), the matrix values.
!
!    UA(N), the index of the diagonal element of each row.
!
!    L(NZ_NUM), the ILU factorization of A.
!
      implicit none

      integer (kind=si), intent (in)    :: n
      integer (kind=si), intent (in)    :: nz_num
      integer (kind=si), intent (in)    :: ia(n+1)
      integer (kind=si), intent (in)    :: ja(nz_num)
      real    (kind=dp), intent (in)    :: a(nz_num)
      integer (kind=si), intent (inout) :: ua(n)
      real    (kind=dp), intent (out)   :: l(nz_num)

      integer (kind=si) :: i
      integer (kind=si) :: iw(n)
      integer (kind=si) :: j
      integer (kind=si) :: jj
      integer (kind=si) :: jrow
      integer (kind=si) :: jw
      integer (kind=si) :: k
      real    (kind=dp) :: tl
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
      end subroutine ilu_cr

      subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )
!*****************************************************************************80
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Parameters:
!
!    Input N, the order of the system.
!
!    Input NZ_NUM, the number of nonzeros.
!
!    Input IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input L(NZ_NUM), the matrix values.
!
!    Input UA(N), the index of the diagonal element
!    of each row.
!
!    Input R(N), the right hand side.
!
!    Output Z(N), the solution of the system M * Z = R.
!
      implicit none

      integer (kind=si), intent (in)  :: n
      integer (kind=si), intent (in)  :: nz_num
      integer (kind=si), intent (in)  :: ia(n+1)
      integer (kind=si), intent (in)  :: ja(nz_num)
      real    (kind=dp), intent (in)  :: l(nz_num)
      integer (kind=si), intent (in)  :: ua(n)
      real    (kind=dp), intent (in)  :: r(n)
      real    (kind=dp), intent (out) :: z(n)

      integer (kind=si) :: i
      integer (kind=si) :: j
      real    (kind=dp) :: w(n)
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
      end subroutine lus_cr

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
!    C, S, the cosine and sine of a Givens
!    rotation.
!
!    K, indicates the location of the first
!    vector entry.
!
!    G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
      implicit none

      real    (kind=dp), intent (in)  :: c
      real    (kind=dp), intent (in)  :: s
      integer (kind=si), intent (in)  :: k
      real    (kind=dp), intent (out) :: g(1:k+1)

      real    (kind=dp) :: g1
      real    (kind=dp) :: g2

      g1 = c * g(k) - s * g(k+1)
      g2 = s * g(k) + c * g(k+1)

      g(k)   = g1
      g(k+1) = g2

      return
      end subroutine mult_givens

      subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs,          &
     &  itr_max, mr, tol_abs, tol_rel )
!*****************************************************************************80
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Parameters:
!
!    N, the order of the linear system.
!
!    NZ_NUM, the number of nonzero matrix values.
!
!    IA(N+1), JA(NZ_NUM), the row and column indices of the matrix values.
!    The row vector has been compressed.
!
!    A(NZ_NUM), the matrix values.
!
!    X(N); on input, an approximation to the solution.
!    On output, an improved approximation.
!
!    RHS(N), the right hand side of the linear system.
!
!    ITR_MAX, the maximum number of (outer) iterations to take.
!
!    MR, the maximum number of (inner) iterations to take.
!    MR must be less than N.
!
!    TOL_ABS, an absolute tolerance applied to the current residual.
!
!    TOL_REL, a relative tolerance comparing the current residual
!    to the initial residual.
!
      implicit none

      integer (kind=si), intent (in)    :: n
      integer (kind=si), intent (in)    :: nz_num
      integer (kind=si), intent (in)    :: ia(n+1)
      integer (kind=si), intent (inout) :: ja(nz_num)
      real    (kind=dp), intent (inout) :: a(nz_num)
      real    (kind=dp), intent (inout) :: x(n)
      real    (kind=dp), intent (in)    :: rhs(n)
      integer (kind=si), intent (in)    :: itr_max
      integer (kind=si), intent (in)    :: mr
      real    (kind=dp), intent (in)    :: tol_abs
      real    (kind=dp), intent (in)    :: tol_rel

      real    (kind=dp) :: av
      real    (kind=dp) :: c(mr+1)
      real    (kind=dp), parameter :: delta = 1.0D-03
      real    (kind=dp) :: g(mr+1)
      real    (kind=dp) :: h(mr+1,mr)
      real    (kind=dp) :: htmp
      integer (kind=si) :: i
      integer (kind=si) :: itr
      integer (kind=si) :: itr_used
      integer (kind=si) :: j
      integer (kind=si) :: k
      integer (kind=si) :: k_copy
      real    (kind=dp) :: l(ia(n+1)+1)
      real    (kind=dp) :: mu
      real    (kind=dp) :: r(n)
      real    (kind=dp) :: rho
      real    (kind=dp) :: rho_tol
      real    (kind=dp) :: s(mr+1)
      integer (kind=si) :: ua(n)
      real    (kind=dp) :: v(n,mr+1);
      logical, parameter :: verbose = .true.
      real    (kind=dp) :: y(mr+1)

      itr_used = 0

      call rearrange_cr ( n, nz_num, ia, ja, a )

      call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

      call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

      if ( verbose ) then
        write (*, 1000) n
1000    format (' ', /, 'PMGMRES_ILU_CR', /,                             &
     &    '  Number of unknowns = ', i4)
      end if

      do itr = 1, itr_max

        call ax_cr ( n, nz_num, ia, ja, a, x, r )

        call lus_cr ( n, nz_num, ia, ja, l, ua, rhs - r, r )

        rho = sqrt ( dot_product ( r, r ) )

        if ( verbose ) then
          write (*, 1010) itr, rho
1010      format ('  ITR = ', i4, '  Residual = ', g14.6)
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

          call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1),            &
     &                                            v(1:n,k+1))

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
            write (*, 1020) k, rho
1020        format ('  K = ', i4, '  Residual = ', g14.6)
          end if

          if ( rho <= rho_tol .and. rho <= tol_abs ) then
            exit
          end if

        end do

        k = k_copy - 1

        y(k+1) = g(k+1) / h(k+1,k+1)

        do i = k, 1, -1
          y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) )     &
     &           / h(i,i)
        end do

        do i = 1, n
          x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
        end do

        if ( rho <= rho_tol .and. rho <= tol_abs ) then
          exit
        end if

      end do

      if ( verbose ) then
        write (*, 1030) itr_used, rho
1030    format (' ', /, 'PMGMRES_ILU_CR:', /, '  Iterations = ', i6, /,  &
     &     '  Final residual = ', g14.6)
      end if

      return
      end subroutine pmgmres_ilu_cr

      subroutine rearrange_cr ( n, nz_num, ia, ja, a )
!*****************************************************************************80
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Parameters:
!
!    N, the order of the system.
!
!    NZ_NUM, the number of nonzeros.
!
!    IA(N+1), the compressed row indices.
!
!    JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
      implicit none

      integer (kind=si), intent (in)    :: n
      integer (kind=si), intent (in)    :: nz_num
      integer (kind=si), intent (in)    :: ia(n+1)
      integer (kind=si), intent (inout) :: ja(nz_num)
      real    (kind=dp), intent (inout) :: a(nz_num)

      integer (kind=si) :: i
      integer (kind=si) :: i4temp
      integer (kind=si) :: k
      integer (kind=si) :: l
      real    (kind=dp) :: r8temp

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
      end subroutine rearrange_cr

      end module mgmres
