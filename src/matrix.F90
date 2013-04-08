module matrix

!-module options

  implicit none

contains

  subroutine test_matrix()

!---local variables

    integer :: N=3    ! order of matrix
    integer :: NNZ=6  ! number of nonzeros in matrix
    integer :: row(6) ! row vector of matrix row elements
    integer :: col(6) ! col vector of matrix column elements
    real(8) :: val(6) ! val vector of matrix value elements
    real(8) :: b(3)   ! right hand side vector
    real(8) :: x(3)   ! solution vector
    real(8) :: tol = 1.e-10_8 ! linear tolerance of GMRES solver
    integer :: iter ! number of iterations that GMRES took
    real(8) :: err ! final error
    integer :: ierr ! error code
    real(8), allocatable :: rwork(:) ! working vector for real
    integer :: lenw ! length of rwork
    integer, allocatable :: iwork(:) ! work vector for integer
    integer :: leniw ! length of iwork

!---begin execution

    ! (1,1) element
    row(1) = 1
    col(1) = 1
    val(1) = 2.0_8

    ! (1,3) element
    row(2) = 1
    col(2) = 3
    val(2) = 4.0_8

    ! (2,2) element
    row(3) = 2
    col(3) = 2
    val(3) = 5.0_8

    ! (2,3) element
    row(4) = 2
    col(4) = 3
    val(4) = 1.0_8

    ! (3,1) element
    row(5) = 3
    col(5) = 1
    val(5) = 3.0_8

    ! (3,3) element
    row(6) = 3
    col(6) = 3
    val(6) = 4.0_8

    ! RHS vector
    b = (/1.0_8,2.0_8,3.0_8/)

    ! guess answer
    x = 1.0_8

    ! create temp working vectors for GMRES (see dslugm.f for details)
    lenw = 131 + 18*N + NNZ
    leniw = 32 + 5*N + NNZ
    allocate(rwork(lenw))
    allocate(iwork(leniw))

    ! solve matrix with GMRES
    call DSLUGM(N, b, x, NNZ, row, col, val, 0, 3, 0, tol, 100, iter, &
                err, ierr, 0, rwork, lenw, iwork, leniw)
    print *, iter,err,ierr
    ! print answer
    print *, 'Solution is:'
    print *, x

    ! return allocated memory
    deallocate(rwork)
    deallocate(iwork)

  end subroutine test_matrix

end module matrix
