module matrix

!-module options

  implicit none

contains

  subroutine test_matrix()

!---references

    use linearsolver,  only: gmres

!---local variables

    integer :: N=3    ! order of matrix
    integer :: NNZ=6  ! number of nonzeros in matrix
    integer :: row(6) ! row vector of matrix row elements
    integer :: col(6) ! col vector of matrix column elements
    real(8) :: val(6) ! val vector of matrix value elements
    real(8) :: b(3)   ! right hand side vector
    real(8) :: x(3)   ! solution vector
    real(8) :: tol = 1.e-10_8 ! linear tolerance of GMRES solver
    integer :: iter=100 ! number of iterations that GMRES took

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

    ! solve matrix
    call gmres(N, NNZ, row, col, val, x, b, tol, iter)

    ! print answer
    print *, 'Solution is:'
    print *, x

  end subroutine test_matrix

end module matrix
