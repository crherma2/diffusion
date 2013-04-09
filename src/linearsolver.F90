module linearsolver 

!-module options

  implicit none

contains

  subroutine gmres(n, nnz, row, col, val, x, b, err, iter)

!---arguments

    integer, intent(in)    :: n        ! order of matrix
    integer, intent(in)    :: nnz      ! number of nonzeros in sparse matrix
    integer, intent(in)    :: row(nnz) ! vector of rows
    integer, intent(in)    :: col(nnz) ! vector of columns
    integer, intent(inout) :: iter     ! max iterations on input, actual on output
    real(8), intent(in)    :: val(nnz) ! vector of matrix values
    real(8), intent(inout) :: x(n)     ! unknown vector
    real(8), intent(inout) :: b(n)     ! right hand side vector
    real(8), intent(inout) :: err      ! tolerance on input, actual error on output

!---local variables

    integer :: ierr                  ! error code
    integer :: lenw                  ! length of rwork
    integer :: leniw                 ! length of iwork
    integer :: iterout               ! acutal iterations
    integer, allocatable :: iwork(:) ! work vector for integer
    real(8) :: errout                ! final converged error
    real(8), allocatable :: rwork(:) ! working vector for real

!---begin execution

    ! create temp working vectors for GMRES (see dslugm.f for details)
    lenw = 131 + 18*n + nnz
    leniw = 32 + 5*n + nnz
    allocate(rwork(lenw))
    allocate(iwork(leniw))

    ! solve matrix with GMRES
    call DSLUGM(n, b, x, nnz, row, col, val, 0, n, 0, err, iter, iterout, &
                errout, ierr, 0, rwork, lenw, iwork, leniw)

    ! return allocated memory
    deallocate(rwork)
    deallocate(iwork)

    ! send out iter and err
    iter = iterout
    err = errout

  end subroutine gmres 

end module linearsolver 
