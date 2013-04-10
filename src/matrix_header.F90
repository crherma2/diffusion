module matrix_header

!-module options

  implicit none
  private
  public :: allocate_matrix

!-module variables

  type, public :: matrix_type

    ! matrix variables
    integer :: n                   ! order of matrix
    integer :: nz                  ! number of non-zero numbers
    integer, allocatable :: row(:) ! rows
    integer, allocatable :: col(:) ! columns
    real(8), allocatable :: val(:) ! values that go into matrix

  end type matrix_type

contains

  subroutine allocate_matrix(self)

!---arguments

    type(matrix_type) :: self

!---begin execution

    ! allocate
    allocate(self % row (self % nz))
    allocate(self % col (self % nz))
    allocate(self % val (self % nz))

  end subroutine allocate_matrix

  subroutine deallocate_matrix(self)
  
!---arguments

    type(matrix_type) :: self

!---begin execution

    ! deallocate
    deallocate(self % row)
    deallocate(self % col)
    deallocate(self % val)

  end subroutine deallocate_matrix

end module matrix_header
