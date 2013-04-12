module loss

!-module options

  implicit none
  private
  public :: create_lossmatrix

contains

!-create loss matrix

  subroutine create_lossmatrix()

    ! start preallocate
    call preallocate_lossmatrix()

  end subroutine create_lossmatrix

  subroutine preallocate_lossmatrix()

    use global, only: geometry
    use global, only: loss_matrix
    use matrix_header, only: allocate_matrix
!---arguments
    
    integer :: n_i  ! number of interior cells
    integer :: n_c  ! number of corner cells
    integer :: n_s  ! number side cells
    integer :: nz_c ! number of non-zero corner cells
    integer :: nz_s ! number of non-zero side cells
    integer :: nz_i ! number of non-zero interior cells

!---begin execution

    n_c = 4                   ! defined number of corners
    n_s = 2*(geometry % nxc+geometry % nyc)-8       ! defined number of sides
    n_i = geometry % nxc*geometry % nyc-(n_c+n_S) ! define # of interiors
    
    nz_c = geometry % ng*2*n_c*(3+geometry % ng-1) ! define # nonzero corners
    nz_s = geometry % ng*3*n_s*(4+geometry % ng-1) ! define # nonzero sides
    nz_i = geometry % ng*4*n_i*(5+geometry % ng-1) ! define # nonzero interiors
    loss_matrix % nz = nz_c+nz_s+nz_i ! total number of non-zero entries

    PRINT*, loss_matrix % nz

    ! allocate matrix dimensions
    call allocate_matrix(loss_matrix)

  end subroutine preallocate_lossmatrix

end module loss
