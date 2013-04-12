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
    
    ! build loss matrix
    call build_lossmatrix()

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
    n_s = 2*(geometry % nx+geometry % ny)-8       ! defined number of sides
    n_i = geometry % nx*geometry % ny-(n_c+n_S) ! define # of interiors
    
    nz_c = geometry % ng*2*n_c*(3+geometry % ng-1) ! define # nonzero corners
    nz_s = geometry % ng*3*n_s*(4+geometry % ng-1) ! define # nonzero sides
    nz_i = geometry % ng*4*n_i*(5+geometry % ng-1) ! define # nonzero interiors
    loss_matrix % nz = nz_c+nz_s+nz_i ! total number of non-zero entries
    loss_matrix % n = geometry % ng*geometry % nx*geometry % ny ! order of matrix
 
    PRINT*, loss_matrix % nz

    ! allocate matrix dimensions
    call allocate_matrix(loss_matrix)

  end subroutine preallocate_lossmatrix

  subroutine build_lossmatrix()

    use global, only: loss_matrix

!---arguments

    integer :: irow ! loop around rows
    integer :: counter ! counter indice of matrix
    integer :: i ! iteration counter for x
    integer :: j ! iteration coutner for y
    integer :: g ! iteration counter for groups
   
!---begin execution
    
    ! set counter
    counter = 1
   
    ! Go row by row to fill in values of matrix
    do irow = 1,loss_matrix % n 
 
      ! call matrix to indice code below           
      call matrix_to_indices(irow,i,j,g)
      PRINT*, irow, i, j, g

    end do 
    
  end subroutine build_lossmatrix

  subroutine matrix_to_indices(irow,i,j,g)

    use global, only: geometry

    integer :: i ! iteration counter for x
    integer :: j ! iteration counter for y
    integer :: g ! iteration counter for groups
    integer :: irow ! iteration counter over row (0 reference)
    integer :: nx ! x total number of fine meshes from geometry
    integer :: ny ! y total number of fine meshes from geometry
    integer :: ng ! number of energy groups from geometry

    ! read in variables from geometry
    nx = geometry % nx
    ny = geometry % ny
    ng = geometry % ng

    ! compute indices
    i = mod(irow-1,nx) + 1
    j = mod(irow-1,ny*nx)/nx + 1
    g = mod(irow-1,ng*ny*nx)/(ny*nx) + 1

  end subroutine matrix_to_indices

  subroutine indices_to_matrix(i,j,g,matidx)

    use global, only: geometry

    integer :: matidx ! the index location in matrix
    integer :: i ! current x index
    integer :: j ! current y index
    integer :: g ! current group index
    integer :: nx ! x total number of fine meshes from geometry
    integer :: ny ! y total number of fine meshes from geometry
    integer :: ng ! number of energy groups from geometry

    ! read in variables from geometry
    nx = geometry % nx
    ny = geometry % ny
    ng = geometry % ng

    ! compute index
    matidx = i + nx*(j - 1) + nx*ny*(g - 1)

  end subroutine indices_to_matrix

end module loss
