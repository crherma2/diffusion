module prod

!-module options

  implicit none
  private
  public :: create_prodmatrix

contains

!-create prod matrix

  subroutine create_prodmatrix()

    ! start preallocate
    call preallocate_prodmatrix()
    
    ! build prod matrix
    call build_prodmatrix()

  end subroutine create_prodmatrix

  subroutine preallocate_prodmatrix()

    use global, only: geometry
    use global, only: prod_matrix
    use matrix_header, only: allocate_matrix

!---begin execution

    prod_matrix % n = geometry % ng*geometry % nx*geometry % ny ! order of matrix
    prod_matrix % nz = prod_matrix % n*geometry % ng ! total number of nonzeros

    PRINT*, prod_matrix % nz

    ! allocate matrix dimensions
    call allocate_matrix(prod_matrix)

  end subroutine preallocate_prodmatrix

  subroutine build_prodmatrix()

    use global, only: material
    use global, only: prod_matrix
    use global, only: geometry
    use geometry_header, only: get_coarsemesh   
    use material_header, only: material_type  
 
!---local variables

    integer :: irow ! loop around rows
    integer :: counter ! counter indice of matrix
    integer :: i ! index for x
    integer :: j ! index for y
    integer :: g ! index for groups
    integer :: h ! do loop for groups
    integer :: colidx ! index for column
    integer :: matid ! material identification
    integer :: ic ! coarse mesh index x
    integer :: jc ! coarse mesh index y
    type(material_type), pointer :: m => null() ! creates a pointer shortens word  

!---begin execution
    
    ! set counter
    counter = 1
   
    ! go row by row to fill in values of matrix
    ROWS: do irow = 1,prod_matrix % n 
 
      ! call matrix to indice code below           
      call matrix_to_indices(irow,i,j,g)

      ! call to get coarse mesh indices
      call get_coarsemesh(geometry,i,j,ic,jc)
      matid = geometry % cmat(ic,jc)
      m => material(matid) 
      
      ! get fission term
      FISSION: do h = 1, geometry % ng

        ! getting column identification
        call indices_to_matrix(i,j,h,colidx)
        
        !put fission terms in matrix
        prod_matrix % row(counter) = irow
        prod_matrix % col(counter) = colidx
        prod_matrix % val(counter) = m % chi(g)*m % xs_nufis(h)   
        counter = counter + 1
                         
      end do FISSION
  
    end do ROWS 
    PRINT*,'counter', counter 
  end subroutine build_prodmatrix

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

end module prod
