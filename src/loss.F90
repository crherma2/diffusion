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
    
    nz_c = geometry % ng*n_c*(3+geometry % ng-1) ! define # nonzero corners
    nz_s = geometry % ng*n_s*(4+geometry % ng-1) ! define # nonzero sides
    nz_i = geometry % ng*n_i*(5+geometry % ng-1) ! define # nonzero interiors
    loss_matrix % nz = nz_c+nz_s+nz_i ! total number of non-zero entries
    loss_matrix % n = geometry % ng*geometry % nx*geometry % ny ! order of matrix
 
    PRINT*, loss_matrix % nz

    ! allocate matrix dimensions
    call allocate_matrix(loss_matrix)

  end subroutine preallocate_lossmatrix

  subroutine build_lossmatrix()

    use global, only: material
    use global, only: loss_matrix
    use global, only: geometry
    use geometry_header, only: get_coarsemesh   
    use material_header, only: material_type  
 
!---arguments

    integer :: irow ! loop around rows
    integer :: counter ! counter indice of matrix
    integer :: i ! index for x
    integer :: j ! index for y
    integer :: g ! index for groups
    integer :: h ! do loop for groups
    integer :: colidx ! index for column
    integer :: matid ! material identification
    integer :: nmatid ! neighbor material indentification
    integer :: s ! loop around surfaces
    integer :: ni ! neighbor index x
    integer :: nj ! neighbor index y
    integer :: ic ! coarse mesh index x
    integer :: jc ! coarse mesh index y
    integer :: nic ! neighbor coarse mesh index x
    integer :: njc ! neighbor coarse mesh index y
    real(8) :: dtilda ! diffusion coupling term
    real(8) :: pm ! plus minus multiplier
    real(8) :: diag ! banking all materials for diagonal
    real(8) :: du ! cell width
    real(8) :: ndu ! neighbor cell width 
    logical :: bound ! check if at boundary
    type(material_type), pointer :: m => null() ! creates a pointer shortens word  
    type(material_type), pointer :: nm => null() ! creates a pointer neihbor

!---begin execution
    
    ! set counter
    counter = 1
   
    ! go row by row to fill in values of matrix
    ROWS: do irow = 1,loss_matrix % n 
     
      ! reset diag value
      diag = 0.0_8
 
      ! call matrix to indice code below           
      call matrix_to_indices(irow,i,j,g)

      ! call to get coarse mesh indices
      call get_coarsemesh(geometry,i,j,ic,jc)
      matid = geometry % cmat(ic,jc)
      m => material(matid) 
      
      ! start do loop by looping around rows
      LEAK: do s = 1,4

        if (s == 1) then ! minus x direction
          if (i == 1) then ! check if you are on left boundary
            bound = .true.
            pm = -1.0_8
          else
            bound = .false.
            ni = i - 1 ! neighbor defined index x direction to left
            nj = j ! neighbor defined index y direction to left
            call get_coarsemesh(geometry,ni,nj,nic,njc)
            ndu = geometry % dx(nic)
          end if
          du = geometry % dx(ic)
        else if (s == 2) then ! plus x direction
          if (i == geometry % nx) then ! check if you are on right boundary
            bound = .true.
            pm = 1.0_8
          else
            bound = .false.
            ni = i + 1 ! neighbor defined index x direction to the right
            nj = j ! neighbor defined index y direction to right
            call get_coarsemesh(geometry,ni,nj,nic,njc)
            ndu = geometry % dx(nic)
          end if
          du = geometry % dx(ic)
        else if (s == 3) then ! minus y direction
          if (j == 1) then ! check if you are on the bottom boundary
            bound = .true.
            pm = -1.0_8
          else
            bound = .false.
            ni = i ! neighbor defiend index x direction to bottom
            nj = j - 1 ! neighbor defined index y direction to bottom
            call get_coarsemesh(geometry,ni,nj,nic,njc)
            ndu = geometry % dy(njc)
          end if
          du = geometry % dy(jc)
        else if (s == 4) then ! plus y direction
          if (j == geometry % ny) then ! check if you are on top boundary
            bound = .true.
            pm = 1.0_8
          else
            bound = .false.
            ni = i ! neighbor defined index x direction to top
            nj = j + 1 ! neighbor defined index y direction to top
            call get_coarsemesh(geometry,ni,nj,nic,njc)
            ndu = geometry % dy(njc)
          end if
          du = geometry % dy(jc)
        end if
        

        ! check if at boundary
        if (bound) then
          ! calculate diffusion coupling coefficient for boundary
          dtilda = (2.0_8*m % difco(g)*(1.0_8 - geometry % bc(s)))/ &
                   (du*(1.0_8 - geometry % bc(s)) + &
                   4.0_8*m % difco(g)*(1.0_8 + geometry % bc(s)))
          
          ! bank in temporary diagonal term
          diag = diag + pm*(dtilda/du)
        else
          
          ! get neighbor material
          nmatid = geometry % cmat(nic,njc)
          nm => material(nmatid)
          
          ! calculate diffusion coupling coefficient for neighbor
          dtilda = -(2.0_8*m % difco(g)*nm % difco(g))/ &
                    (du*m % difco(g) + ndu*nm % difco(g)) 
          
          ! bank in temporary diagonal term
          diag = diag + (dtilda/du)
          
          ! get column index of neighbor
          call indices_to_matrix(ni,nj,g,colidx)

          ! bank negative into off diagonal
          loss_matrix % row(counter) = irow
          loss_matrix % col(counter) = colidx
          loss_matrix % val(counter) = -dtilda/du
          counter = counter + 1

        end if           

      end do LEAK  

      ! get rest of diagonal terms
      diag = diag + m % xs_tot(g) - m % xs_scat(g,g)   

      ! put diagonal term into matrix
      loss_matrix % row(counter) = irow
      loss_matrix % col(counter) = irow
      loss_matrix % val(counter) = diag
      counter = counter + 1
 
      ! get the inscattering term
      INSCATTER: do h = 1, geometry % ng

        if (h == g) cycle    

        ! getting column identification
        call indices_to_matrix(i,j,h,colidx)
        
        !put the inscattering terms in matrix
        loss_matrix % row(counter) = irow
        loss_matrix % col(counter) = colidx
        loss_matrix % val(counter) = -m % xs_scat(g,h)    
        counter = counter + 1
                         
      end do INSCATTER
  
    end do ROWS 
    PRINT*,'counter', counter 
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
