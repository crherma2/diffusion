module output

!-module options

  implicit none
  private
  public :: post_process

contains

!-outer routine

  subroutine post_process()
 
    use global,  only: geometry
 
    ! call the normalize flux subroutine
    call normalize_flux()

    ! call surface currents
    call surface_currents()

    ! write out mesh file for gnuplot if ny > 2
    call write_results()
    if (geometry % ny > 2) then
      call write_gnuplot_mesh()
    else
      call write_gnuplot_1D()
    end if

  end subroutine post_process

!-normalizing the flux

  subroutine normalize_flux()

  use global, only: flux
  use global, only: geometry
  use geometry_header, only: get_coarsemesh
 
!---local variables
    integer :: i ! iteration counter for x
    integer :: j ! iteration counter for y
    integer :: g ! iteration counter for groups
    integer :: irow ! loop around rows
    integer :: ic ! coarse mesh index x
    integer :: jc ! coarse mesh index y
    real(8) :: flux_sum ! summation of flux     

    ! set counter equal to zero
    flux_sum = 0.0_8    

    ! loop around rows
    ROWS: do irow = 1, size(flux)
   
      ! call matrix to indice code below
      call matrix_to_indices(irow,i,j,g)

      ! call to get coarse mesh indices
      call get_coarsemesh(geometry,i,j,ic,jc)
      
      ! append the counter
      flux_sum = flux_sum + flux(irow)*geometry % dx(ic)*geometry % dy(jc)
      
    end do ROWS
    
    ! normalize the flux
    flux = flux/flux_sum


  end subroutine normalize_flux

  subroutine surface_currents()

    use global, only: material
    use global, only: flux
    use global, only: curr
    use global, only: geometry
    use geometry_header, only: get_coarsemesh
    use material_header, only: material_type
 
!---arguments

    integer :: irow ! loop around rows
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
    real(8) :: du ! cell width
    real(8) :: ndu ! neighbor cell width
    logical :: bound ! check if at boundary
    type(material_type), pointer :: m => null() ! creates a pointer shortens word
    type(material_type), pointer :: nm => null() ! creates a pointer neihbor

!---begin execution
    
    ! allocate the current array
    allocate(curr(size(flux),4))   
    
    ! go row by row to fill in values of matrix
    ROWS: do irow = 1, size(flux)
     
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
            pm = -1.0_8
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
            pm = 1.0_8
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
            pm = -1.0_8
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
            pm = 1.0_8
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
          
          curr(irow,s) = pm*dtilda*flux(irow)
        else
          
          ! get neighbor material
          nmatid = geometry % cmat(nic,njc)
          nm => material(nmatid)
          
          ! calculate diffusion coupling coefficient for neighbor
          dtilda = (2.0_8*m % difco(g)*nm % difco(g))/ &
                    (du*m % difco(g) + ndu*nm % difco(g))
         
          ! get column index of neighbor
          call indices_to_matrix(ni,nj,g,colidx)

          ! calculate the current
          curr(irow,s) = -dtilda*(pm*flux(colidx)-pm*flux(irow))

        end if

      end do LEAK 

    end do ROWS
  
  end subroutine surface_currents

  subroutine write_results()

    use global, only : flux
    use global, only : curr

!---local variables
    integer :: i ! loop index

!---open files

    ! open file for the normalized flux
    open(file='flux', unit=15)

    ! open file for the surface currents
    open(file='current', unit=18)

    ! loop around size of rows
    do i = 1, size(flux)
      write(15,'(1PE13.6)') flux(i)
      write(18,'(4(1PE13.6,1X))') curr(i,1), curr(i,2), curr(i,3), curr(i,4)
    end do

    close(15)
    close(18)

  end subroutine write_results

  subroutine write_gnuplot_1D()

  use global, only: flux, curr
  use global, only: geometry

!---local variables
    integer :: i ! iteration counter for x
    integer :: j ! iteration counter for y
    real(8), allocatable :: flux_mat(:,:,:)
    real(8), allocatable :: curr_mat(:,:,:)
    real(8), allocatable :: curr2temp(:,:,:)

    ! open up files
    open(file='flux1.dat',unit=35)
    open(file='flux2.dat',unit=36)
    open(file='curr1.dat',unit=37)
    open(file='curr2.dat',unit=38)
 

    ! allocate temp data vars 
    allocate(flux_mat(geometry % nx, geometry % ny, geometry % ng))
    allocate(curr_mat(geometry % nx + 1, geometry % ny, geometry % ng))
    allocate(curr2temp(geometry % nx, geometry % ny, geometry % ng))

    ! reshape all
    flux_mat = reshape(flux,(/geometry % nx, geometry % ny, geometry % ng/))
    curr_mat(1:geometry % nx,:,:) = reshape(curr(:,1), &
                         (/geometry % nx, geometry % ny, geometry % ng/))
    curr2temp = reshape(curr(:,2),(/geometry % nx, geometry % ny, geometry % ng/))

    ! append last surface on right
    curr_mat(geometry % nx + 1,:,1) = curr2temp(geometry % nx,:,1)
    curr_mat(geometry % nx + 1,:,2) = curr2temp(geometry % nx,:,2)

    ! loop around rows
    do i = 1, geometry % nx
      write(35,'(1PE13.6)') flux_mat(i,:,1)
      write(36,'(1PE13.6)') flux_mat(i,:,2)
      write(37,'(4(1PE13.6,1X))') curr_mat(i,:,1)
      write(38,'(4(1PE13.6,1X))') curr_mat(i,:,2)
    end do
    write(37,'(4(1PE13.6,1X))') curr_mat(geometry % nx + 1,:,1)
    write(38,'(4(1PE13.6,1X))') curr_mat(geometry % nx + 1,:,2)

    ! deallocate variables
    deallocate(flux_mat)
    deallocate(curr_mat)
    deallocate(curr2temp)

    ! close files
    close(35)
    close(36)
    close(37)
    close(38)

  end subroutine write_gnuplot_1D

  subroutine write_gnuplot_mesh()

  use global, only: flux
  use global, only: geometry

!---local variables
    integer :: i ! iteration counter for x
    integer :: j ! iteration counter for y
    real(8), allocatable :: flux_mat(:,:,:)

    ! open up files
    open(file='flux1.dat',unit=30)
    open(file='flux2.dat',unit=31)

    ! reshape flux into matrix
    allocate(flux_mat(geometry % nx, geometry % ny, geometry % ng))
    flux_mat = reshape(flux,(/geometry % nx, geometry % ny, geometry % ng/))

    ! loop around rows
    do i = 1, size(flux_mat,1)
      write(30,*) (flux_mat(i,j,1),j=1,size(flux_mat,2))
      write(31,*) (flux_mat(i,j,2),j=1,size(flux_mat,2))
    end do

    ! deallocate variables
    deallocate(flux_mat)

    ! close files
    close(30)
    close(31)

  end subroutine write_gnuplot_mesh 

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
  
end module output
