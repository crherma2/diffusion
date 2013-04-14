module initialize

!-module options

  implicit none
  private
  public :: initialize_run

contains

!-initialize run

  subroutine initialize_run()

    ! read input
    call read_input()

  end subroutine initialize_run

!-read input

  subroutine read_input()
    
    use geometry_header, only: allocate_geom
    use global,          only: geometry, material
    use global,          only: ktol
    use global,          only: ftol
    use global,          only: itol
    use material_header, only: allocate_mat
    use xml_data_input_t
  
!---define variables

    character(9) :: filename
    integer      :: i        ! loop counter
    integer      :: size_mat ! size of material array
!---begin execution
    
    ! read in input file
    filename = "input.xml"
    call read_xml_file_input_t(filename)

    ! read in coarse mesh
    geometry % nxc = geometry_ % nx
    geometry % nyc = geometry_ % ny
    geometry % ng  = geometry_ % ng

    ! allocate geometry properties
    call allocate_geom(geometry)
    
    ! read in fine mesh
    geometry % nnx = geometry_ % nnx
    geometry % nny = geometry_ % nny

    ! calculate the total number of fine mesh
    geometry % nx = sum(geometry_ % nnx)
    geometry % ny = sum(geometry_ % nny)
   
    ! read in coarse mesh regions
    geometry % gridx = geometry_ % xgrid
    geometry % gridy = geometry_ % ygrid
    
    ! read in coarse material map
    geometry % cmat = reshape(geometry_ % mat, (/geometry % nxc, &
                                                 geometry % nyc/))
    
    ! calculate fine mesh with
    geometry % dx = geometry % gridx / geometry % nnx
    geometry % dy = geometry % gridy / geometry % nny

    ! read in boundary conditions
    geometry % bc = geometry_ % bc
  
    ! read in size of material array
    size_mat = size(material_)
    
    ! allocate the material array
    allocate(material(size_mat))

    ! read in material with loop for different materials
    do i = 1,size_mat
      call allocate_mat(material(i),geometry % ng)
      material(i) % xs_tot   = material_(i) % totalxs
      material(i) % xs_scat  = reshape(material_(i) % scattxs, (/geometry % &
                                                                 ng, geometry &
                                                                 % ng/))
      material(i) % xs_nufis = material_(i) % nfissxs
      material(i) % chi      = material_(i) % chi
      material(i) % difco    = material_(i) % diffcof
      PRINT*, material(i) % xs_tot
   end do

   ! read in tolerances
   ktol = ktol_
   ftol = ftol_
   itol = itol_

  end subroutine read_input

end module initialize
