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
    use global,          only: geometry
    use xml_data_input_t

!---define variables

    character(9) :: filename
    
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

    ! read in coarse mesh regions
    geometry % gridx = geometry_ % xgrid
    geometry % gridy = geometry_ % ygrid
    
    ! read in coarse material map
    geometry % cmat = reshape(geometry_ % mat, (/geometry % nxc, &
                                                 geometry % nyc/))
    
    ! calculate fine mesh with
    geometry % dx = geometry % gridx / geometry % nnx
    geometry % dy = geometry % gridy / geometry % nny
  
  end subroutine read_input

end module initialize