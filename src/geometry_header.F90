module geometry_header

!-module options

  implicit none
  private
  public :: allocate_geom, deallocate_geom, get_coarsemesh

!-module variables

  type, public :: geometry_type
  
    ! coarse indices
    integer :: nxc  ! number of coarse meshes in x-direction
    integer :: nyc  ! number of coarse meshes in y-direction
    integer :: ng   ! number of energy groups

    ! number of fine mesh per coarse mesh
    integer, allocatable :: nnx(:)  ! x-direction
    integer, allocatable :: nny(:)  ! y-direction
   
    ! total number fine meshes
    integer :: nx ! x-direction
    integer :: ny ! y-direction
    
    ! width of course mesh regions
    real(8), allocatable :: gridx(:)
    real(8), allocatable :: gridy(:)
    
    ! coarse material map
    integer, allocatable :: cmat(:,:)

    ! fine mesh width
    real(8), allocatable :: dx(:)  ! x-direction
    real(8), allocatable :: dy(:)  ! y-direction

    ! boundary conditions
    real(8) :: bc(4) = 1.0_8 ! default is reflective

  end type geometry_type

contains 

  subroutine allocate_geom(self)
    
!---arguments

    type(geometry_type) :: self

!---begin execution

    ! allocate
    allocate(self % nnx (self % nxc))
    allocate(self % nny (self % nyc))
    allocate(self % gridx (self % nxc))
    allocate(self % gridy (self % nyc))
    allocate(self % cmat (self % nxc, self % nyc))
    allocate(self % dx (self % nxc))
    allocate(self % dy (self % nyc))

  end subroutine allocate_geom

  subroutine deallocate_geom(self)

!---arguments

    type(geometry_type) :: self

!---begin execution

    ! deallocate
    deallocate(self % nnx)
    deallocate(self % nny)
    deallocate(self % gridx)
    deallocate(self % gridy)
    deallocate(self % cmat)
    deallocate(self % dx)
    deallocate(self % dy)

  end subroutine deallocate_geom   

  subroutine get_coarsemesh(self,i,j,ic,jc)

!---arguments

    integer :: i ! read in x indice
    integer :: j ! read in y indce
    integer :: ic ! determine which mesh is in for x
    integer :: jc ! determine which mesh is in for y
    type(geometry_type) :: self ! for use of variables in above file
    
!---local variables
   
    integer :: icor ! loop counter for x
    integer :: jcor ! loop counter for y
   
!---begin execution

    ! loop around the fine mesh
    do icor = 1,self % nxc
      if (i <= sum(self % nnx(1:icor))) then
        ic = icor
        exit
      end if
    end do 

    do jcor = 1,self % nyc
      if (j <= sum(self % nny(1:jcor))) then
        jc = jcor
        exit
      end if
    end do

  end subroutine get_coarsemesh

end module geometry_header
