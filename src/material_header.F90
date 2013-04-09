module material_header

!-module options

  implicit none
  private
  public :: allocate_mat, deallocate_mat

!-module variables

  type, public :: material_type

    ! cross sections
    real(8), allocatable :: xs_tot(:)    ! total cross section
    real(8), allocatable :: xs_scat(:,:) ! scattering cross section
    real(8), allocatable :: xs_nufis(:)  ! nu fission cross section

    ! fission spectrum
    real(8), allocatable :: chi(:)
    
    ! diffusion coefficient
    real(8), allocatable :: difco(:)

  end type material_type

contains

  subroutine allocate_mat(self, ng)
  
!---arguments

    type(material_type) :: self
    integer :: ng

!---begin execution

    ! allocate
    allocate(self % xs_tot (ng))
    allocate(self % xs_scat (ng, ng))
    allocate(self % xs_nufis (ng))
    allocate(self % chi (ng))
    allocate(self % difco (ng))

  end subroutine allocate_mat

  subroutine deallocate_mat(self)

!---arguments

    type(material_type) :: self

!---begin execution

    ! deallocate
    deallocate(self % xs_tot)
    deallocate(self % xs_scat)
    deallocate(self % xs_nufis)
    deallocate(self % chi)
    deallocate(self % difco)

  end subroutine deallocate_mat

end module material_header


