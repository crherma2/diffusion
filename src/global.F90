module global

!-module references

  use geometry_header, only: geometry_type
  use material_header, only: material_type
  use matrix_header,   only: matrix_type
!-module options

  implicit none
  save
  
!-variables

  ! object variables
  type(geometry_type) :: geometry 
  type(material_type), allocatable, target :: material(:)
  type(matrix_type) :: loss_matrix
  type(matrix_type) :: prod_matrix
 
  ! solver tolerances
  real(8) :: ktol ! eigenvalue tolerance
  real(8) :: ftol ! flux tolerance
  real(8) :: itol ! inner linear solver tolerance
  
  ! results
  real(8) :: keff
  real(8), allocatable :: flux(:)
  real(8), allocatable :: curr(:,:)

end module global
