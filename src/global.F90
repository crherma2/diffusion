module global

!-module references

  use geometry_header, only: geometry_type
  use material_header, only: material_type
  use matrix_header,   only: matrix_type
!-module options

  implicit none
  save
  
!-variables

  type(geometry_type) :: geometry 
  type(material_type), allocatable, target :: material(:)
  type(matrix_type) :: loss_matrix
  type(matrix_type) :: prod_matrix

end module global
