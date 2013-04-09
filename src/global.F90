module global

!-module references

  use geometry_header, only: geometry_type
  use material_header, only: material_type

!-module options

  implicit none
  save
  
!-variables

  type(geometry_type) :: geometry 
  type(material_type), allocatable :: materials(:)

end module global
