module global

!-module references

  use geometry_header, only: geometry_type

!-module options

  implicit none
  save
  
!-variables

  type(geometry_type) :: geometry 

end module global
