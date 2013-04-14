program main

!-references
  
  use initialize, only: initialize_run
  use loss,       only: create_lossmatrix
  use matrix,     only: test_matrix
  use prod,       only: create_prodmatrix

!-program options

  implicit none

!-begin execution

  call test_matrix()
  
  ! initialize input
  call initialize_run()  

  ! create loss matrix
  call create_lossmatrix()

  ! create prod matrix
  call create_prodmatrix()

end program main
