program main

!-references
  
  use initialize, only: initialize_run
  use loss,       only: create_lossmatrix
  use matrix,     only: test_matrix

!-program options

  implicit none

!-begin execution

  call test_matrix()
  
  ! initialize input
  call initialize_run()  

  ! create loss matrix
  call create_lossmatrix()

end program main
