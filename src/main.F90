program main

!-references
  
  use initialize, only: initialize_run
  use matrix,     only: test_matrix

!-program options

  implicit none

!-begin execution

  call test_matrix()
  ! initialize input
  call initialize_run()  

end program main
