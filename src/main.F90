program main

!-references
  
  use eigenval,   only: solve_eigenval
  use initialize, only: initialize_run
  use loss,       only: create_lossmatrix
  use prod,       only: create_prodmatrix

!-program options

  implicit none

!-begin execution

  ! initialize input
  call initialize_run()  

  ! create loss matrix
  call create_lossmatrix()

  ! create prod matrix
  call create_prodmatrix()

  ! power iteration
  call solve_eigenval() 

end program main
