program main

!-references
  
  use eigenval,   only: solve_eigenval
  use initialize, only: initialize_run
  use loss,       only: create_lossmatrix
  use output,     only: post_process
  use prod,       only: create_prodmatrix
  use timing,     only: timer_start
  use timing,     only: timer_stop
  use global,     only: runtime

!-program options

  implicit none

!-begin execution

  ! start timer
  call timer_start(runtime)

  ! initialize input
  call initialize_run()  

  ! create loss matrix
  call create_lossmatrix()

  ! create prod matrix
  call create_prodmatrix()

  ! power iteration
  call solve_eigenval() 

  ! post processing
  call post_process()

  ! stop timer
  call timer_stop(runtime)

  ! print out timer
  write(*,'(A,2X,ES11.4,1X,A)') 'Total run time:', runtime % elapsed, 'seconds.'

end program main
