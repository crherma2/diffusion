module eigenval

!-module options

  implicit none
  private
  public :: solve_eigenval

contains

  subroutine solve_eigenval()

    call power_iter()

  end subroutine solve_eigenval

  subroutine power_iter()
    
    use global, only: itol
    use global, only: ktol
    use global, only: ftol
    use global, only: keff
    use global, only: flux
    use global, only: loss_matrix
    use global, only: prod_matrix
    use linearsolver, only: gmres
    
!---local variables

    integer :: i ! iteration index
    integer :: maxiter
    integer :: iter
    real(8) :: errork
    real(8) :: errorf
    real(8) :: err
    real(8) :: knew ! new eigenvalue
    real(8) :: kold ! old eigenvalue
    real(8), allocatable :: phinew(:) ! new flux
    real(8), allocatable :: phiold(:) ! old flux
    real(8), allocatable :: sold(:) ! old sourcie
    real(8), allocatable :: snew(:) ! new source

!---begin execution

    ! convert sparse matrix notation for linear solver
    call DS2Y(loss_matrix % n, loss_matrix % nz, loss_matrix % row, &
                         loss_matrix % col, loss_matrix % val, 0)

    call DS2Y(prod_matrix % n, prod_matrix % nz, prod_matrix % row, &
                         prod_matrix % col, prod_matrix % val, 0)

    ! allocate temporary variables
    allocate(phinew(loss_matrix % n))
    allocate(phiold(loss_matrix % n))
    allocate(sold(prod_matrix % n))
    allocate(snew(prod_matrix % n))

    ! configure random number generator
    call random_seed()

    ! guess values
    call random_number(knew)
    call random_number(phinew)
   
    ! output write values to file
    open(file='initial',unit=12)
    write(12,*) 'initial k_eff: ',knew
    write(12,*) 'initial fluxes: '
    do i=1, loss_matrix % n
      write(12,*) phinew(i)
    end do
    close(12)
 
    ! calcluate fission source
    call DSMV(prod_matrix % n, phinew, snew, prod_matrix % nz, &
                prod_matrix % row, prod_matrix % col, &
                prod_matrix % val, 0) ! F*phi sparse matrix multiplication


    ! linear solver parameters
    maxiter = 1000

    ! open file for output of keff, iter and error
    open(file='keff',unit=14)
    write(14,*) 'iteration', 'keff', 'error_keff', 'error_flux'

    ! begin power iterations
    do i = 1, 10000

      ! bank new to old
      kold = knew
      phiold = phinew
      sold = snew  

      ! normalized fission source by eigenvalue
      sold = sold/kold
      
      ! solve for new flux estimate using gmres linear solver
      call gmres(loss_matrix % n, loss_matrix % nz, loss_matrix % row, &
                 loss_matrix % col, loss_matrix % val, phinew, sold, &
                 itol, err, maxiter, iter)      

      ! calculate new fission source
      call DSMV(prod_matrix % n, phinew, snew, prod_matrix % nz, &
                prod_matrix % row, prod_matrix % col, &
                prod_matrix % val, 0)

      ! calculate new k
      knew = sum(snew)/sum(sold)   

      ! check convergence criteria     
      errork = abs(1.0_8-(kold/knew))
      errorf = maxval(abs(1.0_8-(phiold/phinew)))

      ! print out information
      print*, i, knew, errork, errorf, iter 
     
      write(14,*) i, knew, errork, errorf
      
      ! check tolerances
      if(errork < ktol .and. errorf < ftol) then
        print*, 'Eigenvalue converged normally.'
        keff = knew
        allocate(flux(prod_matrix % n))
        flux = phinew
        write(14,*) 'End k-effective value: ', keff
        return
      end if
               
    end do    
  
    ! close the opend file
    close(14)

    ! if this is reached max power iterations hit 
    print*, 'Run did not converge.'
    stop

  end subroutine power_iter

end module eigenval 
