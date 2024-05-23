module solve_generalized_eigenvalue_problem
   use global_variables
   use compute_lambda_and_xi
   use compute_AB

   !use lapacks generalized eigenvalue solver dggevx, with A = ham%A and B = ham%B and m=ham%l

   implicit none


   contains


   subroutine solve_generalized_eigenvalue_AB()
    
    integer :: ilo,ihi
    double precision,allocatable,dimension(:) :: alphar,alphai,beta
    double precision,allocatable,dimension(:,:) :: vl,vr
    double precision,allocatable,dimension(:) :: lscale,rscale
    double precision,allocatable,dimension(:) :: rconde,rcondv
    double precision,allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    logical,allocatable,dimension(:) :: bwork
    integer :: lwork,info
    integer :: i,j
    double precision :: abnrm,bbnrm

  
    lwork = 8*ham%l
    allocate(alphar(ham%l),alphai(ham%l),beta(ham%l))
    allocate(vl(ham%l,ham%l),vr(ham%l,ham%l))
    allocate(lscale(ham%l),rscale(ham%l))
    allocate(rconde(ham%l),rcondv(ham%l))
    allocate(work(lwork))
    allocate(iwork(ham%l+6))
    allocate(bwork(ham%l))
    alphar = 0.0
    alphai = 0.0
    beta = 0.0
    vl = 0.0
    vr = 0.0
    BWORK = .false.
    lwork = 8*ham%l
    call dggevx('N', 'N', 'V', 'N', ham%l, ham%A, ham%l, ham%B, ham%l, &
    alphar, alphai, beta, vl, ham%l, vr, ham%l, &
    ilo,ihi, LSCALE, RSCALE, abnrm, bbnrm, &
    rconde, rcondv, work, lwork, iwork, bwork, info)


    do i =  1, ham%l
      !print*,alphar(i),alphai(i),beta(i)
      print*,((alphar(i)/beta(i))*27.2114079527) 
      
    end do
    
            
   end subroutine solve_generalized_eigenvalue_AB

end module solve_generalized_eigenvalue_problem