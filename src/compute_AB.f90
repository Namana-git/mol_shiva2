module compute_AB
    use global_variables
    use compute_lambda_and_xi
    use input_read
    use, intrinsic :: iso_fortran_env
    !use blas_interface
    implicit none



    ! calculates A and B from \lambda and \xi matrices
    !(E_{M} + E_{N} -\epsilon_{m})c_{MN} + \sum_{IJ} \xi^{tot} (M,N,I,J) = 0
    !\xi^{tot}(M,N,I,J) = \xi(M,N,I,J) + \xi^{in}(M,N,I,J) + (E_{I}+E_{J} - \epsilon_{m})\lambda (M,N,I,J)
    !A(x(MN),y(IJ)) = (E_{M} + E_{N})\delta_{MN,IJ) +  \xi(M,N,I,J) + \xi^{in}(M,N,I,J) + (E_{I}+E_{J})\lambda (M,N,I,J) 
    !B(x(MN,Y(IJ))) = \delta_{MN,IJ) + \lambda (M,N,I,J)     
contains

   subroutine compute_xi_and_lambda()

       implicit none
       integer :: I,J,M,N
       integer :: c1,c2
       double precision :: lambda1,xi_ee1,xi_ee2,xi_hh1,xi_hh2,xi_in_ehd1,xi_out_ehd1
       double precision,allocatable,dimension(:,:,:,:) :: bsemat
    
       ham%lambda = 0.0 
       lambda1 = 0.0    
       sys%nb = sys%nc+sys%nv
       allocate(ham%lambda(sys%nex,sys%nex,sys%nex,sys%nex))
       allocate(ham%lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
       allocate(ham%xi(sys%nex,sys%nex,sys%nex,sys%nex))
       ham%lambda = 0.0
       ham%lambda_e = 0.0
       ham%xi = 0.0
    
       allocate(bsemat(sys%nb,sys%nb,sys%nb,sys%nb))
       call read_bsemat(bsemat)
      
        c1 = 0 
        c2 = 0
    
        do M = 1,sys%nex
          do N = 1,sys%nex
              do I = 1,sys%nex
                  do J = 1,sys%nex
                    
                      call compute_lambda(I,J,M,N,lambda1)
                      ham%lambda(I,J,M,N) = lambda1
                      call compute_lambda(I,J,N,M,lambda1)
                      ham%lambda_e(I,J,M,N) = lambda1
                      call compute_xi_ee(I,J,M,N,bsemat,xi_ee1)
                      print*,"xi_ee",I,J,M,N,xi_ee1
                      call compute_xi_hh(I,J,M,N,bsemat,xi_hh1)
                      print*,"xi_hh",I,J,M,N,xi_hh1
                      call compute_xi_in_ehd(I,J,M,N,bsemat,xi_in_ehd1)
                      print*,"xi_in_eh",I,J,M,N,xi_in_ehd1
                      call compute_xi_out_ehd(I,J,M,N,bsemat,xi_out_ehd1)
                      print*,"xi_out_eh",I,J,M,N,xi_out_ehd1
                      ham%xi(I,J,M,N) = xi_ee1 + xi_hh1 - xi_in_ehd1 - xi_out_ehd1

                      !call compute_xi_ee(I,J,M,N,bsemat,xi_ee1,xi_ee2)
                    
                  end do
              end do
          end do
        end do
    
    
    end subroutine compute_xi_and_lambda

   

    
    subroutine compute_A_and_B()

        implicit none
        integer :: I,J,M,N
        integer :: c1,c2

        allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
        ham%xi_in = 0.0
        call compute_xi_in()
        c1 = 0
        c2 = 0
        allocate(ham%A(ham%l,ham%l))
        allocate(ham%B(ham%l,ham%l))
        ham%A = 0.0
        ham%B = 0.0
        ham%xi = 0.0
        ham%xi_in=0.0
            
        do M = 1,sys%nex
            do N = M,sys%nex
                c2 = 0
                c1 = c1 + 1
                do I = 1,sys%nex
                    do J = I,sys%nex
                        c2 = c2 + 1
                        if(I==M .and. J==N) then
                            if(I==J) then
                               ham%A(c1,c2) = 2*(exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
                               ham%B(c1,c2) = 2
                            end if
                            if(I.ne.J) then
                               ham%A(c1,c2) = (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
                               ham%B(c1,c2) = 1
                            end if
                        end if
                        print*,I,J,M,N,ham%lambda(I,J,M,N),ham%lambda(I,J,N,M)
                        ham%A(c1,c2) = ham%A(c1,c2) + ham%xi(I,J,M,N) &
                                       - ham%xi_in(I,J,M,N) - ham%xi_in(I,J,N,M) &
                                       -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N)+ham%lambda(I,J,N,M)))                    
                        ham%B(c1,c2) = ham%B(c1,c2) - ham%lambda(I,J,M,N) -ham%lambda(I,J,N,M)
                        !print*,c1,c2,ham%A(c1,c2),ham%B(c1,c2)
                        
                        
                    end do
                end do
            end do
        end do
        
    end subroutine compute_A_and_B

    subroutine compute_A_and_B_double_basis()

        implicit none
        integer :: I,J,M,N
        integer :: c1,c2
        

        allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%xi_dir_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%xi_in_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%xi_out(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%xi_out_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%E_MNRS_lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
        ham%xi_in = 0.0
        ham%xi_out = 0.0
        ham%xi_dir_lambda_e = 0.0
        ham%xi_in_lambda_e = 0.0
        ham%xi_out_lambda_e = 0.0
        ham%E_MNRS_lambda_lambda_e = 0.0
        ham%lambda_lambda_e = 0.0
        call compute_xi_in()
        c1 = 0
        c2 = 0
        allocate(ham%A(ham%l,ham%l))
        allocate(ham%B(ham%l,ham%l))
        ham%A = 0.0
        ham%B = 0.0
        !ham%xi = 0.0
        !ham%xi_in=0.0
            
        do M = 1,sys%nex
            do N = 1,sys%nex
                c2 = 0
                c1 = c1 + 1
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        c2 = c2 + 1
                        if(I==M .and. J==N) then
                            ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
                            ham%B(c1,c2) = ham%B(c1,c2) + 1
                        end if
                        if(I==N .and. J==M) then
                            ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
                            ham%B(c1,c2) = ham%B(c1,c2) + 1
                        end if
                        !print*,I,J,M,N,ham%lambda(I,J,M,N),ham%lambda(I,J,N,M)
                        ham%A(c1,c2) = ham%A(c1,c2) &
                                       + ham%xi(I,J,M,N) + ham%xi(I,J,N,M) &
                                       - ham%xi_in(I,J,M,N) - ham%xi_in(I,J,N,M) &
                                       -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N)+ham%lambda(I,J,N,M)))                    
                        ham%B(c1,c2) = ham%B(c1,c2) - ham%lambda(I,J,M,N) -ham%lambda(I,J,N,M)
                        print*,c1,c2,ham%A(c1,c2),ham%B(c1,c2)
                        
                        
                    end do
                end do
            end do
        end do
        
    end subroutine compute_A_and_B_double_basis


    

    subroutine compute_A_and_B_vectorized()

        implicit none
        integer :: I,J,M,N
        integer :: x,c,s,t,z,y
        double precision :: lambda,sum, vec_start,vec_end,loop_start,loop_end
        double precision,allocatable,dimension(:,:) :: left_eigp,right_eigp 
        double precision, allocatable,dimension(:,:) :: A_I,A_J
        double precision,allocatable,dimension(:):: B_I,B_J
        
        ! \lambda(M,N,I,J) = \sum_{c,c',v,v'} left_eigp(c,c',v,v')* right_eigp(c,c',v,v') 
        != inner product(left_eigp,right_eigp)
        !right_eigp = outer prduct(A^{I}_cv,A^{J}_c'v')
        
        allocate(right_eigp(sys%nc*sys%nv,sys%nc*sys%nv))
        !allocate(right_eigp(sys%nc,sys%nv,sys%nc,sys%nv))
       !allocate(A_I(sys%nc,sys%nv))
        !allocate(A_J(sys%nc,sys%nv))
        allocate(B_I(sys%nc*sys%nv))
        !allocate(B_J(sys%nc*sys%nv))
        !debug
        z =sys%nv*sys%nc
        sys%nex =1
        do M = 1,sys%nex
            do N = M,sys%nex
                do I = 1,sys%nex
                    do J = I,sys%nex
                        !A_I(:,:) = exciton_sys%eigenvectors(I,:,:)
                        !A_J(:,:) = exciton_sys%eigenvectors(J,:,:)
                        !print*,A_I,"A_I"
                        !B_I = reshape(exciton_sys%eigenvectors(I,:,:),[sys%nv*sys%nc])
                        !B_J = reshape(exciton_sys%eigenvectors(J,:,:),[sys%nv*sys%nc])
                        !print*,B_I,"B_I"
                        
                        !result = matmul(reshape(a, [n,1]), reshape(b, [1,n]))
                        !print*,size(matmul(reshape(reshape(exciton_sys%eigenvectors(I,:,:),[sys%nv*sys%nc]),[x,1]),reshape(B_J,[1,x])))
                        call cpu_time(vec_start)
                        right_eigp = matmul(reshape(reshape(exciton_sys%eigenvectors(I,:,:),[sys%nv*sys%nc]),[z,1]), &
                          reshape(reshape(exciton_sys%eigenvectors(J,:,:),[sys%nv*sys%nc]),[1,z]))
                        print*,size(right_eigp)
                        !print*,reshape(B_I,[x,1]),"right_eigp"
                        
                                            !reshape(exciton_sys%eigenvectors(J,:,:),(/1,sys%nv*sys%nc/)) )
                        !left_eigp = matmul(reshape(exciton_sys%eigenvectors(M,:,:),(/sys%nv*sys%nc,1/)), &
                                      !      reshape(exciton_sys%eigenvectors(N,:,:),(/1,sys%nv*sys%nc/)))
                        !right_eigp = reshape(right_eigp,(/sys%nc*sys%nv*sys%nc*sys%nv/))                    
                        !left_eigp = reshape(left_eigp,(/sys%nc*sys%nv*sys%nc*sys%nv/))
                        lambda = dot_product(reshape(reshape(right_eigp,[sys%nv,sys%nc,sys%nv,sys%nc],order=[1,2,3,4]),[sys%nv*sys%nc*sys%nv*sys%nc]), &
                                     reshape(reshape(right_eigp,[sys%nv,sys%nc,sys%nv,sys%nc],order=[3,2,1,4]),[sys%nv*sys%nc*sys%nv*sys%nc]))
                        call cpu_time(vec_end)
                        print*,lambda,"lambda"
                        B_I = reshape(exciton_sys%eigenvectors(J,:,:),[sys%nv*sys%nc]) 
                        sum =0
                        call cpu_time(loop_start)
                        do t=1,sys%nc ! c'
                            do s = 1,sys%nv !v'
                                 do y = 1,sys%nc !c 
                                   do x = 1, sys%nv !v
                                        sum = sum + (exciton_sys%eigenvectors(I,x,y)*exciton_sys%eigenvectors(J,s,t) &
                                                    *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t)) 
                                        !print*,exciton_sys%eigenvectors(I,x,t)*exciton_sys%eigenvectors(J,s,y),right_eigp(s,t,x,y)
                                        !print*,s,t,exciton_sys%eigenvectors(I,s,t),B_I(c)
                                    end do
                                end do
                            end do
                        end do
                        call cpu_time(loop_end)
                        print*,sum,"sum"
                    end do
                end do
            end do
        end do
        
        print*,loop_end - loop_start,"loop"
        print*,vec_end - vec_start,"vec"
    end subroutine compute_A_and_B_vectorized
    
end module compute_AB