module compute_lambda_and_xi
    use global_variables
    implicit none
    
contains

    subroutine compute_lambda(I,J,M,N,lambda1)

        integer,intent(in) :: I,J,M,N
        integer :: c,cp,v,vp
        double precision,intent(inout) :: lambda1
        double precision :: sum1
        double precision :: AI_cv,AJ_cpvp,AM_cvp,AN_cpv
        
        sum1 = 0
    
        !call cpu_time(loop_start)
        do cp=1,sys%nc ! c'
            do vp = 1,sys%nv !v'
                do c = 1,sys%nc !c 
                   do v = 1, sys%nv !v
                        AI_cv = exciton_sys%eigenvectors(I,v,c)
                        AJ_cpvp = exciton_sys%eigenvectors(J,vp,cp)
                        AM_cvp = exciton_sys%eigenvectors(M,vp,c)
                        AN_cpv = exciton_sys%eigenvectors(N,v,cp)
                        sum1 = sum1 + (AI_cv*AJ_cpvp *AM_cvp*AN_cpv)
                         
                                    ! use complex conjugate for complex eigen vectors
                       ! sum2 = sum2 + (exciton_sys%eigenvectors(J,x,y)*exciton_sys%eigenvectors(I,s,t) &
                        !                *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t))
                        !print*,I,J,M,N,c,cp,v,vp,exciton_sys%eigenvectors(I,v,c),exciton_sys%eigenvectors(J,vp,cp),exciton_sys%eigenvectors(M,vp,c),exciton_sys%eigenvectors(N,v,cp)
                                        
                    end do
                end do
            end do
        end do
        lambda1 = sum1
        !print*,I,J,M,N,lambda1
        

    end subroutine compute_lambda


    subroutine compute_xi_ee(Iex,Jex,Mex,Nex,bsemat,xi_ee1)
        !\xi_ee(I,J,M,N) = \sum_{i,j,c,c',v,v'}  A_{M}^{i,v}A_{N}^{j,v'}bsemat(cp,j,i,c)A_{I)^{c,v}A_{J}^{c'v'}
        integer,intent(in) :: Iex,Jex,Mex,Nex
        integer :: c,v,cp,vp,i,j
        integer :: nv
        double precision,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
        double precision,intent(inout) :: xi_ee1
        double precision :: sum1,sum2,AI_cv,AJ_cpvp,AM_iv,AN_jvp
        nv = sys%nv
        AI_cv = 0.0
        AJ_cpvp = 0.0
        AM_iv = 0.0
        AN_jvp = 0.0
        sum1 = 0.0
        sum2 = 0.0
        do i = 1,sys%nc
            do j = 1,sys%nc
                do c = 1,sys%nc
                    do cp = 1,sys%nc
                        do v = 1,sys%nv
                            do vp = 1,sys%nv
                                !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
                                !bsemat()   = <ij|W|ccp>
                                AI_cv= exciton_sys%eigenvectors(Iex,v,c)
                                AJ_cpvp = exciton_sys%eigenvectors(Jex,vp,cp)
                                AM_iv = exciton_sys%eigenvectors(Mex,v,i)
                                AN_jvp = exciton_sys%eigenvectors(Nex,vp,j)
                                !AI_cpvp = exciton_sys%eigenvectors(Iex,vp,cp)
                                !AJ_cv = exciton_sys%eigenvectors(Jex,v,c)
                                sum1 = sum1 + AM_iv*AN_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AI_cv*AJ_cpvp 
                                !sum2 = sum2 + AM_iv*AN_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AJ_cv*AI_cpvp
                                !print*,sum1,"sum1",bsemat(cp+nv,j+nv,i+nv,c+nv),"bsemat",AI_cv,"AI_cv",AJ_cpvp,"AJ_cpvp",AM_iv,"AM_iv",AN_jvp,"AN_jvp"
                                !print*,"i,j,c,cp,v,vp",i,j,c,cp,v,vp
                            end do
                        end do
                    end do
                end do
            end do
        end do
        xi_ee1 = sum1
        !xi_ee2 = sum2

    end subroutine compute_xi_ee


    subroutine compute_xi_hh(Iex,Jex,Mex,Nex,bsemat,xi_hh1)
        !\xi_hh(I,J,M,N) = \sum_{\alpha,\beta,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{c',beta}bsemat(cp,j,i,c)A_{I)^{c,v}A_{J}^{c'v'}
        ! bsemat(eta,beta,alpha,gama)) = <alpha,beta|W|gama,eta>
        ! bsemat(beta,vp,v,alpha)) = <vvp|W|alpha,beta> 
        
        integer,intent(in) :: Iex,Jex,Mex,Nex
        integer :: c,v,cp,vp,alpha,beta
        double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
        double precision,intent(inout) :: xi_hh1
        double precision :: sum1,sum2,AI_cv,AJ_cpvp,AM_calpha,AN_cpbeta,AJ_cv,AI_cpvp
        sum1 = 0
        sum2 = 0
        do alpha = 1,sys%nv
            do beta = 1,sys%nv
                do c = 1,sys%nc
                    do cp = 1,sys%nc
                        do v = 1,sys%nv
                            do vp = 1,sys%nv
                                !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
                                !bsemat()   = <ij|W|ccp>
                                AI_cv= exciton_sys%eigenvectors(Iex,v,c)
                                AJ_cpvp = exciton_sys%eigenvectors(Jex,vp,cp)
                                AM_calpha = exciton_sys%eigenvectors(Mex,alpha,c)
                                AN_cpbeta = exciton_sys%eigenvectors(Nex,beta,cp)
                                sum1 = sum1 + AM_calpha*AN_cpbeta*bsemat(beta,vp,v,alpha)*AI_cv*AJ_cpvp 
                                
                            end do
                        end do
                    end do
                end do
            end do 
        end do
        xi_hh1 = sum1
                                  
    end subroutine compute_xi_hh

    subroutine compute_xi_in_ehd(Iex,Jex,Mex,Nex,bsemat,xi_in_ehd1)
        !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
        !\xi_out_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{i,v'}bsemat(alpha,v,i+nv,c'+nv)A_{I)^{c,v}A_{J}^{c'v'}
        !bsemat(alpha,v',i+nv,c+nv) = <iv'|W|c,alpha>
        !bsemat(alpha,v,i+nv,c'+nv)= <iv|W|c',alpha>

        integer,intent(in) :: Iex,Jex,Mex,Nex
        integer :: c,v,cp,vp,i,alpha,nv
        double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
        double precision,intent(out) :: xi_in_ehd1 !,xi_in_ehd2,xi_out_ehd1,xi_out_ehd2
        double precision :: sum1,sum2,sum3,sum4,AI_cv,AJ_cpvp,AM_iv,AN_cpalpha,AJ_cv,AI_cpvp,AM_calpha,AN_ivp
        sum1 = 0
        sum2 = 0
        sum3 = 0
        sum4 = 0
        nv = sys%nv
        do i = 1,sys%nc  !
            do alpha = 1,sys%nv
                do c = 1,sys%nc
                    do cp = 1,sys%nc
                        do v = 1,sys%nv
                            do vp = 1,sys%nv
                                !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
                                !bsemat()   = <ij|W|ccp>
                                AI_cv= exciton_sys%eigenvectors(Iex,v,c)
                                AJ_cpvp = exciton_sys%eigenvectors(Jex,vp,cp)
                                AM_iv = exciton_sys%eigenvectors(Mex,v,i)
                                AN_cpalpha = exciton_sys%eigenvectors(Nex,alpha,cp)
                                !AI_cpvp = exciton_sys%eigenvectors(Iex,vp,cp)
                                !AJ_cv = exciton_sys%eigenvectors(Jex,v,c)
                                !AM_calpha = exciton_sys%eigenvectors(Mex,alpha,c)
                                !AN_ivp = exciton_sys%eigenvectors(Nex,vp,i)
                                sum1 = sum1 + AM_iv*AN_cpalpha*bsemat(alpha,vp,i+nv,c+nv)*AI_cv*AJ_cpvp
                                !sum2 = sum2 + AM_iv*AN_calpha*bsemat(alpha,vp,i+nv,c+nv)*AJ_cv*AI_cpvp
                                !sum3 = sum3 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
                                !sum4 = sum4 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AJ_cv*AI_cpvp
                                print*,sum1,"sum1"       
                            end do
                        end do
                    end do
                end do
            end do
        end do
        xi_in_ehd1 = sum1
        !xi_in_ehd2 = sum2
        !xi_out_ehd1 = sum3
        !xi_out_ehd2 = sum4
    end subroutine compute_xi_in_ehd
    
    subroutine compute_xi_out_ehd(Iex,Jex,Mex,Nex,bsemat,xi_out_ehd1)
        !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
        !\xi_out_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{i,v'}bsemat(alpha,v,i+nv,c'+nv)A_{I)^{c,v}A_{J}^{c'v'}
        !bsemat(alpha,v',i+nv,c+nv) = <iv'|W|c,alpha>
        !bsemat(alpha,v,i+nv,c'+nv)= <iv|W|c',alpha>

        integer,intent(in) :: Iex,Jex,Mex,Nex
        integer :: c,v,cp,vp,i,alpha,nv
        double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
        double precision,intent(out) :: xi_out_ehd1 !,xi_in_ehd2,xi_out_ehd1,xi_out_ehd2
        double precision :: sum1,sum2,sum3,sum4,AI_cv,AJ_cpvp,AM_iv,AN_calpha,AJ_cv,AI_cpvp,AM_calpha,AN_ivp
        sum1 = 0
        sum2 = 0
        sum3 = 0
        sum4 = 0
        nv = sys%nv
        do i = 1,sys%nc  !
            do alpha = 1,sys%nv
                do c = 1,sys%nc
                    do cp = 1,sys%nc
                        do v = 1,sys%nv
                            do vp = 1,sys%nv
                                !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
                                !bsemat()   = <ij|W|ccp>
                                AI_cv= exciton_sys%eigenvectors(Iex,v,c)
                                AJ_cpvp = exciton_sys%eigenvectors(Jex,vp,cp)
                                !AM_iv = exciton_sys%eigenvectors(Mex,v,i)
                                !AN_calpha = exciton_sys%eigenvectors(Nex,alpha,c)
                                !AI_cpvp = exciton_sys%eigenvectors(Iex,vp,cp)
                                !AJ_cv = exciton_sys%eigenvectors(Jex,v,c)
                                AM_calpha = exciton_sys%eigenvectors(Mex,alpha,c)
                                AN_ivp = exciton_sys%eigenvectors(Nex,vp,i)
                                !sum1 = sum1 + AM_iv*AN_calpha*bsemat(alpha,vp,i+nv,c+nv)*AI_cv*AJ_cpvp
                                !sum2 = sum2 + AM_iv*AN_calpha*bsemat(alpha,vp,i+nv,c+nv)*AJ_cv*AI_cpvp
                                sum3 = sum3 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
                                !sum4 = sum4 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AJ_cv*AI_cpvp         
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !xi_in_ehd1 = sum1
        !xi_in_ehd2 = sum2
        xi_out_ehd1 = sum3
        !xi_out_ehd2 = sum4
    end subroutine compute_xi_out_ehd


    subroutine compute_xi_in()

        implicit none
        integer :: I,J,M,N,R,S

        do M = 1,sys%nex
            do N = 1,sys%nex   
                do I= 1,sys%nex
                    do J = 1,sys%nex 
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                                ham%xi_in(I,J,M,N) = ham%xi_in(I,J,M,N) + ham%lambda(R,S,M,N)*ham%xi(I,J,R,S)
                                !print*,"xi_in",I,J,M,N,ham%xi_in(I,J,M,N)
                            end do
                        end do
                    end do
                end do
            end do
        end do



    end subroutine compute_xi_in


    subroutine compute_xi_out()
        implicit none
        integer :: I,J,M,N,R,S
        do M = 1,sys%nex
            do N = 1,sys%nex
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                                ham%xi_out(I,J,M,N) = ham%xi_out(I,J,M,N) + ham%xi(R,S,M,N)*ham%lambda(I,J,R,S)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    subroutine compute_xi_dir_lambda_e()
        implicit none
        integer :: I,J,M,N,R,S
        do M = 1,sys%nex
            do N = 1,sys%nex
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                                ham%xi_dir_lambda_e(I,J,M,N) = ham%xi_dir_lambda_e(I,J,M,N) &
                                                              + ham%xi(R,S,M,N)*ham%lambda_e(I,J,R,S)
                            end do
                        end do
                    end do 
                end do
            end do 
        end do
    end subroutine compute_xi_dir_lambda_e

    subroutine compute_xi_in_lambda_e()
        implicit none
        integer :: I,J,M,N,R,S
        do M = 1,sys%nex
            do N = 1,sys%nex
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                               ham%xi_in_lambda_e(I,J,M,N) = ham%xi_in_lambda_e(I,J,M,N) &
                                                             + ham%xi_in(R,S,M,N)*ham%lambda_e(I,J,R,S) 
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine compute_xi_in_lambda_e()

    subroutine compute_xi_out_lambda_e()
        implicit none
        integer :: I,J,M,N,R,S
        do M = 1,sys%nex
            do N = 1,sys%nex
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                                ham%xi_out_lambda_e(I,J,M,N) = ham%xi_out_lambda_e(I,J,M,N) &
                                                               + ham%xi_out(R,S,M,N)*ham%lambda_e(I,J,R,S)  
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine compute_xi_out_lambda_e

    subroutine compute_E_MNRS_lambda_lambda_e()
     integer :: I,J,M,N,R,S
        do M = 1,sys%nex    
            do N = 1,sys%nex
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                               ham%E_MNRS_lambda_lambda_e(I,J,M,N) = ham%E_MNRS_lambda_lambda_e(I,J,M,N) &
                                    + ((exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N) + exciton_sys%eigenvalues(R) + exciton_sys%eigenvalues(S))&
                                    * ham%lambda(R,S,M,N)*ham%lambda_e(I,J,R,S))
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine compute_E_MNRS_lambda_lambda_e


    subroutine compute_lambda_lambda_e

         integer :: I,J,M,N,R,S

         do M = 1,sys%nex
            do N = 1,sys%nex
                do I = 1,sys%nex
                    do J = 1,sys%nex
                        do R = 1,sys%nex
                            do S = 1,sys%nex
                                ham%lambda_lambda_e(I,J,M,N) = ham%lambda_lambda_e(I,J,M,N) &
                                                              + ham%lambda(R,S,M,N)*ham%lambda_e(I,J,R,S)
                            end do
                        end do
                    end do
                end do
            end do
        end do


    end subroutine compute_lambda_lambda_e


    
end module compute_lambda_and_xi