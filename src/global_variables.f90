module global_variables
    implicit none
    

    type :: input
        integer :: nc,nv,nex,nb  !no of conduction, valence, exciton states and c+v total no of states 
    end type input
    type(input) :: sys
    

    type :: exciton
        double precision, allocatable :: eigenvalues(:)
        double precision, allocatable :: eigenvectors(:,:,:) !I,v,c
    end type exciton
    type(exciton) :: exciton_sys


    type :: biexciton
        double precision, allocatable :: eigenvalues(:)
        double precision, allocatable :: eigenvectors(:,:)
    end type biexciton
    type(biexciton) :: biexciton_sys
    
    type :: Hamiltonian
          integer :: l
          double precision, allocatable :: A(:,:),B(:,:)
          double precision, allocatable :: lambda(:,:,:,:),xi(:,:,:,:),xi_in(:,:,:,:),xi_out(:,:,:,:)
          double precision, allocatable :: lambda_e(:,:,:,:)
          double precision, allocatable :: xi_dir_lambda_e(:,:,:,:)
          double precision, allocatable :: xi_in_lambda_e(:,:,:,:)
          double precision , allocatable :: xi_out_lambda_e(:,:,:,:)
          double precision, allocatable :: E_MNRS_lambda_lambda_e(:,:,:,:)
          double precision, allocatable :: lambda_lambda_e(:,:,:,:)
    end  type Hamiltonian
    type(Hamiltonian) :: ham



          
end module global_variables