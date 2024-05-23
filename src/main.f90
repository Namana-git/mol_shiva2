program mol_shiva
    use compute_AB
    use global_variables
    use input_read
    use solve_generalized_eigenvalue_problem
    implicit none
     
    ! read no of valence, conduction,and exciton states included.
    call read_input()
    ! read eigenvalues and eigenvectors
    call read_exciton_eigenvalues()
    call read_exciton_eigenvectors()
    sys%nex =4
    !ham%l = sys%nex*(sys%nex + 1)/2
    ham%l = sys%nex*sys%nex
    !ham%l = 4
    !sys%nv = 20
    !sys%nc = 40
    
    !compute A and B Hamiltonian from lambda and xi matrices
    call compute_xi_and_lambda()
    call compute_A_and_B_double_basis() 

    !solve Ax = \lambdaBx 
    call solve_generalized_eigenvalue_AB()
    
    
    
end program  

