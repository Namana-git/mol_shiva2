#MKLROOT='/home/shinjan/intel/oneapi/mkl/2023.2.0'
LINK=' -L${MKLROOT}/lib -L/home/namanav/hdf5/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl  -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5'
COMP='-I/home/namanav/hdf5/include  -I${MKLROOT}/include'
cd /home/namanav/mol_shiva/src
ifort -c  global_variables.f90 $LINK $COMP 
ifort  -c  input_read.f90 $LINK $COMP
ifort -c compute_lambda_and_xi.f90 $LINK $COMP 
ifort -c compute_AB.f90 $LINK $COMP
ifort -c generalized_eigen_solver.f90 $LINK $COMP
ifort  -c  main.f90  $LINK $COMP 
ifort -g -traceback  -fp-stack-check  main.o generalized_eigen_solver.o compute_AB.o compute_lambda_and_xi.o input_read.o global_variables.o $LINK $COMP

cp a.out ../


