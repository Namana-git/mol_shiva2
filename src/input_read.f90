module input_read
  use global_variables
  use hdf5
  implicit none
    
  contains


  subroutine read_input()

    !reads the no of valence and conduction states from "bsemat.h5"

    
    character(len=9), parameter :: filename = "bsemat.h5"
    character(len=21), parameter :: dsetname_1 = "/bse_header/bands/nvb"
    character(len=21), parameter :: dsetname_2 = "/bse_header/bands/ncb"
    !character(len=25), parameter :: dsetname_3 = "/mf_header/crystal/celvol"
  
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id,dset2_id,dset3_id,dset4_id
    integer     ::   error
    integer(hsize_t), DIMENSION(0) :: data1_dims,data2_dims,data3_dims,data4_dims
  
  
  
    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
  
    call h5dopen_f(file_id, dsetname_1, dset1_id, error)
    call h5dread_f(dset1_id, H5T_NATIVE_INTEGER, sys%nv, data1_dims, error)
    call h5dclose_f(dset1_id, error)
  
    call h5dopen_f(file_id, dsetname_2, dset2_id, error)
    call h5dread_f(dset2_id, H5T_NATIVE_INTEGER, sys%nc, data2_dims, error)
    call h5dclose_f(dset2_id, error)
 
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  
    !debug
    print *,"conduction states,valence states",sys%nc,sys%nv
  end subroutine read_input

  subroutine read_exciton_eigenvalues()
    !reads eigenvalues from eig.dat file
    !and reads the no of excitons
    integer :: i
    integer :: error
    double precision :: x
    OPEN(UNIT=10, FILE="eig.dat", STATUS="OLD", ACTION="READ", IOSTAT=error)
    IF(error/=0) STOP "Error opening file eig.dat"
    READ(10, *, IOSTAT=error) sys%nex
    IF(error/=0) STOP "Error reading number of excitons"
    
     
    !debug
     PRINT *, "Number of excitons: ", sys%nex

    allocate(exciton_sys%eigenvalues(sys%nex))
    
    do i = 1, sys%nex
      !READ(10, *, IOSTAT=error)
      !if (error /= 0) stop "Error reading eigenvalue"
      READ(10, *, IOSTAT=error),x
      exciton_sys%eigenvalues(i) = (x/27.2114079527) 
      if (error /= 0) stop "Error reading eigenvalue"
    end do
    CLOSE(10)



  end subroutine read_exciton_eigenvalues

  subroutine read_exciton_eigenvectors()
    !reads eigenvectors from eigenvec.dat file
    integer :: i,j,k
    integer :: error

    open(UNIT=10, FILE="eigvec.dat", STATUS="OLD", ACTION="READ", IOSTAT=error)
    if(error/=0) stop "Error opening file eigvec.dat"

    allocate(exciton_sys%eigenvectors(sys%nex,sys%nv,sys%nc))

    do i = 1,sys%nex
      do j = 1,(sys%nv)
         do k = 1,(sys%nc)
            read(10, *, IOSTAT=error) exciton_sys%eigenvectors(i,j,k)
            !print*,i,j,k,exciton_sys%eigenvectors(i,j,k)
            if (error /= 0) stop "Error reading eigenvector"
         end do
      end do
    end do

    !print*,exciton_sys%eigenvectors(5,20)
    close(10)



  end subroutine read_exciton_eigenvectors
 
  subroutine read_bsemat(bsemat)
    double precision, dimension(sys%nb,sys%nb,sys%nb,sys%nb),intent(inout) :: bsemat
    character(len=9), parameter :: filename = "bsemat.h5"
    character(len=10), parameter :: name_1 = "/mats/head"
    character(len=10), parameter :: name_2 = "/mats/wing"
    character(len=10), parameter :: name_3 = "/mats/body"
    character(len=14), parameter :: name_4 = "/mats/exchange"
    double precision, DIMENSION(:,:,:,:,:,:,:), allocatable :: data1_out
  
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id,dset2_id
    integer     ::   error,iv,ivp,ic,icp,imatrix
    integer(hsize_t), DIMENSION(:), allocatable :: data1_dims
  
    allocate(data1_dims(7))
    data1_dims(1) = 2
    data1_dims(2) = sys%nb
    data1_dims(3) = sys%nb
    data1_dims(4) = sys%nb
    data1_dims(5) = sys%nb
    data1_dims(6) = 1
    data1_dims(7) = 1
    imatrix = 3
    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    if (imatrix== 1) then
       call h5dopen_f(file_id, name_1, dset1_id, error)
    elseif (imatrix== 2) then
       call h5dopen_f(file_id, name_2, dset1_id, error)
    elseif (imatrix== 3) then
       call h5dopen_f(file_id, name_3, dset1_id, error)
    elseif (imatrix== 4) then
       call h5dopen_f(file_id, name_4, dset1_id, error)
    end if
    allocate(data1_out(2,sys%nb,sys%nb,sys%nb,sys%nb,1,1))
    print*,sys%nb
    call h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, data1_out, data1_dims, error)
    call h5dclose_f(dset1_id, error)
    do icp=1,sys%nb
       do ic=1,sys%nb
          do ivp=1,sys%nb
              do iv=1,sys%nb
                  bsemat(iv,ivp,ic,icp) = data1_out(1,iv,ivp,ic,icp,1,1)
                  !print*,bse_mat(iv,ivp,ic,icp),iv,ivp,ic,icp
              enddo
          enddo
       enddo
    enddo
    !print*,bse_mat(1,1,2,1),imatrix
    call h5fclose_f(file_id, error)
    call h5close_f(error)
    deallocate(data1_out)
  
  end subroutine read_bsemat

  
    
end module input_read



  
    
