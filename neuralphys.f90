!***************************************************************************
!
! Name: neuralphys
!
! Language: FORTRAN                           Type - MODULE
!
! Version: 1.0          Date: 03-25-20         
!     
!
! **************************************************************
!
! Module contains all subroutines require to initialize and 
! calculate ensemble of NNs for GFS model physics.
! 
! Originally created by Alex Belochitski for physics emulation
! Modified by Tse-Chun Chen for bias correction
! **************************************************************
!
       module neuralphys

        use sim_nc_mod,        only: open_ncfile,  &
                                     close_ncfile, &
                                     get_var2_real,&
                                     get_var1_real,&
                                     handle_err
        use module_radlw_parameters, only: sfcflw_type,  &
                                           topflw_type
        use module_radsw_parameters, only: sfcfsw_type,  &
                                           topfsw_type
        use machine,           only: kind_phys

        use omp_lib
                 
        implicit none 
#include <netcdf.inc>

        private
 
        public :: init_nn, eval_nn
                  
! Number of members in the NN ensemble
        integer, parameter     :: nn_num_of_members = 4
        real(kind=4), parameter     :: pi = 3.1415927

! Files containing NN weights and biases
        character(*), parameter::nn_file_name(nn_num_of_members)= (/& 
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_t.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_u.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_v.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_q.nc' /)
             
        integer,    allocatable :: nn_sizes(:,:) ! nn_num_of_members x num_layers
        integer num_layers,num_layersm1
 ! Internal types and variables
        type nndata
            real(kind=4), allocatable :: w(:,:)
            real(kind=4), allocatable :: b(:)
        end type nndata

        type nndata2
            real(kind=4) w0(4096,523),b0(4096),w1(4096,4096),b1(4096),w2(127,4096),b2(127)
        end type nndata2
        !type nndata_1d
        !   real, allocatable :: a(:)
        !end type nndata_1d

        !type nndata_2d
        !   real, allocatable :: a(:,:)
        !end type nndata_2d
        
        !type nndata_3d
        !   real, allocatable :: a(:,:,:)
        !end type nndata_3d

! Get number of NN layers from nn_t.nc attributes
        !call open_ncfile(nn_file_name(1), 104)
        !STATUS = NF_INQ_ATTLEN (104, NF_GLOBAL, 'nn_sizes', num_layers) ! get number of NN layers
        !IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
        !call close_ncfile(104)
        
        !allocate(nn_sizes(nn_num_of_members,num_layers))        

! Assume all NNs to have the same number of layers (number of neurons can be different).
        !type(nndata), allocatable :: nns(:,:)
        !type(nndata) :: nns(4,3)
        type(nndata2) nns(4)
! NN hidden weights
        !type(nndata_3d) :: nn_whid(nn_num_of_members)
! NN input and output weights, hidden biases
        !type(nndata_2d) :: nn_win(nn_num_of_members),nn_wout(nn_num_of_members),nn_bhid(nn_num_of_members)
! NN output biases
        !type(nndata_1d) :: nn_bout(nn_num_of_members)

      contains
!
! Initialize TMP NNs
        
        subroutine init_nn(me) 
!
! --- This subroutine initializes NN ensemble, i.e. reads NNs coefficients
!
!   
          integer, intent(in) ::  me
          integer iin,ihid,iout,member,layer,l0,l1,STATUS,ncid !,num_layers,num_layersm1
          character(len=2) :: varname 
          
          if (me == 0) print*,'Module NEURALPHYS: Number of NNs:', &
                               nn_num_of_members

! Get number of NN layers from nn_t.nc attributes
          call open_ncfile(nn_file_name(1), ncid)
          STATUS = NF_INQ_ATTLEN (ncid, NF_GLOBAL, 'nn_sizes', num_layers) ! get number of NN layers
          IF (STATUS .NE. NF_NOERR) CALL handle_err("inq_attlen",STATUS)
          call close_ncfile(ncid)

          num_layersm1 = num_layers - 1

          allocate(nn_sizes(nn_num_of_members,num_layers))

! Assume all NNs to have the same number of layers (number of neurons can be
! different).
          !allocate(nns(nn_num_of_members,num_layers - 1))

! Load NNs weights and biases

          do member=1,nn_num_of_members
              call open_ncfile(nn_file_name(member), ncid)
        
              STATUS = NF_GET_ATT_INT (ncid, NF_GLOBAL, 'nn_sizes', nn_sizes(member,:))
              IF (STATUS .NE. NF_NOERR) CALL handle_err("get_att_int",STATUS)
              
              !do layer=1,num_layers-1
              !    l0 = nn_sizes(member,layer)
              !    l1 = nn_sizes(member,layer+1)
                  !allocate(nns(member,layer)%w(l1,l0))
                  !allocate(nns(member,layer)%b(l1))
                  !nns(member,layer)%w(:,:) = 0.0
                  !nns(member,layer)%b(:)   = 0.0
                  nns(member)%w0 = 0.
                  nns(member)%w1 = 0.
                  nns(member)%w2 = 0.
                  nns(member)%b0 = 0.
                  nns(member)%b1 = 0.
                  nns(member)%b2 = 0.
              !    write(varname, "(A1,I1)") "w",layer-1
              !    call get_var2_real(ncid, varname, l1, l0, nns(member,layer)%w(:,:))
                  call get_var2_real(ncid, 'w0', nn_sizes(member,2), nn_sizes(member,1), nns(member)%w0) 
                  call get_var2_real(ncid, 'w1', nn_sizes(member,3), nn_sizes(member,2), nns(member)%w1)
                  call get_var2_real(ncid, 'w2', nn_sizes(member,4), nn_sizes(member,3),  nns(member)%w2)
                  !call get_var2_double(ncid, varname, l1, l0,nns(member,layer)%w(:,:))
              !    write(varname, "(A1,I1)") "b",layer-1
                  !call get_var1_real(ncid, varname, l1, nns(member,layer)%b(:))
                  call get_var1_real(ncid, 'b0', nn_sizes(member,2), nns(member)%b0)
                  call get_var1_real(ncid, 'b1', nn_sizes(member,3), nns(member)%b1)
                  call get_var1_real(ncid, 'b2', nn_sizes(member,4),  nns(member)%b2)
                  !call get_var1_double(ncid, varname, l1, nns(member,layer)%b(:))
                  if (me == 0)  print*,'Module NEURALPHYS: NN w:',member,shape(nns(member)%w0),shape(nns(member)%w1),shape(nns(member)%w2)
                  if (me == 0)  print*,'Module NEURALPHYS: NN b:',member,shape(nns(member)%b0),shape(nns(member)%b1),shape(nns(member)%b2)
              !end do
              call close_ncfile(ncid)
              
              if (me == 0)  print*,'Module NEURALPHYS: NN File Loaded: ', & 
                                   nn_file_name(member)
          end do
        if (me == 0)  print*,'Module NEURALPHYS: NN SIZES: ', nn_sizes, num_layers
        end subroutine init_nn

        subroutine  eval_nn(                                  & 
!  Inputs:
! Prognostic variables:
     &     pgr,ugrs,vgrs,tgrs,shm,  &
! Surface variables:
     &     sfcflw,sfcfsw,topflw,topfsw,slmsk,        &
! Metavariables:
     &     hour,doy,glon,glat,dt,  &
! Outputs:
     &     gu0,gv0,gt0,oshm)

! Inputs
          real(kind=4), intent(in) :: slmsk,glat,glon,hour,doy,dt
          type (sfcflw_type), intent(in) :: sfcflw
          type (sfcfsw_type), intent(in) :: sfcfsw
          type (topflw_type), intent(in) :: topflw
          type (topfsw_type), intent(in) :: topfsw
          real(kind=8), intent(in):: ugrs(127),vgrs(127),tgrs(127),shm(127),pgr
! Outputs     
          real(kind=8), intent(out):: gu0(127),gv0(127),gt0(127),oshm(127)
!
! Local variables
          real(kind=4)  nn_input_vector(nn_sizes(1,1)), nn_output_vector(nn_sizes(1,4)*4) 
! Create NN input vector:
          nn_input_vector(1:127)  = tgrs(127:1:-1)        ! temperature (inverted in vertical)
          nn_input_vector(128)    = log(pgr)              ! log(surface pressure)
          nn_input_vector(129:255)= ugrs(127:1:-1)        ! u-wind (inverted in vertical)
          nn_input_vector(256:382)= vgrs(127:1:-1)        ! v-wind (inverted in vertical)
          nn_input_vector(383:509)= shm(127:1:-1)         ! specific humidity (inverted in vertical)
          nn_input_vector(510)    = real(sfcflw%dnfx0,4)  ! surface longwave downward flux
          nn_input_vector(511)    = real(sfcfsw%dnfx0,4)  ! surface shortwave downward flux
          nn_input_vector(512)    = real(sfcflw%upfx0,4)  ! surface longwave upward flux
          nn_input_vector(513)    = real(topflw%upfx0,4)  ! TOA longwave upward flux
          nn_input_vector(514)    = real(sfcfsw%upfx0,4)  ! surface shortwave upward flux
          nn_input_vector(515)    = real(topfsw%upfx0,4)  ! TOA shortwave upward flux
          nn_input_vector(516)    = slmsk                 ! land-sea-ice mask
          nn_input_vector(517)    = glat*180.0/pi         ! latitude in degree
          nn_input_vector(518)    = sin(glon)             ! sin(longitude)
          nn_input_vector(519)    = cos(glon)             ! cos(longitude)
          nn_input_vector(520)    = sin(2.0* pi * hour/24.0) ! sin(hour of the day)
          nn_input_vector(521)    = sin(2.0* pi * doy/365.0) ! sin(day of the year)
          nn_input_vector(522)    = cos(2.0* pi * hour/24.0) ! cos(hour of the day)
          nn_input_vector(523)    = cos(2.0* pi * doy/365.0) ! cos(day of the year)

          nn_output_vector(:)     = 0.0

! Call NN computation
          call compute_nn(nn_input_vector,nn_output_vector) !,nn_num_of_members,& 

! Unpack NN output vector
          gu0(:)       = ugrs(:) + nn_output_vector(254:128:-1)*dt ! u component of layer wind (inverted in vertical)
          gv0(:)       = vgrs(:) + nn_output_vector(381:255:-1)*dt ! v component of layer wind (inverted in vertical)
          gt0(:)       = tgrs(:) + nn_output_vector(127:1:-1)*dt   ! layer mean temperature (inverted in vertical) 
          oshm(:)      = shm(:)  + nn_output_vector(508:382:-1)*dt ! specific humidity (inverted in vertical)
        end subroutine eval_nn

        subroutine relum(x,y)
          integer i
          real(kind=4), intent(in) :: x(4096)
          real(kind=4), intent(out) :: y(4096)
          do i =1,4096
             y(i) = max(x(i),0.0)
          end do
        end subroutine relum
 
        subroutine  compute_nn(X,Y) !,num_of_members,w1,w2,b1,b2,nhid)
 !  Input:
 !            X(IN) NN input vector 
          real(kind=4), intent(in):: X(nn_sizes(1,1))

 !   Ouput:
 !            Y(OUT) NN output vector (composition coefficients for SNN)

          real(kind=4), intent(out):: Y(nn_sizes(1,4)*4)

! Local variables 
          integer i, nout
          real(kind=4) :: x_tmp1(nn_sizes(1,2)), x_tmp2(nn_sizes(1,3))! x2(:),x3(:)
          !real(kind=kind_phys), allocatable :: x_tmp1(:), x_tmp2(:)! x2(:),x3(:)
          integer member,layer, l0, l1
          Y(:) = 0.0

          do member = 1,4 !nn_num_of_members
              !call relum( matmul( nns(1,1)%w,X )+nns(1,1)%b, x_tmp2 )
              call relum( matmul( nns(member)%w0,X )+nns(member)%b0, x_tmp2 )
! Internal layers 
                  x_tmp1(:) = x_tmp2(:)   
              !call relum( matmul( nns(1,2)%w,x_tmp1 )+nns(1,2)%b, x_tmp2 )
              call relum( matmul( nns(member)%w1,x_tmp1 )+nns(member)%b1, x_tmp2 )
! Output layer     
                  x_tmp1(:) = x_tmp2(:)
              !Y((member-1)*127+1:member*127) = matmul( nns(1,3)%w,x_tmp1 )+nns(1,3)%b
              !Y(1:127) =  matmul( nns(1)%w2,x_tmp1 )+nns(1)%b2
              Y((member-1)*127+1:member*127) = matmul( nns(member)%w2,x_tmp1 )+nns(member)%b2
          end do    
    end  subroutine  compute_nn

  end module neuralphys
