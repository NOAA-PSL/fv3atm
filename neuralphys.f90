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

        !use module_iounitdef, only : nonetf
        use sim_nc_mod,        only: open_ncfile,  &
                                     close_ncfile, &
                                     get_var2_real,&
                                     get_var1_real,&
                                     get_var2_double,&
                                     get_var1_double,&
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
!
        real(kind=kind_phys), parameter     :: pi = 3.1415927
! Files containing NN weights and biases

        character(*), parameter::nn_file_name(nn_num_of_members)= (/& 
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_t.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_u.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_v.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_q.nc' /)
             
        !integer,      parameter::nn_ncid(nn_num_of_members)= (/101,102,103,104/)
        
        integer,    allocatable :: nn_sizes(:,:) ! nn_num_of_members x num_layers
        integer num_layers,num_layersm1

 ! Internal types and variables
        type nndata
            real(kind=kind_phys), allocatable :: w(:,:)
            real(kind=kind_phys), allocatable :: b(:)
        end type nndata

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
        type(nndata), allocatable :: nns(:,:)

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
          allocate(nns(nn_num_of_members,num_layersm1))

! Load NNs weights and biases

          do member=1,nn_num_of_members
              call open_ncfile(nn_file_name(member), ncid)
        
              STATUS = NF_GET_ATT_INT (ncid, NF_GLOBAL, 'nn_sizes', nn_sizes(member,:))
              IF (STATUS .NE. NF_NOERR) CALL handle_err("get_att_int",STATUS)
              
              do layer=1,num_layers-1
                  l0 = nn_sizes(member,layer)
                  l1 = nn_sizes(member,layer+1)
                  allocate(nns(member,layer)%w(l0,l1))
                  allocate(nns(member,layer)%b(l1))
                  
                  write(varname, "(A1,I1)") "w",layer-1
                  !call get_var2_real(ncid, varname, l0, l1, nns(member,layer)%w(:,:))
                  call get_var2_double(ncid, varname, l0, l1,nns(member,layer)%w(:,:))
                  write(varname, "(A1,I1)") "b",layer-1
                  !call get_var1_real(ncid, varname, l1, nns(member,layer)%b(:))
                  call get_var1_double(ncid, varname, l1, nns(member,layer)%b(:))
              end do
              call close_ncfile(ncid)
              
              if (me == 0)  print*,'Module NEURALPHYS: NN File Loaded: ', & 
                                   nn_file_name(member)
          end do
        if (me == 0)  print*,'Module NEURALPHYS: NN SIZES: ', nn_sizes
        end subroutine init_nn

! TMP emulator

        subroutine  eval_nn(                                  & 
!  Inputs:
! Prognostic variables:
     &     pgr,ugrs,vgrs,tgrs,shm,  &
! Surface variables:
     &     cmm,evcw,evbs,sbsno,snohf,snowc,srunoff,  &
     &     trans,tsfc,tisfc,q2m,epi,zorl,alboldxy,   &
     &     sfcflw,sfcfsw,topflw,topfsw,slmsk,        &
! Metavariables:
     &     hour,doy,glon,glat,dt,  &
! Outputs:
     &     gu0,gv0,gt0,oshm)

! Inputs
          !integer, intent(in) ::  me
          real(kind=kind_phys), intent(in) :: pgr,cmm,evcw,evbs,sbsno,snohf,snowc,srunoff,trans,tsfc,tisfc,q2m,epi,zorl,alboldxy,slmsk,glat,glon,hour,doy,dt
          type (sfcflw_type), intent(in) :: sfcflw
          type (sfcfsw_type), intent(in) :: sfcfsw
          type (topflw_type), intent(in) :: topflw
          type (topfsw_type), intent(in) :: topfsw

          real(kind=kind_phys), intent(in):: ugrs(127),vgrs(127),tgrs(127),shm(127)
! Outputs     
          real(kind=kind_phys), intent(out):: gu0(127),gv0(127),gt0(127),oshm(127)
!
! Local variables
          real(kind=kind_phys)  nn_input_vector(nn_sizes(1,1)),  nn_output_vector(nn_sizes(1,num_layers)) 

!             
! Create NN input vector:
!        
          nn_input_vector(1:127)  = tgrs
          nn_input_vector(128)    = log(pgr)
          nn_input_vector(129:255)= ugrs
          nn_input_vector(256:382)= vgrs
          nn_input_vector(383:509)= shm
          nn_input_vector(510)    = 0. !cmm ! or uustar; may not be the same
          nn_input_vector(511)    = 0. !evcw
          nn_input_vector(512)    = 0. !evbs
          nn_input_vector(513)    = 0. !sbsno
          nn_input_vector(514)    = 0. !snohf
          nn_input_vector(515)    = 0. !snowc !may not be the same
          nn_input_vector(516)    = 0. !srunoff !may not be the same
          nn_input_vector(517)    = 0. !trans !may not be the same
          nn_input_vector(518)    = tsfc
          nn_input_vector(519)    = tisfc
          nn_input_vector(520)    = q2m
          nn_input_vector(521)    = 0. !epi  !may not be the same
          nn_input_vector(522)    = 0. !zorl !may not be the same
          nn_input_vector(523)    = 0. !alboldxy !may not be the same
          nn_input_vector(524)    = sfcflw%dnfx0
          nn_input_vector(525)    = sfcfsw%dnfx0
          nn_input_vector(526)    = sfcflw%upfx0
          nn_input_vector(527)    = topflw%upfx0
          nn_input_vector(528)    = sfcfsw%upfx0
          nn_input_vector(529)    = topfsw%upfx0
          nn_input_vector(530)    = slmsk
          nn_input_vector(531)    = glat*180./pi
          nn_input_vector(532)    = sin(glon)
          nn_input_vector(533)    = cos(glon)
          nn_input_vector(534)    = sin(2.* pi * hour/24.)
          nn_input_vector(535)    = sin(2.* pi * doy/365.)
          nn_input_vector(536)    = cos(2.* pi * hour/24.)
          nn_input_vector(537)    = cos(2.* pi * doy/365.)


!             
! Call NN computation
          call compute_nn(nn_input_vector,nn_output_vector) !,nn_num_of_members,& 
!
! Unpack NN output vector
          !nn_output_vector(:) = 0.
          gu0(:)       = ugrs + nn_output_vector(128:254)*dt !+ ugrs   ! u component of layer wind
          gv0(:)       = vgrs + nn_output_vector(255:381)*dt !+ vgrs ! v component of layer wind
          gt0(:)       = tgrs + nn_output_vector(1:127)*dt   !+ tgrs! layer mean temperature 
          oshm(:)      = shm  + nn_output_vector(382:508)*dt !+ shm                      ! specific humidity
!
          !if (me == 0)  print*,'Module NEURALPHYS: Before end eval_nn'
        end subroutine eval_nn
        
        subroutine  compute_nn(X,Y) !,num_of_members,w1,w2,b1,b2,nhid)
 !  Input:
 !            X(IN) NN input vector 
          real(kind=kind_phys), intent(in)::X(nn_sizes(1,1))
          !integer, intent(in):: me

 !   Ouput:
 !            Y(OUT) NN output vector (composition coefficients for SNN)

          real(kind=kind_phys), intent(out):: Y(508) !(sum(nn_sizes(:,num_layers)))

! Local variables 
          integer i, nout
          real(kind=kind_phys), allocatable :: x_tmp1(:), x_tmp2(:)! x2(:),x3(:)
          integer member,layer, l0, l1
          !if (me==0) print*,"NEURAL: size:", size(X), size(Y)
          do member = 1, nn_num_of_members
          !member=1
              if (allocated(x_tmp1)) deallocate(x_tmp1)
              if (allocated(x_tmp2)) deallocate(x_tmp2) 
              do layer = 1,num_layers-2 ! loop from 1st to final-1 layer
                  l0 = nn_sizes(member,layer)
                  l1 = nn_sizes(member,layer+1)
                  !if (me == 0) print*,"NEURAL: member,layer:", member, layer,nn_num_of_members,num_layers
                  !if (me == 0) print*,"NEURAL: allocate 2", allocated(x_tmp2), l1
                  if (allocated(x_tmp2)) deallocate(x_tmp2)
                  allocate(x_tmp2(l1))
                  
                  if (layer == 1) then
                      !if (me == 0) print*,"NEURAL: allocate 1", allocated(x_tmp1), l0
                      allocate(x_tmp1(l0))
                      x_tmp1(:)=X(:)
                  endif
! Internal layers                  
!$OMP PARALLEL default (shared) private (i)
!$OMP DO
                  do i = 1,l1
                      x_tmp2(i) = tanh(sum(x_tmp1*nns(member,layer)%w(:,i))+nns(member,layer)%b(i))
                  end do
!$OMP END DO
!$OMP END PARALLEL 
                  !if (me == 0) print*,"NEURAL: deallocate 1", allocated(x_tmp1), size(x_tmp1)
                  if (allocated(x_tmp1)) deallocate(x_tmp1)
                  !if (me == 0) print*,"NEURAL: allocate 1", allocated(x_tmp1), l1
                  allocate(x_tmp1(l1))
                  x_tmp1(:) = x_tmp2(:)
                  !if (me == 0) print*,"NEURAL: deallocate 2", allocated(x_tmp2),size(x_tmp2)
                  if (allocated(x_tmp2)) deallocate(x_tmp2)
              end do
              
! Output layer              
!$OMP PARALLEL default (shared) private (i) 
!$OMP DO              
              do i = 1,127 !l1
                  Y((member-1)*127+i) = sum(x_tmp1*nns(member,layer)%w(:,i))+nns(member,layer)%b(i)
              end do
!$OMP END DO 
!$OMP END PARALLEL 
              !if (me == 1) print*,"NEURAL: member done:",member,nn_num_of_members,((member-1)*127+i)

              !print*,'Module NEURALPHYS: End compute' 
              if (allocated(x_tmp1)) deallocate(x_tmp1)
          end do    
         
    end  subroutine  compute_nn

  end module neuralphys
