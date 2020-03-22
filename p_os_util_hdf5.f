! this module provides parallel hdf5 writing functionality 
! for 2d and 3d datasets with decomposition in any dimension
! also included an append mode which is used in pipeline code
! updated : Apr 18, 2010  
! Chengkun Huang

      module m_pdiagnostic_utilities

        use m_h5_diagnostic_utilities, only : detect_precision,         &
     & add_h5_axes, add_h5_attrs, open_h5_dataset, add_h5_atribute, start_hdf5, stop_hdf5
        use HDF5

        implicit none

        include 'mpif.h'
        private

        integer :: io_comm, io_group
        integer :: nio, ntpg, idproc, nvp, ioLeader


        public :: pwrite_hdf_sds, pwrite_hdf_sds_direct
        public :: pwrite_hdf_raw
        public :: create_io_comm, destroy_io_comm
        
          
        interface pwrite_hdf_sds
          module procedure pwrite_hdf_3d
          module procedure pwrite_hdf_2d
        end interface

        interface pwrite_hdf_sds_direct
          module procedure pwrite_hdf_2d_direct
        end interface

        interface pwrite_hdf_raw
          module procedure pwrite_hdf_raw
          module procedure pwrite_hdf_raw_tag
          module procedure pwrite_hdf_raws
          module procedure pwrite_hdf_plasma_raws
        end interface

       contains
       
!---------------------------------------------------
        subroutine create_io_comm()
            integer :: world_group
            integer :: ierr
            integer :: ii
            integer, dimension(:), allocatable :: group_ranks

            ! common block for parallel processing
               integer :: nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
            ! nproc = number of real or virtual processors obtained
            ! lgrp = current communicator
            ! lworld = MPI_COMM_WORLD communicator
               common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

! determine the rank of the calling process in the communicator 
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)

      ! determine the number of tasks doing I/O
      ! I found it best to use only 1 when doing field files
      ! may need to be re-evaluated
      nio = 1
!      if ( nvp <= 16 ) then
!         nio = nvp
!      else
!         nio = max((nvp+64-1)/64,16)
!         nio = min(nio,256)
!      endif

        ! Now group nodes for I/O.  
        ntpg = (nvp+nio-1)/nio      ! max number of tasks per group
        ioLeader = (idproc/ntpg)*ntpg

         ! calculate the members of the IO group
         ! get the world group
         call MPI_COMM_GROUP(lworld, world_group, ierr)
   
         allocate(group_ranks(nio))

         do ii=1,nio
            group_ranks(ii) = (ii-1)*ntpg
         enddo
      
         ! create the I/O group
         call MPI_GROUP_INCL(world_group, nio, group_ranks, io_group, ierr)

         ! create the I/O communicator
         call MPI_COMM_CREATE(lworld, io_group, io_comm, ierr)

        end subroutine create_io_comm

!---------------------------------------------------
      subroutine destroy_io_comm()
         integer :: ierr

         if (io_comm .ne. MPI_COMM_NULL) then
            call MPI_COMM_FREE(io_comm, ierr)
         endif

         call MPI_GROUP_FREE(io_group, ierr)

      end subroutine destroy_io_comm

!---------------------------------------------------
        subroutine set_phdf5_properties(rank, flplID, dcplID, xferID)

        integer, intent(in) :: rank
        integer(hid_t), intent(inout) :: flplID, dcplID, xferID

        integer, dimension(rank) :: dimsChunk
        integer :: ierr

! crete property list for file creation                    
        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_CREATE_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

! set chunk size
!        dimsChunk(1) = chunk_size(1)
!        dimsChunk(2) = chunk_size(2)
!        dimsChunk(3) = chunk_size(3)
      
!        call h5pset_chunk_f(dcplID, rank, dimsChunk, ierr)

! setting tranfer/data conversion buffer size
!        call h5pset_buffer_f( xferID, h5_tune_buffer_size, ierr )

        end subroutine set_phdf5_properties 
!---------------------------------------------------       


!---------------------------------------------------
        subroutine pwrite_hdf_3d(path, file,                            &
     &                    dataset, name, label, units,                  &
     &                     n, t, dt, xmin, xname,                       &
     &                    xlabel, xunits, nvpy, nvpz, in_stride)
!---------------------------------------------------   
! this subroutine collects distributed 3d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:,:) :: dataset
        character(len=*), intent(in) :: path, file
        character(len=*), intent(in) :: name, label, units
        integer, intent(in) :: n, nvpy, nvpz
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        character(len=*), dimension(:), intent(in) :: xname, xlabel, xunits

!       local constants
        integer, parameter :: rank = 3                 ! 3D Data
        integer, parameter :: X = 1, Y = 2, Z = 3

        integer, dimension(rank), intent(in) :: in_stride

!       local variables
        ! real, dimension(:,:,:), allocatable :: ldataset
        integer :: j, ii, nx, ny, nz  
        character(len_trim(path)+1+len_trim(file)) :: file_name
!        character(len = 80) :: tmpstr1, tmpstr2 
        integer sd_id, sds_id, sds_index
        integer(hsize_t), dimension(rank) :: dim_sizes, start, nodestart, edges, stride
        integer(hsize_t), dimension(rank) :: dim_sizes_grid, dimsChunk
        integer :: sdfileID, sdsID, dimID
        integer :: status
        real, dimension(rank) :: xmax
        integer :: info, ierr         
        integer, dimension(1) :: readyflag = 1

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld


! hdf5 variables
        integer(HID_T), dimension(rank) :: bc, mc 
        integer(hid_t) :: d_float=0
        integer(hid_t) :: flplID=0, xferID=0, dcplID=0, memspaceID=0 
        integer(hid_t) :: file_id=0, rootID=0, dset_id=0, dspace_id=0

! extras BEN
	integer, dimension(MPI_STATUS_SIZE) :: stat
	integer :: readystatus
	integer :: nnodes
    real, dimension(:,:,:), allocatable :: tmpdataset1
    real, dimension(:,:,:), allocatable :: tmpdataset2

! local data
      integer istatus, lrec, np, ioff, id, nrec0, i

	readystatus=0
	istatus=0
      call start_hdf5(ierr)

! hard wire boundary condition and moving window properties for the moment
      bc = 1
      mc = 0

	call MPI_COMM_SIZE( lworld, nnodes, ierr)

   ! setting leader node as 0 -- I found previous blocking sends worked best
   ! with only one leader... I have not tested speed again for multiple groups 
   !if (idproc==ioLeader) then
   if (idproc==0) then
      file_name = trim(path) // trim(file)

      ! I was going to write the data on the leader nodes to a scratch file then read it back in,
      ! but then I realized that after dumping it, the data in sfield is not used, so there's no point.

      call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

      info = MPI_INFO_NULL
      !call MPI_INFO_CREATE( info, ierr )

      call h5pset_fapl_mpio_f(flplID, io_comm, info, ierr)

      ! create new file
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
      & access_prp=flplID) 
      

      ! close property list
      call h5pclose_f(flplID, ierr)


        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_CREATE_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif


! get data size
!      idiv = division

      nx = size(dataset,1)
      ny = size(dataset,2)
      nz = size(dataset,3)
      
! set up one tmpdataset for initial output to the current h5 file
	allocate( tmpdataset1(nx,ny,nz) )
	tmpdataset1 = dataset
	
! extend data size
!      dim_sizes(idiv) = dim_sizes(idiv)*nproc
      dim_sizes_grid(1) = (nx + in_stride(1) - 1)/in_stride(1)
      dim_sizes_grid(2) = (ny*nvpy + in_stride(2) - 1)/in_stride(2)
      dim_sizes_grid(3) = (nz*nvpz + in_stride(3) - 1)/in_stride(3)

! chunk stuff
      dimsChunk(1) = dim_sizes_grid(1)
      dimsChunk(2) = (ny + in_stride(2) - 1)/in_stride(2)
      dimsChunk(3) = (nz + in_stride(3) - 1)/in_stride(3)

!        Calcualte Axis Data
         do j=1, rank
            xmax(j) = xmin(j) + (dim_sizes_grid(j)*in_stride(j) - (in_stride(j) - 1))
         enddo

!     Determine the precision of real variable
      d_float = detect_precision()

! try chunking for faster operations
      call h5pset_chunk_f(dcplID, rank, dimsChunk, ierr)
! create file dataspace
      call h5screate_simple_f(rank, dim_sizes_grid, dspace_id, ierr)

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc )

! add data axes 
         call add_h5_axes(file_id, rank, xunits, xname, xlabel, xmin,   &
     & xmax, xferID)

      stride = 1

! create dataset in dataspace
      call h5gopen_f( file_id, '/', rootID, ierr )
      call h5dcreate_f (rootID, name, d_float, dspace_id, dset_id, ierr, dcplID)

      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    


      ! loop over the nodes in the group
      ! may be worth testing again for multiple groups
      !do ii=idproc,min(ntpg+idproc-1,nvp-1)
      do ii=0,nnodes-1

   ! before other calculations and writing of current data are done
   ! 	issue a ready-flag and do non-blocking receive of data for next ii-loop
   ! alternates between tmpdataset1 and tmpdataset2, depending on which needs to be present
   ! 	for later writing and should therefore not be the object of a current non-blocking receive 
         if (ii < nnodes - 1) then ! don't have to recieve for last time through loop
              ! send is necessary to stop buffer overflow on uBGL
            if ( mod( ii, 2 ) == 0 ) then
               allocate( tmpdataset2(nx,ny,nz) )
               call MPI_ISEND(readyflag, size(readyflag), mint, ii+1, 99, lworld, readystatus,ierr)
               call MPI_IRECV(tmpdataset2,size(tmpdataset2),mreal,ii+1,99,lworld,istatus,ierr)
            else
               allocate( tmpdataset1(nx,ny,nz) )
               call MPI_ISEND(readyflag, size(readyflag), mint, ii+1, 99, lworld, readystatus,ierr)
               call MPI_IRECV(tmpdataset1,size(tmpdataset1),mreal,ii+1,99,lworld,istatus,ierr)
	    endif
         endif


   ! calculate start of the data this node writes..
   !..gets a bit complicated b/c I allow for arbitrary strides.
   ! hyperslab coordinates are 0 indexed
!print *,'iter starting: ',ii
         start = 0
         start(2) = ( ((ii - ii / nvpy * nvpy) * ny + in_stride(2) - 1) / in_stride(2) )
         start(3) = ( (ii / nvpy * nz + in_stride(3) - 1) / in_stride(3) )

   ! calculate where to start writing from this dataset
         nodestart(2) = start(2) * in_stride(2) + 1 - (ii - ii / nvpy * nvpy) * ny
         nodestart(3) = start(3) * in_stride(3) + 1 - ii / nvpy * nz

   ! the size of the data *this* node will write
         dim_sizes(1) = dim_sizes_grid(1)
         dim_sizes(2) = (ny - (nodestart(2) - 1) + in_stride(2) - 1)/in_stride(2)
         dim_sizes(3) = (nz - (nodestart(3) - 1) + in_stride(3) - 1)/in_stride(3)


   ! create memory dataspace
         call h5screate_simple_f(rank, dim_sizes, memspaceID, ierr )

   ! select hyperslab in the file

   !      call h5dget_space_f(dset_id, dspace_id, ierr)
         call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F, start,   &
        & dim_sizes, ierr, stride)


   ! write the dataset collectively  
   ! depending on which tmpdataset has previously been collected for writing
         if ( mod( ii, 2 ) == 0 ) then
            call h5dwrite_f(dset_id, d_float, &
        &     tmpdataset1(1:nx:in_stride(1), nodestart(2):ny:in_stride(2), nodestart(3):nz:in_stride(3) ), &
        &     dim_sizes,ierr, memspaceID, dspace_id, xfer_prp=xferID)
         else
             call h5dwrite_f(dset_id, d_float, &
        &      tmpdataset2(1:nx:in_stride(1), nodestart(2):ny:in_stride(2), nodestart(3):nz:in_stride(3) ), &
        &      dim_sizes,ierr, memspaceID, dspace_id, xfer_prp=xferID)
         endif

   !      call h5dwrite_f(dset_id, d_float, dataset(1:nx:in_stride(1), nodestart(2):ny:in_stride(2), nodestart(3):nz:in_stride(3) ), dim_sizes,ierr, memspaceID, dspace_id, xfer_prp=xferID)
         
         call h5sclose_f(memspaceID, ierr)

   ! waiting for non-blocking calls to complete
         if (ii < nnodes - 1) then
              call MPI_WAIT(readystatus,stat,ierr)
              call MPI_WAIT(istatus,stat,ierr)
         endif

         if ( mod( ii, 2 ) == 0 ) then
            deallocate( tmpdataset1 )
         else
            deallocate( tmpdataset2 )
         endif

      enddo

   ! write dataset arrtibutes
      call add_h5_atribute(dset_id,'UNITS', units) 
      call add_h5_atribute(dset_id,'LONG_NAME', label) 

   ! close resources
      call h5sclose_f(dspace_id, ierr)
         
   ! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5dclose_f(dset_id, ierr)

      call h5pclose_f(xferID, ierr)
      call h5pclose_f(dcplID, ierr)

      !flush the file....not sure if this helps, but it doesn't seem to hurt
      call h5fflush_f(file_id, H5F_SCOPE_LOCAL_F, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error flushing file", idproc
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

      call h5fclose_f(file_id, ierr)

   else ! send the data to the leader node

      ! leader has been set to idproc = 0

      call MPI_IRECV(readyflag, size(readyflag), mint, 0, 99, lworld, readystatus, ierr)
      call MPI_WAIT(readystatus,stat,ierr)

      call MPI_ISEND(dataset,size(dataset),mreal,0,99,lworld,istatus,ierr)
      call MPI_WAIT(istatus,stat,ierr)

   endif

      call stop_hdf5(ierr)
     
      end subroutine pwrite_hdf_3d
!---------------------------------------------------       

!---------------------------------------------------
        subroutine pwrite_hdf_2d(path, file,                            &
     &                    dataset ,division,                            &
     &                   name, label, units,                            &
     &                       n, t, dt, xmin, dx,                        &
     &                xname, xlabel, xunits)
!---------------------------------------------------   
! this subroutine collects distributed 2d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:) :: dataset
        character(len=*), intent(in) :: path, file
        character(len=*), intent(in) :: name, label, units
        integer, intent(in) :: n
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(:), intent(in) :: dx
        character(len=*), dimension(:), intent(in) :: xname, xlabel, xun&
     &its
        integer, intent(in) :: division
        integer :: idiv

!       local constants
        integer, parameter :: rank = 2                 ! 2D Data
        integer, parameter :: X = 1, Y = 2

!       local variables
        real, dimension(:,:), allocatable :: ldataset
        integer :: j   
        character(len_trim(path)+1+len_trim(file)) :: file_name
!        character(len = 80) :: tmpstr1, tmpstr2 
        integer sd_id, sds_id, sds_index
        integer(hsize_t), dimension(rank) :: dim_sizes, start, edges,   &
     & stride
        integer(hsize_t), dimension(rank) :: dim_sizes_append
        integer :: sdfileID, sdsID, dimID
        integer :: status
        real, dimension(rank) :: xmax
        integer(hid_t) :: d_float
        logical :: filexisted

        integer(HID_T), dimension(rank) :: bc, mc 

! hdf5 variables
        integer(hid_t) :: flplID, xferID, dcplID, memspaceID, driver 
        integer(hid_t) :: file_id, rootID, dset_id, dspace_id
        integer :: info, ierr 

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i
      dimension istatus(lstat)
! determine the rank of the calling process in the communicator 
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
!      print *,"ierr,idproc==",ierr,idproc  

! hard wire boundary condition and moving window properties for the moment
      bc = 0
      mc = 0

      file_name = trim(path) // trim(file)

! init phdf5 interface
      !call h5open_f(ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error starting HDF5"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
!     Determine the precision of real variable
      d_float = detect_precision()
      call set_phdf5_properties(rank, flplID, dcplID, xferID)
      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, lgrp, info, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pset_fapl_mpio_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
! double checking file driver is MPI-IO      
!      CALL h5pget_driver_f(flplID, driver, ierr)
!      write(message,*), "driver, H5FD_MPIO_F = ", driver, H5FD_MPIO_F
!      call Write_Log(message, 0)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pset_dxpl_mpio_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

! get data size
      idiv = division
      if ((division.ne.X) .and. (division.ne.Y)) idiv=Y 
      dim_sizes(1) = size(dataset,1)   
      dim_sizes(2) = size(dataset,2)    
      edges(1:2) = dim_sizes(1:2)

! extend data size
!      dim_sizes(idiv) = dim_sizes(idiv)*nproc
      dim_sizes_append = dim_sizes
      dim_sizes_append(idiv) = dim_sizes(idiv)*nproc
      dim_sizes_append(rank) = dim_sizes_append(rank)
      
! preprocess data
      select case (idiv)
        case (X)
          allocate(ldataset(dim_sizes(1),dim_sizes(2)))
          ldataset = dataset
!          allocate(ldataset(dim_sizes(2),dim_sizes(1)))
!          ldataset = reshape(dataset, (/dim_sizes(2),dim_sizes(1)/))
        case (Y)
          allocate(ldataset(dim_sizes(1),dim_sizes(2)))
          ldataset = dataset
        case default
         print*, "p_os_util_hdf5: no match for idiv, idiv=", idiv
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      end select

! open file and dataset
! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error creating file ", file_name
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc )

!        Write Axis Data
         do j=1, rank
            xmax(j) = xmin(j) + (dim_sizes_append(j)-1)*dx(j)
         enddo

! add data axes 
        call add_h5_axes(file_id, rank, xunits, xname, xlabel, xmin,    &
     & xmax, xferID)

! synchronization
      call MPI_BARRIER(lgrp,ierr)
! create file dataspace
      call h5screate_simple_f(rank, dim_sizes_append, dspace_id, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5screate_simple_f creating dataspace"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
! create memory dataspace
      call h5screate_simple_f(rank, dim_sizes, memspaceID, ierr )
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5screate_simple_f creating memory dataspace"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
! create dataset in dataspace
      call h5gopen_f( file_id, '/', rootID, ierr )
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5gopen_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

      call h5dcreate_f(rootID, name, d_float        , dspace_id,     &
     & dset_id, ierr, dcplID)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5dcreate_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = 0
      start(idiv) = idproc*dim_sizes(idiv) 
      !if (ldoappend) then
      !   start(rank) = dim_sizes(rank)*append_pos
      !endif
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F, start,   &
     & dim_sizes, ierr, stride)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, ldataset, dim_sizes,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
!     &     memspaceID, dspace_id)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5dwrite_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif

! write dataset arrtibutes
      call add_h5_atribute(dset_id,'UNITS', units) 
      call add_h5_atribute(dset_id,'LONG_NAME', label) 


! close resources
      call h5sclose_f(memspaceID, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5sclose_f, memspace"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      call h5sclose_f(dspace_id, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5sclose_f, dspace"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      call h5pclose_f(xferID, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pclose_f, xferID"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      call h5pclose_f(dcplID, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pclose_f, dcplID"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      call h5pclose_f(flplID, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5pclose_f, flplID"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
! close dataset and file
      call h5gclose_f(rootID, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5gclose_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      call h5dclose_f(dset_id, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5dclose_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      call h5fclose_f(file_id, ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error in h5fclose_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
      !call h5close_f(ierr)
      !if (ierr.ne.0) then 
      !   print*, "p_os_util_hdf5: Error shutting down HDF5"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      !endif

      deallocate(ldataset)
     
      end subroutine pwrite_hdf_2d
!---------------------------------------------------       

!---------------------------------------------------
        subroutine pwrite_hdf_2d_direct(path, file,                     &
     &                       dataset , idsend,                          &
     &                       name, label, units,                        &
     &                       n, t, dt, xmin, dx,                        &
     &                       xname, xlabel, xunits,                     &
     &                       doappend, nappend, append_pos, docompress)
!---------------------------------------------------   
! this subroutine collect local 2d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:) :: dataset
        character(len=*), intent(in) :: path, file
        character(len=*), intent(in) :: name, label, units
        integer, intent(in) :: n, nappend, append_pos
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(:), intent(in) :: dx
        character(len=*), dimension(:), intent(in) :: xname, xlabel, xun&
     &its
        integer :: idsend
        logical :: doappend, docompress

!       local constants
        integer, parameter :: rank = 2                 ! 2D Data
        integer, parameter :: X = 1, Y = 2      

!       local variables
!        real, dimension(:,:), allocatable :: ldataset
        integer :: j   
        character(len_trim(path)+1+len_trim(file)) :: file_name
!        character(len = 80) :: tmpstr1, tmpstr2 
        integer sd_id, sds_id, sds_index
        integer(hsize_t), dimension(rank) :: dim_sizes, start, edges,   &
     & stride
        integer(hsize_t), dimension(rank) :: dim_sizes_append
        integer :: sdfileID, sdsID, dimID
        integer :: status
        real, dimension(rank) :: xmax
        integer(hid_t) :: d_float
        logical :: ldoappend, filexisted

        integer(HID_T), dimension(rank) :: bc, mc 
        
! common block for parallel processing
        integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
        parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
        common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble,lworld
! local data
        integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, ierr
        dimension istatus(lstat)
        logical :: fileexisted

! hdf5 variables
        integer(HID_T) :: file_id, dset_id, dspace_id, memspaceID
        integer(hid_t) :: rootID
        integer(HSIZE_T), dimension(2) :: dims
      
! determine the rank of the calling process in the communicator 
        call MPI_COMM_RANK(lgrp,idproc,ierr)
!      print *,"ierr,idproc,lworld==",ierr,idproc,lworld
! determine the size of the group associated with a communicator
        call MPI_COMM_SIZE(lgrp,nvp,ierr)
!      print *,"ierr,idproc==",ierr,idproc  

! hard wire boundary condition and moving window properties for the moment
        bc = 0
        mc = 0

        file_name = trim(path) // trim(file)
        ldoappend = doappend
!        if (ldoappend) then
! check if file existed         
!           inquire(file=file_name,exist=filexisted)
!           if (.not. filexisted) ldoappend=.false.
!        endif
!
! synchronization
!      call MPI_BARRIER(lgrp,ierr)
!     
      dim_sizes(1) = size(dataset,1)   
      dim_sizes(2) = size(dataset,2)    
      edges(1:rank) = dim_sizes(1:rank)

! extend data size in append mode 
      dim_sizes_append = dim_sizes
      dim_sizes_append(rank) = dim_sizes_append(rank) * nappend

! only the node that contains the data will write 
      if (idproc == idsend) then
!         allocate(ldataset(dim_sizes(1),dim_sizes(2)))
                
! initialize fortran predifined datatype            
         call h5open_f(ierr)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error starting HDF5"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
! detect the precision of the system
         d_float = detect_precision()

! open file and dataset
! create new file
            call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr) 
            if (ierr.ne.0) then 
               print*, "p_os_util_hdf5: Error creating file ", file_name
               !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
            endif
!         set data label and units
            call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p',xmin&
     &, xmax, bc, mc )

!        Write Axis Data
            do j=1, rank
               xmax(j) = xmin(j) + (dim_sizes_append(j)-1)*dx(j)
            enddo

! add data axes 
            call add_h5_axes(file_id, rank, xunits, xname, xlabel, xmin,&
     & xmax)

! create file dataspace
         call h5screate_simple_f(rank, dim_sizes_append, dspace_id,ierr)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5screate_simple_f, dataspace"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
! create memory dataspace
         call h5screate_simple_f(rank, dim_sizes, memspaceID, ierr)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5screate_simple_f, memory dataspace"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
! create dataset in dataspace
         call h5gopen_f( file_id, '/', rootID, ierr )
         call h5dcreate_f (rootID, name, d_float, dspace_id, dset_id,&
     & ierr)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5dcreate_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
         call h5gclose_f(rootID, ierr)

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
         start = 0
         stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
         call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F, start,&
     & dim_sizes, ierr, stride)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f, H5S_SELECT_SET_F"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

! write the dataset collectively  
         call h5dwrite_f(dset_id, d_float, dataset, dim_sizes, ierr,    &
     &     memspaceID, dspace_id)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error writing dataset"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

!       Compress SDS
         if (docompress) then
!       place holder        
         endif

! write dataset arrtibutes
         call add_h5_atribute(dset_id,'UNITS', units) 
         call add_h5_atribute(dset_id,'LONG_NAME', label) 
        
! close resources
         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
! close dataset and file
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         !call h5close_f(ierr)

!         deallocate(ldataset)

      endif
                
     
      end subroutine pwrite_hdf_2d_direct
!-----------------------------------------------------------------------              
      subroutine pwrite_hdf_raw(path, file, tnps,tnpp, soffset,x0, dx,&
     &               name, n, t, dt, xmin, dspl,                      &
     & dataset,doappend, docompress)
!---------------------------------------------------   
! this subroutine collects distributed 1d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:,:), pointer, intent(in) :: dataset
        real, intent(in) :: soffset
        real, dimension(:) :: x0, dx
        character(len=*), intent(in) :: path, file
        integer :: tnps, tnpp, dspl
        integer :: pgrp, nvp, idproc, tp, tpo, color, i, j
        
        integer(HSIZE_T), dimension(1) :: maxdim
        integer, dimension(:,:), allocatable:: dims
        integer(hsize_t), dimension(1) :: ldim
        integer, dimension( : ), allocatable :: np
        real, dimension(:), allocatable :: write_buffer
        
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(3) :: xmax
        logical, intent(in) :: doappend,docompress

        character(len_trim(path)+1+len_trim(file)) :: file_name

        integer(hsize_t), dimension(1) :: start, stride
        integer(hid_t) :: d_float
        logical :: ldoappend

        integer(HID_T), dimension(1) :: bc, mc 

! hdf5 variables
        integer(hid_t) :: flplID, xferID, dcplID, memspaceID, attrID 
        integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
        integer :: info, ierr 

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, ioff, id, nrec0
      dimension istatus(lstat)
      
      tpo = 0
      ldim(1) = 1

! hard wire boundary condition and moving window properties for the moment
      bc = 0
      mc = 0

      file_name = trim(path)//trim(file)
      ldoappend = doappend

      !write(message,*) "tnpp", tnps,tnpp,dspl
      !call Write_Log(message,2)
      
      tnpp = int(tnpp/dspl)
      tnps = 1
      !write(message,*) "tnpp", tnps,tnpp,dspl
      !call Write_Log(message,2)

! If total number of particles to output is 0 then node 0 will create 
! an empty file
      call MPI_ALLREDUCE( tnpp, tp, 1, MPI_INTEGER, MPI_SUM, lgrp, ierr )

      if ( tp == 0 ) then
      
! init phdf5 interface
      call h5open_f(ierr)
      if (ierr.ne.0) then 
         print*, "p_os_util_hdf5: Error starting HDF5"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
      endif
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  

         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, lgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)

! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
! create dataset in dataspace

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
          call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)
! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      
      return
      
      else

! create a communicator that includes only nodes that have particles to output
         if ( tnpp > 0 ) then 
            color = 1
         else
            color = MPI_UNDEFINED
         endif

! this is a global call
         call MPI_COMM_SPLIT( lgrp, color, 0, pgrp, ierr )

         if (tnpp > 0) then

         call MPI_COMM_RANK( pgrp, idproc, ierr )

! get number of nodes with particles to output
		 call MPI_COMM_SIZE( pgrp, nvp, ierr )

! gather number of particles on all 
		 allocate( np( nvp ), stat = ierr )
		 allocate(dims( 2, nvp ), stat = ierr)
		 
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error allocating memory for RAW diagnostic"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
		 
		 call MPI_ALLGATHER( tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER, &
		 &pgrp, ierr)
	 
		 ! get total number of particles and positions on file
		 dims( 1, 1 ) = 1
		 dims( 2, 1 ) = np(1) 
		 do i = 2, nvp
		    dims( 1, i ) = dims( 2, i-1 ) + 1
		    dims( 2, i ) = dims( 1, i ) + np(i) - 1
		 enddo
	 
		 ! get temporary memory for write buffer
		 allocate( write_buffer( tnpp ) , stat = ierr )
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error allocating memory for RAW diagnostic"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif



! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         

         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  

         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! open file
! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
! create dataset in dataspace

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
          call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! synchronization
!      call MPI_BARRIER(lgrp,ierr)

! write data

		 do i = 1, 3
		 
	  if (i == 3) then
!	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) + soffset - x0(i)) * dx(i)
	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp*dspl-1):dspl&
	  &,1) + soffset - x0(i)) * dx(i)
      else	     
!	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) - x0(i)) * dx(i)
	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp*dspl-1):dspl&
	  &,1) - x0(i)) * dx(i)
      endif
      
      
      maxdim = (/H5S_UNLIMITED_F/)
      
      ldim(1) = 1

      call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      ldim(1) = tp

      call h5pset_chunk_f(dcplID, 1, ldim, ierr)

      call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr, dcplID)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5dcreate_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

      ldim(1) = tp

      call h5dextend_f(dset_id, ldim, ierr)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5dextend_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = tpo + dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
         if (ierr.ne.0) then 
            print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'x_'//char(iachar('0')+i)) 

      call h5pclose_f(dcplID, ierr)
      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo


	  do i = 1, 3
		 
!			write_buffer(1:tnpp) = dataset((i+3),tnps:(tnps+tnpp-1),1) 
			write_buffer(1:tnpp) = dataset((i+3),tnps:(tnps+tnpp*dspl-1):dspl,1) 
      
      maxdim = (/H5S_UNLIMITED_F/)

      ldim(1) = 1

      call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      ldim(1) = tp

      call h5pset_chunk_f(dcplID, 1, ldim, ierr)

      call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr, dcplID)
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5dcreate_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5dextend_f(dset_id, ldim, ierr)

		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5extend_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = tpo + dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F, start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim, ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error writing dataset"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'p_'//char(iachar('0')+i)) 

      call h5pclose_f(dcplID, ierr)
	
      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)

	  enddo

!       Compress SDS

      if (docompress) then 
        ! place holder        
      endif


! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      

      deallocate(np)
      deallocate(dims)
      deallocate(write_buffer)

      endif
      
      endif

      if ( pgrp /= MPI_COMM_NULL ) then
         call MPI_COMM_FREE( pgrp, ierr )
      endif

     
      end subroutine pwrite_hdf_raw
!-----------------------------------------------------------------------              
      subroutine pwrite_hdf_raws(path,file, tnps,tnpp, soffset,x0, dx,&
     &               name, n, t, dt, xmin,                        &
     & dataset, docompress)
!---------------------------------------------------   
! this subroutine collects distributed 1d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:,:), pointer, intent(in) :: dataset
        real, intent(in) :: soffset
        real, dimension(:) :: x0, dx
        character(len=*), intent(in) :: path, file
        integer :: tnps, tnpp
        integer :: pgrp, nvp, idproc, tp, color, i, j
        
        integer(HSIZE_T), dimension(1) :: maxdim
        integer, dimension(:,:), allocatable:: dims
        integer(hsize_t), dimension(1) :: ldim
        integer, dimension( : ), allocatable :: np
        real, dimension(:), allocatable :: write_buffer
        
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(3) :: xmax
        logical, intent(in) :: docompress

        character(len_trim(path)+1+len_trim(file)) :: file_name

        integer(hsize_t), dimension(1) :: start, stride
        integer(hid_t) :: d_float
        logical :: ldoappend

        integer(HID_T), dimension(1) :: bc, mc 

! hdf5 variables
        integer(hid_t) :: flplID, xferID, dcplID, memspaceID, attrID 
        integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
        integer :: info, ierr 

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, ioff, id, nrec0
      dimension istatus(lstat)
      
      ldim(1) = 1

! hard wire boundary condition and moving window properties for the moment
      bc = 0
      mc = 0

      file_name = trim(path)//trim(file)

! If total number of particles to output is 0 then node 0 will create 
! an empty file
      call MPI_ALLREDUCE( tnpp, tp, 1, MPI_INTEGER, MPI_SUM, lgrp, ierr )

      if ( tp == 0 ) then

! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, lgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! create new file
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
!         set data label and units
      call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

      call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
      call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )

          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)


      else

! create a communicator that includes only nodes that have particles to output
         if ( tnpp > 0 ) then 
            color = 1
         else
            color = MPI_UNDEFINED
         endif

! this is a global call
         call MPI_COMM_SPLIT( lgrp, color, 0, pgrp, ierr )

         if (tnpp > 0) then

         call MPI_COMM_RANK( pgrp, idproc, ierr )

! get number of nodes with particles to output
		 call MPI_COMM_SIZE( pgrp, nvp, ierr )

! gather number of particles on all 
		 allocate( np( nvp ), stat = ierr )
		 allocate(dims( 2, nvp ), stat = ierr)
		 
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error allocating memory for RAW write"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
		 
		 call MPI_ALLGATHER( tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER, &
		 &pgrp, ierr)
	 
		 ! get total number of particles and positions on file
		 dims( 1, 1 ) = 1
		 dims( 2, 1 ) = np(1) 
		 do i = 2, nvp
		    dims( 1, i ) = dims( 2, i-1 ) + 1
		    dims( 2, i ) = dims( 1, i ) + np(i) - 1
		 enddo
	 
		 ! get temporary memory for write buffer
		 allocate( write_buffer( tnpp ) , stat = ierr )
		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error allocating memory for RAW write"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif



! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         

		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  

		 if (ierr /= 0) then 
         print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! create new file
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
!         set data label and units
      call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

      call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
      call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )

          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! synchronization
!      call MPI_BARRIER(lgrp,ierr)

! write data

		 do i = 1, 3
		 
	  if (i == 3) then
	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) + soffset - x0(i)) * dx(i)
      else	     
	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) - x0(i)) * dx(i)
      endif
            
      ldim(1) = tp

      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'x_'//char(iachar('0')+i)) 

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo


	  do i = 1, 3
		 
			write_buffer(1:tnpp) = dataset((i+3),tnps:(tnps+tnpp-1),1) 
            
      ldim(1) = tp

      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dwrite_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'p_'//char(iachar('0')+i)) 

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo

!       Compress SDS

      if (docompress) then 
        ! place holder        
      endif


! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      

      deallocate(np)
      deallocate(dims)
      deallocate(write_buffer)

      endif

      if ( pgrp /= MPI_COMM_NULL ) then
         call MPI_COMM_FREE( pgrp, ierr )
      endif

      endif


     
      end subroutine pwrite_hdf_raws

      subroutine pwrite_hdf_plasma_raws(path,file, tnps,tnpp,x0, dx,&
     &               name, n, t, dt, xmin,                        &
     & dataset, docompress,tag)
!---------------------------------------------------   
! this subroutine collects distributed 1d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:,:), pointer, intent(in) :: dataset
        real, dimension(:) :: x0, dx
        character(len=*), intent(in) :: path, file
        integer :: tnps, tnpp
        integer :: pgrp, nvp, idproc, tp, color, i, j
        
        integer, dimension(:,:), allocatable:: dims
        integer(hsize_t), dimension(1) :: ldim
        integer(hsize_t), dimension(2) :: tdim
        integer, dimension( : ), allocatable :: np
        real, dimension(:), allocatable :: write_buffer
        integer, dimension(:,:), allocatable :: tag_buffer
        
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(3) :: xmax
        logical, intent(in) :: docompress
        integer, dimension(:,:,:), pointer, intent(in), optional :: tag

        character(len_trim(path)+1+len_trim(file)) :: file_name

        integer(hsize_t), dimension(1) :: start, stride
        integer(hsize_t), dimension(2) :: tstart, tstride
        integer(hid_t) :: d_float
        logical :: ldoappend

        integer(HID_T), dimension(1) :: bc, mc 

! hdf5 variables
        integer(hid_t) :: flplID, xferID, dcplID, memspaceID, attrID 
        integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
        integer :: info, ierr 

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, ioff, id, nrec0
      dimension istatus(lstat)
      
      ldim(1) = 1

! hard wire boundary condition and moving window properties for the moment
      bc = 0
      mc = 0

      file_name = trim(path)//trim(file)

! If total number of particles to output is 0 then node 0 will create 
! an empty file
      call MPI_ALLREDUCE( tnpp, tp, 1, MPI_INTEGER, MPI_SUM, lgrp, ierr )

      if ( tp == 0 ) then

! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, lgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! create new file
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
!         set data label and units
      call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

      call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
      call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )

          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)


      else

! create a communicator that includes only nodes that have particles to output
         if ( tnpp > 0 ) then 
            color = 1
         else
            color = MPI_UNDEFINED
         endif

! this is a global call
         call MPI_COMM_SPLIT( lgrp, color, 0, pgrp, ierr )

         if (tnpp > 0) then

         call MPI_COMM_RANK( pgrp, idproc, ierr )

! get number of nodes with particles to output
		 call MPI_COMM_SIZE( pgrp, nvp, ierr )

! gather number of particles on all 
		 allocate( np( nvp ), stat = ierr )
		 allocate(dims( 2, nvp ), stat = ierr)
		 
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW diagnostic"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
		 
		 call MPI_ALLGATHER( tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER, &
		 &pgrp, ierr)
	 
		 ! get total number of particles and positions on file
		 dims( 1, 1 ) = 1
		 dims( 2, 1 ) = np(1) 
		 do i = 2, nvp
		    dims( 1, i ) = dims( 2, i-1 ) + 1
		    dims( 2, i ) = dims( 1, i ) + np(i) - 1
		 enddo
	 
		 ! get temporary memory for write buffer
		 allocate( write_buffer( tnpp ) , stat = ierr )
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW diagnostic"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
		 
		 if (present(tag)) then
		    allocate( tag_buffer( 2,tnpp ) , stat = ierr )
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW diagnostic"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
         endif


! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! create new file
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
!         set data label and units
      call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

      call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
      call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )

          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! synchronization
!      call MPI_BARRIER(lgrp,ierr)

! write data

         if (present(tag)) then
		 
	     tag_buffer(1,1:tnpp) = tag(1,tnps:(tnps+tnpp-1),1)
	     tag_buffer(2,1:tnpp) = tag(2,tnps:(tnps+tnpp-1),1)



            
	  tdim(1) = 2
	  tdim(2) = tp

      call h5screate_simple_f(2, tdim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'tag', H5T_NATIVE_INTEGER, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      tdim(2) = tnpp

! create memory dataspace
      call h5screate_simple_f(2, tdim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      tstart(1) = 0
      tstart(2) = dims(1,idproc+1) - 1
      tstride(1) = 1
      tstride(2) = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,tstart,  &
     & tstride, ierr, tstride, tdim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, tag_buffer, tdim,ierr,&
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
         
         endif
      

		 do i = 1, 3
		 
		 if (i == 3) then
	     write_buffer(:) = dx(i)		 
		 else
	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) - x0(i)) * dx(i)
         endif   
      ldim(1) = tp

       call h5screate_simple_f(1, ldim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing h5dwrite_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'x_'//char(iachar('0')+i)) 

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo


	  do i = 1, 3
		 
!			write_buffer(1:tnpp) = dataset((i+2),tnps:(tnps+tnpp-1),1)/&
!	 &dataset(6,tnps:(tnps+tnpp-1),1) * dx(1) 
			write_buffer(1:tnpp) = dataset((i+2),tnps:(tnps+tnpp-1),1)*&
	 &dx(1) 
            
      ldim(1) = tp

      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'p_'//char(iachar('0')+i)) 

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo

			write_buffer(1:tnpp) = dataset(6,tnps:(tnps+tnpp-1),1)
            
      ldim(1) = tp

      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'gamma', d_float, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error it h5dwrite_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error it h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'gamma') 

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)

			write_buffer(1:tnpp) = dataset(15,tnps:(tnps+tnpp-1),1)
            
      ldim(1) = tp

      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'onepluspsi', d_float, dspace_id, &
                           dset_id, ierr)
                           
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p')
         call add_h5_atribute(dset_id,'LONG_NAME', 'onepluspsi') 

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)

!       Compress SDS

      if (docompress) then 
        ! place holder        
      endif


! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      

      deallocate(np)
      deallocate(dims)
      deallocate(write_buffer)
      deallocate(tag_buffer)

      endif

      if ( pgrp /= MPI_COMM_NULL ) then
         call MPI_COMM_FREE( pgrp, ierr )
      endif

      endif


     
      end subroutine pwrite_hdf_plasma_raws
!-----------------------------------------------------------------------              
      subroutine ptrite_hdf_raw_tag(path, file, tnps,tnpp, soffset,x0, d&
     &x,               name, n, t, dt, xmin, dspl,                      &
     & dataset,doappend, docompress,tag)
! beam raw data diagnostic, for single stage
!---------------------------------------------------   
! this subroutine collects distributed 1d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:,:), pointer, intent(in) :: dataset
        integer, dimension(:,:,:), pointer, intent(in) :: tag
        integer, dimension(:), pointer :: ptrack
        real, intent(in) :: soffset
        real, dimension(:) :: x0, dx
        character(len=*), intent(in) :: path, file
        integer :: tnps, tnpp, dspl, tnp
        integer :: pgrp, nvp, idproc, tp, tpo, color, i, j
        
        integer(HSIZE_T), dimension(1) :: maxdim
        integer(HSIZE_T), dimension(2) :: maxtdim
        integer, dimension(:,:), allocatable:: dims
        integer(hsize_t), dimension(1) :: ldim
        integer(hsize_t), dimension(2) :: tdim
        integer, dimension( : ), allocatable :: np
        real, dimension(:), allocatable :: write_buffer
        integer, dimension(:,:), allocatable :: tag_buffer
        
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(3) :: xmax
        logical, intent(in) :: doappend,docompress

        character(len_trim(path)+1+len_trim(file)) :: file_name

        integer(hsize_t), dimension(1) :: start, stride
        integer(hsize_t), dimension(2) :: tstart, tstride
        integer(hid_t) :: d_float
        logical :: ldoappend

        integer(HID_T), dimension(1) :: bc, mc 

! hdf5 variables
        integer(hid_t) :: flplID, xferID, dcplID, memspaceID, attrID 
        integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
        integer :: info, ierr 

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, ioff, id, nrec0
      dimension istatus(lstat)
      
      tpo = 0
      tnp = 0
      ldim(1) = 1

! hard wire boundary condition and moving window properties for the moment
      bc = 0
      mc = 0

      file_name = trim(path)//trim(file)
      ldoappend = doappend

      
      allocate(ptrack(tnpp))
      
         do i=1, tnpp
            if (tag(3,i,1) == 1) then
               tnp = tnp + 1
               ptrack(tnp) = i
            endif
         enddo
      
      tnpp = tnp
      tnps = 1

! If total number of particles to output is 0 then node 0 will create 
! an empty file
      call MPI_ALLREDUCE( tnpp, tp, 1, MPI_INTEGER, MPI_SUM, lgrp, ierr )

      if ( tp == 0 ) then

! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, lgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)

! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
! create dataset in dataspace

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
          call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)
! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      
      return
      
      else

! create a communicator that includes only nodes that have particles to output
         if ( tnpp > 0 ) then 
            color = 1
         else
            color = MPI_UNDEFINED
         endif

! this is a global call
         call MPI_COMM_SPLIT( lgrp, color, 0, pgrp, ierr )

         if (tnpp > 0) then

         call MPI_COMM_RANK( pgrp, idproc, ierr )

! get number of nodes with particles to output
		 call MPI_COMM_SIZE( pgrp, nvp, ierr )
		 
! gather number of particles on all 
		 allocate( np( nvp ), stat = ierr )
		 allocate(dims( 2, nvp ), stat = ierr)
		 
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW write"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
		 
		 call MPI_ALLGATHER( tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER, &
		 &pgrp, ierr)
	 
		 ! get total number of particles and positions on file
		 dims( 1, 1 ) = 1
		 dims( 2, 1 ) = np(1) 
		 do i = 2, nvp
		    dims( 1, i ) = dims( 2, i-1 ) + 1
		    dims( 2, i ) = dims( 1, i ) + np(i) - 1
		 enddo
	 
		 ! get temporary memory for write buffer
		 allocate( write_buffer( tnpp ) , stat = ierr )
		 allocate( tag_buffer(2,tnpp ) , stat = ierr )
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW write"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif



! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! open file
! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
! create dataset in dataspace

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
          call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! synchronization
!      call MPI_BARRIER(lgrp,ierr)

! write data

		 do i = 1, 3
		 
         do j = 1, tnpp
      if (i == 3) then
!	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) + soffset - x0(i)) * dx(i)
	     write_buffer(j) = (dataset(i,ptrack(j),1) + soffset - x0(i)) * dx(i)
      else	     
!	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) - x0(i)) * dx(i)
	     write_buffer(j) = (dataset(i,ptrack(j),1) - x0(i)) * dx(i)
      endif
         enddo
      
      maxdim = (/H5S_UNLIMITED_F/)
      
      ldim(1) = 1

      call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      ldim(1) = tp

      call h5pset_chunk_f(dcplID, 1, ldim, ierr)

      call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr, dcplID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      ldim(1) = tp

      call h5dextend_f(dset_id, ldim, ierr)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dextend_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = tpo + dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'x_'//char(iachar('0')+i)) 

      call h5pclose_f(dcplID, ierr)

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo


	  do i = 1, 3
		 
			
         do j = 1, tnpp
            write_buffer(j) = dataset(i+3,ptrack(j),1)
         enddo			
			
      
      maxdim = (/H5S_UNLIMITED_F/)

      ldim(1) = 1

      call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      ldim(1) = tp

      call h5pset_chunk_f(dcplID, 1, ldim, ierr)

      call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr, dcplID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5dextend_f(dset_id, ldim, ierr)

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dextend_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = tpo + dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F, start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim, ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'p_'//char(iachar('0')+i)) 

      call h5pclose_f(dcplID, ierr)
	
      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)

	  enddo


         do j = 1, tnpp
            tag_buffer(1,j) = tag(1,ptrack(j),1)
            tag_buffer(2,j) = tag(2,ptrack(j),1)
         enddo
      
	  tdim(1) = 2
	  tdim(2) = tp

      call h5screate_simple_f(2, tdim, dspace_id, ierr)

      call h5dcreate_f(rootID, 'tag', H5T_NATIVE_INTEGER, dspace_id, &
                           dset_id, ierr)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      
      tdim(2) = tnpp

! create memory dataspace
      call h5screate_simple_f(2, tdim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      tstart(1) = 0
      tstart(2) = dims(1,idproc+1) - 1
      tstride(1) = 1
      tstride(2) = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,tstart,  &
     & tstride, ierr, tstride, tdim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, tag_buffer, tdim,ierr,&
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	

!       Compress SDS

      if (docompress) then 
        ! place holder        
      endif


! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      

      deallocate(np)
      deallocate(dims)
      deallocate(write_buffer)
      deallocate(tag_buffer)

      endif
      
      endif

      if ( pgrp /= MPI_COMM_NULL ) then
         call MPI_COMM_FREE( pgrp, ierr )
      endif

      deallocate(ptrack)
     
      end subroutine ptrite_hdf_raw_tag

!-----------------------------------------------------------------------              
      subroutine pwrite_hdf_raw_tag(path, file, tnps,tnpp, soffset,x0, d&
     &x,               name, n, t, dt, xmin, dspl,                      &
     & dataset,doappend, docompress,tag)
! beam raw data diagnostic, for multi stages
!---------------------------------------------------   
! this subroutine collects distributed 1d datasets and 
! writes to a hdf file
!---------------------------------------------------     
        implicit none
      
        real, dimension(:,:,:), pointer, intent(in) :: dataset
        integer, dimension(:,:,:), pointer, intent(in) :: tag
        integer, dimension(:), pointer :: ptrack
        real, intent(in) :: soffset
        real, dimension(:) :: x0, dx
        character(len=*), intent(in) :: path, file
        integer :: tnps, tnpp, dspl, tnp
        integer :: pgrp, nvp, idproc, tp, tpo, color, i, j
        
        integer(HSIZE_T), dimension(1) :: maxdim
        integer(HSIZE_T), dimension(2) :: maxtdim
        integer, dimension(:,:), allocatable:: dims
        integer(hsize_t), dimension(1) :: ldim
        integer(hsize_t), dimension(2) :: tdim
        integer, dimension( : ), allocatable :: np
        real, dimension(:), allocatable :: write_buffer
        integer, dimension(:,:), allocatable :: tag_buffer
        
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real, intent(in) :: t, dt
!        integer, dimension(:), intent(in) :: nx
        real, dimension(:), intent(in) :: xmin
        real, dimension(3) :: xmax
        logical, intent(in) :: doappend,docompress

        character(len_trim(path)+1+len_trim(file)) :: file_name

        integer(hsize_t), dimension(1) :: start, stride
        integer(hsize_t), dimension(2) :: tstart, tstride
        integer(hid_t) :: d_float
        logical :: ldoappend

        integer(HID_T), dimension(1) :: bc, mc 

! hdf5 variables
        integer(hid_t) :: flplID, xferID, dcplID, memspaceID, attrID 
        integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
        integer :: info, ierr 

! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=8)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, ioff, id, nrec0
      dimension istatus(lstat)
      
      tpo = 0
      tnp = 0
      ldim(1) = 1

! hard wire boundary condition and moving window properties for the moment
      bc = 0
      mc = 0

      file_name = trim(path)//trim(file)
      ldoappend = doappend
      
      allocate(ptrack(tnpp))
      
         do i=1, tnpp
            if (tag(3,i,1) == 1) then
               tnp = tnp + 1
               ptrack(tnp) = i
            endif
         enddo
      
      tnpp = tnp
      tnps = 1

! If total number of particles to output is 0 then node 0 will create 
! an empty file
      call MPI_ALLREDUCE( tnpp, tp, 1, MPI_INTEGER, MPI_SUM, lgrp, ierr )

      if ( tp == 0 ) then

! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, lgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)

! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
! create dataset in dataspace

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
          call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)
! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
          
      
      return
      
      else

! create a communicator that includes only nodes that have particles to output
         if ( tnpp > 0 ) then 
            color = 1
         else
            color = MPI_UNDEFINED
         endif

! this is a global call
         call MPI_COMM_SPLIT( lgrp, color, 0, pgrp, ierr )

         if (tnpp > 0) then

         call MPI_COMM_RANK( pgrp, idproc, ierr )


! get number of nodes with particles to output
		 call MPI_COMM_SIZE( pgrp, nvp, ierr )

		 
! gather number of particles on all 
		 allocate( np( nvp ), stat = ierr )
		 allocate(dims( 2, nvp ), stat = ierr)

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW write"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
		 
		 call MPI_ALLGATHER( tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER, &
		 &pgrp, ierr)
	 
		 ! get total number of particles and positions on file
		 dims( 1, 1 ) = 1
		 dims( 2, 1 ) = np(1) 
		 do i = 2, nvp
		    dims( 1, i ) = dims( 2, i-1 ) + 1
		    dims( 2, i ) = dims( 1, i ) + np(i) - 1
		 enddo
	 
		 ! get temporary memory for write buffer
		 allocate( write_buffer( tnpp ) , stat = ierr )
		 allocate( tag_buffer(2,tnpp ) , stat = ierr )
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error allocating memory for RAW write"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif



! init phdf5 interface
      call h5open_f(ierr)
!     Determine the precision of real variable
      d_float = detect_precision()

        call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_FILE_ACCESS_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

!        call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         

!        write(message,*) "h5pcreate_f ierr, dcplID=", ierr, dcplID
!        call Write_Log(message, 2)

! create property list for dataset transfer                  
        call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5pcreate_f, H5P_DATASET_XFER_F"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      info = MPI_INFO_NULL
      call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
      call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

! open file
! create new file
         call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr,    &
     & access_prp=flplID) 
! create dataset in dataspace

!         set data label and units
         call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin, &
     & xmax, bc, mc, 'particles' )

      call h5gopen_f( file_id, '/', rootID, ierr )

         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         
          call h5acreate_f( rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, tp, ldim, ierr)

          call h5aclose_f(attrID, ierr)

          call h5sclose_f(aspace_id, ierr)

! synchronization
!      call MPI_BARRIER(lgrp,ierr)

! write data

		 do i = 1, 3
		 
         do j = 1, tnpp
      if (i == 3) then
!	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) + soffset - x0(i)) * dx(i)
	     write_buffer(j) = (dataset(i,ptrack(j),1) + soffset - x0(i)) * dx(i)
      else	     
!	     write_buffer(1:tnpp) = (dataset(i,tnps:(tnps+tnpp-1),1) - x0(i)) * dx(i)
	     write_buffer(j) = (dataset(i,ptrack(j),1) - x0(i)) * dx(i)
      endif
         enddo
      

      
      maxdim = (/H5S_UNLIMITED_F/)
      
      ldim(1) = 1

      call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      ldim(1) = tp

      call h5pset_chunk_f(dcplID, 1, ldim, ierr)

      call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr, dcplID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      ldim(1) = tp

      call h5dextend_f(dset_id, ldim, ierr)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dextend_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = tpo + dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim,ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dwrite_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'x_'//char(iachar('0')+i)) 

      call h5pclose_f(dcplID, ierr)

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
	
	  enddo


	  do i = 1, 3
		 
			
         do j = 1, tnpp
            write_buffer(j) = dataset(i+3,ptrack(j),1)
         enddo			
			
      
      maxdim = (/H5S_UNLIMITED_F/)

      ldim(1) = 1

      call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      ldim(1) = tp

      call h5pset_chunk_f(dcplID, 1, ldim, ierr)

      call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), d_float, dspace_id, &
                           dset_id, ierr, dcplID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5dextend_f(dset_id, ldim, ierr)

		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dextend_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(1, ldim, dspace_id, ierr)

      ldim(1) = tnpp

! create memory dataspace
      call h5screate_simple_f(1, ldim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      start = tpo + dims(1,idproc+1) - 1
      stride = 1

!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F, start,   &
     & stride, ierr, stride, ldim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id, d_float, write_buffer, ldim, ierr,       &
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error writing dataset"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

         call add_h5_atribute(dset_id,'UNITS', 'c/\omega_p') 
         call add_h5_atribute(dset_id,'LONG_NAME', 'p_'//char(iachar('0')+i)) 

      call h5pclose_f(dcplID, ierr)
	
      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)

	  enddo


         do j = 1, tnpp
            tag_buffer(1,j) = tag(1,ptrack(j),1)
            tag_buffer(2,j) = tag(2,ptrack(j),1)
         enddo

      maxtdim = (/H5S_UNLIMITED_F,H5S_UNLIMITED_F/)
      
      tdim(1) = 2
      tdim(2) = 1

      call h5screate_simple_f(2, tdim, dspace_id, ierr, maxtdim)

!Modify dataset creation properties, i.e. enable chunking

      call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)

      tdim(2) = tp

      call h5pset_chunk_f(dcplID, 2, tdim, ierr)

      call h5dcreate_f(rootID, 'tag', H5T_NATIVE_INTEGER, dspace_id, &
                           dset_id, ierr, dcplID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dcreate_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      tdim(2) = tp

      call h5dextend_f(dset_id, tdim, ierr)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dextend_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif
      call h5sclose_f(dspace_id, ierr)
      
! create file dataspace
      call h5screate_simple_f(2, tdim, dspace_id, ierr)

      tdim(2) = tnpp

! create memory dataspace
      call h5screate_simple_f(2, tdim, memspaceID, ierr )

! select hyperslab in the file    
! hyperslab coordinates are 0 indexed
      tstart(1) = 0
      tstart(2) = tpo + dims(1,idproc+1) - 1
      tstride(1) = 1
      tstride(2) = 1


!      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f( dspace_id, H5S_SELECT_SET_F,tstart,  &
     & tstride, ierr, tstride, tdim)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5sselect_hyperslab_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

! write the dataset collectively  
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER, tag_buffer, tdim,ierr,&
     &     memspaceID, dspace_id, xfer_prp=xferID)
		 if (ierr /= 0) then 
          print*, "p_os_util_hdf5: Error in h5dwrite_f"
          !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
		 endif

      call h5pclose_f(dcplID, ierr)

      call h5sclose_f(memspaceID, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      


!       Compress SDS

      if (docompress) then 
        ! place holder        
      endif


! close resources
      call h5pclose_f(xferID, ierr)
      call h5pclose_f(flplID, ierr)
! close dataset and file
      call h5gclose_f(rootID, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      

      deallocate(np)
      deallocate(dims)
      deallocate(write_buffer)
      deallocate(tag_buffer)

      endif
      
      endif

      if ( pgrp /= MPI_COMM_NULL ) then
         call MPI_COMM_FREE( pgrp, ierr )
      endif

      deallocate(ptrack)
     
      end subroutine pwrite_hdf_raw_tag

      end module m_pdiagnostic_utilities
