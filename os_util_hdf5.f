      module m_h5_diagnostic_utilities
         
        use hdf5

        implicit none
        private

        public :: start_hdf5, stop_hdf5
        public :: write_hdf_sds, read_hdf_sds
        public :: detect_precision
        public :: add_h5_axes, add_h5_attrs, open_h5_dataset
        public :: add_h5_atribute

        interface write_hdf_sds
          module procedure write_hdf_sds1D
          module procedure write_hdf_sds2D
          module procedure write_hdf_sds3D
        end interface

        interface read_hdf_sds
          module procedure read_hdf_sds2D
          module procedure read_hdf_sds3D
        end interface

       interface add_h5_axes
          module procedure add_h5_axes_serial
       end interface

       interface add_h5_attrs
          module procedure add_h5_attrs_serial
       end interface

       interface open_h5_dataset 
          module procedure open_h5_file_dataset
       end interface

! osiris hdf5 data attribute interface
! copy over from osiris os-dutil-hdf5.f90 to ensure compatible data format with visXD
       interface add_h5_atribute
         module procedure add_h5_atribute_str
         module procedure add_h5_atribute_str_v1
         module procedure add_h5_atribute_logical
         module procedure add_h5_atribute_logical_v1
       
         module procedure add_h5_atribute_single
         module procedure add_h5_atribute_v1_single
       
!         module procedure add_h5_atribute_double
!         module procedure add_h5_atribute_v1_double
       
         module procedure add_h5_atribute_int
         module procedure add_h5_atribute_v1_int
       
       end interface

       contains

!---------------------------------------------------
        subroutine start_hdf5(ierr)
             implicit none
   
            integer, intent(out) :: ierr
   
            call h5open_f(ierr) 
        end subroutine start_hdf5

!---------------------------------------------------
        subroutine stop_hdf5(ierr)
             implicit none
   
            integer, intent(out) :: ierr
   
            call h5close_f(ierr) 
        end subroutine stop_hdf5

!---------------------------------------------------
        subroutine write_hdf_sds1D( path, file,                         &
     &                              dataset, name, label, units,        &
     &                              n, t, dt, nx, xmin, dx,             &
     &                              xname, xlabel, xunits,              &
     &                              docompress)
!---------------------------------------------------
!       this routine writes a visXD hdf5 file for 1D dataset
!---------------------------------------------------
          implicit none
  
          real, dimension(:), intent(in) :: dataset
          character(len=*),   intent(in) :: path, file
          character(len=*),   intent(in) :: name, label, units
          integer,               intent(in) :: n
          real,                  intent(in) :: t, dt
          integer, dimension(:), intent(in) :: nx
          real, dimension(:),    intent(in) :: xmin
          real, dimension(:),    intent(in) :: dx
          character(len=*),dimension(:), intent(in) :: xname, xlabel, xu&
     &nits
          logical, intent(in) :: docompress

!       local variables
          integer, parameter :: rank = 1

          character(len_trim(path)+1+len_trim(file)) :: filename
          integer(HID_T) :: file_id, dspace_id, dset_id, rootID
          real, dimension(rank) :: xmax 
          integer(HID_T), dimension(rank) :: bc, mc 

          integer(HID_T) :: d_float
          integer :: ierror
          integer(hsize_t), dimension(rank) :: dims

!         generate filename
          filename = trim(path) // trim(file)

!         generate xmax
          xmax = xmin + dx*nx

! hard wire boundary condition and moving window properties for the moment
          bc = 0
          mc = 0

          ! Create HDF file
         
          call h5open_f(ierror)
!       Determine the precision of real variable
          d_float = detect_precision()
          call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierror)

          ! Create dataspace 
          dims = nx
          call h5screate_simple_f(rank, dims, dspace_id, ierror) 
   
          ! create dataset
          call h5gopen_f( file_id, '/', rootID, ierror )
          call h5dcreate_f(rootID, name, d_float, dspace_id,            &
     & dset_id, ierror)
          call h5gclose_f(rootID, ierror)

          ! Compress SDS
          if (docompress) then
          ! place holder
          endif

!         write data and attributes
          call h5dwrite_f(dset_id, d_float, dataset, dims, ierror)
          call add_h5_atribute(dset_id,'UNITS', units) 
          call add_h5_atribute(dset_id,'LONG_NAME', label) 

!         set file attributes 
          call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin,&
     & xmax, bc, mc )

! add data axes 
          call add_h5_axes(file_id, rank, xunits, xname, xlabel, xmin,  &
     & xmax)

          ! Close dataset, dataspace, file, fortran interface 
          call h5dclose_f(dset_id, ierror)
          call h5sclose_f(dspace_id, ierror)
          call h5fclose_f(file_id, ierror)
          call h5close_f(ierror)

        end subroutine write_hdf_sds1D
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_hdf_sds2D( path, file,                         &
     &                              dataset, name, label, units,        &
     &                              n, t, dt, nx, xmin, dx,             &
     &                              xname, xlabel, xunits,              &
     &                              docompress)
!---------------------------------------------------
!       this routine writes a visXD hdf5 file for 2D dataset
!---------------------------------------------------
          implicit none
  
          real, dimension(:,:), intent(in) :: dataset
          character(len=*),   intent(in) :: path, file
          character(len=*),   intent(in) :: name, label, units
          integer,               intent(in) :: n
          real,                  intent(in) :: t, dt
          integer, dimension(:), intent(in) :: nx
          real, dimension(:),    intent(in) :: xmin
          real, dimension(:),    intent(in) :: dx
          character(len=*),dimension(:), intent(in) :: xname, xlabel, xu&
     &nits
          logical, intent(in) :: docompress

!       local variables
          integer, parameter :: rank = 2 

          character(len_trim(path)+1+len_trim(file)) :: filename
          integer(HID_T) :: file_id, dspace_id, dset_id, rootID
          real, dimension(rank) :: xmax 
          integer(HID_T), dimension(rank) :: bc, mc 

          integer(HID_T) :: d_float
          integer :: ierror, ierror2
          integer(hsize_t), dimension(rank) :: dims

!         generate filename
          filename = trim(path) // trim(file)

!         generate xmax
          xmax = xmin + dx*nx

! hard wire boundary condition and moving window properties for the moment
          bc = 0
          mc = 0

          ! Create HDF file
         
          call h5open_f(ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error starting HDF5"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
!       Determine the precision of real variable
          d_float = detect_precision()
          call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error creating file ", filename
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

          ! Create dataspace 
          dims = nx
          call h5screate_simple_f(rank, dims, dspace_id, ierror) 
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error creating dataspace"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
   
          ! create dataset
          call h5gopen_f( file_id, '/', rootID, ierror )
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error in h5gopen_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call h5dcreate_f(rootID, name, d_float, dspace_id,            &
     & dset_id, ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error in h5dcreate_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call h5gclose_f(rootID, ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error in h5gclose_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

          ! Compress SDS
          if (docompress) then
          ! place holder
          endif

!         write data and attributes
          call h5dwrite_f(dset_id, d_float, dataset, dims, ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error writing dataset"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call add_h5_atribute(dset_id,'UNITS', units) 
          call add_h5_atribute(dset_id,'LONG_NAME', label) 

!         set data label and units
          call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin,&
     & xmax, bc, mc )

! add data axes 
          call add_h5_axes(file_id, rank, xunits, xname, xlabel, xmin,  &
     & xmax)

          ! Close dataset, dataspace, file, fortran interface 
          call h5dclose_f(dset_id, ierror)
          call h5sclose_f(dspace_id, ierror)
          call h5fclose_f(file_id, ierror)
          call h5close_f(ierror)

        end subroutine write_hdf_sds2D
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_hdf_sds3D( path, file,                         &
     &                              dataset, name, label, units,        &
     &                              n, t, dt, nx, xmin, dx,             &
     &                              xname, xlabel, xunits,              &
     &                              docompress)
!---------------------------------------------------
!       this routine writes a visXD hdf5 file for 3D dataset
!---------------------------------------------------
          implicit none
  
          real, dimension(:,:,:), intent(in) :: dataset
          character(len=*),   intent(in) :: path, file
          character(len=*),   intent(in) :: name, label, units
          integer,               intent(in) :: n
          real,                  intent(in) :: t, dt
          integer, dimension(:), intent(in) :: nx
          real, dimension(:),    intent(in) :: xmin
          real, dimension(:),    intent(in) :: dx
          character(len=*),dimension(:), intent(in) :: xname, xlabel, xu&
     &nits
          logical, intent(in) :: docompress

!       local variables
          integer, parameter :: rank = 3 

          character(len_trim(path)+1+len_trim(file)) :: filename
          integer(HID_T) :: file_id, dspace_id, dset_id, rootID
          real, dimension(rank) :: xmax 
          integer, dimension(rank) :: bc, mc 

          integer(HID_T) :: d_float
          integer(hsize_t), dimension(rank) :: dims
          integer :: ierror

!         generate filename
          filename = trim(path) // trim(file)

!         generate xmax
          xmax = xmin + dx*nx

! hard wire boundary condition and moving window properties for the moment
          bc = 0
          mc = 0

          ! Create HDF file
         
          call h5open_f(ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error starting HDF5"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
!       Determine the precision of real variable
          d_float = detect_precision()
          call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierror)
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error creating file ", filename
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

          ! Create dataspace 
          dims = nx
          call h5screate_simple_f(rank, dims, dspace_id, ierror) 
         if (ierror.ne.0) then 
            print*, "os_util_hdf5: Error creating dataspace"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
   
          ! create dataset
          call h5gopen_f( file_id, '/', rootID, ierror )
          call h5dcreate_f(rootID, name, d_float, dspace_id,       &
     & dset_id, ierror)
          call h5gclose_f(rootID, ierror)

          ! Compress SDS
          if (docompress) then
          ! place holder
          endif

!         write data and attributes
          call h5dwrite_f(dset_id, d_float, dataset, dims, ierror)
          call add_h5_atribute(dset_id,'UNITS', units) 
          call add_h5_atribute(dset_id,'LONG_NAME', label) 

!         set data label and units
          call add_h5_attrs(file_id, name, t, n, dt, '1/\omega_p', xmin,&
     & xmax, bc, mc )

! add data axes 
          call add_h5_axes(file_id, rank, xunits, xname, xlabel, xmin,  &
     & xmax)

          ! Close dataset, dataspace, file, fortran interface 
          call h5dclose_f(dset_id, ierror)
          call h5sclose_f(dspace_id, ierror)
          call h5fclose_f(file_id, ierror)
          call h5close_f(ierror)

        end subroutine write_hdf_sds3D
!---------------------------------------------------


!---------------------------------------------------
        subroutine read_hdf_sds2D( path, file,                          &
     &                              dataset, timestep, t, ierr,         &
     &                             info_only,dim_sizes,data_type,n_atts)
!---------------------------------------------------
!       this routine reads a HDF SD File into a 2D dataset
!---------------------------------------------------

!         <use-statements for this subroutine>
          
          implicit none
  
          character(len=*),              intent(in) :: path, file
          real, dimension(:,:), intent(inout) :: dataset
          integer,               intent(inout) :: timestep
          real,               intent(inout) :: t
          integer, intent(out) :: ierr 
          logical, intent(in) :: info_only
          integer, dimension(:), intent(inout), optional :: dim_sizes 
          integer, intent(inout), optional :: data_type, n_atts

!       local constants       
          integer, parameter                        :: rank = 2        

!       local variables

          character(len_trim(path)+1+len_trim(file)) :: filename
          integer :: i,attr_index
          integer :: status, irank, idata_type, in_atts
          character(len=80) :: name
          integer(hid_t) :: file_id, rootID, attrID, stringType
          integer(hid_t) :: attr_step_id, attr_time_id 
          integer(hid_t) :: dset_id, dspace_id 
          character(len=80) :: string
          integer(hsize_t) :: str_len
          integer(hsize_t), dimension(1) :: dims 
          integer(HID_T) :: d_float
          integer(hsize_t), dimension(rank) :: idim_sizes, max_dims

!         executable statements
          ierr = 0

          filename = trim(path) // trim(file)

! open file
          call h5open_f(ierr)
          d_float = detect_precision()
          call h5fopen_f(filename,H5F_ACC_RDONLY_F, file_id, ierr)
!          call h5gopen_f(file_id, "/", rootID, ierr)
          
! open dataset in the file
          call open_h5_dataset(file_id, dset_id)

! get the properties of the dataset 
          call H5Dget_space_f(dset_id, dspace_id, ierr)
          call H5Sget_simple_extent_ndims_f(dspace_id, irank, ierr);
          call H5Sget_simple_extent_dims_f(dspace_id, idim_sizes,       &
     & max_dims, ierr);

! read dataset
          if (.not.(info_only)) then
             call H5Dread_f(dset_id, d_float, dataset, idim_sizes,      &
     & ierr)
             call h5dclose_f(dset_id, ierr)
          endif

!         read attribute ITER and TIME
          dims = 1
          call h5aopen_by_name_f(file_id, "/","ITER",attr_step_ID, ierr)
          call h5aread_f(attr_step_ID, H5T_NATIVE_INTEGER, timestep,    &
     & dims, ierr)  
          call h5aopen_by_name_f(file_id, "/","TIME",attr_time_ID, ierr)
          call h5aread_f(attr_time_ID, d_float , t, dims, ierr)  

! close dataset and file
          if (present(dim_sizes)) dim_sizes=idim_sizes
!          if (present(data_type)) data_type = idata_type 
!          if (present(n_atts))  n_atts=in_atts 
          
! 	  print *, "ITER = ", timestep
! 	  print *, "TIME = ", t
! 	  print *, "status = ", status	  

! close file
        call h5sclose_f(dspace_id, ierr)
        call h5fclose_f(file_id, ierr)
        call h5close_f(ierr)

        end subroutine read_hdf_sds2D
!---------------------------------------------------


!---------------------------------------------------
        subroutine read_hdf_sds3D( path, file,                          &
     &                              dataset, timestep, t, ierr,         &
     &                             info_only,dim_sizes,data_type,n_atts)
!---------------------------------------------------
!       this routine reads a HDF SD File into a 2D dataset
!---------------------------------------------------

!         <use-statements for this subroutine>
          
          implicit none
  
          character(len=*),              intent(in) :: path, file
          real, dimension(:,:,:), intent(inout) :: dataset
          integer,               intent(inout) :: timestep
          real,               intent(inout) :: t
          integer, intent(out) :: ierr 
          logical, intent(in) :: info_only
          integer, dimension(:), intent(inout),optional :: dim_sizes 
          integer, optional :: data_type, n_atts

!       local constants       
          integer, parameter                        :: rank = 3        

!       local variables

          character(len_trim(path)+1+len_trim(file)) :: filename
          integer :: i,attr_index
          integer :: status, irank, idata_type, in_atts
          character(len=80) :: name
          integer(hid_t) :: file_id, rootID, attrID, stringType
          integer(hid_t) :: attr_step_id, attr_time_id
          integer(hid_t) :: dset_id, dspace_id
          character(len=80) :: string
          integer(hsize_t) :: str_len
          integer(hsize_t), dimension(1) :: dims 
          integer(HID_T) :: d_float
          integer(hsize_t), dimension(rank) :: idim_sizes, max_dims

!         executable statements
          ierr = 0

          filename = trim(path) // trim(file)

! open file
          call h5open_f(ierr)
          d_float = detect_precision()
          call h5fopen_f(filename,H5F_ACC_RDONLY_F, file_id, ierr)
!          call h5gopen_f(file_id, "/", rootID, ierr)
          
! open dataset in the file
          call open_h5_dataset(file_id, dset_id)

! get the properties of the dataset 
          call H5Dget_space_f(dset_id, dspace_id, ierr)
          call H5Sget_simple_extent_ndims_f(dspace_id, irank, ierr);
          call H5Sget_simple_extent_dims_f(dspace_id, idim_sizes,       &
     & max_dims, ierr);

! read dataset
          if (.not.(info_only)) then
             call H5Dread_f(dset_id, d_float, dataset, idim_sizes,      &
     & ierr)
             call h5dclose_f(dset_id, ierr)
          endif

!         read attribute ITER and TIME
          dims = 1
          call h5aopen_by_name_f(file_id, "/","ITER",attr_step_ID, ierr)
          call h5aread_f(attr_step_ID, H5T_NATIVE_INTEGER, timestep,    &
     & dims, ierr)  
          call h5aopen_by_name_f(file_id, "/","TIME",attr_time_ID, ierr)
          call h5aread_f(attr_time_ID, d_float , t, dims, ierr)  

! close dataset and file
          if (present(dim_sizes)) dim_sizes=idim_sizes
!          if (present(data_type)) data_type = idata_type 
!          if (present(n_atts))  n_atts=in_atts 
          
! 	  print *, "ITER = ", timestep
! 	  print *, "TIME = ", t
! 	  print *, "status = ", status	  

! close file
        call h5sclose_f(dspace_id, ierr)
        call h5fclose_f(file_id, ierr)
        call h5close_f(ierr)

        end subroutine read_hdf_sds3D
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine open_h5_file_dataset(file_id, dset_id)
! open default dataset in hdf5 file        
!---------------------------------------------------------------------------------------------------
          integer(hid_t), intent(in) :: file_id
          integer(hid_t), intent(inout) :: dset_id

          integer(hid_t) :: stringType, attrID, ierr 
          character(len=80) :: string
          integer(size_t) :: str_len
!          integer :: str_len
          integer(hsize_t), dimension(1) :: dims 

! first read "NAME" attribute
          call h5aopen_by_name_f(file_id, "/", "NAME", attrID, ierr)
          str_len = len(string)
          call h5tcopy_f(H5T_NATIVE_CHARACTER, stringType, ierr)
          call h5tset_size_f(stringType, str_len, ierr)
          dims = 1
          call h5aread_f(attrID, stringType, string, dims, ierr)  

! open the dataset with the "NAME" attribute
          call h5dopen_f(file_id, string, dset_id, ierr, H5P_DEFAULT_F)
          call h5aclose_f(attrID, ierr)
          call h5tclose_f(stringType, ierr)

        end subroutine open_h5_file_dataset 
!---------------------------------------------------



!---------------------------------------------------------------------------------------------------
        subroutine add_h5_axes_serial( file_id, ndims, xunits, xname,   &
     & xlabel, xmin, xmax, xfer_prop)
!---------------------------------------------------------------------------------------------------
! add visXD data axes via hdf5 serial interface
! modified from osiris os-dutil.f90i open_diag_file function
  
          integer, intent(in) :: ndims
          integer(hid_t) :: file_id 
          real, dimension(:), intent(in) :: xmax, xmin
          character(len=*), dimension(:), intent(in) :: xname, xlabel,  &
     &xunits

          integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID
          integer(hid_t), optional :: xfer_prop
          integer(hid_t) :: plistID
          integer(hsize_t), dimension(1) :: dims
          integer, parameter :: izero = ichar('0')
!          integer, parameter :: double=selected_real_kind(15, 307)
!          real(kind=double), dimension(2) :: axis_range
          real, dimension(2) :: axis_range
          integer(hid_t) :: d_float 
          integer :: i, ierr

          d_float = detect_precision()

          if (present(xfer_prop)) then 
             plistID = xfer_prop
          else
! default serial write property
             plistID = H5P_DEFAULT_F
          endif

          call h5gopen_f( file_id, '/', rootID, ierr )
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: Error in h5gopen_f"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

          ! add axis information
          call h5gcreate_f( rootID, 'AXIS', axisGroupID, ierr) 

          dims(1) = 2 
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          do i = 1, ndims
!             call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+i),      &
!     & H5T_NATIVE_DOUBLE, dataspaceID,   datasetID, ierr )
             call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+i),      &
     & d_float, dataspaceID,   datasetID, ierr )
            if (ierr.ne.0) then 
               print*, "os_util_hdf5: h5dcreate_f encountered an error creating an axis"
               !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
            endif

             call add_h5_atribute( datasetID, 'TYPE', 'linear' ) 
             call add_h5_atribute( datasetID, 'UNITS', xunits(i) ) 
             call add_h5_atribute( datasetID, 'NAME', xname(i) ) 
             call add_h5_atribute( datasetID, 'LONG_NAME', xlabel(i) ) 

             axis_range(1) = xmin( i )
             axis_range(2) = xmax( i )

!             call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, real(axis_range), &
!     & dims, ierr, xfer_prp = plistID )
             call h5dwrite_f( datasetID, d_float, axis_range,           &
     & dims, ierr, xfer_prp = plistID )
            if (ierr.ne.0) then 
               print*, "os_util_hdf5: h5dwrite_f encountered an error writing an axis"
               !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
            endif

             call h5dclose_f( datasetID, ierr )
          enddo
  
          call h5sclose_f( dataspaceID, ierr ) 
          call h5gclose_f( rootID, ierr ) 
          call h5gclose_f( axisGroupID, ierr ) 
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5gclose_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

        end subroutine add_h5_axes_serial 
!---------------------------------------------------



!---------------------------------------------------------------------------------------------------
        subroutine add_h5_attrs_serial( file_id, name, t, n, dt,        &
     & timeUnits, xmin, xmax, bc, mc, ty )
!---------------------------------------------------------------------------------------------------
! add visXD data attributes via hdf5 serial interface
  
          integer(hid_t), intent(in) :: file_id 
          character(len=*), intent(in) :: name, timeUnits
          integer, intent(in) :: n
          real, intent(in) :: t, dt
          real, dimension(:), intent(in) :: xmin, xmax
          integer, dimension(:), intent(in) :: bc, mc
          character(len=*), intent(in), optional :: ty
          

          integer(hid_t) :: rootID
          integer, parameter :: izero = ichar('0')
          integer :: ierr

          call h5gopen_f( file_id, '/', rootID, ierr )
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5gopen_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

          ! add name property
          call add_h5_atribute( rootID, 'NAME', name ) 
          if(.not. present(ty)) then
          call add_h5_atribute( rootID, 'TYPE', 'grid' ) 
          else
          call add_h5_atribute( rootID, 'TYPE', ty ) 
          endif
          call add_h5_atribute( rootID, 'TIME', t )
          call add_h5_atribute( rootID, 'ITER', n ) 
          call add_h5_atribute( rootID, 'DT', dt )
          call add_h5_atribute( rootID, 'TIME UNITS', timeUnits ) 
          call add_h5_atribute( rootID, 'XMIN', xmin ) 
          call add_h5_atribute( rootID, 'XMAX', xmax ) 
          call add_h5_atribute( rootID, 'PERIODIC', bc ) 
          call add_h5_atribute( rootID, 'MOVE C', mc ) 
  
          call h5gclose_f( rootID, ierr )
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5gclose_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif

        end subroutine add_h5_attrs_serial 
!---------------------------------------------------



!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_str( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          character( len = * ), intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, typeID, attrID
          integer(hsize_t), dimension(1) :: dims
          integer(size_t) :: size
          
          integer :: ierr
          
          dims(1) = 1
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
          
          size = len(attribute)
          call h5tset_size_f(typeID, size, ierr)
          
          call h5acreate_f( objID, name, typeID, dataspaceID, attrID,   &
     &ierr )
          call h5awrite_f( attrID, typeID, attribute, dims, ierr)
          call h5aclose_f( attrID, ierr )
          call h5tclose_f( typeID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
        end subroutine add_h5_atribute_str
!---------------------------------------------------------------------------------------------------
        
!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_str_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          character( len = * ), dimension(:), intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, typeID, attrID
          integer(hsize_t), dimension(1) :: dims
          integer(size_t) :: maxlen
          integer :: i, ierr
          
          dims(1) = size(attribute)
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          
          maxlen = 0
          do i = 1, size(attribute)-1
            if (len(attribute(i)) > maxlen) maxlen = len(attribute(i))
          enddo
          
          call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
          call h5tset_size_f(typeID, maxlen, ierr)
         
          call h5acreate_f( objID, name, typeID, dataspaceID, attrID,   &
     & ierr )
          call h5awrite_f( attrID, typeID, attribute, dims, ierr)
          call h5aclose_f( attrID, ierr )
          call h5tclose_f( typeID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
        end subroutine add_h5_atribute_str_v1
!---------------------------------------------------------------------------------------------------
        
!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_logical( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          logical, intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, attrID
          integer(hsize_t), dimension(1) :: dims
          integer :: bool, ierr
          
          dims(1) = 1
          if ( attribute ) then
            bool = 1
          else
            bool = 0
          endif
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
          call h5aclose_f( attrID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
        end subroutine add_h5_atribute_logical
!---------------------------------------------------------------------------------------------------
        
!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_logical_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          logical, dimension(:), intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, attrID
          integer :: i, ierr
          integer(hsize_t), dimension(1) :: dims
          integer, dimension(:), allocatable :: bool
          
          dims(1) = size(attribute)
          allocate( bool(dims(1)) )
          do i = 1, dims(1)
             if ( attribute(i) ) then
                bool(i) = 1
             else
                bool(i) = 0
             endif
          enddo
          
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
          call h5aclose_f( attrID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
          deallocate(bool)
          
        end subroutine add_h5_atribute_logical_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_single( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          real, intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, attrID
          integer(hid_t) :: d_float 
          integer(hsize_t), dimension(1) :: dims
          integer :: ierr
          
          d_float = detect_precision()
          dims(1) = 1
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5screate_simple_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call h5acreate_f( objID, name, d_float, dataspaceID,          & 
     &attrID, ierr )
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5acreate_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call h5awrite_f( attrID, d_float, attribute, dims,            &
     &ierr)
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5awrite_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call h5aclose_f( attrID, ierr )
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5aclose_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
          call h5sclose_f( dataspaceID, ierr )
         if (ierr.ne.0) then 
            print*, "os_util_hdf5: h5sclose_f encountered an error"
            !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
         endif
        
        end subroutine add_h5_atribute_single
!---------------------------------------------------------------------------------------------------
        
!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_v1_single( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          real, dimension(:), intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, attrID
          integer(hid_t) :: d_float
          integer(hsize_t), dimension(1) :: dims
          integer :: ierr
          
          d_float = detect_precision()
          dims(1) = size(attribute)
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          call h5acreate_f( objID, name, d_float, dataspaceID,          &
     & attrID, ierr )
          call h5awrite_f( attrID, d_float, attribute, dims,            &
     & ierr)
          call h5aclose_f( attrID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
        end subroutine add_h5_atribute_v1_single
!---------------------------------------------------------------------------------------------------


!!---------------------------------------------------------------------------------------------------
!        subroutine add_h5_atribute_double( objID, name, attribute )
!!---------------------------------------------------------------------------------------------------
!        
!          implicit none
!          
!          integer(hid_t), intent(in) :: objID
!          character( len = * ), intent(in) :: name
!          double precision, intent(in) :: attribute
!          
!          integer(hid_t) :: dataspaceID, attrID
!          integer(hsize_t), dimension(1) :: dims
!          integer :: ierr
!          
!          dims(1) = 1
!          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
!          call h5acreate_f( objID, name, H5T_NATIVE_DOUBLE, dataspaceID,&
!     &attrID, ierr )
!          call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, real(attribute), dims,  &
!     &ierr)
!          call h5aclose_f( attrID, ierr )
!          call h5sclose_f( dataspaceID, ierr )
!        
!        end subroutine add_h5_atribute_double
!!---------------------------------------------------------------------------------------------------
!        
!!---------------------------------------------------------------------------------------------------
!        subroutine add_h5_atribute_v1_double( objID, name, attribute )
!!---------------------------------------------------------------------------------------------------
!        
!          implicit none
!          
!          integer(hid_t), intent(in) :: objID
!          character( len = * ), intent(in) :: name
!          double precision,  dimension(:), intent(in) :: attribute
!          
!          integer(hid_t) :: dataspaceID, attrID
!          integer(hsize_t), dimension(1) :: dims
!          integer :: ierr
!          
!          dims(1) = size(attribute)
!          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
!          call h5acreate_f( objID, name, H5T_NATIVE_DOUBLE, dataspaceID,&
!     & attrID, ierr )
!          call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, real(attribute), dims,  &
!     & ierr)
!          call h5aclose_f( attrID, ierr )
!          call h5sclose_f( dataspaceID, ierr )
!        
!        end subroutine add_h5_atribute_v1_double
!!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_int( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          integer, intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, attrID
          integer(hsize_t), dimension(1) :: dims
          integer :: ierr
          
          dims(1) = 1
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, &
     &ierr)
          call h5aclose_f( attrID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
        end subroutine add_h5_atribute_int
!---------------------------------------------------------------------------------------------------
        
!---------------------------------------------------------------------------------------------------
        subroutine add_h5_atribute_v1_int( objID, name, attribute )
!---------------------------------------------------------------------------------------------------
        
          implicit none
          
          integer(hid_t), intent(in) :: objID
          character( len = * ), intent(in) :: name
          integer, dimension(:), intent(in) :: attribute
          
          integer(hid_t) :: dataspaceID, attrID
          integer(hsize_t), dimension(1) :: dims
          integer :: ierr
          
          dims(1) = size(attribute)
          call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID&
     &, attrID, ierr )
          call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, &
     & ierr)
          call h5aclose_f( attrID, ierr )
          call h5sclose_f( dataspaceID, ierr )
        
        end subroutine add_h5_atribute_v1_int
!---------------------------------------------------------------------------------------------------


      function detect_precision() result(prec)
      real :: small
      integer(HID_T) :: prec
      integer :: ierror

!      logical, save :: h5_inited = .false. 

!      if (.not. h5_inited) then
!         call h5open_f(ierror)
!         write(message,*) "h5open_f ierr=", ierror
!         call write_log(message,0)
!         h5_inited = .true.
!      endif 
! detect the precision of the system
      small = 1e-12
      small = 1.0 + small
      if (small>1.0) then 
          prec = H5T_NATIVE_DOUBLE 
      else
          prec = H5T_NATIVE_REAL
      endif
!      call h5close_f(ierror)
!      write(message,*) "h5close_f ierr=", ierror
!      call write_log(message,0)

      end function detect_precision

!---------------------------------------------------
      end module m_h5_diagnostic_utilities
