! serial HDF5 routines, used for restarts in QuickPIC      

      module m_data_utils
      
      use  m_h5_diagnostic_utilities
      use  m_string

      include 'mpif.h'

      interface write_data
        module procedure write_hdf3d
        module procedure write_hdf2d
      end interface

      interface read_data
        module procedure read_hdf2d
        module procedure read_hdf3d
      end interface
      
      contains

! ----------------------------------------------------------------------      
      subroutine write_hdf3d(datain,dataname, timestep,dt,path,label,uni&
     &ts,xlabel,xunits,xmin,dx,docompress)
! ----------------------------------------------------------------------      

      implicit none
      
      real, dimension(:,:,:), intent(in) :: datain
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      character(20)  :: chindx
      real :: dt
      character(len=*), optional :: path 
      character(len=80), optional ::  label, units
      character(len=80),dimension(3),optional:: xlabel, xunits
      real, dimension(3), optional :: xmin
      real, dimension(3), optional :: dx  
      logical, optional :: docompress

      character(len=len(label)) ::  l_label
      character(len=len(units)) :: l_units
      character(len=80),dimension(3) :: l_xlabel, l_xunits
      real, dimension(3) :: l_xmin
      real, dimension(3) :: l_dx  
      logical :: l_docompress

      integer, dimension(3) :: ssize
      character(len=80) :: l_path, file
      character(len=80) :: name
      integer :: n
      real :: t
      character(len=80),dimension(3) :: xname
        
      call int2str(timestep,chindx,4)
     
      file=trim(dataname)//'_'//trim(chindx)//'.h5'
      name=trim(dataname)//'-'//trim(chindx)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
      if (.not.present(docompress)) then 
         l_docompress=.true.
      else 
         l_docompress=docompress
      endif
      if (.not.present(label)) then 
         l_label=name
      else
         l_label=label
      endif 
      if (.not.present(units)) then 
         l_units='arb.'
      else
         l_units=units
      endif      
      n=timestep
      t=dt * n
        
      if (.not.present(xmin)) then 
         l_xmin=0.
      else
         l_xmin=xmin
      endif
      if (.not.present(dx)) then 
         l_dx=1.
      else
         l_dx=dx
      endif
     
      xname(1)='x1'
      xname(2)='x2'
      xname(3)='x3'
      if (.not.present(xlabel)) then 
         l_xlabel=xname
      else
         l_xlabel=xlabel
      endif
 
      if (.not.present(xunits)) then 
         l_xunits='arb.' 
      else
         l_xunits=xunits
      endif
      !write(*,*) 'Size  : ', size(datain)
      !write(*,*) 'Size1 : ', size(datain,1)
      !write(*,*) 'Size2 : ', size(datain,2)
      !write(*,*) 'Size3 : ', size(datain,3)
      
        
      ssize(1) = size(datain,1)
      ssize(2) = size(datain,2)
      ssize(3) = size(datain,3)
      
      if (.not.present(path)) then
         l_path = ''
      else
         l_path = trim(path)
      endif
      
      
        call write_hdf_sds( l_path, file, &
      &                     datain, name, l_label, l_units,&
      &                     n, t, dt, ssize, l_xmin, l_dx, &
      &                     xname, l_xlabel, l_xunits, l_docompress)

      end subroutine write_hdf3d
     
     
! ----------------------------------------------------------------------      
      subroutine write_hdf2d(datain, dataname,timestep,dt,path, label,un&
     &its,xlabel,xunits,xmin,dx, docompress)
! ----------------------------------------------------------------------      

      implicit none
      
      real, dimension(:,:), intent(in) :: datain
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real :: dt
      character(len=*), optional::  path 
      character(len=80), optional::  label, units
      character(len=80),dimension(2), optional :: xlabel, xunits
      real, dimension(2), optional :: xmin
      real, dimension(2), optional :: dx
      logical, optional :: docompress

      character(len=len(label)) ::  l_label
      character(len=len(units)) ::  l_units
      character(len=80),dimension(2) :: l_xlabel, l_xunits
      real, dimension(2) :: l_xmin
      real, dimension(2) :: l_dx
      logical :: l_docompress

      integer, dimension(2) :: ssize
      character(len=80) :: l_path, file
      character(len=80) :: name
      character(len=80),dimension(2) :: xname
      integer :: n
      real :: t
      character(20)  :: chindx
        

      call int2str(timestep,chindx,4)
     
      file=trim(dataname)//'_'//trim(chindx)//'.h5'
      name=trim(dataname)//'-'//trim(chindx)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
      if (.not.present(docompress)) then 
          l_docompress=.true.
      else
          l_docompress=docompress
      endif
      if (.not.present(label)) then 
          l_label=name
      else
          l_label=label
      endif
      if (.not.present(units)) then 
          l_units='arb.'
      else
          l_units=units
      endif
      n=timestep
      t=dt * n
        
      if (.not.present(xmin)) then 
          l_xmin=0.
      else
          l_xmin=xmin
      endif
      if (.not.present(dx)) then 
          l_dx=1.
      else
          l_dx=dx
      endif   
      xname(1)='x1'
      xname(2)='x2'
      if (.not.present(xlabel)) then 
          l_xlabel=xname
      else
          l_xlabel=xlabel
      endif
      if (.not.present(xunits)) then 
          l_xunits='arb.'
      else
          l_xunits=xunits
      endif  

      !write(*,*) 'Size  : ', size(datain)
      !write(*,*) 'Size1 : ', size(datain,1)
      !write(*,*) 'Size2 : ', size(datain,2)

      ssize(1) = size(datain,1)
      ssize(2) = size(datain,2)
      
      if (.not.present(path)) then
         l_path = ''
      else
         l_path = trim(path)
      endif
      
        call write_hdf_sds( l_path, file, &
      &                     datain, name, l_label, l_units,&
      &                     n, t, dt, ssize, l_xmin, l_dx, &
      &                     xname, l_xlabel, l_xunits, l_docompress)

      end subroutine write_hdf2d
     
! ----------------------------------------------------------------------      
      subroutine read_hdf2d(dataset,dataname,timestep,t,dimsizes,info_on&
     &ly,ierror,path)
! ---------------------------------------------------------------------- 
      implicit none
      
      real, dimension(:,:), intent(inout) :: dataset
      character(len=*), intent(in), optional :: path
      character(len=*), intent(in) :: dataname
      integer, intent(inout) :: timestep
      real, intent(inout) :: t  
      integer, dimension(:), intent(inout) :: dimsizes
      logical, intent(in) :: info_only
      integer, intent(out) :: ierror

      integer, dimension(2) :: ssize
      character(len=80) :: l_path, file
      character(len=80) :: name
      character(len=80),dimension(2) :: xname
      integer  :: ierr
      integer, dimension(2) :: idim_sizes
      character(20)  :: chindx
      logical :: filexisted

      !for MPI
      integer :: ierr2
      
      call int2str(timestep,chindx,4)
     
      file=trim(dataname)//'_'//trim(chindx)//'.h5'
      name=trim(dataname)//'-'//trim(chindx)
 
      ssize(1) = size(dataset,1)
      ssize(2) = size(dataset,2)

      if (.not.present(path)) then
         l_path = ''
      else
         l_path = trim(path)
      endif

      inquire(file = trim(l_path)//file, exist = filexisted) 
      if (filexisted) then 
         call read_hdf_sds(l_path, file, dataset,                       &
     &timestep, t, ierr,info_only,dim_sizes=idim_sizes)
         if (ierr.ne.0) then 
            print*, "data_utils: Error reading file ", file
            call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr2)
         endif
         ierror = ierr
         dimsizes=idim_sizes
      else
! file not existed
         print*, "data_utils: can't find ", file
         call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr2)
         ierror = -1
         dimsizes=0
      endif
      
      end subroutine read_hdf2d

! ----------------------------------------------------------------------      
      subroutine read_hdf3d(dataset,dataname,timestep,t,dimsizes,info_on&
     &ly,ierror,path)
! ---------------------------------------------------------------------- 
      implicit none
      
      real, dimension(:,:,:), intent(inout) :: dataset
      character(len=*), intent(in), optional :: path
      character(len=*), intent(in) :: dataname
      integer, intent(inout) :: timestep
      real, intent(inout) :: t  
      integer, dimension(:), intent(inout) :: dimsizes
      logical, intent(in) :: info_only
      integer, intent(out) :: ierror

      integer, dimension(3) :: ssize
      character(len=80) :: l_path, file
      character(len=80) :: name
      character(len=80),dimension(3) :: xname
      integer :: sdsIndex = 0
      integer  ::  ierr
      integer, dimension(3) :: idim_sizes
      character(20)  :: chindx
      logical :: filexisted

      !for MPI
      integer :: ierr2
      
      call int2str(timestep,chindx,4)
     
      file=trim(dataname)//'_'//trim(chindx)//'.h5'
      name=trim(dataname)//'-'//trim(chindx)
 
      ssize(1) = size(dataset,1)
      ssize(2) = size(dataset,2)
      ssize(3) = size(dataset,3)

      if (.not.present(path)) then
         l_path = ''
      else
         l_path = trim(path)
      endif
       
      inquire(file = trim(l_path)//file, exist = filexisted) 
      if (filexisted) then 
         call read_hdf_sds(l_path, file, dataset,                       &
     &timestep, t, ierr,info_only=info_only,dim_sizes=idim_sizes)
         if (ierr.ne.0) then 
            print*, "data_utils: Error reading file ", file
            call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr2)
         endif
         ierror = ierr
         dimsizes=idim_sizes
       else
! file not existed
         print*, "data_utils: can't find ", file
         call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr2)
         ierror = -1
         dimsizes=0
       endif
      
      end subroutine read_hdf3d

      end module m_data_utils
