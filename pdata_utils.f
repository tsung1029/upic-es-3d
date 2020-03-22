      module m_pdata_utils
      
      use  m_pdiagnostic_utilities
      use  m_string

      interface pwrite_data
        module procedure pwrite_hdf3d
        module procedure pwrite_hdf2d
      end interface
      
      interface pwrite_data_direct
        module procedure pwrite_hdf2d_direct
      end interface      

      interface pwrite_raw_data
        module procedure pwrite_raw_data_tag
        module procedure pwrite_raw_data
        module procedure pwrite_raws_data
        module procedure pwrite_plasma_raws_data
      end interface
      
      contains
! ----------------------------------------------------------------------      
      subroutine pwrite_plasma_raws_data(datain,timestep, tnps, tnpp,  &
                & dt, deltax, path, dataname, docompress,tagin)
! ----------------------------------------------------------------------      
     
      implicit none
      
      real, dimension(:,:,:), pointer, intent(in) :: datain
      integer, dimension(:,:,:), pointer, intent(in),optional :: tagin
      integer, intent(in) :: tnps,tnpp
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real :: dt,deltax
      logical :: docompress
      character(len=*), intent(in) ::  path

      character(20)  :: chindx
      logical :: l_docompress
      logical :: l_doappend

      character(len=80) :: file
      character(len=80) :: name
      integer :: n
      real :: t
        
      real, dimension(3) :: x0  
      real, dimension(3) :: xmin
      real, dimension(3) :: dx  

      x0 = (/0.0,0.0,0.0/)
      xmin = (/0.0,0.0,0.0/)
      dx = (/deltax,deltax,dt*timestep/)
      
      call int2str(timestep,chindx,6)

      file=trim(dataname)//'-'//trim(chindx)//'.h5'     
      name=trim(dataname)//'-'//trim(chindx)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
      l_docompress=docompress

      n=timestep
      t=dt * n     
 
                         
      call pwrite_hdf_raw(path, file, tnps,tnpp,x0, dx,&
     &               name, n, t, dt, xmin,                        &
     & datain,l_docompress,tagin)
      

     end subroutine pwrite_plasma_raws_data
! ----------------------------------------------------------------------      
      subroutine pwrite_raws_data(datain,timestep, tnps, tnpp, x0,  dt, &
                          &path, dataname, soffset, xmin, dx,&
                          &docompress)
! ----------------------------------------------------------------------      
     
      implicit none
      
      real, dimension(:,:,:), pointer, intent(in) :: datain
      integer, intent(in) :: tnps,tnpp
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real, dimension(:), intent(in):: x0  
      real :: dt, soffset
      real, dimension(:) :: xmin
      real, dimension(:) :: dx  
      logical :: docompress
      character(len=*), intent(in) ::  path

      character(20)  :: chindx
      logical :: l_docompress
      logical :: l_doappend

      character(len=80) :: file
      character(len=80) :: name
      integer :: n
      real :: t
        
      call int2str(timestep,chindx,6)
     
      file=trim(dataname)//'.h5'
      name=trim(dataname)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
      l_docompress=docompress

      n=timestep
      t=dt * n     
 
                         
      call pwrite_hdf_raw(path, file, tnps,tnpp, soffset,x0, dx,&
     &               name, n, t, dt, xmin,                        &
     & datain,l_docompress)
      

     end subroutine pwrite_raws_data

! ----------------------------------------------------------------------      
      subroutine pwrite_raw_data(datain,timestep, tnps, tnpp, x0,  dt, &
                          &path, dataname, soffset, xmin, dx, dspl,&
                          &docompress, doappend)
! ----------------------------------------------------------------------      
     
      implicit none
      
      real, dimension(:,:,:), pointer, intent(in) :: datain
      integer, intent(in) :: tnps,tnpp,dspl
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real, dimension(:), intent(in):: x0  
      real :: dt, soffset
      logical :: doappend
      real, dimension(:) :: xmin
      real, dimension(:) :: dx  
      logical :: docompress
      character(len=*), intent(in) ::  path

      character(20)  :: chindx
      logical :: l_docompress
      logical :: l_doappend

      character(len=80) :: file
      character(len=80) :: name
      integer :: n
      real :: t
        
      call int2str(timestep,chindx,6)
     
      file=trim(dataname)//'.h5'
      name=trim(dataname)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
      l_docompress=docompress

      n=timestep
      t=dt * n     
 
      l_doappend = doappend          
                         
      call pwrite_hdf_raw(path, file, tnps,tnpp, soffset,x0, dx,&
     &               name, n, t, dt, xmin, dspl,                &
     & datain,l_doappend, l_docompress)
      

     end subroutine pwrite_raw_data
! ----------------------------------------------------------------------      
      subroutine pwrite_raw_data_tag(datain,timestep, tnps, tnpp, x0,  dt, &
                          &path, dataname, soffset, xmin, dx, dspl,&
                          &docompress, doappend, tag)
! ----------------------------------------------------------------------      
     
      implicit none
      
      real, dimension(:,:,:), pointer, intent(in) :: datain
      integer, dimension(:,:,:), pointer, intent(in) :: tag
      integer, intent(in) :: tnps,tnpp,dspl
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real, dimension(:), intent(in):: x0  
      real :: dt, soffset
      logical :: doappend
      real, dimension(:) :: xmin
      real, dimension(:) :: dx  
      logical :: docompress
      character(len=*), intent(in) ::  path

      character(20)  :: chindx
      logical :: l_docompress
      logical :: l_doappend

      character(len=80) :: file
      character(len=80) :: name
      integer :: n
      real :: t
        
      call int2str(timestep,chindx,6)
     
      file=trim(dataname)//'.h5'
      name=trim(dataname)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
      l_docompress=docompress

      n=timestep
      t=dt * n     
 
      l_doappend = doappend          
                         
      call pwrite_hdf_raw(path, file, tnps,tnpp, soffset,x0, dx,&
     &               name, n, t, dt, xmin, dspl,                &
     & datain,l_doappend, l_docompress, tag)
      

     end subroutine pwrite_raw_data_tag
! ----------------------------------------------------------------------      
      subroutine pwrite_hdf3d(datain, dataname, timestep, dt, &
                          & nvpy, nvpz, in_stride, label, units, &
                          & xlabel, xunits, xmin, path)
! ----------------------------------------------------------------------      

      implicit none
      
      real, dimension(:,:,:) :: datain  
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep, nvpy, nvpz
      integer, dimension(:), intent(in) :: in_stride
      real :: dt
      character(len=*), intent(in), optional ::  label, units
      character(len=*),dimension(:), intent(in), optional:: xlabel, xunits
      real, dimension(:), optional :: xmin
      character(len=*), intent(in), optional ::  path

      character(20)  :: chindx
      character(len=80) ::  l_label
      character(len=80) :: l_units
      character(len=80),dimension(3) :: l_xlabel
      character(len=80),dimension(3) :: l_xunits
      real, dimension(3) :: l_xmin

      integer, dimension(3) :: ssize, stride
      character(len=80) :: l_path, file
      character(len=80) :: name
      character(len=80),dimension(3) :: xname
      integer :: n
      real :: t
        
      call int2str(timestep,chindx,6)
     
      file=trim(dataname)//'-'//trim(chindx)//'.h5'
      name=trim(dataname)//'-'//trim(chindx)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
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
      
      
        call pwrite_hdf_sds(l_path, file, &
      &                     datain, name, l_label, l_units,&
      &                     n, t, dt, l_xmin, xname, &
      &                     l_xlabel, l_xunits, nvpy, nvpz, in_stride)

     end subroutine pwrite_hdf3d
     
 
! ----------------------------------------------------------------------      
      subroutine pwrite_hdf2d(datain, division, dataname, timestep,dt,&
                          &label,units,xlabel,xunits,xmin,dx,&
                          &path)
! ----------------------------------------------------------------------      

      implicit none
      
      real, dimension(:,:), intent(in) :: datain
      integer, intent(in) :: division      
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real :: dt
      character(len=*), optional::  label, units
      character(len=*),dimension(:), optional :: xlabel, xunits
      real, dimension(:), optional :: xmin
      real, dimension(:), optional :: dx
      character(len=*), intent(in), optional::  path
      
      character(len=80) ::  l_label
      character(len=80) ::  l_units
      character(len=80),dimension(2) :: l_xlabel
      character(len=80),dimension(2) :: l_xunits
      real, dimension(2) :: l_xmin
      real, dimension(2) :: l_dx

      integer, dimension(2) :: ssize
      character(len=80) :: l_path, file
      character(len=80) :: name
      character(len=80),dimension(2) :: xname
      integer :: n
      real :: t
      character(20)  :: chindx
        
      call int2str(timestep,chindx,6)
     
      file=trim(dataname)//'_'//trim(chindx)//'.h5'
      name=trim(dataname)//'-'//trim(chindx)
      
      !write(*,*) 'Filename :',file
      !write(*,*) 'Data Name:',name
 
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

!      write(*,*) 'Size  : ', size(datain)
!      write(*,*) 'Size1 : ', size(datain,1)
!      write(*,*) 'Size2 : ', size(datain,2)

      ssize(1) = size(datain,1)
      ssize(2) = size(datain,2)
      
      if (.not.present(path)) then
         l_path = ''
      else
         l_path = trim(path)
      endif

        call pwrite_hdf_sds(l_path, file, &
      &                     datain, division, name, l_label, l_units,&
      &                     n, t, dt, l_xmin, l_dx, &
      &                     xname, l_xlabel, l_xunits)

     end subroutine pwrite_hdf2d
    

! ----------------------------------------------------------------------      
      subroutine pwrite_hdf2d_direct(datain, idsend, dataname, timestep,&
                          &dt,doappend, nappend, append_pos, label,units&
                          &,xlabel,xunits,xmin,dx,                      &
                          &docompress,path)
! ----------------------------------------------------------------------      

      implicit none
      
      real, dimension(:,:), intent(in) :: datain
      integer, intent(in) :: idsend      
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: timestep
      real :: dt
      logical, optional :: doappend
      integer, optional :: nappend, append_pos
      character(len=*), optional::  label, units
      character(len=*),dimension(:), optional :: xlabel, xunits
      real, dimension(:), optional :: xmin
      real, dimension(:), optional :: dx
      logical, optional :: docompress
      character(len=*), optional, intent(in) ::  path

      character(len=80) ::  l_label
      character(len=80) ::  l_units
      character(len=80),dimension(2) :: l_xlabel
      character(len=80),dimension(2) :: l_xunits
      real, dimension(2) :: l_xmin
      real, dimension(2) :: l_dx
      logical :: l_docompress
      logical :: l_doappend
      integer :: l_nappend, l_append_pos

      integer, dimension(2) :: ssize
      character(len=80) :: l_path, file
      character(len=80) :: name
      character(len=80),dimension(2) :: xname
      integer :: n
      real :: t
      character(20)  :: chindx
        
     
      call int2str(timestep,chindx,6)

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

      if (.not.present(doappend)) then
         l_doappend = .false.
         l_nappend = 1
         l_append_pos = 0 
      else
         l_doappend = doappend 
         l_nappend = nappend
         l_append_pos = append_pos
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
     
        call pwrite_hdf_sds_direct(l_path, file, &
     &                     datain, idsend, name, l_label, l_units,&
     &                     n, t, dt, l_xmin, l_dx, &
     &                     xname, l_xlabel, l_xunits, l_doappend,   &
     &                     l_nappend, l_append_pos, l_docompress)

      end subroutine pwrite_hdf2d_direct     
     
      end module m_pdata_utils
