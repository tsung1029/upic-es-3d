      module hdf_write32_jf
      
      use globals, only: LINEAR,QUADRATIC
      use pinit32d_jf
      implicit none
      private
      public :: writef, write_jf
      public :: EX,EY,JX,JY,DEN,VDOTEX,VDOTEY,PHSL_XX,PHSL_XY,PHSL_YX,PHSL_YY
      public :: ESPOYNT_X,ESPOYNT_Y,FIN_VS_INIT_ENE,FIN_VS_INIT_VX,FIN_VS_INIT_VY
      public :: VDOTEX_PART,VDOTEY_PART,VDOTEX_INT,VDOTEY_INT,FVX_XY_LABEL,FVY_XY_LABEL
      public :: VDOTEX_FOLLOW_PART,VDOTEY_FOLLOW_PART,ESPOYNT_INT_X,ESPOYNT_INT_Y
      public :: EX_ENE_INT,EY_ENE_INT,DIV_ESPOYNTINT,VXVY,PH_VX_X,PH_VY_X,PH_VX_Y,PH_VY_Y
      public :: DIVESPOYNT,BIN_EX_LABEL,BIN_EY_LABEL,POT,EZ
		
		save
      integer,parameter :: EX = 1,EY = 2,JX =3,JY = 4,DEN = 5,VDOTEX = 6,VDOTEY = 7
      integer,parameter :: PHSL_XX = 8,PHSL_XY=9,PHSL_YX=10,PHSL_YY=11
      integer,parameter :: ESPOYNT_X=12,ESPOYNT_Y=13,FIN_VS_INIT_ENE=14,FIN_VS_INIT_VX=15
      integer,parameter :: FIN_VS_INIT_VY=16,VDOTEX_PART=17,VDOTEY_PART=18
      integer,parameter :: VDOTEX_INT=19,VDOTEY_INT=20,FVX_XY_LABEL=21,FVY_XY_LABEL=22
      integer,parameter :: VDOTEX_FOLLOW_PART = 23,VDOTEY_FOLLOW_PART=24
      integer,parameter :: ESPOYNT_INT_X=25,ESPOYNT_INT_Y=26,EX_ENE_INT=27,EY_ENE_INT=28
      integer,parameter :: DIV_ESPOYNTINT=29,VXVY=30,PH_VX_X=31,PH_VY_X=32,PH_VX_Y=33,PH_VY_Y=34
      integer,parameter :: DIVESPOYNT=35,BIN_EX_LABEL=36,BIN_EY_LABEL=37,POT=38,EZ=39

		interface writef
			 module procedure ipwrite2_hdf
			 module procedure ipwrite3_hdf
!         module procedure ipwrite2
		end interface

		interface write_jf
			module procedure write_jf2d_hdf
			module procedure write_jf3d_hdf
		end interface
      
      contains
!
! FST
! modified by JF to write a single array without collecting from other nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine write_jf2d_hdf(f,nx,ny,iunit,nrec,fname,label_code,time,it,dvx,dvy,xoff,yoff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine writes real data f to a file
! f = input data to be written, modified on node 0
! nx/ny = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! iunit = fortran unit number
! nrec = current record number for write (if negative, open file with
! recl=-nren)
! name = file name (used only if nren < 0)
! dvx and dvy are for changing the space between points on the x and y axis, ie, 
!		dvx = xaxis(2)-xaxis(1)
! xoff and yoff are the positions of the xaxis(1) and yaxis(1)
! input: f, nx, ny, iunit, nrec, fname
! output: nrec
      implicit none
      integer nx,ny, iunit, nrec,it
      real f
      character*(*) fname
      dimension f(nx,ny,1)
      real time, dvx,dvy,xoff,yoff
      integer label_code
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k
      integer ierr
      dimension istatus(lstat)
      character(len=5) :: nlabel
      character(len=100) :: name
! HDF stuff
!
!     Function declaration.
!
      integer sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      integer sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      integer sfcreate
      external sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      external sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      external sfcreate
      integer DFACC_CREATE,DFNT_INT32,DFNT_FLOAT32
      integer DFNT_FLOAT64
      parameter(DFACC_CREATE = 4, DFNT_INT32 = 24)
      parameter(DFNT_FLOAT32 = 5)
      parameter(DFNT_FLOAT64 = 6)
      integer SD_UNLIMITED
      parameter(SD_UNLIMITED = 0)

! 
!**** Variable declaration *******************************************
!

! file handles for HDF
      integer sd_id, sds_id, sds_index, status
      integer dimid1,dimid2,dimid3
! array informations needed for HDF
      integer start(2), edges(2), stride(2), dims(2)
      integer ledge(2)
      
      character*40 title
      character*40 xtitle,xunits
      character*40 ytitle,yunits
      character*40 ztitle,zunits
      
! arrays needed for HDF-write      
      real*4 xaxis(nx),yaxis(ny)
      real*4 ftemp(nx,ny)
      integer ix,iy
!      
! temporary attributes, needed for Ricardo's IDL routines
!
      real*8 time_double
      
! DEBUG
!     write(*,*)'in pwrite_hdf'
!     write(*,*)'(pwrite_hdf) nproc= ',nproc
!     write(*,*)'(pwrite_hdf) nx = ',nx
!     write(*,*)'(pwrite_hdf) ny = ',ny
!     write(*,*)'(pwrite_hdf)'
! DEBUG
! define HDF "chunk" here
      start(1)=0
      start(2)=0
      stride(1)=1
      stride(2)=1
      ledge(1)=nx
      ledge(2)=ny
! for another implemenation, not used
!      gedge(1)=nx
!      gedge(2)=ny*nproc
!
      dims(1)=(ledge(1)-start(1))/stride(1)
      dims(2)=(ledge(2)-start(2))/stride(2)
      
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
			call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
			write (nlabel,'(i5)') it + 10000
			name = trim(fname)//'_'//trim(adjustl(nlabel))//'.hdf'

      if (idproc.eq.0) then
         sd_id = sfstart(trim(name), DFACC_CREATE)
! new hdf write
         sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT64, 2, dims)
! HDF
         status = sfwdata(sds_id,start,stride,ledge,f)
         nrec = nrec + 1
         ! now write HDF header needed for Ricardo's package
         time_double=time
         status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
         status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,it)

         do ix=1,nx
             xaxis(ix)=(ix-1)*dvx + xoff
         enddo
         dimid1 = sfdimid(sds_id,0)
         status = sfsdmname(dimid1, 'x1')
         status = sfsdscale(dimid1, nx, DFNT_FLOAT64, xaxis)
	 			 status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         
         do iy=1,(ny)
             yaxis(iy)=(iy-1)*dvy + yoff
         enddo
         
         if ((label_code .eq. PHSL_XX) .OR. (label_code .eq. PHSL_XY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvxmax / (real(ny)/2.)
				 endif
         if ((label_code .eq. PHSL_YX) .OR. (label_code .eq. PHSL_YY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvymax / (real(ny)/2.)
         endif
         
         dimid2 = sfdimid(sds_id,1)
         status = sfsdmname(dimid2, 'x2')
         status = sfsdscale(dimid2, (ny), DFNT_FLOAT64, yaxis)
	 			 status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')
	 		
				 status = sfsdtstr(sds_id,trim(ztitle),trim(zunits),'F5.2','cartesian')
	 		
! HDF equivalent of the "close" statement
         status=sfendacc(sds_id)
         status=sfend(sd_id)
!

      endif
      return
      end subroutine write_jf2d_hdf

!
! FST
! modified by JF to write a single array without collecting from other nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine write_jf3d_hdf(f,nx,ny,nz,iunit,nrec,fname,label_code,time,it,&
      	&dvx,dvy,dvz,xoff,yoff,zoff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine writes real data f to a file
! f = input data to be written, modified on node 0
! nx/ny = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! iunit = fortran unit number
! nrec = current record number for write (if negative, open file with
! recl=-nren)
! name = file name (used only if nren < 0)
! dvx and dvy are for changing the space between points on the x and y axis, ie, 
!		dvx = xaxis(2)-xaxis(1)
! xoff and yoff are the positions of the xaxis(1) and yaxis(1)
! input: f, nx, ny, iunit, nrec, fname
! output: nrec
      implicit none
      integer nx,ny,nz,iunit, nrec,it
      real f
      character*(*) fname
      dimension f(nx,ny,nz)
      real time, dvx,dvy,dvz,xoff,yoff,zoff
      integer label_code
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k
      integer ierr
      dimension istatus(lstat)
      character(len=5) :: nlabel
      character(len=100) :: name
! HDF stuff
!
!     Function declaration.
!
      integer sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      integer sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      integer sfcreate
      external sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      external sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      external sfcreate
      integer DFACC_CREATE,DFNT_INT32,DFNT_FLOAT32
      integer DFNT_FLOAT64
      parameter(DFACC_CREATE = 4, DFNT_INT32 = 24)
      parameter(DFNT_FLOAT32 = 5)
      parameter(DFNT_FLOAT64 = 6)
      integer SD_UNLIMITED
      parameter(SD_UNLIMITED = 0)

! 
!**** Variable declaration *******************************************
!

! file handles for HDF
      integer sd_id, sds_id, sds_index, status
      integer dimid1,dimid2,dimid3
! array informations needed for HDF
      integer start(3), edges(3), stride(3), dims(3)
      integer ledge(3)
      
      character*40 title
      character*40 xtitle,xunits
      character*40 ytitle,yunits
      character*40 ztitle,zunits
      
! arrays needed for HDF-write      
      real*4 xaxis(nx),yaxis(ny),zaxis(nz)
      real*4 ftemp(nx,ny)
      integer ix,iy,iz
!      
! temporary attributes, needed for Ricardo's IDL routines
!
      real*8 time_double
      
! DEBUG
!     write(*,*)'in pwrite_hdf'
!     write(*,*)'(pwrite_hdf) nproc= ',nproc
!     write(*,*)'(pwrite_hdf) nx = ',nx
!     write(*,*)'(pwrite_hdf) ny = ',ny
!     write(*,*)'(pwrite_hdf)'
! DEBUG
! define HDF "chunk" here
      start(1)=0
      start(2)=0
      start(3)=0
      stride(1)=1
      stride(2)=1
      stride(3)=1
      ledge(1)=nx
      ledge(2)=ny
      ledge(3)=nz
! for another implemenation, not used
!      gedge(1)=nx
!      gedge(2)=ny*nproc
!
      dims(1)=(ledge(1)-start(1))/stride(1)
      dims(2)=(ledge(2)-start(2))/stride(2)
      dims(3)=(ledge(3)-start(3))/stride(3)
      
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
			call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
			write (nlabel,'(i5)') it + 10000
			name = trim(fname)//'_'//trim(adjustl(nlabel))//'.hdf'

      if (idproc.eq.0) then
         sd_id = sfstart(trim(name), DFACC_CREATE)
! new hdf write
!         sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT32, 3, dims)
         sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT64, 3, dims)
! HDF
         status = sfwdata(sds_id,start,stride,ledge,f)
         nrec = nrec + 1
         ! now write HDF header needed for Ricardo's package
         time_double=time
         status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
         status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,it)

         do ix=1,nx
             xaxis(ix)=(ix-1)*dvx + xoff
         enddo
         dimid1 = sfdimid(sds_id,0)
         status = sfsdmname(dimid1, 'x1')
         status = sfsdscale(dimid1, nx, DFNT_FLOAT64, xaxis)
	 			 status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         
         do iy=1,(ny)
             yaxis(iy)=(iy-1)*dvy + yoff
         enddo
         do iz=1,(nz)
             zaxis(iz)=(iz-1)*dvz + zoff
         enddo
         
         if ((label_code .eq. PHSL_XX) .OR. (label_code .eq. PHSL_XY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvxmax / (real(ny)/2.)
				 endif
         if ((label_code .eq. PHSL_YX) .OR. (label_code .eq. PHSL_YY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvymax / (real(ny)/2.)
         endif
         
         dimid2 = sfdimid(sds_id,1)
         status = sfsdmname(dimid2, 'x2')
         status = sfsdscale(dimid2, (ny), DFNT_FLOAT64, yaxis)
	 			 status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')

         dimid3 = sfdimid(sds_id,2)
         status = sfsdmname(dimid3, 'x3')
         status = sfsdscale(dimid3, (nz), DFNT_FLOAT64, zaxis)
	 			 status = sfsdmstr(dimid3, trim(ztitle),trim(zunits),'F5.2')
	 		
!				 status = sfsdtstr(sds_id,trim(ztitle),trim(zunits),'F5.2','cartesian')
	 		
! HDF equivalent of the "close" statement
         status=sfendacc(sds_id)
         status=sfend(sd_id)
!

      endif
      return
      end subroutine write_jf3d_hdf


! FST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine PWRITE2_HDF(f,nx,nxv,kyp,iunit,nrec,name,label_code,iter,time,iorder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine collects distributed real data f and writes to a file
! f = input data to be written, modified on node 0
! nx/kyp = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! iunit = fortran unit number
! nrec = current record number for write (if negative, open file with
! recl=-nren)
! name = file name (used only if nren < 0)
! input: f, nx, kyp, nxv, iunit, nrec, fname
! output: nrec
      implicit none
      integer nx,kyp, nxv, iunit, nrec,iter
      real f
      character*(*) name
      dimension f(nxv,kyp,1)
      real time
      integer iorder, label_code
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k,l,temp
      integer ierr,ny
      integer per
      dimension per(2)
      dimension istatus(lstat)
! HDF stuff
!
!     Function declaration.
!
      integer sfstart, sfselect, sfwdata, sfendacc, sfend
      integer sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt, sfsdtstr
      integer sfcreate
      external sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      external sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      external sfcreate
      integer DFACC_CREATE,DFNT_INT32,DFNT_FLOAT32
      integer DFNT_FLOAT64,DFNT_CHAR
      parameter(DFACC_CREATE = 4, DFNT_INT32 = 24)
      parameter(DFNT_FLOAT32 = 5)
      parameter(DFNT_FLOAT64 = 6)
      integer SD_UNLIMITED
      parameter(SD_UNLIMITED = 0)
      parameter(DFNT_CHAR=4)
      

! 
!**** Variable declaration *******************************************
!

! file handles for HDF
      integer sd_id, sds_id, sds_index, status
      integer dimid1,dimid2,dimid3
! array informations needed for HDF
      integer start(2), edges(2), stride(2), dims(2)
      integer ledge(2)
      
      character*70 title
      character*70 xtitle,xunits
      character*70 ytitle,yunits
      character*70 ztitle,zunits
      
! arrays needed for HDF-write      
!      real*4 xaxis(nx),yaxis((kyp-3)*nproc)
!      real*4 ftemp(nx,kyp-3)
!      real*4 ftemprcv(nxv,kyp)
			real,dimension(:), allocatable :: xaxis,yaxis
			real,dimension(:,:), allocatable :: ftemp,ftemprcv
			
      integer ix,iy
!      
! temporary attributes, needed for Ricardo's IDL routines
!
      real*8 time_double
      
      allocate(xaxis(nx),yaxis((kyp-3)*nproc))
      allocate(ftemp(nx,kyp-3))
      allocate(ftemprcv(nxv,kyp))
      

! DEBUG
!     write(*,*)'in pwrite_hdf'
!     write(*,*)'(pwrite_hdf) nproc= ',nproc
!     write(*,*)'(pwrite_hdf) nx = ',nx
!     write(*,*)'(pwrite_hdf) nxv = ',nxv
!     write(*,*)'(pwrite_hdf) kyp = ',kyp
!     write(*,*)'(pwrite_hdf)'
! DEBUG
			ny = kyp-3
! define HDF "chunk" here
      start(1)=0
      start(2)=0
      stride(1)=1
      stride(2)=1
      ledge(1)=nx
      ledge(2)=ny
! for another implemenation, not used
!      gedge(1)=nx
!      gedge(2)=kyp*nproc
!
      dims(1)=(ledge(1)-start(1))/stride(1)
      dims(2)=((ledge(2)-start(2))/stride(2))*nproc
!
! this segment is used for shared memory computers
!     if (nrec.lt.0) then
!        lrec = -nrec
!        open(unit=iunit,file=name,form='unformatted',access='direct',re
!    1cl=lrec,status='unknown')
!        nrec = 1
!     endif
!     write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,kyp)
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)

		call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)

! node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            lrec = -nrec
!            open(unit=iunit,file=name,form='unformatted',access='direct'
!     1,recl=lrec,status='unknown')
! HDF
! the following lines replaces the generic FORTRAN open
             sd_id = sfstart(name, DFACC_CREATE)
! add an attribute called "nproc", which tells me how many processors there are.
!             status = sfsnatt(sd_id,"nproc",DFNT_INT32,1,nproc)
            nrec = 1
! 
! new hdf write
            sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT64, 2, dims)

! HDF
! DEBUG
!	    write(*,*)'sds_id = ',sds_id
! DEBUG
         endif
! no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
! special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
            id = 1
            call MPI_RECV(ftemprcv,nx*kyp,mreal,id,99,lworld,istatus,ierr)
         endif
! first write data for node 0
         nrec0 = nrec
         start(2)=(nrec-1)*(ny)
! if no special node present, then copy data to ftemp instead of MPI_RECV
         if(nvp.eq.nproc) then
            if (iorder .eq. QUADRATIC) then
            do ix=1,nx
                do iy=1,ny
                    ftemp(ix,iy)=f(ix+1,iy+1,1)
                enddo
            enddo            
            else
            do ix=1,nx
                do iy=1,kyp
                    ftemp(ix,iy)=f(ix,iy,1)
                enddo
            enddo                     
            endif 
         endif       
! OLD Binary write
!        write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,kyp)
! HDF replacement
         status = sfwdata(sds_id,start,stride,ledge,ftemp)
!        write(*,*)'status (sfwdata_) = ', status
         nrec = nrec + 1
! then write data from remaining nodes
         do i = 2, np
            id = i - ioff
! calculate the starting point for             
            start(2)=(nrec-1)*(ny)
            call MPI_RECV(ftemprcv,nxv*kyp,mreal,id,99,lworld,istatus,ierr)
            do ix=1,nx
            	do iy=1,ny
            		ftemp(ix,iy) = ftemprcv(ix+1,iy+1)
            	enddo
            enddo
! old binary write            
!           write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,kyp)
! 
! new hdf write
            status = sfwdata(sds_id,start,stride,ledge,ftemp)
!           write(*,*)'status (sfwdata_) = ', status
            nrec = nrec + 1
         enddo
         per(1) = 1
         per(2) = nvp
         
         
         ! now write HDF header needed for Ricardo's package
         time_double=time
         status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
         status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,iter)
! read data back for node 0
!         read (unit=iunit,rec=nrec0) ((f(j,k),j=1,nx),k=1,kyp)

         do ix=1,nx
             xaxis(ix)=(ix-1)
         enddo
         dimid1 = sfdimid(sds_id,0)
! 	     write(*,*)'in dim 1',dimid1
         status = sfsdmname(dimid1, 'x1')
!         status = sfsdscale(dimid1, nx, DFNT_FLOAT32, xaxis)
         status = sfsdscale(dimid1, nx, DFNT_FLOAT64, xaxis)
	 		status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         
         do iy=1,(ny)*nproc
             yaxis(iy)=(iy-1)
         enddo
         dimid2 = sfdimid(sds_id,1)
         status = sfsdmname(dimid2, 'x2')
         status = sfsdscale(dimid2, (ny)*nproc, DFNT_FLOAT64, yaxis)
!         status = sfsdscale(dimid2, (ny)*nproc, DFNT_FLOAT32, yaxis)
	 		status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')
	 		
			status = sfsdtstr(sds_id,trim(ztitle),trim(zunits),'F5.2','cartesian')
			
! HDF equivalent of the "close" statement
         status=sfendacc(sds_id)
         status=sfend(sd_id)
!

! other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
!            if (iorder .eq. QUADRATIC) then
!            do ix=1,nx
!                do iy=1,kyp
!                    ftemp(ix,iy)=f(ix+1,iy+1,1)
!                enddo
!            enddo            
!            else
!            do ix=1,nx
!                do iy=1,kyp
!                    ftemp(ix,iy)=f(ix,iy,1)
!                end do
!            enddo                     
!            endif
         call MPI_SEND(f,nxv*kyp,mreal,0,99,lworld,ierr)
      endif
      
      deallocate(xaxis,yaxis)
      deallocate(ftemp)
      deallocate(ftemprcv)

      
      return
      end subroutine pwrite2_hdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
	subroutine PWRITE3_HDF(f,nx,nxv,kyp,kzp,nvpy,nvpz,name,label_code,&
		& iter,time,iorder,in_stride)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! this subroutine collects distributed real data f and writes to a file
	! f = input data to be written, modified on node 0
	! nx/kyp = length of data f in x/y on each processor to write
	! nxv = first dimension of data array f, must be >= nx
	! name = file name (used only if nren < 0)
	! input: f, nx, kyp, nxv, fname
	implicit none
	integer nx,kyp,kzp, nxv,iter,nvpy,nvpz,in_stride
	real f
	character*(*) name
	dimension f(nxv,kyp,kzp,1)
	dimension in_stride(3)
	real time
	integer iorder, label_code
	! common block for parallel processing
	integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
	! lstat = length of status array
	parameter(lstat=10)
	! nproc = number of real or virtual processors obtained
	! lgrp = current communicator
	! lworld = MPI_COMM_WORLD communicator
	common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
	! local data
	integer istatus, lrec, nvp, idproc, np, ioff, id, i, j, k,l,temp
	integer ierr,ny,nz
	integer per
	dimension per(2)
	dimension istatus(lstat)
	! HDF stuff
	integer sfstart, sfselect, sfwdata, sfendacc, sfend
	integer sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt, sfsdtstr
	integer sfcreate
	external sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
	external sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
	external sfcreate
	integer DFACC_CREATE,DFNT_INT32,DFNT_FLOAT32
	integer DFNT_FLOAT64,DFNT_CHAR
	parameter(DFACC_CREATE = 4, DFNT_INT32 = 24)
	parameter(DFNT_FLOAT32 = 5)
	parameter(DFNT_FLOAT64 = 6)
	integer SD_UNLIMITED
	parameter(SD_UNLIMITED = 0)
	parameter(DFNT_CHAR=4)

	!**** Variable declaration *******************************************
	! file handles for HDF
	integer sd_id, sds_id, sds_index, status
	integer dimid1,dimid2,dimid3
	! array informations needed for HDF
	integer start(3), edges(3), stride(3), dims(3)
	integer ledge(3)
	
	character*70 title
	character*70 xtitle,xunits
	character*70 ytitle,yunits
	character*70 ztitle,zunits
	
	! arrays needed for HDF-write      
	real,dimension(:), allocatable :: xaxis,yaxis, zaxis
	real,dimension(:,:,:), allocatable :: ftemp,ftemprcv
	
	integer ix,iy,iz
	real*8 time_double

	stride = 1
	ny = kyp-3; nz = kzp-3
	allocate(xaxis(nx/in_stride(1)),yaxis((kyp-3)*nvpy/in_stride(2)),zaxis((kzp-3)*nvpz/in_stride(3)))
!	allocate(ftemp(nx,kyp-3,nz))
	allocate(ftemp(nx/in_stride(1),ny/in_stride(2),nz/in_stride(3)))
	allocate(ftemprcv(nxv,kyp,kzp))
	! define HDF "chunk" here
	start = 0
	ledge(1) = nx/in_stride(1); ledge(2) = ny/in_stride(2); ledge(3) = nz/in_stride(3);
!	dims(1)=nx
!	dims(2)=ny*nvpy/stride(2)
!	dims(3)=nz*nvpz/stride(3)
!	dims(1)=(ledge(1)-start(1))
!	dims(2)=((ledge(2)-start(2)))*nvpy
!	dims(3)=((ledge(3)-start(3)))*nvpz
	dims(1)=(ledge(1)-start(1))
	dims(2)=((ledge(2)-start(2)))*nvpy
	dims(3)=((ledge(3)-start(3)))*nvpz
	
	! this segment is used for mpi computers
	! determine the rank of the calling process in the communicator
	call MPI_COMM_RANK(lworld,idproc,ierr)
	! determine the size of the group associated with a communicator
	call MPI_COMM_SIZE(lworld,nvp,ierr)
	
	call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
	
	! node 0 receives messages from other nodes
	if (idproc.eq.0) then
	! HDF
	! the following lines replaces the generic FORTRAN open
	sd_id = sfstart(name, DFACC_CREATE)
	! add an attribute called "nproc", which tells me how many processors there are.
	!             status = sfsnatt(sd_id,"nproc",DFNT_INT32,1,nproc)
	! new hdf write
	sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT64, 3, dims)

	np = nvpy*nvpz
		! first write data for node 0
	! if no special node present, then copy data to ftemp instead of MPI_RECV
	if (iorder .eq. QUADRATIC) then
	 do ix=1,nx/in_stride(1)
		 do iy=1,ny/in_stride(2)
			 do iz=1,nz/in_stride(3)
				 ftemp(ix,iy,iz)=f((ix-1)*in_stride(1)+3,(iy-1)*in_stride(2)+3,(iz-1)*in_stride(3)+3,1)
			 if (ix ==5 .and. iy == 5) then
			 endif
			 enddo
		 enddo
	 enddo            
	else
	 do ix=1,nx
		 do iy=1,kyp
			 do iz=1,kzp
				 ftemp(ix,iy,iz)=f(ix,iy,iz,1)
			 enddo
		 enddo
	 enddo                     
	endif 

	status = sfwdata(sds_id,start,stride,ledge,ftemp)
	!        write(*,*)'status (sfwdata_) = ', status
	! then write data from remaining nodes
	do i = 0, nvpz-1
		do j = 0, nvpy-1
			if (i + j .ne. 0) then	!Don't do the first one, cause we did it above
				id = i*nvpy + j
				! calculate the starting point for             
				start(2) = j * ny / in_stride(2)
				start(3) = i * nz / in_stride(3)
				call MPI_RECV(ftemprcv,nxv*kyp*kzp,mreal,id,99,lworld,istatus,ierr)
				do ix=1,nx/in_stride(1)
					do iy=1,ny/in_stride(2)
						do iz=1,nz/in_stride(3)
							ftemp(ix,iy,iz)=ftemprcv((ix-1)*in_stride(1)+3,(iy-1)*in_stride(2)+3,(iz-1)*in_stride(3)+3)
						enddo
					enddo
				enddo
				! new hdf write
				status = sfwdata(sds_id,start,stride,ledge,ftemp)
				!           write(*,*)'status (sfwdata_) = ', status
			endif
		enddo
	enddo

	per(1) = 1
	per(2) = nvp
		
		! now write HDF header needed for Ricardo's package
		time_double=time
		status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
		status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,iter)
	! read data back for node 0
	!         read (unit=iunit,rec=nrec0) ((f(j,k),j=1,nx),k=1,kyp)
	
		do ix=1,nx/in_stride(1)
			 xaxis(ix)=(ix-1)*in_stride(1)
		enddo
		dimid1 = sfdimid(sds_id,0)
		! 	     write(*,*)'in dim 1',dimid1
		status = sfsdmname(dimid1, 'x1')
		!         status = sfsdscale(dimid1, nx, DFNT_FLOAT32, xaxis)
		status = sfsdscale(dimid1, nx/in_stride(1), DFNT_FLOAT64, xaxis)
		status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         

		do iy=1,(ny)*nvpy/in_stride(2)
			yaxis(iy)=(iy-1)*in_stride(2)
		enddo
		dimid2 = sfdimid(sds_id,1)
		status = sfsdmname(dimid2, 'x2')
		status = sfsdscale(dimid2, (ny)*nvpy/in_stride(2), DFNT_FLOAT64, yaxis)
		! status = sfsdscale(dimid2, (ny)*nproc, DFNT_FLOAT32, yaxis)
		status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')
		
		do iz=1,nz*nvpz/in_stride(3)
			zaxis(iz)= (iz - 1)*in_stride(3)
		enddo
		dimid3 = sfdimid(sds_id,2)
		status = sfsdmname(dimid3, 'x3')
		status = sfsdscale(dimid3, nz*nvpz/in_stride(3), DFNT_FLOAT64, zaxis)
		status = sfsdmstr(dimid3, trim(ztitle),trim(zunits),'F5.2')
		
		! HDF equivalent of the "close" statement
		status=sfendacc(sds_id)
		status=sfend(sd_id)

	! other nodes send data to node 0
		else
			call MPI_SEND(f,nxv*kyp*kzp,mreal,0,99,lworld,ierr)
		endif
	
		deallocate(xaxis,yaxis,zaxis)
		deallocate(ftemp)
		deallocate(ftemprcv)
		
		
		return
	end subroutine pwrite3_hdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! HDF
! this is the wrapper function for pwrite2_hdf
			subroutine ipwrite2_hdf(f,nxv,nypmx,iunit,nrec,it,time,label_code,name,inorder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine collects distributed real 2d data f 
! and it writes to a file.  Only some of the data can be written.
! this = ufield2d descriptor of data
! f = input data to be written
! iunit = fortran unit number
! nrec = current record number for write for old file
!        or record length for new file
! name = file name (used only if new file is to be opened)
! it = integer to append to file name
! inorder = interpolation order, 1 = LINEAR, 2 = QUADRATIC
			implicit none
			!         type (ufield2d), intent(in) :: this
			real, dimension(:,:,:) :: f
			integer :: iunit, nrec, nxv, nypmx
			character(len=*), optional :: name
			integer :: it, label_code
			integer, optional :: inorder
			real :: time
			! local data
			integer :: mx, lrec, order,nblok,i,temp
			character(len=5) :: nlabel
			character(len=60) :: fname
			character(len=10), save :: sname = ':ipwrite2:'
! check for errors
!         if (monitor > 0) then
!         if (monitor==2) call werrfl(class//sname//' started')
!         if (checkinit(sname,this) /= 0) return
!         if (this%layout /= XLOCAL) then
!            erstr = ' invalid layout'
!            UFIELD2D_ERR = 2; EXCEPTION = EXCEPTION + 1
!            call ehandler(EDEFAULT,class//sname//erstr); return
!         endif
!         endif
! unpack arguments
!         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
			nblok = 1	!This is only appropriate for when there is no shared memory
!         mx = 2*this%nd1; kyp = this%nd2p; nxv = size(f,1)
			lrec = nrec
			fname = ' '
			if (present(name)) then
				fname = name
				lrec = -lrec
			!            if (present(it)) then
			write (nlabel,'(i5)') it + 10000
			fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
			!            endif
			endif
			order = QUADRATIC
			if (present(inorder)) then
				if ((inorder >= 1) .and. (inorder <= 2)) then
					 order = inorder
				endif
			endif
!         if (.not.present(time)) then
!             time = 0.0
!         endif
! write guard cells if present

			if (order==QUADRATIC) then
 !           if ((mx+2) <= nxv) mx = mx + 1
				lrec=-1
				call PWRITE2_HDF(f,nxv-4,nxv,nypmx,iunit,lrec,&
     &      trim(fname),label_code,it,time,QUADRATIC)
			else
				call PWRITE2_HDF(f,nxv-1,nxv,nypmx,iunit,lrec,&
     &      trim(fname),label_code,it,time,LINEAR)
			endif
		end subroutine ipwrite2_hdf

! this is the wrapper function for pwrite3_hdf
			subroutine ipwrite3_hdf(f,nxv,nypmx,nzpmx,nvpy,nvpz,it,time,label_code,name,inorder,stride)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine collects distributed real 3d data f 
! and it writes to a file.  Only some of the data can be written.
! this = ufield2d descriptor of data
! f = input data to be written
! name = file name (used only if new file is to be opened)
! it = integer to append to file name
! inorder = interpolation order, 1 = LINEAR, 2 = QUADRATIC
			implicit none
			real, dimension(:,:,:,:) :: f
			integer :: nxv, nypmx, nzpmx,nvpy, nvpz
			character(len=*), optional :: name
			integer :: it, label_code
			integer, optional :: inorder
			real :: time
			integer, dimension(3) :: stride
			! local data
			integer :: mx, lrec, order,nblok,i,temp
			character(len=5) :: nlabel
			character(len=60) :: fname
			character(len=10), save :: sname = ':ipwrite3:'
			nblok = 1	!This is only appropriate for when there is no shared memory
			fname = ' '
			if (present(name)) then
				fname = name
			write (nlabel,'(i5)') it + 10000
			fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
			!            endif
			endif
			order = QUADRATIC
			if (present(inorder)) then
				if ((inorder >= 1) .and. (inorder <= 2)) then
					 order = inorder
				endif
			endif
			
			if (order==QUADRATIC) then
				call PWRITE3_HDF(f,nxv-4,nxv,nypmx,nzpmx,nvpy,nvpz,trim(fname), &
				& label_code,it,time,QUADRATIC,stride)
			else
				call PWRITE3_HDF(f,nxv-1,nxv,nypmx,nzpmx,nvpy,nvpz,trim(fname), &
				& label_code,it,time,LINEAR,stride)
			endif
		end subroutine ipwrite3_hdf

!takes the label_code and generates titles  and units for all the different diagnostics
      subroutine do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
      	implicit none
      	integer :: label_code
      	character(len=*) :: title,xtitle,xunits,ytitle,yunits,ztitle,zunits
      	
      	select case(label_code)
      	
      	case(EX)
      		title = 'Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'z'
      		zunits = 'k!N!DD!N'
      	case(EY)
      		title = 'Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(EZ)
      		title = 'Ez'
      		xtitle = 'z'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(JX)
      		title = 'jx'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'j!Dx!N'
      		zunits = 'env!Dth!N'
      	case(JY)
      		title = 'jy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'j!Dy!N'
      		zunits = 'env!Dth!N'
      	case(VDOTEX)
      		title = 'vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY)
      		title = 'vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(DEN)
      		title = 'Electron Density'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Density'
      		zunits = 'n!D0!N'
      	case(PHSL_XX)
      		title = 'Phase Space Slice - vx vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_XY)
      		title = 'Phase Space Slice - vx vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_YX)
      		title = 'Phase Space Slice - vy vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_YY)
      		title = 'Phase Space Slice - vy vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(ESPOYNT_X)
      		title = 'ES Comp. Poynting Vector Px'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dx!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(ESPOYNT_Y)
      		title = 'ES Comp. Poynting Vector Py'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dy!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(FIN_VS_INIT_ENE)
      		title = 'Initial vs. Final Energy'
      		xtitle = 'Initial Ene'
      		xunits = 'mv!Dth!N!U2!N'
      		ytitle = 'Final Ene'
      		yunits = 'mv!Dth!N!U2!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(FIN_VS_INIT_VX)
      		title = 'Initial vs. Final vx'
      		xtitle = 'Initial v!Dx!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dx!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(FIN_VS_INIT_VY)
      		title = 'Initial vs. Final vy'
      		xtitle = 'Initial v!Dy!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dy!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(VDOTEX_PART)
      		title = 'vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_PART)
      		title = 'vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEX_INT)
      		title = 'Time Integrated vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_INT)
      		title = 'Time Integrated vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(FVX_XY_LABEL)
      		title = 'Phase Space vx vs x and y'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N'
      		zunits = '1/v!Dth!N'
      	case(FVY_XY_LABEL)
      		title = 'Phase Space vy vs x and y'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N'
      		zunits = '1/v!Dth!N'
      	case(VXVY)
      		title = 'Phase Space vx vs vy'
      		xtitle = 'v!Dx!N'
      		xunits = '1/v!Dth!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(VDOTEX_FOLLOW_PART)
      		title = 'Summed by following particle vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_FOLLOW_PART)
      		title = 'Summed by following particle vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(ESPOYNT_INT_X)
      		title = 'Time Integrated ES Comp. Poynting Vector Px'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dx!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(ESPOYNT_INT_Y)
      		title = 'Time Integrated ES Comp. Poynting Vector Py'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dy!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(EX_ENE_INT)
      		title = 'x Component of ES Field Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Ex!U2!N'
      		zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	case(EY_ENE_INT)
      		title = 'y Component of ES Field Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Ey!U2!N'
      		zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	case(DIV_ESPOYNTINT)
      		title = 'Time Integrated div ES Poynting Vector'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'div P!D!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(DIVESPOYNT)
      		title = 'div ES Poynting Vector'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'div P!D!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(PH_VX_X)
      		title = 'Phase Space - vx vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VY_X)
      		title = 'Phase Space - vy vs. x'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VX_Y)
      		title = 'Phase Space - vx vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VY_Y)
      		title = 'Phase Space - vy vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(BIN_EX_LABEL)
      		title = 'Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(BIN_EY_LABEL)
      		title = 'Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dy!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(POT)
      		title = 'Potential'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Potential'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N!U2!N'


      	end select
      end subroutine do_labels


		end module hdf_write32_jf
