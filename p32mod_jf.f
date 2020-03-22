!-----------------------------------------------------------------------
! This module adds some functions written by Jay Fahlen to fix problems
! with the original p2mod.f file by Viktor Decyk.  Feb. 2008
		module p32d_jf
		implicit none
		private
		public :: write_jf,plsum

		interface write_jf
			module procedure write_jf3d
			module procedure write_jf1d
		end interface

		interface plsum
			 module procedure ipsum2
			 module procedure ipsum3
		end interface
		
		contains

		subroutine ipsum2(f)
		! perform global sum of 2d real array
		implicit none
		real, dimension(:,:) :: f
		! local data
		integer :: nxyp, nblok
		real, dimension(size(f,1),size(f,2)) :: g
		nxyp = size(f,1)*size(f,2); nblok = 1
		call PSUM(f,g,nxyp,nblok)
		end subroutine ipsum2

		subroutine ipsum3(f)
! perform global sum of 2d real array
			implicit none
			real, dimension(:,:,:) :: f
! local data
			integer :: nxyp, nblok
!			real, dimension(size(f,1),size(f,2),size(f,3)) :: g
			real,dimension(:,:,:),allocatable,save :: g
			
			if(.not. allocated(g)) then
				allocate( g(size(f,1),size(f,2),size(f,3)))
			endif
			
			nxyp = size(f,1)*size(f,2)*size(f,3); nblok = 1
			call PSUM(f,g,nxyp,nblok)
		end subroutine ipsum3

		subroutine write_jf1d(f,nx,iunit,nrec,name,order)
! write a single 1 dim array to file, it does not collect from other nodes
         implicit none
         integer :: nx, iunit, nrec,j,idproc,ierr
         real, dimension(:) :: f
         character(len=*), optional :: name
         integer, optional :: order
	      integer :: nproc, lgrp, lstat=10, mreal, mint, mcplx, mdouble,&
	      	lworld
	      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
         integer :: lrec, nxv, kypmx, nblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1)
            lrec = nx*lrec
         endif
         nxv = size(f,1)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
	      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
!   	   call MPI_COMM_SIZE(lworld,nvp,ierr)
! node 0 receives messages from other nodes
			if (idproc.eq.0) then

				if (nrec.lt.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='replace')
					nrec = 1
	! open old file
				else if (nrec.eq.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='old')
				endif
				write (unit=iunit,rec=nrec) (f(j),j=1,nx)
				nrec = nrec + 1
			endif
		end subroutine write_jf1d
		
		subroutine write_jf3d(f,nx,kyp,iunit,nrec,name,order)
! write a single 3 dim array to file, it does not collect from other nodes
         implicit none
         integer :: nx, kyp, iunit, nrec,k,l,j,idproc,ierr
         real, dimension(:,:,:) :: f
         character(len=*), optional :: name
         integer, optional :: order
	      integer :: nproc, lgrp, lstat=10, mreal, mint, mcplx, mdouble,&
	      	lworld
	      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
         integer :: lrec, nxv, kypmx, nblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = nx*kyp*lrec
         endif
         nxv = size(f,1); kypmx = size(f,2); nblok = size(f,3)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
	      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
!   	   call MPI_COMM_SIZE(lworld,nvp,ierr)
! node 0 receives messages from other nodes
			if (idproc.eq.0) then
	
				if (nrec.lt.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='replace')
					nrec = 1
	! open old file
				else if (nrec.eq.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='old')
				endif
				write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,&
					nblok)
				nrec = nrec + 1
			endif
		end subroutine write_jf3d
		
		end module p32d_jf
