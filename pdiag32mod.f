!-----------------------------------------------------------------------
!
      module pdiag32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library pdiag32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: march 28, 2008
!
      use globals, only: LINEAR, QUADRATIC
      use pinit32d, only: idrun, indx, indy, indz, ntp, ntd, nta, nte,  &
     &psolve, tend, dt, omx, omy, omz, ci, t0, ceng, indian, rlprec,    &
     &inprec, pden32d, modesxd, modesyd, modeszd, ndrec, fdname, ppot32d&
     &, modesxp, modesyp, modeszp, nprec, fpname, pvpot32d, modesxa,    &
     &modesya, modesza, narec, faname, pem32d, modesxe, modesye, modesze&
     &, nerec, fename
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: idrun, indx, indy, indz, ntp, ntd, nta, nte, psolve
      public :: tend, dt, omx, omy, omz, ci, t0, ceng, indian, rlprec
      public :: inprec
      public :: pden32d, modesxd, modesyd, modeszd, ndrec, fdname
      public :: ppot32d, modesxp, modesyp, modeszp, nprec, fpname
      public :: pvpot32d, modesxa, modesya, modesza, narec, faname
      public :: pem32d, modesxe, modesye, modesze, nerec, fename
      public :: vdist, bfopen
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PVDIST32(part,fv,fvm,npp,idimp,npmax,mnblok,nmv,nmvf&
     &)
         implicit none
         integer :: idimp, npmax, mnblok, nmv, nmvf
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nmvf,3,mnblok) :: fv
         real, dimension(2,3,mnblok) :: fvm
         integer, dimension(mnblok) :: npp
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface vdist
         module procedure ipvdist32
      end interface
!
      interface bfopen
         module procedure bfopen32
         module procedure bfcopen32
         module procedure bfvcopen32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipvdist32(part,fv,fvm,npp,nmv)
! calculates 3d velocity distribution and velocity moments
         integer :: nmv
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fv, fvm
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, mnblok, nmvf
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); nmvf = size(fv,1)
         call PVDIST32(part,fv,fvm,npp,idimp,npmax,mnblok,nmv,nmvf)
        end subroutine ipvdist32
!
         subroutine bfopen32(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for real 3d scalar data with 2d domain decomposition
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, iunit, nrec
         real, dimension(:,:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1,1,1,1)
         lrec = nx*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfopen32
!
         subroutine bfcopen32(f,nz,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 3d scalar data with 2d domain decomposition
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nz, iunit, nrec
         complex, dimension(:,:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1,1,1,1)
         lrec = nz*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfcopen32
!
         subroutine bfvcopen32(f,nz,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 3d vector data with 2d domain decomposition
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nz, iunit, nrec
         complex, dimension(:,:,:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, nnz, ierr
         if (nrec > 0) return
         nnz = size(f,1)*nz
         inquire(iolength=lrec) f(1,1,1,1,1)
         lrec = nnz*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfvcopen32
!
      end module pdiag32d
