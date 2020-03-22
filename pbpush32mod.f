!-----------------------------------------------------------------------
!
      module pbpush32d
!
! Fortran90 interface to 3d PIC Fortran77 library pbpush32lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: january 16, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: djpost, push, retard
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PJDOST32(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,np&
     &max,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: cux, cuy, cuz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGJPOST32(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npma&
     &x,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSJPOST32(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npm&
     &ax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PSJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,nx,idimp,n&
     &pmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n81
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n81,npd,mnblok) :: nn
         real, dimension(n81,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PGSJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n81, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n81,npd,mnblok) :: nn
         real, dimension(n81,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PJDOST32L(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,n&
     &pmax,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: cux, cuy, cuz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGJPOST32L(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npm&
     &ax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSJPOST32L(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,np&
     &max,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PSJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,nx,idimp,&
     &npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n24
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n24,npd,mnblok) :: nn
         real, dimension(n24,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PGSJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n24, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n24,npd,mnblok) :: nn
         real, dimension(n24,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PBPUSH32(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,n&
     &x,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz, bx, by, &
     &bz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PBPUSH32L(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,&
     &nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz, bx, by, &
     &bz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx&
     &,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PBPUSH32C(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,&
     &nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz, bx, by, &
     &bz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH32C(part,fxyz,bxyz,npp,noff,qbm,dt,ek,nx,ny,n&
     &z,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PBPUSH32CL(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek&
     &,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz, bx, by, &
     &bz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH32CL(part,fxyz,bxyz,npp,noff,qbm,dt,ek,nx,ny,&
     &nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PRETARD32(part,npp,dtc,nx,ny,nz,idimp,npmax,mnblok,i&
     &pbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, ipbc
         real :: dtc
         real, dimension(idimp,npmax,mnblok) :: part
         integer, dimension(mnblok) :: npp
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure ipgjpost32
      end interface
!
      interface push
         module procedure ipgbpush32
      end interface
!
      interface retard
         module procedure ipretard32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipgjpost32(part,cu,npp,noff,qm,dt,tdjpost,nx,ny,nz,i&
     &pbc,inorder,djopt)
! deposit current, 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: order, opt
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, n24 = 24, n81 = 81
         integer, dimension(n81,npd,size(part,3)) :: nn
         real, dimension(n81,npd,size(part,3)) :: amxyz
         real :: tj
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSJPOST32L(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,np&
     &max,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else if (opt==VECTOR) then
               call PGSJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n24,ipbc)
            else
               call PGJPOST32L(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npm&
     &ax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSJPOST32(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npm&
     &ax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else if (opt==VECTOR) then
               call PGSJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n81,ipbc)
            else
               call PGJPOST32(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npma&
     &x,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         endif
! record time
         call wtimer(tj,dtime)
         tdjpost = tdjpost + tj
         end subroutine ipgjpost32
!
         subroutine ipgbpush32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,tpu&
     &sh,nx,ny,nz,ipbc,inorder,popt)
! push particles with 3d electromagnetic fields, 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: fxyz, bxyz
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); idds = size(noff,1)
         nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx&
     &,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgbpush32
!
         subroutine ipretard32(part,npp,dtc,nx,ny,nz,ipbc)
! retards particle positions half time-step
         implicit none
         integer :: nx, ny, nz, ipbc
         real :: dtc
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, mnblok
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         call PRETARD32(part,npp,dtc,nx,ny,nz,idimp,npmax,mnblok,ipbc)
         end subroutine ipretard32
!
      end module pbpush32d
