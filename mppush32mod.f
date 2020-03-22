!-----------------------------------------------------------------------
!
      module mppush32d
!
! Fortran90 interface to parallel 3d PIC F77 library mppush32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: april 4, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use ppush32d, only: wtimer, countp, prmove, initmomt3, premoment3,&
     & primoment3
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: dpost, push, sortp, countp, prmove
      public :: initmomt3, premoment3, primoment3
!
! define interface to Fortran77 procedures
      interface
          subroutine MPDOST32(part,q,npp,nps,noff,qm,nx,idimp,npmax,mnbl&
     &ok,nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxv,nypmx,nzpmx,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPOST32(part,q,npp,nps,noff,qm,idimp,npmax,mnblok,&
     &nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxv,nypmx,nzpmx,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPOST32(part,q,npp,nps,noff,qm,idimp,npmax,mnblok&
     &,nxv,nypmx,nxyzp,idds,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxyzp,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxyzp,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSOST32X(part,q,npp,nps,noff,nn,amxyz,qm,nx,idimp,n&
     &pmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n27,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, n27, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n27,npd,mnblok,nmt+1) :: nn
         real, dimension(n27,npd,mnblok,nmt+1) :: amxyz
         real, dimension(nxvyzp,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSOST32X(part,q,npp,nps,noff,nn,amxyz,qm,idimp,npm&
     &ax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n27,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, n27, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n27,npd,mnblok,nmt+1) :: nn
         real, dimension(n27,npd,mnblok,nmt+1) :: amxyz
         real, dimension(nxvyzp,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPDOST32L(part,q,npp,nps,noff,qm,nx,idimp,npmax,mnbl&
     &ok,nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxv,nypmx,nzpmx,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPOST32L(part,q,npp,nps,noff,qm,idimp,npmax,mnblok&
     &,nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxv,nypmx,nzpmx,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPOST32L(part,q,npp,nps,noff,qm,idimp,npmax,mnblo&
     &k,nxv,nypmx,nxyzp,idds,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxyzp,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxyzp,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSOST32XL(part,q,npp,nps,noff,nn,amxyz,qm,nx,idimp,&
     &npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,ieight,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, ieight, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(ieight,npd,mnblok,nmt+1) :: nn
         real, dimension(ieight,npd,mnblok,nmt+1) :: amxyz
         real, dimension(nxvyzp,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSOST32XL(part,q,npp,nps,noff,nn,amxyz,qm,idimp,np&
     &max,mnblok,nxv,nypmx,nxvyzp,idds,npd,ieight,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, ieight, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(ieight,npd,mnblok,nmt+1) :: nn
         real, dimension(ieight,npd,mnblok,nmt+1) :: amxyz
         real, dimension(nxvyzp,mnblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPPUSH32(part,fx,fy,fz,npp,nps,noff,qbm,dt,ek,nx,idi&
     &mp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPPUSH32L(part,fx,fy,fz,npp,nps,noff,qbm,dt,ek,nx,id&
     &imp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,n&
     &z,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr&
     &)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSORTP32YZ(part,pt,ip,npic,npp,nps,noff,nyzp,idimp,&
     &npmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds, nmt, ierr
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(npmax,mnblok) :: pt
         integer, dimension(npmax,mnblok) :: ip
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp, nps, noff
         integer, dimension(idds,mnblok) :: nyzp
         integer, dimension(nyzpm1,mnblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSORTP32YZL(part,pt,ip,npic,npp,nps,noff,nyzp,idimp&
     &,npmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds, nmt, ierr
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(npmax,mnblok) :: pt
         integer, dimension(npmax,mnblok) :: ip
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp, nps, noff
         integer, dimension(idds,mnblok) :: nyzp
         integer, dimension(nyzpm1,mnblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPDSORTP32YZ(parta,partb,npic,npp,nps,noff,nyzp,idim&
     &p,npmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds, nmt, ierr
         real, dimension(idimp,npmax,mnblok) :: parta, partb
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp, nps, noff
         integer, dimension(idds,mnblok) :: nyzp
         integer, dimension(nyzpm1,mnblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPDSORTP32YZL(parta,partb,npic,npp,nps,noff,nyzp,idi&
     &mp,npmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds, nmt, ierr
         real, dimension(idimp,npmax,mnblok) :: parta, partb
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp, nps, noff
         integer, dimension(idds,mnblok) :: nyzp
         integer, dimension(nyzpm1,mnblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure impgpost32
      end interface
!
      interface push
         module procedure impgpush32
      end interface
!
      interface sortp
         module procedure impsortp32yz
         module procedure imdpsortp32yz
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgpost32(part,q,qm,npp,noff,tdpost,inorder,dopt)
! multi-tasking charge deposit, 2d partition
         implicit none
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, ieight = 8, n27 = 27
         integer, dimension(ntasks) :: idtask
         real, dimension(size(q,1),size(q,2),size(q,3),size(q,4),ntasks)&
     & :: qp
         integer, dimension(n27,npd,size(part,3),ntasks+1) :: nn
         real, dimension(n27,npd,size(part,3),ntasks+1) :: amxyz
         real :: td
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         call wtimer(td,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSPOST32L(part,q,npp,nps,noff,qm,idimp,npmax,mnblo&
     &k,nxv,nypmx,nxyzp,idds,qp,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSOST32XL(part,q,npp,nps,noff,nn,amxyz,qm,idimp,np&
     &max,mnblok,nxv,nypmx,nxyzp,idds,npd,ieight,qp,idtask,nmt,ierr)
            else
               call MPGPOST32L(part,q,npp,nps,noff,qm,idimp,npmax,mnblok&
     &,nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSPOST32(part,q,npp,nps,noff,qm,idimp,npmax,mnblok&
     &,nxv,nypmx,nxyzp,idds,qp,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSOST32X(part,q,npp,nps,noff,nn,amxyz,qm,idimp,npm&
     &ax,mnblok,nxv,nypmx,nxyzp,idds,npd,n27,qp,idtask,nmt,ierr)
            else
               call MPGPOST32(part,q,npp,nps,noff,qm,idimp,npmax,mnblok,&
     &nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(td,dtime)
         tdpost = tdpost + td
         end subroutine impgpost32
!
         subroutine impgpush32(part,fxyz,npp,noff,qbm,dt,ek,tpush,nx,ny,&
     &nz,ipbc,inorder,popt)
! multi-tasking particle push with 3d electrostatic fields, 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: fxyz
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,n&
     &z,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr&
     &)
            else
               call MPGPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine impgpush32
!
         subroutine impsortp32yz(part,pt,ip,npp,noff,nyzp,npic,tsort,ino&
     &rder)
! multi-tasking particle sort by y,z grid,
! using memory-conserving bin sort
         implicit none
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: pt
         integer, dimension(:,:), pointer :: ip, npic
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         integer :: idimp, npmax, mnblok, nyzpm1, idds, nmt, order, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(size(npic,1),size(npic,2),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         real :: ts
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); idds = size(nyzp,1)
         nyzpm1 = size(npic,1)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call MPSORTP32YZL(part,pt,ip,npic,npp,nps,noff,nyzp,idimp,np&
     &max,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         else
            call MPSORTP32YZ(part,pt,ip,npic,npp,nps,noff,nyzp,idimp,npm&
     &ax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts
         end subroutine impsortp32yz
!
         subroutine imdpsortp32yz(parta,partb,npp,noff,nyzp,npic,tsort,i&
     &norder)
! multi-tasking particle sort by y,z grid using optimized bin sort
         implicit none
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff, nyzp
         integer, dimension(:,:), pointer :: npic
! local data
         integer :: idimp, npmax, mnblok, nyzpm1, idds, nmt, order, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(size(npic,1),size(npic,2),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         real, dimension(:,:,:), pointer :: part
         real :: ts
         double precision :: dtime
         idimp = size(parta,1); npmax = size(parta,2)
         mnblok = size(parta,3); idds = size(nyzp,1)
         nyzpm1 = size(npic,1)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call MPDSORTP32YZL(parta,partb,npic,npp,nps,noff,nyzp,idimp,&
     &npmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         else
            call MPDSORTP32YZ(parta,partb,npic,npp,nps,noff,nyzp,idimp,n&
     &pmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts
         end subroutine imdpsortp32yz
!
      end module mppush32d
