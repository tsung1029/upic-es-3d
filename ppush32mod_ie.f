!-----------------------------------------------------------------------
!
      module ppush32d_ie
!
! Fortran90 interface to parallel 3d PIC Fortran77 library ppush32lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: june 6, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: dpost, push, sortp, countp, prmove
      public :: initmomt3, premoment3, primoment3
!
! define interface to Fortran77 procedures
      interface
         subroutine PDOST32(part,q,npp,noff,qm,nx,idimp,npmax,mnblok,nxv&
     &,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGPOST32(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,n&
     &ypmx,nzpmx,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSPOST32(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,&
     &nypmx,nxyzp,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxyzp,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PSOST32X(part,q,npp,noff,nn,amxyz,qm,nx,idimp,npmax,&
     &mnblok,nxv,nypmx,nxvyzp,idds,npd,n27)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, n27
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n27,npd,mnblok) :: nn
         real, dimension(n27,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PGSOST32X(part,q,npp,noff,nn,amxyz,qm,idimp,npmax,mn&
     &blok,nxv,nypmx,nxvyzp,idds,npd,n27)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, n27
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n27,npd,mnblok) :: nn
         real, dimension(n27,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PDOST32L(part,q,npp,noff,qm,nx,idimp,npmax,mnblok,nx&
     &v,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGPOST32L(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,&
     &nypmx,nzpmx,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSPOST32L(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv&
     &,nypmx,nxyzp,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxyzp,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PSOST32XL(part,q,npp,noff,nn,amxyz,qm,nx,idimp,npmax&
     &,mnblok,nxv,nypmx,nxvyzp,idds,npd,ieight)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, ieight
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(ieight,npd,mnblok) :: nn
         real, dimension(ieight,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PGSOST32XL(part,q,npp,noff,nn,amxyz,qm,idimp,npmax,m&
     &nblok,nxv,nypmx,nxvyzp,idds,npd,ieight)
         implicit none
         integer :: idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, ieight
         real :: qm
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxvyzp,mnblok) :: q
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(ieight,npd,mnblok) :: nn
         real, dimension(ieight,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PPUSH32(part,fx,fy,fz,npp,noff,qbm,dt,ek,nx,idimp,np&
     &max,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGPUSH32(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSPUSH32(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSPUSH32F(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,&
     &nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PPUSH32L(part,fx,fy,fz,npp,noff,qbm,dt,ek,nx,idimp,n&
     &pmax,mnblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGPUSH32L(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSPUSH32L(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idi&
     &mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PSORTP32YZ(part,pt,ip,npic,npp,noff,nyzp,idimp,npmax&
     &,mnblok,nyzpm1,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(npmax,mnblok) :: pt
         integer, dimension(npmax,mnblok) :: ip
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PSORTP32YZL(part,pt,ip,npic,npp,noff,nyzp,idimp,npma&
     &x,mnblok,nyzpm1,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(npmax,mnblok) :: pt
         integer, dimension(npmax,mnblok) :: ip
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PDSORTP32YZ(parta,partb,npic,npp,noff,nyzp,idimp,npm&
     &ax,mnblok,nyzpm1,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds
         real, dimension(idimp,npmax,mnblok) :: parta, partb
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PDSORTP32YZL(parta,partb,npic,npp,noff,nyzp,idimp,np&
     &max,mnblok,nyzpm1,idds)
         implicit none
         integer :: idimp, npmax, mnblok, nyzpm1, idds
         real, dimension(idimp,npmax,mnblok) :: parta, partb
         integer, dimension(nyzpm1,mnblok) :: npic
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PCOUNT32YZL(part,isign,npicyz,npp,noff,nyzp,idimp,np&
     &max,mnblok,myzpm1,idds)
         implicit none
         integer :: isign, idimp, npmax, mnblok, myzpm1, idds
         real, dimension(idimp,npmax,mnblok) :: part
         integer, dimension(myzpm1,idds,mnblok) :: npicyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PRMOVE32(part,npp,ihole,jss,nx,ny,nz,idimp,npmax,mnb&
     &lok,idds,ntmax,ipbc)
         implicit none
         integer nx, ny, nz, idimp, npmax, mnblok, idds, ntmax, ipbc
         real, dimension(idimp,npmax,mnblok) :: part
         integer, dimension(mnblok) :: npp
         integer, dimension(ntmax,mnblok) :: ihole
         integer, dimension(idds,mnblok) :: jss
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure ipgpost32
      end interface
!
      interface push
         module procedure ipgpush32
      end interface
!
      interface sortp
         module procedure ipsortp32yz
         module procedure ipdsortp32yz
      end interface
!
      interface countp
         module procedure ipcount32yz
      end interface
!
      interface prmove
         module procedure iprmove32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine ipgpost32(part,q,qm,npp,noff,tdpost,inorder,dopt)
! deposit charge, 2d partition
         implicit none
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: order, opt
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, ieight = 8, n27 = 27
         integer, dimension(n27,npd,size(part,3)) :: nn
         real, dimension(n27,npd,size(part,3)) :: amxyz
         real :: td
         double precision :: dtime 
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         call wtimer(td,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSPOST32L(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv&
     &,nypmx,nxyzp,idds)
            else if (opt==VECTOR) then
               call PGSOST32XL(part,q,npp,noff,nn,amxyz,qm,idimp,npmax,m&
     &nblok,nxv,nypmx,nxyzp,idds,npd,ieight)
            else
               call PGPOST32L(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,&
     &nypmx,nzpmx,idds)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSPOST32(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,&
     &nypmx,nxyzp,idds)
            else if (opt==VECTOR) then
               call PGSOST32X(part,q,npp,noff,nn,amxyz,qm,idimp,npmax,mn&
     &blok,nxv,nypmx,nxyzp,idds,npd,n27)
            else
               call PGPOST32(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,n&
     &ypmx,nzpmx,idds)
            endif
         endif
! record time
         call wtimer(td,dtime)
         tdpost = tdpost + td
         end subroutine ipgpost32
!
         subroutine ipgpush32(part,fxyz,npp,noff,qbm,dt,ek,tpush,nx,ny,n&
     &z,ipbc,inorder,popt)
! push particles with 3d electrostatic fields, 2d partition
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
         integer :: order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSPUSH32L(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idi&
     &mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGPUSH32L(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               if (idimp > 7) then  ! if tracking forces for test particles
                  call PGSPUSH32F(part,fxyz,npp,noff,qbm,dt,ek,nx&
        &,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
               else
                  call PGSPUSH32(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,i&
        &dimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
               endif
            else
               call PGPUSH32(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgpush32
!
         subroutine ipsortp32yz(part,pt,ip,npp,noff,nyzp,npic,tsort,inor&
     &der)
! sort particles by y-z grid using memory-conserving bin sort
         implicit none
         real :: tsort
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: pt
         integer, dimension(:,:), pointer :: ip, npic
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff, nyzp
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, mnblok, nyzpm1, idds, order
         real :: ts
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); idds = size(noff,1)
         nyzpm1 = size(npic,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call PSORTP32YZL(part,pt,ip,npic,npp,noff,nyzp,idimp,npmax,m&
     &nblok,nyzpm1,idds)
         else
            call PSORTP32YZ(part,pt,ip,npic,npp,noff,nyzp,idimp,npmax,mn&
     &blok,nyzpm1,idds)
         endif
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts
         end subroutine ipsortp32yz
!
         subroutine ipdsortp32yz(parta,partb,npp,noff,nyzp,npic,tsort,in&
     &order)
! sort particles by y-z grid using optimized bin sort
         implicit none
         real :: tsort
         real, dimension(:,:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff, nyzp, npic
         integer, optional :: inorder
! local data
         real, dimension(:,:,:), pointer :: part
         integer :: idimp, npmax, mnblok, nyzpm1, idds, order
         real :: ts
         double precision :: dtime
         idimp = size(parta,1); npmax = size(parta,2)
         mnblok = size(parta,3); idds = size(noff,1)
         nyzpm1 = size(npic,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call PDSORTP32YZL(parta,partb,npic,npp,noff,nyzp,idimp,npmax&
     &,mnblok,nyzpm1,idds)
         else
            call PDSORTP32YZ(parta,partb,npic,npp,noff,nyzp,idimp,npmax,&
     &mnblok,nyzpm1,idds)
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts
         end subroutine ipdsortp32yz
!
         subroutine ipcount32yz(part,npicyz,npp,noff,nyzp,kstrt,nvpy,nvp&
     &z,mblok)
! counts number of particles per cell in y, z grid
         implicit none
         integer :: kstrt, nvpy, nvpz, mblok
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:,:), pointer :: npicyz
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         integer, dimension(size(npicyz,1),1,size(npicyz,3)) :: npicyzt
         integer :: isign = 1
         integer :: idimp, npmax, nblok, mnblok, myzpm1, idds
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok; mnblok = mblok*nblok
         myzpm1 = size(npicyz,1); idds = size(noff,1)
         call PCOUNT32YZL(part,isign,npicyz,npp,noff,nyzp,idimp,npmax,mn&
     &blok,myzpm1,idds)
! integrate y distribution in z direction
         call PISUM2(npicyz(1,1,1),npicyzt,nyzp(1,1),kstrt,nvpy,nvpz,2,m&
     &blok,nblok)
! integrate z distribution in y direction
         call PISUM2(npicyz(1,2,1),npicyzt,nyzp(2,1),kstrt,nvpy,nvpz,1,m&
     &blok,nblok)
         end subroutine ipcount32yz
!
         subroutine iprmove32(part,npp,nx,ny,nz,nbmax,idds,ipbc)
! removes particles which would normally be reflected
         implicit none
         integer :: nx, ny, nz, nbmax, idds, ipbc
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(2*nbmax,size(part,3)) :: ihole
         integer, dimension(idds,size(part,3)) :: jss
         integer :: idimp, npmax, mnblok, ntmax
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         ntmax = 2*nbmax
         call PRMOVE32(part,npp,ihole,jss,nx,ny,nz,idimp,npmax,mnblok,id&
     &ds,ntmax,ipbc)
         end subroutine iprmove32
!
         subroutine initmomt3(part,npp,px,py,pz)
! calculate local initial momentum for each processor, for 3d code
         real :: px, py, pz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: j, m
         double precision :: sum1, sum2, sum3
! calculate momentum at t=t-dt/2
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         do m = 1, size(part,3)
         do j = 1, npp(m)
         sum1 = sum1 + part(4,j,m)
         sum2 = sum2 + part(5,j,m)
         sum3 = sum3 + part(6,j,m)
         enddo
         enddo
         px = sum1
         py = sum2
         pz = sum3
         end subroutine initmomt3
!
         subroutine premoment3(part,itime,npp,id0,iunit,px,py,pz,sx,sy,s&
     &z,wx,wy,wz,nprint)
! print out electron and field momentum, calculate total momentum.
! assumes values of px, py, pz, sx, sy, sz are local to this processor
! for 3d code
         integer :: itime, id0, iunit
         real :: px, py, pz, sx, sy, sz, wx, wy, wz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         integer, optional :: nprint
! local data
         integer :: j, m, np
         real :: tx, ty, tz
         double precision :: sum1, sum2, sum3
         double precision, dimension(6) :: sum6, work6
  991    format (' T = ',i7)
  994    format (' electron momentum = ',3e14.7)
  996    format (' field momentum = ',3e14.7)
         np = 1
         if (present(nprint)) np = nprint
         if (np < 0) return
         if (np==0) then
            px = 0.0; py = 0.0; pz = 0.0
         else if (np==1) then
            if (id0==0) write (iunit,991) itime
         endif
! calculate and print electron momentum
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         do m = 1, size(part,3)
         do j = 1, npp(m)
         sum1 = sum1 + part(4,j,m)
         sum2 = sum2 + part(5,j,m)
         sum3 = sum3 + part(6,j,m)
         enddo
         enddo
         sum6(1) = 0.5*(px + sum1)
         sum6(2) = 0.5*(py + sum2)
         sum6(3) = 0.5*(pz + sum3)
         sum6(4) = sx
         sum6(5) = sy
         sum6(6) = sz
         call PDSUM(sum6,work6,6,1)
         px = sum6(1)
         py = sum6(2)
         pz = sum6(3)
         tx = sum6(4)
         ty = sum6(5)
         tz = sum6(6)
         if (np==1) then
            if (id0==0) then
! print electron momentum
               write (iunit,994) px, py, pz
! print field momentum
               write (iunit,996) tx, ty, tz
            endif
! calculate total momentum
            wx = px + tx
            wy = py + ty
            wz = pz + tz
         endif
         px = sum1
         py = sum2
         pz = sum3
         end subroutine premoment3
!
         subroutine primoment3(parti,nppi,id0,iunit,rmass,px,py,pz,wx,wy&
     &,wz,nprint)
! print out ion momentum, adds total momentum, for 3d code
! assumes values of px, py, pz, are local to this processor
         integer :: id0, iunit
         real :: rmass, px, py, pz, wx, wy, wz
         real, dimension(:,:,:), pointer :: parti
         integer, dimension(:), pointer :: nppi
         integer, optional :: nprint
! local data
         integer :: j, m, np
         real :: at1
         double precision :: sum0, sum1, sum2
         double precision, dimension(3) :: sum3, work3
  995    format (' ion momentum = ',3e14.7)
         np = 1
         if (present(nprint)) np = nprint
         if (np < 0) return
         at1 = 0.5*rmass
         if (np==0) then
            px = 0.0; py = 0.0; pz = 0.0
         endif
! calculate and print ion momentum
         sum0 = 0.0d0
         sum1 = 0.0d0
         sum2 = 0.0d0
         do m = 1, size(parti,3)
         do j = 1, nppi(m)
         sum0 = sum0 + parti(4,j,m)
         sum1 = sum1 + parti(5,j,m)
         sum2 = sum2 + parti(6,j,m)
         enddo
         enddo
         sum3(1) = at1*(px + sum0)
         sum3(2) = at1*(py + sum1)
         sum3(3) = at1*(pz + sum2)
         call PDSUM(sum3,work3,3,1)
         px = sum3(1)
         py = sum3(2)
         pz = sum3(3)
         if (np==1) then
! print ion momentum
            if (id0==0) write (iunit,995)  px, py, pz
! add to total momentum
            wx = wx + px
            wy = wy + py
            wz = wz + pz
         endif
         px = sum0
         py = sum1
         pz = sum2
         end subroutine primoment3
!
      end module ppush32d_ie
