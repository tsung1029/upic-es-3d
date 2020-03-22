!-----------------------------------------------------------------------
!
      module p32d_ie
!
! Fortran90 interface to 3d parallel PIC Fortran77 library p32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: june 24, 2008
!
      use p0d
      use par_track_ie

      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PPINIT, PPID, PPEXIT, HARTBEAT
      public :: get_funit, pwtimer, plsum, plmax, plbcast
      public :: writebf, readbf, wtimer, wrdata, rddata
      public :: fcomp, dcomp, pmove, pcguard, pcguardp, pnlcguard
      public :: paguard, paguardp, pncguardp, pnlaguard, pnaguardp
      public :: pfmove, repart
      public :: fnoff, trpsin, haftrp, dblsin, hafdbl, ztrp, pmoves
!
! buffer data for particle managers
      real, dimension(:,:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(:,:), allocatable :: ihole
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine DCOMP32(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,i&
     &dds,mblok,nblok)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idps, idds, mblok, nblok
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(idds,mblok*nblok) :: nyzp, noff
         end subroutine
      end interface
      interface
         subroutine DCOMP32L(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,&
     &idds,mblok,nblok)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idps, idds, mblok, nblok
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(idds,mblok*nblok) :: nyzp, noff
         end subroutine
      end interface
      interface
         subroutine PCGUARD32(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblo&
     &k,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PNCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,&
     &nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mbl&
     &ok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzp&
     &mx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: ngds, idds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,ngds,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzp&
     &mx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,ngds,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nyp&
     &mx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PTRPSIN32C(cu,cu3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,n&
     &xv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp
         integer :: kypd, kzpd, kyp2, kzp2, kblok, lblok, k2blok, l2blok
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,2*nxv,2*kypd,kzp2,k2blok*l2blok) :: cu3
         real, dimension(3,nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(3,nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PTRPSIN32D(q,q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nxv&
     &,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp
         integer :: kypd, kzpd, kyp2, kzp2, kblok, lblok, k2blok, l2blok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,2*kypd,kzp2,k2blok*l2blok) :: q3
         real, dimension(nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PHAFTRP32C(fxyz,fxyz3,scb,scd,nx,ny,nz,kstrt,nvpy,nx&
     &v,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd
         integer :: kyp2, kzp2, kblok, lblok, k2blok, l2blok
!        real, dimension(*) :: fxyz
         real :: fxyz
         real, dimension(3,2*nxv,2*kypd,kzp2,k2blok*l2blok) :: fxyz3
         real, dimension(3,nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(3,nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PHAFTRP32D(q,q3,scb,scd,nx,ny,nz,kstrt,nvpy,nxv,kyp,&
     &kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd
         integer :: kyp2, kzp2, kblok, lblok, k2blok, l2blok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,2*kypd,kzp2,k2blok*l2blok) :: q3
         real, dimension(nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PLCGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PLDGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PNLCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNLDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNLCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: ngds, idds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,ngds,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PLAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,ngds,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNLACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNLAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLACGUARDS32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm&
     &x,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PLAGUARDS32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx&
     &,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PNLACGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm&
     &x,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
          subroutine PNLAGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm&
     &x,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,&
     &nzpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, kyp, kzp
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PLAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, kyp, kzp
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNLACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,&
     &nypmx,nzpmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNLAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PMCGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PMDGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PNMCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNMDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PMCGUARD32L(f,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,&
     &nblok,kzp)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: kzp
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         end subroutine
      end interface

      interface
         subroutine PNMCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok
         integer :: ngds, idds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PMACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,ngds,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PMAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,ngds,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNMACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNMAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter, nter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx*ngds,2,mblok*nblok) :: scs
         real, dimension(nxv,nypmx*ngds,2,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PMACGUARDS32(f,scs,kstrt,nvpy,nx,nxv,nypmx,nzpmx,mbl&
     &ok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nx, nxv, nypmx, nzpmx, mblok, nblok
         integer :: kyp, kzp, ngds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PMAGUARDS32(f,scs,kstrt,nvpy,nx,nxv,nypmx,nzpmx,mblo&
     &k,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nx, nxv, nypmx, nzpmx, mblok, nblok
         integer :: kyp, kzp, ngds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,ngds,mblok*nblok) :: scs
         end subroutine
      end interface
      interface
         subroutine PNMACGUARDS32(f,scs,nyzp,kstrt,nvpy,nx,nxv,nypmx,nzp&
     &mx,mblok,nblok,ngds,idds,mter)
         implicit none
         integer :: kstrt, nvpy, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNMAGUARDS32(f,scs,nyzp,kstrt,nvpy,nx,nxv,nypmx,nzpm&
     &x,mblok,nblok,ngds,idds,mter)
         implicit none
         integer :: kstrt, nvpy, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds, mter
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,ngds,mblok*nblok) :: scs
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PMACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,&
     &nzpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, kyp, kzp
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PMAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, kyp, kzp
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,mblok*nblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNMACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,&
     &nypmx,nzpmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds
         real, dimension(3,nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(3,nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(3,nxv,nypmx,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PNMAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds)
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer :: mblok, nblok, ngds, idds
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: f
         real, dimension(nxv,nzpmx,2*ngds,mblok*nblok) :: scs
         real, dimension(nxv,nypmx,mblok*nblok) :: scr
         integer, dimension(idds,mblok*nblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PDBLSIN32C(cu,cu2,scb,scd,nx,ny,kstrt,nvpy,nxv,kyp,k&
     &zp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd, kyp2
         integer :: kblok, lblok, k2blok
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,2*nxv,2*kypd,kzp,k2blok*lblok) :: cu2
         real, dimension(3,nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(3,nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PDBLSIN32D(q,q2,scb,scd,nx,ny,kstrt,nvpy,nxv,kyp,kzp&
     &,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd, kyp2
         integer :: kblok, lblok, k2blok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,2*kypd,kzp,k2blok*lblok) :: q2
         real, dimension(nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PHAFDBL32C(fxyz,fxyz2,scb,scd,nx,ny,kstrt,nvpy,nvpz,&
     &nxv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
         integer :: kyp2, kblok, lblok, k2blok
!        real, dimension(*) :: fxyz
         real :: fxyz
         real, dimension(3,2*nxv,2*kypd,kzp,k2blok*lblok) :: fxyz2
         real, dimension(3,nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(3,nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PHAFDBL32D(q,q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxv,ky&
     &p,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
         integer :: kyp2, kblok, lblok, k2blok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,2*kypd,kzp,k2blok*lblok) :: q2
         real, dimension(nxv,kypd,kzpd,kblok*lblok) :: scb
         real, dimension(nxv,kzpd,2,kblok*lblok) :: scd
         end subroutine
      end interface
      interface
         subroutine PZTRP32D(q,q3,scb,nx,ny,nz,kstrt,nvpy,nxv,kyp,kzp,ky&
     &pd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nxv, kyp, kzp
         integer :: kypd, kzpd, kyp2, kzp2, kblok, lblok, k2blok, l2blok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,2*kypd,kzp2,k2blok*l2blok) :: q3
         real, dimension(nxv,kypd,kzpd,kblok*lblok) :: scb
         end subroutine
      end interface
      interface
         subroutine PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
     &,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nb&
     &max,idds,ntmax,info)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine PMOVE32_tracks(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
     &,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nb&
     &max,idds,ntmax,info,tracks)
         use par_track_ie
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         type(t_track_set) :: tracks
         end subroutine
      end interface
      interface
         subroutine PXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
     &,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nb&
     &max,idds,ntmax,maskp,info)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(npmax,mblok*nblok) :: maskp
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine WPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idp&
     &s,nbmax,idds,ntmax,info)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine WPXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idp&
     &s,nbmax,idds,ntmax,maskp,info)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(npmax,mblok*nblok) :: maskp
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine WPMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,info)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine WPMOVES32_tracks(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,info,tracks)
         use par_track_ie
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         type(t_track_set) :: tracks
         end subroutine
      end interface
      interface
         subroutine WPXMOVS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,maskp,info)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(npmax,mblok*nblok) :: maskp
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine PMOVEH32(part,edges,npp,ihole,jss,idimp,npmax,mblok,&
     &nblok,idps,idds,ntmax,n)
         implicit none
         integer :: idimp, npmax, mblok, nblok, idps, idds, ntmax, n
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         integer, dimension(idds,mblok*nblok) :: jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         end subroutine
      end interface
      interface
         subroutine PMOVEHX32(part,edges,npp,ihole,jss,idimp,npmax,mblok&
     &,nblok,idps,idds,ntmax,maskp,n)
         implicit none
         integer :: idimp, npmax, mblok, nblok, idps, idds, ntmax, n
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(npmax,mblok*nblok) :: maskp
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         integer, dimension(idds,mblok*nblok) :: jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         end subroutine
      end interface
      interface
         subroutine PMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n&
     &bmax,idds,ntmax,info,n)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax, n
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine PMOVESS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,&
     &nbmax,idds,ntmax,info,n)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax, n
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         end subroutine
      end interface
      interface
         subroutine PFMOVE32(f,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,jsr&
     &,jsl,isign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nblok,idd&
     &s,mter,nter,ierr)
         implicit none
         integer :: isign, kyp, kzp, kstrt, nvpy, nvpz, nxv
         integer :: nypmx, nzpmx, mblok, nblok, idds, mter, nter, ierr
!        real, dimension(*) :: f
         real :: f
         real, dimension(nxv,nypmx*nzpmx,mblok*nblok) :: g, h
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         integer, dimension(idds,mblok*nblok) :: noffs, nyzps
         integer, dimension(idds,mblok*nblok) :: noffd, nyzpd
         integer, dimension(idds,mblok*nblok) :: jsl, jsr
         end subroutine
      end interface
      interface
         subroutine REPARTD32(edges,edg,eds,eg,es,et2,npicyz,noff,nyzp,a&
     &npav,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,nvpz,mblok,nblok,idps,&
     &idds,myzpm1)
         implicit none
         integer :: nypmin, nypmax, nzpmin, nzpmax, kstrt, nvpy, nvpz
         integer :: mblok, nblok, idps, idds, myzpm1
         real :: anpav
         real, dimension(idps,mblok*nblok) :: edges
         real, dimension(myzpm1,mblok*nblok) :: edg, eds
         real, dimension(idds,mblok*nblok) :: eg, es
         real, dimension(2*idds,mblok*nblok) :: et2
         integer, dimension(myzpm1,idds,mblok*nblok) :: npicyz
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine FNOFF32(edges,noff,nyzp,nypmin,nypmax,nzpmin,nzpmax,&
     &mnblok,idps,idds)
         implicit none
         integer  :: nypmin, nypmax, nzpmin, nzpmax, mnblok, idps, idds
         real, dimension(idps,mnblok) :: edges
         integer, dimension(idds,mnblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PTPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,&
     &kxypd,kypd,kzpd,jblok,kblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok
         complex, dimension(nxv,kypd,kzpd,kblok*lblok) :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(kxyp,kyp,kzp,kblok*lblok) :: s
         complex, dimension(kxyp,kyp,kzp,jblok*lblok) :: t
         end subroutine
      end interface
      interface
         subroutine PTPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp&
     &,kxypd,kyzpd,kzpd,jblok,mblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(kyzp,kxyp,kzp,jblok*lblok) :: s
         complex, dimension(kyzp,kxyp,kzp,jblok*mblok) :: t
         end subroutine
      end interface
      interface
         subroutine P3TPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp&
     &,kxypd,kypd,kzpd,jblok,kblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok
         complex, dimension(3,nxv,kypd,kzpd,kblok*lblok) :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,kxyp,kyp,kzp,kblok*lblok) :: s
         complex, dimension(3,kxyp,kyp,kzp,jblok*lblok) :: t
         end subroutine
      end interface
      interface
         subroutine P3TPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kz&
     &p,kxypd,kyzpd,kzpd,jblok,mblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(3,kyzp,kxyp,kzp,jblok*lblok) :: s
         complex, dimension(3,kyzp,kxyp,kzp,jblok*mblok) :: t
         end subroutine
      end interface
      interface
         subroutine PTPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxy&
     &pd,kypd,kzpd,jblok,kblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok
         complex, dimension(nxv,kypd,kzpd,kblok*lblok) :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         end subroutine
      end interface
      interface
         subroutine PTPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kx&
     &ypd,kyzpd,kzpd,jblok,mblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         end subroutine
      end interface
      interface
         subroutine P3TPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kx&
     &ypd,kypd,kzpd,jblok,kblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok
         complex, dimension(3,nxv,kypd,kzpd,kblok*lblok) :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         end subroutine
      end interface
      interface
         subroutine P3TPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k&
     &xypd,kyzpd,kzpd,jblok,mblok,lblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         end subroutine
      end interface
      interface
         subroutine PNTPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp&
     &,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok, ndim
         complex, dimension(ndim,nxv,kypd,kzpd,kblok*lblok) :: f
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,kxyp,kyp,kzp,kblok*lblok) :: s
         complex, dimension(ndim,kxyp,kyp,kzp,jblok*lblok) :: t
         end subroutine
      end interface
      interface
         subroutine PNTPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kz&
     &p,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok, ndim
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(ndim,kyzp,kxyp,kzp,jblok*lblok) :: s
         complex, dimension(ndim,kyzp,kxyp,kzp,jblok*mblok) :: t
         end subroutine
      end interface
      interface
         subroutine PNTPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kx&
     &ypd,kypd,kzpd,jblok,kblok,lblok,ndim)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok, ndim
         complex, dimension(ndim,nxv,kypd,kzpd,kblok*lblok) :: f
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         end subroutine
      end interface
      interface
         subroutine PNTPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k&
     &xypd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok, ndim
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd,jblok*mblok) :: h
         end subroutine
      end interface
      interface
         subroutine PN2TPOS3A(f1,f2,g1,g2,s,t,nx,ny,nz,kstrt,nxv,nyv,kxy&
     &p,kyp,kzp,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim1,ndim2)
         implicit none
         integer :: nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
         integer :: kxypd, kypd, kzpd, jblok, kblok, lblok, ndim1, ndim2
         complex, dimension(ndim1,nxv,kypd,kzpd,kblok*lblok) :: f1
         complex, dimension(ndim2,nxv,kypd,kzpd,kblok*lblok) :: f2
         complex, dimension(ndim1,nyv,kxypd,kzpd,jblok*lblok) :: g1
         complex, dimension(ndim2,nyv,kxypd,kzpd,jblok*lblok) :: g2
         complex, dimension(ndim1+ndim2,kxyp,kyp,kzp,kblok*lblok) :: s
         complex, dimension(ndim1+ndim2,kxyp,kyp,kzp,jblok*lblok) :: t
         end subroutine
      end interface
      interface
         subroutine PN2TPOS3B(g1,g2,h1,h2,s,t,nx,ny,nz,kstrt,nyv,nzv,kxy&
     &p,kyzp,kzp,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim1,ndim2)
         implicit none
         integer :: nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
         integer :: kxypd, kyzpd, kzpd, jblok, mblok, lblok
         integer :: ndim1, ndim2
         complex, dimension(ndim1,nyv,kxypd,kzpd,jblok*lblok) :: g1
         complex, dimension(ndim2,nyv,kxypd,kzpd,jblok*lblok) :: g2
         complex, dimension(ndim1,nzv,kxypd,kyzpd,jblok*mblok) :: h1
         complex, dimension(ndim2,nzv,kxypd,kyzpd,jblok*mblok) :: h2
         complex, dimension(ndim1+ndim2,kyzp,kxyp,kzp,jblok*lblok) :: s
         complex, dimension(ndim1+ndim2,kyzp,kxyp,kzp,jblok*mblok) :: t
         end subroutine
      end interface
      interface
         subroutine PSUM2(f,g,nxp,kstrt,nvpy,nvpz,nd,mblok,nblok)
         implicit none
         integer :: nxp, kstrt, nvpy, nvpz, nd, mblok, nblok
         real, dimension(nxp,mblok*nblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PISUM2(if,ig,nxp,kstrt,nvpy,nvpz,nd,mblok,nblok)
         implicit none
         integer :: nxp, kstrt, nvpy, nvpz, nd, mblok, nblok
         integer, dimension(nxp,mblok*nblok) :: if, ig
         end subroutine
      end interface
      interface
         subroutine PSCAN2(f,g,s,nxp,kstrt,nvpy,nvpz,nd,mblok,nblok)
         implicit none
         integer :: nxp, kstrt, nvpy, nvpz, nd, mblok, nblok
         real, dimension(nxp,mblok*nblok) :: f, g, s
         end subroutine
      end interface
      interface
         subroutine PWRITE32(f,nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok,iunit,n&
     &rec,lrec,name)
         implicit none
         integer :: nx, kyp, kzp, nxv, kypmx, kzpmx, mnblok
         integer :: iunit, nrec, lrec
         character(len=*) :: name
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine PREAD32(f,nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok,iunit,nr&
     &ec,lrec,name,ierror)
         implicit none
         integer :: nx, kyp, kzp, nxv, kypmx, kzpmx, mnblok
         integer :: iunit, nrec, lrec, ierror
         character(len=*) :: name
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine PCWRITE32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jbl&
     &ok,mblok,iunit,nrec,lrec,name)
         implicit none
         integer :: nx, ny, nz, kxyp, kyzp, nzv, kxypd, kyzpd
         integer :: jblok, mblok, iunit, nrec, lrec
         character(len=*) :: name
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: f
         complex, dimension(nzv,kxypd*kyzpd,jblok*mblok) :: g
         end subroutine
      end interface
      interface
         subroutine PCREAD32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblo&
     &k,mblok,iunit,nrec,lrec,name,ierror)
         implicit none
         integer :: nx, ny, nz, kxyp, kyzp, nzv, kxypd,kyzpd
         integer :: jblok, mblok, iunit, nrec, lrec, ierror
         character(len=*) :: name
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: f
         complex, dimension(nzv,kxypd*kyzpd,jblok*mblok) :: g
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dcomp
         module procedure idcomp32
      end interface
!
      interface pmove
         module procedure ipmove32_tracks
         module procedure ipdmove32
         module procedure idmove32
!        module procedure iwpmove32
      end interface
!
      interface pmoves
         module procedure iwpmoves32_tracks
      end interface
!
      interface pcguard
         module procedure ipcguard32
         module procedure ipdguard32
      end interface
!
      interface pcguardp
         module procedure ipcguard32p
         module procedure ipdguard32p
      end interface
!
      interface pnlcguard
         module procedure ipnlcguard32
         module procedure ipnldguard32
      end interface
!
      interface pncguardp
         module procedure ipncguard32p
         module procedure ipndguard32p
      end interface
!
      interface paguard
         module procedure ipacguard32
         module procedure ipaguard32
      end interface
!
      interface paguardp
         module procedure ipacguard32p
         module procedure ipaguard32p
      end interface
!
      interface pnlaguard
         module procedure ipnlacguard32
         module procedure ipnlaguard32
      end interface
!
      interface pnaguardp
         module procedure ipnacguard32p
         module procedure ipnaguard32p
      end interface
!
      interface pfmove
         module procedure ipfmove32
         module procedure ipfcmove32
         module procedure ipnfmove32
      end interface
!
      interface repart
         module procedure irepartd32
      end interface
!
      interface fnoff
         module procedure ifnoff32
      end interface
!
      interface trpsin
         module procedure iptrpsin32c
         module procedure iptrpsin32d
      end interface
!
      interface haftrp
         module procedure iphaftrp32c
         module procedure iphaftrp32d
      end interface
!
      interface dblsin
         module procedure ipdblsin32c
         module procedure ipdblsin32d
      end interface
!
      interface hafdbl
         module procedure iphafdbl32c
         module procedure iphafdbl32d
      end interface
!
      interface ztrp
         module procedure izptrp32d
      end interface
!
      interface plsum
         module procedure ipsum2
      end interface
!
      interface plbcast
         module procedure ipbcast2
      end interface
!
      interface writebf
         module procedure ipwrite32
         module procedure ipcwrite32
      end interface
!
      interface readbf
         module procedure ipread32
         module procedure ipcread32
      end interface
!
      interface wrdata
         module procedure ipwrdata32
         module procedure ipwrrdata32
         module procedure ipwrcdata32
      end interface
!
      interface rddata
         module procedure iprddata32
         module procedure iprdrdata32
         module procedure iprdcdata32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!    
         subroutine fcomp(nvp,nx,ny,nz,nvpy,nvpz,ierr)
! determines optimal partition for nvp processors
! input: nvp, number of processors, nx, ny, nz = number of grids
! output: nvpy, nvpz, processors in y, z direction, ierr = error code
         implicit none
         integer, intent(in) :: nvp, nx, ny, nz
         integer, intent(out) :: nvpy, nvpz, ierr
! local data
         integer :: i, lnvp, ly, lz, ind
         ierr = 0
! set lnvp to largest power of 2 contained in nvp
         lnvp = 1
   10    i = lnvp + lnvp
         if (i.le.nvp) then
            lnvp = i
            go to 10
         endif
! check if nvp is a power of 2
         if (nvp.gt.lnvp) then
            write (2,*) 'nvp is not power of 2: nvp, actual nvp used = '&
     &, nvp, lnvp
            ierr = 1
         endif
         ly = min(nx,ny)
         lz = min(ny,nz)
         if (lnvp > ly*lz) then
            write (2,*) 'too many processors to be efficient: nvp, max='&
     &, lnvp, ly*lz
            ly = ny
            lz = nz
            ierr = 2
            if (lnvp > ly*lz) then
               write (2,*) 'too many processors: nvp, actual nvp used= '&
     &, lnvp, ly*lz
            lnvp = ly*lz
            ierr = 3
            endif
         endif
         ind = 0.5*(alog(real(lnvp))/alog(2.) + .00001)
!        ind = sqrt(alog(real(lnvp))/alog(2.) + .00001)
         nvpy = 2**ind
         if (nvpy > ly) nvpy = ly
         nvpz = lnvp/nvpy
         if (nvpz > lz) nvpz = lz
         nvpy = lnvp/nvpz
         end subroutine fcomp
!
         subroutine idcomp32(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,mblok&
     &,inorder)
! find uniform 2d partition boundaries in 3d code
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, mblok
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:,:), pointer :: nyzp, noff
! local data
         integer :: idps, idds, nblok, order
         idps = size(edges,1); nblok = size(edges,2)/mblok
         idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DCOMP32L(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,idd&
     &s,mblok,nblok)
         else
            call DCOMP32(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,idds&
     &,mblok,nblok)
         endif
         end subroutine idcomp32
!
         subroutine ipcguard32(f,kstrt,nvpy,nvpz,kyp,kzp,ngds,mblok,inor&
     &der)
! copy guard cells in y, z for uniform, periodic 3d vector data
         implicit none
         integer :: kstrt, nvpy, nvpz, kyp, kzp, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1)*size(f,2),size(f,4),2*ngds,size(f,5))&
     & :: scs
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,&
     &nblok,kyp,kzp,ngds)
         else
            call PCGUARD32(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,n&
     &blok,kyp,kzp,ngds)
         endif
         end subroutine ipcguard32
!
         subroutine ipdguard32(f,kstrt,nvpy,nvpz,kyp,kzp,ngds,mblok,inor&
     &der)
! copy guard cells in y, z for uniform, periodic 3d scalar data
         implicit none
         integer :: kstrt, nvpy, nvpz, kyp, kzp, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1),size(f,3),2*ngds,size(f,4)) :: scs
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,&
     &nblok,kyp,kzp,ngds)
         else
            call PCGUARD32(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,n&
     &blok,kyp,kzp,ngds)
         endif
         end subroutine ipdguard32
!
         subroutine ipcguard32p(f,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,m&
     &blok,inorder)
! copy guard cells in y, z for uniform 3d vector data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1)*size(f,2),size(f,4),2*ngds,size(f,5))&
     & :: scs
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            nxv = size(f,1)*size(f,2)
            if (order==LINEAR) then
               call PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mbl&
     &ok,nblok,kyp,kzp,ngds)
            else
               call PCGUARD32(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblo&
     &k,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==2) then
            if (order /= LINEAR) then
               call PLCGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==3) then
            if (order /= LINEAR) then
               call PMCGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
            endif
         endif
         end subroutine ipcguard32p
!
         subroutine ipdguard32p(f,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,m&
     &blok,inorder)
! copy guard cells in y, z for uniform 3d scalar data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1),size(f,3),2*ngds,size(f,4)) :: scs
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mbl&
     &ok,nblok,kyp,kzp,ngds)
            else
               call PCGUARD32(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblo&
     &k,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==2) then
            if (order /= LINEAR) then
               call PLDGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==3) then
            if (order /= LINEAR) then
               call PMDGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
            endif
         endif
         end subroutine ipdguard32p
!
         subroutine ipnlcguard32(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ngd&
     &s,mblok,inorder)
! copy guard cells in z for non-uniform 3d vector data, for ipbc=2
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer, dimension(:,:), pointer :: nyzp
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1)*size(f,2),size(f,4)*ngds,2,size(f,5))&
     & :: scs
         real, dimension(size(f,1)*size(f,2),size(f,3)*ngds,2,size(f,5))&
     & :: scr
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            nxv = size(f,1)*size(f,2)
            call PNLCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx&
     &,mblok,nblok,ngds,idds)
         else
            call PNLCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx&
     &,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         endif
         end subroutine ipnlcguard32
!
         subroutine ipnldguard32(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ngd&
     &s,mblok,inorder)
! copy guard cells in z for non-uniform 3d scalar data, for ipbc=2
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1),size(f,3)*ngds,2,size(f,4)) :: scs
         real, dimension(size(f,1),size(f,2)*ngds,2,size(f,4)) :: scr
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PNLCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx&
     &,mblok,nblok,ngds,idds)
         else
            call PNLDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx&
     &,nzpmx,mblok,nblok,ngds,idds,mter,nter)
         endif
         end subroutine ipnldguard32
!
         subroutine ipncguard32p(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ngd&
     &s,ipbc,mblok,inorder)
! copy guard cells in z for non-uniform 3d vector data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer, dimension(:,:), pointer :: nyzp
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1)*size(f,2),size(f,4)*ngds,2,size(f,5))&
     & :: scs
         real, dimension(size(f,1)*size(f,2),size(f,3)*ngds,2,size(f,5))&
     & :: scr
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            nxv = size(f,1)*size(f,2)
            if (order==LINEAR) then
               call PNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzp&
     &mx,mblok,nblok,ngds,idds)
            else
               call PNCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,&
     &nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               nxv = size(f,1)*size(f,2)
               call PNLCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds)
            else
               call PNLCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               nxv = size(f,1)*size(f,2)
               call PNMCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds)
            else
               call PNMCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         endif
         end subroutine ipncguard32p
!
         subroutine ipndguard32p(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ngd&
     &s,ipbc,mblok,inorder)
! copy guard cells in z for non-uniform 3d scalar data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1),size(f,3)*ngds,2,size(f,4)) :: scs
         real, dimension(size(f,1),size(f,2)*ngds,2,size(f,4)) :: scr
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzp&
     &mx,mblok,nblok,ngds,idds)
            else
               call PNCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,&
     &nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PNLCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds)
            else
               call PNLDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PNMCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds)
            else
               call PNMDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         endif
         end subroutine ipndguard32p
!
         subroutine ipacguard32(f,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,mblok,&
     &inorder)
! add guard cells in y, z for uniform, periodic 3d vector data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1),size(f,2),size(f,4),2*ngds,size(f,5))&
     & :: scs
         real, dimension(size(f,1),size(f,2),size(f,3),ngds,size(f,5)) :&
     &: scr
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm&
     &x,mblok,nblok,kyp,kzp,ngds)
         else
            call PACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx&
     &,mblok,nblok,kyp,kzp,ngds)
         endif
         end subroutine ipacguard32
!
         subroutine ipaguard32(f,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,mblok,i&
     &norder)
! add guard cells in y, z for uniform, periodic 3d scalar data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1),size(f,3),2*ngds,size(f,4)) :: scs
         real, dimension(size(f,1),size(f,2),ngds,size(f,4)) :: scr
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx&
     &,mblok,nblok,kyp,kzp,ngds)
         else
            call PAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,&
     &mblok,nblok,kyp,kzp,ngds)
         endif
         end subroutine ipaguard32
!
         subroutine ipacguard32p(f,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,&
     &mblok,inorder)
! add guard cells in z for uniform 3d vector data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1),size(f,2),size(f,4),2*ngds,size(f,5))&
     & :: scs
         real, dimension(size(f,1),size(f,2),size(f,3),ngds,size(f,5)) :&
     &: scr
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
            else
               call PACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               call PLACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,&
     &nzpmx,mblok,nblok,kyp,kzp,ngds)
            else
               call PLACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
               call PLACGUARDS32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm&
     &x,mblok,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               call PMACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,&
     &nzpmx,mblok,nblok,kyp,kzp,ngds)
            else
               call PMACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
               call PMACGUARDS32(f,scs,kstrt,nvpy,nx,nxv,nypmx,nzpmx,mbl&
     &ok,nblok,kyp,kzp,ngds)
            endif
         endif
         end subroutine ipacguard32p
!
         subroutine ipaguard32p(f,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,m&
     &blok,inorder)
! add guard cells in z for uniform 3d scalar data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, nzpmx, nblok, order
         real, dimension(size(f,1),size(f,3),2*ngds,size(f,4)) :: scs
         real, dimension(size(f,1),size(f,2),ngds,size(f,4)) :: scr
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
            else
               call PAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzp&
     &mx,mblok,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               call PLAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
            else
               call PLAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
               call PLAGUARDS32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx&
     &,mblok,nblok,kyp,kzp,ngds)
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               call PMAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,kyp,kzp,ngds)
            else
               call PMAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,kyp,kzp,ngds)
               call PMAGUARDS32(f,scs,kstrt,nvpy,nx,nxv,nypmx,nzpmx,mblo&
     &k,nblok,kyp,kzp,ngds)
            endif
         endif
         end subroutine ipaguard32p
!
         subroutine ipnlacguard32(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ng&
     &ds,mblok,inorder)
! add guard cells in z for non-uniform 3d vector data, for ipbc=2
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1)*size(f,2),size(f,4)*ngds,2,size(f,5))&
     & :: scs
         real, dimension(size(f,1)*size(f,2),size(f,3)*ngds,2,size(f,5))&
     & :: scr
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ngds.eq.1) then
            call PNLACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nyp&
     &mx,nzpmx,mblok,nblok,ngds,idds)
         else
            call PNLACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm&
     &x,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            call PNLACGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,n&
     &zpmx,mblok,nblok,ngds,idds,mter,nter)
         endif
         end subroutine ipnlacguard32
!
         subroutine ipnlaguard32(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ngd&
     &s,mblok,inorder)
! add guard cells in z for non-uniform 3d scalar data, for ipbc=2
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1),size(f,3)*ngds,2,size(f,4)) :: scs
         real, dimension(size(f,1),size(f,2)*ngds,2,size(f,4)) :: scr
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ngds.eq.1) then
            call PNLAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm&
     &x,nzpmx,mblok,nblok,ngds,idds)
         else
            call PNLAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx&
     &,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            call PNLAGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
     &pmx,mblok,nblok,ngds,idds,mter,nter)
         endif
         end subroutine ipnlaguard32
!
         subroutine ipnacguard32p(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ng&
     &ds,ipbc,mblok,inorder)
! add guard cells in z for non-uniform 3d vector data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1)*size(f,2),size(f,4)*ngds,2,size(f,5))&
     & :: scs
         real, dimension(size(f,1)*size(f,2),size(f,3)*ngds,2,size(f,5))&
     & :: scr
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PNACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds)
            else
               call PNACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               call PNLACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,&
     &nypmx,nzpmx,mblok,nblok,ngds,idds)
            else
               call PNLACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
               call PNLACGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm&
     &x,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               call PNMACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,&
     &nypmx,nzpmx,mblok,nblok,ngds,idds)
            else
               call PNMACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
               call PNMACGUARDS32(f,scs,nyzp,kstrt,nvpy,nx,nxv,nypmx,nzp&
     &mx,mblok,nblok,ngds,idds,mter)
            endif
         endif
         end subroutine ipnacguard32p
!
         subroutine ipnaguard32p(f,nyzp,kstrt,nvpy,nvpz,nx,mter,nter,ngd&
     &s,ipbc,mblok,inorder)
! add guard cells in z for non-uniform 3d scalar data
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, mter, nter, ngds, ipbc, mblok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real, dimension(size(f,1),size(f,3)*ngds,2,size(f,4)) :: scs
         real, dimension(size(f,1),size(f,2)*ngds,2,size(f,4)) :: scr
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds)
            else
               call PNAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nyp&
     &mx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               call PNLAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds)
            else
               call PNLAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
               call PNLAGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx&
     &,nzpmx,mblok,nblok,ngds,idds,mter,nter)
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               call PNMAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,n&
     &ypmx,nzpmx,mblok,nblok,ngds,idds)
            else
               call PNMAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,ny&
     &pmx,nzpmx,mblok,nblok,ngds,idds,mter,nter)
               call PNMAGUARDS32(f,scs,nyzp,kstrt,nvpy,nx,nxv,nypmx,nzpm&
     &x,mblok,nblok,ngds,idds,mter)
            endif
         endif
         end subroutine ipnaguard32p
!
         subroutine iptrpsin32c(cu,cu3,nx,ny,nz,kstrt,nvpy,nvpz,kyp,kyp2&
     &,kzp,kzp2,kblok,k2blok,inorder)
! double array in each dimension for 3d vector data
! for dirichlet boundary conditions
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nvpz, kyp, kyp2, kzp, kzp2
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: cu
         real, dimension(:,:,:,:,:), pointer :: cu3
! local data
         integer :: nxv, kypd, kzpd, lblok, l2blok, order
         real, dimension(size(cu,1),size(cu,2),size(cu,3),size(cu,4),siz&
     &e(cu,5)) :: scb
         real, dimension(size(cu,1),size(cu,2),size(cu,4),2,size(cu,5)) &
     &:: scd
         nxv = size(cu,2);  kypd = size(cu,3);  kzpd = size(cu,4)
         lblok = size(cu,5)/kblok; l2blok = size(cu3,5)/k2blok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PTRPSIN32C(cu(1,1,1,1,1),cu3,scb,scd,nx,ny,nz,kstrt,nvp&
     &y,nvpz,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         else
            call PTRPSIN32C(cu(1,2,2,2,1),cu3,scb,scd,nx,ny,nz,kstrt,nvp&
     &y,nvpz,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         endif
         end subroutine iptrpsin32c
!
         subroutine iptrpsin32d(q,q3,nx,ny,nz,kstrt,nvpy,nvpz,kyp,kyp2,k&
     &zp,kzp2,kblok,k2blok,inorder)
! double array in each dimension for 3d scalar data
! for dirichlet boundary conditions
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, nvpz, kyp, kyp2, kzp, kzp2
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: q3
! local data
         integer :: nxv, kypd, kzpd, lblok, l2blok, order
         real, dimension(size(q,1),size(q,2),size(q,3),size(q,4)) :: scb
         real, dimension(size(q,1),size(q,3),2,size(q,4)) :: scd
         nxv = size(q,1);  kypd = size(q,2);  kzpd = size(q,3)
         lblok = size(q,4)/kblok; l2blok = size(q3,4)/k2blok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PTRPSIN32D(q(1,1,1,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nv&
     &pz,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         else
            call PTRPSIN32D(q(2,2,2,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nv&
     &pz,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         endif
         end subroutine iptrpsin32d
!
         subroutine iphaftrp32c(fxyz,fxyz3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,&
     &kzp,kzp2,kblok,k2blok,inorder)
! copy from double to normal array in each dimension for 3d vector data
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, kyp, kyp2, kzp, kzp2
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: fxyz3
! local data
         real, dimension(size(fxyz,1),size(fxyz,2),size(fxyz,3),size(fxy&
     &z,4),size(fxyz,5)) :: scb
         real, dimension(size(fxyz,1),size(fxyz,2),size(fxyz,4),2,size(f&
     &xyz,5)) :: scd
         integer :: nxv, kypd, kzpd, lblok, l2blok, order
         nxv = size(fxyz,2);  kypd = size(fxyz,3);  kzpd = size(fxyz,4)
         lblok = size(fxyz,5)/kblok; l2blok = size(fxyz3,5)/k2blok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PHAFTRP32C(fxyz(1,1,1,1,1),fxyz3,scb,scd,nx,ny,nz,kstrt&
     &,nvpy,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         else
            call PHAFTRP32C(fxyz(1,2,2,2,1),fxyz3,scb,scd,nx,ny,nz,kstrt&
     &,nvpy,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         endif
         end subroutine iphaftrp32c
!
         subroutine iphaftrp32d(q,q3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,kz&
     &p2,kblok,k2blok,inorder)
! copy from double to normal array in each dimension for 3d scalar data
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, kyp, kyp2, kzp, kzp2
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: q3
! local data
         integer :: nxv, kypd, kzpd, lblok, l2blok, order
         real, dimension(size(q,1),size(q,2),size(q,3),size(q,4)) :: scb
         real, dimension(size(q,1),size(q,3),2,size(q,4)) :: scd
         nxv = size(q,1);  kypd = size(q,2);  kzpd = size(q,3)
         lblok = size(q,4)/kblok; l2blok = size(q3,4)/k2blok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PHAFTRP32D(q(1,1,1,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nx&
     &v,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         else
            call PHAFTRP32D(q(2,2,2,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nx&
     &v,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         endif
         end subroutine iphaftrp32d
!
         subroutine ipdblsin32c(cu,cu2,nx,ny,kstrt,nvpy,kyp,kyp2,kzp,kbl&
     &ok,k2blok,inorder)
! double array in x,y dimension for 3d vector data
! for mixed periodic/dirichlet boundary conditions
         implicit none
         integer :: nx, ny, kstrt, nvpy, kyp, kyp2, kzp, kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: cu
         real, dimension(:,:,:,:,:), pointer :: cu2
! local data
         integer :: nxv, kypd, kzpd, lblok, order
         real, dimension(size(cu,1),size(cu,2),size(cu,3),size(cu,4),siz&
     &e(cu,5)) :: scb
         real, dimension(size(cu,1),size(cu,2),size(cu,4),2,size(cu,5)) &
     &:: scd
         nxv = size(cu,2);  kypd = size(cu,3);  kzpd = size(cu,4)
         lblok = size(cu,5)/kblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDBLSIN32C(cu(1,1,1,1,1),cu2,scb,scd,nx,ny,kstrt,nvpy,n&
     &xv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         else
            call PDBLSIN32C(cu(1,2,2,2,1),cu2,scb,scd,nx,ny,kstrt,nvpy,n&
     &xv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         endif
         end subroutine ipdblsin32c
!
         subroutine ipdblsin32d(q,q2,nx,ny,kstrt,nvpy,kyp,kyp2,kzp,kblok&
     &,k2blok,inorder)
! double array in x,y dimension for 3d scalar data
! for mixed periodic/dirichlet boundary conditions
         implicit none
         integer :: nx, ny, kstrt, nvpy, kyp, kyp2, kzp, kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: q2
! local data
         integer :: nxv, kypd, kzpd, lblok, order
         real, dimension(size(q,1),size(q,2),size(q,3),size(q,4)) :: scb
         real, dimension(size(q,1),size(q,3),2,size(q,4)) :: scd
         nxv = size(q,1);  kypd = size(q,2);  kzpd = size(q,3)
         lblok = size(q,4)/kblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDBLSIN32D(q(1,1,1,1),q2,scb,scd,nx,ny,kstrt,nvpy,nxv,k&
     &yp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         else
            call PDBLSIN32D(q(2,2,2,1),q2,scb,scd,nx,ny,kstrt,nvpy,nxv,k&
     &yp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         endif
         end subroutine ipdblsin32d
!
         subroutine iphafdbl32c(fxyz,fxyz2,nx,ny,kstrt,nvpy,nvpz,kyp,kyp&
     &2,kzp,kblok,k2blok,inorder)
! copy from double to normal array in x,y dimension for 3d vector data
         implicit none
         integer :: nx, ny, kstrt, nvpy, nvpz, kyp, kyp2, kzp
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: fxyz2
! local data
         real, dimension(size(fxyz,1),size(fxyz,2),size(fxyz,3),size(fxy&
     &z,4),size(fxyz,5)) :: scb
         real, dimension(size(fxyz,1),size(fxyz,2),size(fxyz,4),2,size(f&
     &xyz,5)) :: scd
         integer :: nxv, kypd, kzpd, lblok, order
         nxv = size(fxyz,2);  kypd = size(fxyz,3);  kzpd = size(fxyz,4)
         lblok = size(fxyz,5)/kblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PHAFDBL32C(fxyz(1,1,1,1,1),fxyz2,scb,scd,nx,ny,kstrt,nv&
     &py,nvpz,nxv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         else
            call PHAFDBL32C(fxyz(1,2,2,2,1),fxyz2,scb,scd,nx,ny,kstrt,nv&
     &py,nvpz,nxv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         endif
         end subroutine iphafdbl32c
!
         subroutine iphafdbl32d(q,q2,nx,ny,kstrt,nvpy,nvpz,kyp,kyp2,kzp,&
     &kblok,k2blok,inorder)
! copy from double to normal array in x,y dimension for 3d scalar data
         implicit none
         integer :: nx, ny, kstrt, nvpy, nvpz, kyp, kyp2, kzp
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: q2
! local data
         integer :: nxv, kypd, kzpd, lblok, order
         real, dimension(size(q,1),size(q,2),size(q,3),size(q,4)) :: scb
         real, dimension(size(q,1),size(q,3),2,size(q,4)) :: scd
         nxv = size(q,1);  kypd = size(q,2);  kzpd = size(q,3)
         lblok = size(q,4)/kblok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PHAFDBL32D(q(1,1,1,1),q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,&
     &nxv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         else
            call PHAFDBL32D(q(2,2,2,1),q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,&
     &nxv,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
         endif
         end subroutine iphafdbl32d
!
         subroutine izptrp32d(q,q3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,kzp2&
     &,kblok,k2blok,inorder)
! double array in each dimension for 3d scalar data, zeroing copies
! for open boundary conditions
         implicit none
         integer :: nx, ny, nz, kstrt, nvpy, kyp, kyp2, kzp, kzp2
         integer :: kblok, k2blok
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: q3
! local data
         integer :: nxv, kypd, kzpd, lblok, l2blok, order
         real, dimension(size(q,1),size(q,2),size(q,3),size(q,4)) :: scb
         nxv = size(q,1);  kypd = size(q,2);  kzpd = size(q,3)
         lblok = size(q,4)/kblok; l2blok = size(q3,4)/k2blok
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PZTRP32D(q(1,1,1,1),q3,scb,nx,ny,nz,kstrt,nvpy,nxv,kyp,&
     &kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         else
            call PZTRP32D(q(2,2,2,1),q3,scb,nx,ny,nz,kstrt,nvpy,nxv,kyp,&
     &kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
         endif
         end subroutine izptrp32d
!
         subroutine ipmove32_tracks(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,&
     &nbmax,idds,mblok,vt,ierr,tracks)
! particle manager: moves particles to appropriate processor
! non-uniform 2d partition boundaries in 3d code
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
         type (t_track_set),optional :: tracks
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(idds,size(part,3)) :: jsl, jsr, jss
         integer, dimension(9) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         idps = size(edges,1)
         ntmax = 2*nbmax
! check if size of buffers has changed
         if (szbuf < idimp*nbmax*size(part,3)) then
            !if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
            if (allocated(sbufl)) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = idimp*nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (.not. present(tracks) ) then
            if (vt.eq.1) then
               call PXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
           &r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax&
           &,idds,ntmax,maskp,info)
            else
               call PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
           &r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax&
           &,idds,ntmax,info)
            endif
         else
            if (vt.eq.1) then
               call PXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
           &r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax&
           &,idds,ntmax,maskp,info)
            else
               call PMOVE32_tracks(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
           &r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax&
           &,idds,ntmax,info,tracks)
            endif
         endif

         ! need to flush or we never see anything in fort.2 on BGL
         ! Note: doing this after each write call seriously degrades performance!
         ! so we do it here instead
         ! flush(2)
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         ierr = info(1)
         end subroutine ipmove32_tracks
!
         subroutine idmove32(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,&
     &nbmax,idds,mblok,info,tinfo,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 2d partition boundaries in 3d code
! debug version
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok
         integer :: ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
         integer, dimension(9) :: info
         real, dimension(13) :: tinfo
! local data
         integer, dimension(idds,size(part,3)) :: jsl, jsr, jss
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         idps = size(edges,1)
         ntmax = 2*nbmax
! check if size of buffers has changed
         if (szbuf < idimp*nbmax*size(part,3)) then
            !if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
            if (allocated(sbufl)) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = idimp*nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         call DMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,j&
     &sl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,id&
     &ds,ntmax,info,tinfo)
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         ierr = info(1)
         end subroutine idmove32
!
         subroutine ipdmove32(part,edges,npp,anpav,pibal,tmove,ny,nz,kst&
     &rt,nvpy,nvpz,nbmax,idds,mblok,vt,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 2d partition boundaries in 3d code
! returns load imbalance
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real :: anpav, pibal
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(idds,size(part,3)) :: jsl, jsr, jss
         integer, dimension(9) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         idps = size(edges,1)
         ntmax = 2*nbmax
! check if size of buffers has changed
         if (szbuf < idimp*nbmax*size(part,3)) then
            !if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
            if (allocated(sbufl)) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = idimp*nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call PXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
     &r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax&
     &,idds,ntmax,maskp,info)
         else
            call PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
     &r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax&
     &,idds,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
! calculate percent imbalance
         anpav = real(info(9))/real(nvpy*nvpz)
         if (anpav > 0.0) then
            pibal = max(real(info(2))-anpav,anpav-real(info(3)))/anpav
         endif
         ierr = info(1)
         end subroutine ipdmove32
!
         subroutine iwpmove32(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz&
     &,nbmax,idds,mblok,vt,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 2d partition boundaries in 3d code
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(idds,size(part,3)) :: jsl, jsr, jss
         integer, dimension(9) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         idps = size(edges,1)
         ntmax = 2*nbmax
         th = 0.0
! check if size of buffers has changed
         if (szbuf < idimp*nbmax*size(part,3)) then
            !if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
            if (allocated(sbufl)) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = idimp*nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call WPXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j&
     &sr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n&
     &bmax,idds,ntmax,maskp,info)
         else
            call WPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j&
     &sr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n&
     &bmax,idds,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine iwpmove32
!
         subroutine iwpdmove32(part,edges,npp,anpav,pibal,tmove,ny,nz,ks&
     &trt,nvpy,nvpz,nbmax,idds,mblok,vt,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 2d partition boundaries in 3d code
! returns load imbalance
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real :: anpav, pibal
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(idds,size(part,3)) :: jsl, jsr, jss
         integer, dimension(9) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         idps = size(edges,1)
         ntmax = 2*nbmax
         th = 0.0
! check if size of buffers has changed
         if (szbuf < idimp*nbmax*size(part,3)) then
            !if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
            if (allocated(sbufl)) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = idimp*nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call WPXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j&
     &sr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n&
     &bmax,idds,ntmax,maskp,info)
         else
            call WPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j&
     &sr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n&
     &bmax,idds,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
! calculate percent imbalance
         anpav = real(info(9))/real(nvpy*nvpz)
         if (anpav > 0.0) then
            pibal = max(real(info(2))-anpav,anpav-real(info(3)))/anpav
         endif
         ierr = info(1)
         end subroutine iwpdmove32
!
         subroutine iwpmoves32_tracks(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvp&
     &z,nbmax,idds,mblok,vt,ierr,tracks)
! particle manager: moves particles to appropriate processor
! non-uniform 2d partition boundaries in 3d code
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
         type (t_track_set),optional :: tracks
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(idds,size(part,3)) :: jsl, jsr, jss
         integer, dimension(9) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         idps = size(edges,1)
         ntmax = 2*nbmax
! set maximum number of passes
         info(6) = 2; info(7) = 2
         th = 0.0
! check if size of buffers has changed
         if (szbuf < idimp*nbmax*size(part,3)) then
            !if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
            if (allocated(sbufl)) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = idimp*nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (.not. present(tracks) ) then
            if (vt.eq.1) then
               call WPXMOVS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&
           &jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,&
           &nbmax,idds,ntmax,maskp,info)
            else
               call WPMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&
           &jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,&
           &nbmax,idds,ntmax,info)
            endif
         else
            if (vt.eq.1) then
               call WPXMOVS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&
           &jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,&
           &nbmax,idds,ntmax,maskp,info)
            else
               call WPMOVES32_tracks(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&
           &jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,&
           &nbmax,idds,ntmax,info,tracks)
            endif
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPABORT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine iwpmoves32_tracks
!
         subroutine ipfcmove32(f,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nv&
     &py,nvpz,mblok,mter,nter,ierr,inorder)
! field manager: moves 3d vector data between uniform and non-uniform
! partitions
         implicit none
         integer :: isign, kyp, kzp, kstrt, nvpy, nvpz, mblok
         integer :: mter, nter, ierr
         real :: tfmove
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         real, dimension(size(f,1),size(f,2),size(f,3),size(f,4),size(f,&
     &5)) :: g, h
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: noffs, nyzps
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: noffd, nyzpd
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: jsl, jsr
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real :: tfm
         double precision :: dtime
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nblok = size(f,5)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfm,dtime,-1)
         if (order==LINEAR) then
            call PFMOVE32(f(1,1,1,1,1),g,h,noff,nyzp,noffs,nyzps,noffd,n&
     &yzpd,jsr,jsl,isign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,n&
     &blok,idds,mter,nter,ierr)
         else
            call PFMOVE32(f(1,1,2,2,1),g,h,noff,nyzp,noffs,nyzps,noffd,n&
     &yzpd,jsr,jsl,isign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,n&
     &blok,idds,mter,nter,ierr)
         endif
! record time
         call wtimer(tfm,dtime)
         tfmove = tfmove + tfm
         end subroutine ipfcmove32
!
         subroutine ipfmove32(f,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvp&
     &y,nvpz,mblok,mter,nter,ierr,inorder)
! field manager: moves 3d scalar data between uniform and non-uniform
! partitions
         implicit none
         integer :: isign, kyp, kzp, kstrt, nvpy, nvpz, mblok
         integer :: mter, nter, ierr
         real :: tfmove
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         real, dimension(size(f,1),size(f,2),size(f,3),size(f,4)) :: g,h
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: noffs, nyzps
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: noffd, nyzpd
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: jsl, jsr
         integer :: nxv, nypmx, nzpmx, nblok, idds, order
         real :: tfm
         double precision :: dtime
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfm,dtime,-1)
         if (order==LINEAR) then
            call PFMOVE32(f(1,1,1,1),g,h,noff,nyzp,noffs,nyzps,noffd,nyz&
     &pd,jsr,jsl,isign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nbl&
     &ok,idds,mter,nter,ierr)
         else
            call PFMOVE32(f(1,2,2,1),g,h,noff,nyzp,noffs,nyzps,noffd,nyz&
     &pd,jsr,jsl,isign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nbl&
     &ok,idds,mter,nter,ierr)
         endif
! record time
         call wtimer(tfm,dtime)
         tfmove = tfmove + tfm
         end subroutine ipfmove32
!
         subroutine ipnfmove32(f,noff,nyzp,noffs,nyzps,isign,tfmove,kstr&
     &t,nvpy,nvpz,mblok,ierr,inorder)
! field manager: moves 3d scalar data between two different non-uniform
! 2d partitions
! noffs and nyzps are modified by this call
         implicit none
         integer :: isign, kstrt, nvpy, nvpz, mblok, ierr
         real :: tfmove
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         real, dimension(size(f,1),size(f,2),size(f,3),size(f,4)) :: g,h
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: noffs, nyzps
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: noffd, nyzpd
         integer, dimension(size(nyzp,1),size(nyzp,2)) :: jsl, jsr
         integer :: mter = 0, nter = 0, kyp = 1, kzp = 1
         integer :: jsign, nxv, nypmx, nzpmx, nblok, idds, order
         real :: tfm
         double precision :: dtime
         nxv = size(f,1); nypmx = size(f,2); nzpmx = size(f,3)
         nblok = size(f,4)/mblok; idds = size(nyzp,1)
         jsign = 0
         if (isign /= 0) jsign = 2
         noffd = noff; nyzpd = nyzp
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfm,dtime,-1)
         if (order==LINEAR) then
            call PFMOVE32(f(1,1,1,1),g,h,noff,nyzp,noffs,nyzps,noffd,nyz&
     &pd,jsr,jsl,jsign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nbl&
     &ok,idds,mter,nter,ierr)
         else
            call PFMOVE32(f(1,2,2,1),g,h,noff,nyzp,noffs,nyzps,noffd,nyz&
     &pd,jsr,jsl,jsign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nbl&
     &ok,idds,mter,nter,ierr)
         endif
! record time
         call wtimer(tfm,dtime)
         tfmove = tfmove + tfm
         end subroutine ipnfmove32
!
         subroutine irepartd32(edges,npicyz,noff,nyzp,anpav,kstrt,nvpy,n&
     &vpz,nypmx,nzpmx,mblok,mterg,nterg,ierr,inorder)
! finds new 2d partitions from old partition and particle information
         implicit none
         integer :: kstrt, nvpy, nvpz, nypmx, nzpmx, mblok, mterg, nterg
         integer :: ierr
         real :: anpav
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:,:,:), pointer :: npicyz
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         real, dimension(size(npicyz,1),size(npicyz,2)) :: edg, eds
         real, dimension(size(nyzp,1),size(nyzp,2)) :: eg, es
         real, dimension(2*size(nyzp,1),size(nyzp,2)) :: et2
         integer :: nypmin, nypmax, nzpmin, nzpmax
         integer :: nblok, idps, myzpm1, idds, order
         idps = size(edges,1); nblok = size(edges,2)/mblok
         myzpm1 = size(npicyz,1); idds = size(nyzp,1)
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
         call REPARTD32(edges,edg,eds,eg,es,et2,npicyz,noff,nyzp,anpav,n&
     &ypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,nvpz,mblok,nblok,idps,idds,m&
     &yzpm1)
         if (order==LINEAR) then
            nypmax = nypmax + 1
            nzpmax = nzpmax + 1
         else
            nypmax = nypmax + 3
            nzpmax = nzpmax + 3
         endif
         if ((nypmin.lt.1).or.(nypmax.gt.nypmx)) then
            write (2,*) 'Field size error: nypmin,nypmax=',nypmin,nypmax
            ierr = 1
         endif
         if ((nzpmin.lt.1).or.(nzpmax.gt.nzpmx)) then
            write (2,*) 'Field size error: nzpmin,nzpmax=',nzpmin,nzpmax
            ierr = ierr + 2
         endif
         mterg = nypmin - 1
         nterg = nzpmin - 1
         end subroutine irepartd32
!
         subroutine ifnoff32(edges,noff,nyzp,nypmx,nzpmx,mterg,nterg,ier&
     &r,inorder)
! finds new 2d partitions arrays from edges
         implicit none
         integer :: nypmx, nzpmx, mterg, nterg, ierr
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         integer :: nypmin, nypmax, nzpmin, nzpmax, mnblok, idps, idds
         integer :: order
         idps = size(edges,1); mnblok = size(edges,2)
         idds = size(nyzp,1)
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
         call FNOFF32(edges,noff,nyzp,nypmin,nypmax,nzpmin,nzpmax,mnblok&
     &,idps,idds)
         if (order==LINEAR) then
            nypmax = nypmax + 1
            nzpmax = nzpmax + 1
         else
            nypmax = nypmax + 3
            nzpmax = nzpmax + 3
         endif
         if ((nypmin.lt.1).or.(nypmax.gt.nypmx)) then
            write (2,*) 'Field size error: nypmin,nypmax=',nypmin,nypmax
            ierr = 1
         endif
         if ((nzpmin.lt.1).or.(nzpmax.gt.nzpmx)) then
            write (2,*) 'Field size error: nzpmin,nzpmax=',nzpmin,nzpmax
            ierr = ierr + 2
         endif
         mterg = nypmin - 1
         nterg = nzpmin - 1
         end subroutine ifnoff32
!
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
!
         subroutine ipbcast2(f)
! broadcast 2d real array
         implicit none
         real, dimension(:,:) :: f
! local data
         integer :: nxp
         nxp = size(f)
         call PBCAST(f,nxp)
         end subroutine ipbcast2
!
         subroutine ipwrite32(f,nx,kyp,kzp,iunit,nrec,name,order)
! collects a subset of a distributed real 3d scalar array and writes it
! to a direct access binary file, for uniform 2d partitions
         implicit none
         integer :: nx, kyp, kzp, iunit, nrec
         integer, optional :: order
         real, dimension(:,:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: lrec, nxv, kypmx, kzpmx, mnblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1,1)
            lrec = nx*kyp*kzp*lrec
         endif
         nxv = size(f,1); kypmx = size(f,2); kzpmx = size(f,3)
         mnblok = size(f,4)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call PWRITE32(f(1,1,1,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblo&
     &k,iunit,nrec,lrec,name)
            else
               call PWRITE32(f(2,2,2,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblo&
     &k,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call PWRITE32(f(1,1,1,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblo&
     &k,iunit,nrec,lrec,noname)
            else
               call PWRITE32(f(2,2,2,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblo&
     &k,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ipwrite32
!
         subroutine ipread32(f,nx,kyp,kzp,iunit,nrec,ierr,name,order)
! reads a subset of a distributed real 3d scalar array from a direct
! access binary file and distributes it, for uniform 2d partitions
         implicit none
         integer :: nx, kyp, kzp, iunit, nrec, ierr
         integer, optional :: order
         real, dimension(:,:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: lrec, nxv, kypmx, kzpmx, mnblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1,1)
            lrec = nx*kyp*kzp*lrec
         endif
         nxv = size(f,1); kypmx = size(f,2); kzpmx = size(f,3)
         mnblok = size(f,4)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call PREAD32(f(1,1,1,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok&
     &,iunit,nrec,lrec,name,ierr)
            else
               call PREAD32(f(2,2,2,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok&
     &,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call PREAD32(f(1,1,1,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok&
     &,iunit,nrec,lrec,noname,ierr)
            else
               call PREAD32(f(2,2,2,1),nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok&
     &,iunit,nrec,lrec,noname,ierr)
            endif
         endif
         end subroutine ipread32
!
         subroutine ipcwrite32(f,nx,ny,nz,kxyp,kyzp,jblok,iunit,nrec,nam&
     &e)
! collects a subset of a distributed complex 3d scalar array and writes
! it to a direct access binary file, for uniform 2d partitions
         implicit none
         integer :: nx, ny, nz, kxyp, kyzp, jblok, iunit, nrec
         complex, dimension(:,:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: lrec, nzv, kxypd, kyzpd, mblok
         character(len=1) :: noname = ' '
         complex, dimension(size(f,1),size(f,2)*size(f,3),size(f,4)) :: &
     &g
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1,1)
            lrec = nz*lrec
         endif
         nzv = size(f,1); kxypd = size(f,2); kyzpd = size(f,3)
         mblok = size(f,4)/jblok
         if (present(name)) then
            call PCWRITE32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblok,&
     &mblok,iunit,nrec,lrec,name)
         else
            call PCWRITE32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblok,&
     &mblok,iunit,nrec,lrec,noname)
         endif
         end subroutine ipcwrite32
!
         subroutine ipcread32(f,nx,ny,nz,kxyp,kyzp,jblok,iunit,nrec,ierr&
     &,name)
! reads a subset of a distributed complex 3d scalar array from a direct
! access binary file and distributes it, for uniform 2d partitions
         implicit none
         integer :: nx, ny, nz, kxyp, kyzp, jblok, iunit, nrec, ierr
         complex, dimension(:,:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: lrec, nzv, kxypd, kyzpd, mblok
         character(len=1) :: noname = ' '
         complex, dimension(size(f,1),size(f,2)*size(f,3),size(f,4)) :: &
     &g
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1,1)
            lrec = nz*lrec
         endif
         nzv = size(f,1); kxypd = size(f,2); kyzpd = size(f,3)
         mblok = size(f,4)/jblok
         if (present(name)) then
            call PCREAD32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblok,m&
     &blok,iunit,nrec,lrec,name,ierr)
         else
            call PCREAD32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblok,m&
     &blok,iunit,nrec,lrec,noname,ierr)
         endif
         end subroutine ipcread32
!
         subroutine ipwrdata32(f,nvp,iunit)
! collects distributed real 3d scalar data and writes it to a
! fortran unformatted sequential file, for uniform 2d partitions
         implicit none
         integer :: nvp, iunit
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kzp, nblok
         nxv = size(f,1); nypmx = size(f,2); kzp = size(f,3)
         nblok = size(f,4)
         call PWRDATA(f,nvp,nxv*nypmx*kzp,nblok,iunit)
         end subroutine ipwrdata32
!
         subroutine iprddata32(f,nvp,iunit,ierror)
! reads real 3d scalar data from a fortran unformatted sequential file
! and distributes it, for uniform 2d partitions
         implicit none
         integer :: nvp, iunit, ierror
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kzp, nblok
         nxv = size(f,1); nypmx = size(f,2); kzp = size(f,3)
         nblok = size(f,4)
         call PRDDATA(f,nvp,nxv*nypmx*kzp,nblok,iunit,ierror)
         end subroutine iprddata32
!
         subroutine ipwrrdata32(f,nvp,iunit)
! collects distributed real 3d vector data and writes it to a
! fortran unformatted sequential file, for uniform 2d partitions
         implicit none
         integer :: nvp, iunit
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kzp, nblok
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); kzp = size(f,4)
         nblok = size(f,5)
         call PWRDATA(f,nvp,nxv*nypmx*kzp,nblok,iunit)
         end subroutine ipwrrdata32
!
         subroutine iprdrdata32(f,nvp,iunit,ierror)
! reads real 3d vector data from a fortran unformatted sequential file
! and distributes it, for uniform 2d partitions
         implicit none
         integer :: nvp, iunit, ierror
         real, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kzp, nblok
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); kzp = size(f,4)
         nblok = size(f,5)
         call PRDDATA(f,nvp,nxv*nypmx*kzp,nblok,iunit,ierror)
         end subroutine iprdrdata32
!
         subroutine ipwrcdata32(f,nvp,iunit)
! collects distributed complex 3d vector data and writes it to a
! fortran unformatted sequential file, for uniform 2d partitions
         implicit none
         integer :: nvp, iunit
         complex, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kzp, nblok
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); kzp = size(f,4)
         nblok = size(f,5)
         call PWRDATA(f,nvp,2*nxv*nypmx*kzp,nblok,iunit)
         end subroutine ipwrcdata32
!
         subroutine iprdcdata32(f,nvp,iunit,ierror)
! reads complex 3d vector data from a fortran unformatted sequential file
! and distributes it, for uniform 2d partitions
         implicit none
         integer :: nvp, iunit, ierror
         complex, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kzp, nblok
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); kzp = size(f,4)
         nblok = size(f,5)
         call PRDDATA(f,nvp,2*nxv*nypmx*kzp,nblok,iunit,ierror)
         end subroutine iprdcdata32
!
      end module p32d_ie
