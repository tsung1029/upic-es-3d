!-----------------------------------------------------------------------
!
      module pbfield32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library
! pbfield32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: october 5, 2005
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PMSCGUARD32X, PMSGUARD32X, PMSCGUARD32XL, PMSGUARD32XL
      public :: poism_init, poism, cuperpm, bpoism, ibpoism
      public :: maxwelm, cmfieldm, emfieldm, cpfieldm, avpotm
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PMSCGUARD32X(cu,kstrt,nvpy,noff,nyzp,xj0,yj0,zj0,nx,&
     &ny,ngx,ngy,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvpy, nx, ny, ngx, ngy, nxe, nypmx, nzpmx
         integer :: idds, mblok, nblok
         real, dimension(3,nxe,nypmx,nzpmx,mblok*nblok) :: cu
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PMSGUARD32X(q,kstrt,nvpy,noff,nyzp,qi0,nx,ny,ngx,ngy&
     &,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvpy, nx, ny, ngx, ngy, nxe, nypmx, nzpmx
         integer :: idds, mblok, nblok
         real, dimension(nxe,nypmx,nzpmx,mblok*nblok) :: q
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PMSCGUARD32XL(cu,kstrt,nvpy,noff,nyzp,xj0,yj0,zj0,nx&
     &,ny,ngx,ngy,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvpy, nx, ny, ngx, ngy, nxe, nypmx, nzpmx
         integer :: idds, mblok, nblok
         real, dimension(3,nxe,nypmx,nzpmx,mblok*nblok) :: cu
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PMSGUARD32XL(q,kstrt,nvpy,noff,nyzp,qi0,nx,ny,ngx,ng&
     &y,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvpy, nx, ny, ngx, ngy, nxe, nypmx, nzpmx
         integer :: idds, mblok, nblok
         real, dimension(nxe,nypmx,nzpmx,mblok*nblok) :: q
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PPOISMX32(q,fx,fy,fz,isign,ffb,ax,ay,az,affp,we,nx,n&
     &y,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         complex, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: q, fx, fy,&
     & fz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPOISMX332(q,fxyz,isign,ffb,ax,ay,az,affp,we,nx,ny,n&
     &z,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         complex, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: q
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPOISM332(q,fxyz,isign,ffb,ax,ay,az,affp,we,nx,ny,nz&
     &,kstrt,nzvh,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         real, dimension(2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: q
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PCUPERPMX32(cu,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok&
     &,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPM32(cu,nx,ny,nz,kstrt,nzvh,kxyp2,kyzp2,j2blok&
     &,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISMX332(cu,bxyz,isign,ffb,ax,ay,az,affp,ci,wm,nx&
     &,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ax, ay, az, affp, ci, wm
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PBPOISM332(cu,bxyz,isign,ffb,ax,ay,az,affp,ci,wm,nx,&
     &ny,nz,kstrt,nzvh,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ax, ay, az, affp, ci, wm
         integer :: isign, nx, ny, nz, kstrt, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine IPBPOISMX332(cu,bxyz,ffb,ci,wm,nx,ny,nz,kstrt,nzv,kx&
     &yp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         integer :: nzhd
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine IPBPOISM332(cu,bxyz,ffb,ci,wm,nx,ny,nz,kstrt,nzvh,kx&
     &yp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, nz, kstrt, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PMAXWELMX32(exyz,bxyz,cu,ffb,affp,ci,dt,wf,wm,nx,ny,&
     &nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         integer :: nzhd
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: exyz, bx&
     &yz, cu
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PMAXWELM32(exyz,bxyz,cu,ffb,affp,ci,dt,wf,wm,nx,ny,n&
     &z,kstrt,nzvh,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, nz, kstrt, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: exyz, bx&
     &yz, cu
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PDMFIELDM32(q2,q,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,kyzp2&
     &,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: q2
         real, dimension(2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: q
         end subroutine
      end interface
      interface
         subroutine PCMFIELDM32(cu2,cu,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,kyz&
     &p2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu2
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PEMFIELDM32(fxyz,exyz,ffb,isign,nx,ny,nz,kstrt,nzv,n&
     &zvh,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         implicit none
         integer :: isign, nx, ny, nz, kstrt, nzv, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzhd
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: exyz
         complex, dimension(nzhd,kxyp2,kyzp2,j2blok*m2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPMFIELDM32(pot2,pot,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,k&
     &yzp2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: pot2
         real, dimension(2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: pot
         end subroutine
      end interface
      interface
         subroutine PCPFIELDM32(fxyz,exyz,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,&
     &kyzp2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: exyz
         end subroutine
      end interface
      interface
         subroutine PAVPOTMX332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2&
     &,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         complex, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: bxyz, ax&
     &yz
         end subroutine
      end interface
      interface
         subroutine PAVPOTM332(bxyz,axyz,nx,ny,nz,kstrt,nzvh,kxyp2,kyzp2&
     &,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzvh, kxyp2, kyzp2
         integer :: j2blok, m2blok
         real, dimension(3,2*nzvh,kxyp2,kyzp2,j2blok*m2blok) :: bxyz, ax&
     &yz
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface poism_init
         module procedure ippoism332init
      end interface
!
      interface poism
         module procedure ippoism32
         module procedure ipspoism32
         module procedure ippoism332
      end interface
!
      interface cuperpm
         module procedure ipcuperpmx32
      end interface
!
      interface bpoism
         module procedure jpbpoism332
      end interface
!
      interface ibpoism
         module procedure jipbpoism332
      end interface
!
      interface maxwelm
         module procedure ipmaxwelm32
      end interface
!
      interface cmfieldm
         module procedure ipcmfieldm32
         module procedure ipdmfieldm32
      end interface
!
      interface emfieldm
         module procedure ipemfieldm32
      end interface
!
      interface cpfieldm
         module procedure ipcpfieldm32
      end interface
!
      interface avpotm
         module procedure ipavpotm332
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ippoism32(q,fx,ffb,we,nx,ny,nz,kstrt,j2blok)
! poisson solver for 3d potential, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, fx, ffb
! local data
         integer :: isign = 1, nzv, kxyp2, kyzp2, m2blok, nzhd
         real :: ax, ay, az, affp
         complex, dimension(1,1,1,1) :: fy, fz
         nzv = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         j2blok = size(q,4)/j2blok; nzhd = size(ffb,1)
         call PPOISMX32(q,fx,fy,fz,isign,ffb,ax,ay,az,affp,we,nx,ny,nz,k&
     &strt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine ippoism32
!
         subroutine ipspoism32(q,fy,ffb,nx,ny,nz,kstrt,j2blok)
! smoother for 3d scalar field, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:), pointer :: q, fy, ffb
! local data
         integer :: isign = 2, nzv, kxyp2, kyzp2, m2blok, nzhd
         real :: ax, ay, az, affp, we
         complex, dimension(1,1,1,1) :: fx, fz
         nzv = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         j2blok = size(q,4)/j2blok; nzhd = size(ffb,1)
         call PPOISMX32(q,fx,fy,fz,isign,ffb,ax,ay,az,affp,we,nx,ny,nz,k&
     &strt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine ipspoism32
!
         subroutine ippoism332init(ffb,ax,ay,az,affp,nx,ny,nz,kstrt,j2bl&
     &ok)
! initialize 3d electric field solver,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: ax, ay, az, affp
         complex, dimension(:,:,:,:), pointer :: ffb
! local data
         integer :: isign = 0, nzv, kxyp2, kyzp2, m2blok, nzhd
         real :: we
         complex, dimension(1,1,1,1) :: q
         complex, dimension(3,1,1,1,1) :: fxyz
         nzhd = size(ffb,1); kxyp2 = size(ffb,2); kyzp2 = size(ffb,3)
         m2blok = size(ffb,4)/j2blok; nzv = size(q,1)
         call PPOISMX332(q,fxyz,isign,ffb,ax,ay,az,affp,we,nx,ny,nz,kstr&
     &t,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine ippoism332init
!
         subroutine ippoism332(q,fxyz,ffb,we,nx,ny,nz,kstrt,j2blok)
! poisson solver for 3d electric field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, ffb
         complex, dimension(:,:,:,:,:), pointer :: fxyz
! local data
         integer :: isign = -1, nzv, kxyp2, kyzp2, m2blok, nzhd
         real :: ax, ay, az, affp
         nzv = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok; nzhd = size(ffb,1)
         call PPOISMX332(q,fxyz,isign,ffb,ax,ay,az,affp,we,nx,ny,nz,kstr&
     &t,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine ippoism332
!
         subroutine ipcuperpmx32(cu,nx,ny,nz,kstrt,j2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok
         nzv = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok
         call PCUPERPMX32(cu,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blo&
     &k)
         end subroutine ipcuperpmx32
!
         subroutine ipcuperpm32(cu,nx,ny,nz,kstrt,j2blok)
! calculates transverse part of 3d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nzvh, kxyp2, kyzp2, m2blok
         nzvh = size(cu,2)/2; kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok
         call PCUPERPM32(cu,nx,ny,nz,kstrt,nzvh,kxyp2,kyzp2,j2blok,m2blo&
     &k)
         end subroutine ipcuperpm32
!
         subroutine jpbpoism332(cu,bxyz,ffb,ci,wm,nx,ny,nz,kstrt,j2blok)
! caculates static vector potential for 3d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: ffb
         complex, dimension(:,:,:,:,:), pointer :: cu, bxyz
! local data
         integer :: isign = 1, nzv, kxyp2, kyzp2, m2blok, nzhd
         real :: ax, ay, az, affp
         nzv = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok; nzhd = size(ffb,1)
         call PBPOISMX332(cu,bxyz,isign,ffb,ax,ay,az,affp,ci,wm,nx,ny,nz&
     &,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine jpbpoism332
!
         subroutine jipbpoism332(cu,bxyz,ffb,ci,wm,nx,ny,nz,kstrt,j2blok&
     &)
! calculates static magnetic field for periodic 3d vector field
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: ffb
         real, dimension(:,:,:,:,:), pointer :: cu, bxyz
! local data
         integer :: nzvh, kxyp2, kyzp2, m2blok, nzhd
         nzvh = size(cu,2)/2; kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok; nzhd = size(ffb,1)
         call IPBPOISM332(cu,bxyz,ffb,ci,wm,nx,ny,nz,kstrt,nzvh,kxyp2,ky&
     &zp2,j2blok,m2blok,nzhd)
         end subroutine jipbpoism332
!
         subroutine ipmaxwelm32(exyz,bxyz,cu,ffb,affp,ci,dt,wf,wm,nx,ny,&
     &nz,kstrt,j2blok)
! calculates maxwell's equation for 3d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: affp, ci, dt, wf, wm
         complex, dimension(:,:,:,:), pointer :: ffb
         real, dimension(:,:,:,:,:), pointer :: exyz, bxyz, cu
! local data
         integer :: nzvh, kxyp2, kyzp2, m2blok, nzhd
         nzvh = size(cu,2)/2; kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok; nzhd = size(ffb,1)
         call PMAXWELM32(exyz,bxyz,cu,ffb,affp,ci,dt,wf,wm,nx,ny,nz,kstr&
     &t,nzvh,kxyp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine ipmaxwelm32
!
         subroutine ipcmfieldm32(cu2,cu,nx,ny,nz,kstrt,j2blok)
! copies from double to normal array in y,z dimension for 3d vector data
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: cu2
         real, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok, nzvh
         nzv = size(cu2,2); kxyp2 = size(cu2,3); kyzp2 = size(cu2,4)
         m2blok = size(cu2,5)/j2blok; nzvh = size(cu,2)/2
         call PCMFIELDM32(cu2,cu,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,kyzp2,j2b&
     &lok,m2blok)
         end subroutine ipcmfieldm32
!
         subroutine ipdmfieldm32(q2,q,nx,ny,nz,kstrt,j2blok)
! copies from double to normal array in y,z dimension for 3d vector data
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:), pointer :: q2
         real, dimension(:,:,:,:), pointer :: q
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok, nzvh
         nzv = size(q2,1); kxyp2 = size(q2,2); kyzp2 = size(q2,3)
         m2blok = size(q2,4)/j2blok; nzvh = size(q,1)/2
         call PDMFIELDM32(q2,q,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,kyzp2,j2blo&
     &k,m2blok)
         end subroutine ipdmfieldm32
!
         subroutine ipemfieldm32(fxyz,exyz,ffb,isign,nx,ny,nz,kstrt,j2bl&
     &ok)
! combines and smooths 3d vector fields,
! mixed conducting/periodic boundaries
         implicit none
         integer :: isign, nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: exyz
         complex, dimension(:,:,:,:), pointer :: ffb
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok, nzvh, nzhd
         nzv = size(fxyz,2); kxyp2 = size(fxyz,3); kyzp2 = size(fxyz,4)
         m2blok = size(fxyz,5)/j2blok; nzvh = size(exyz,2)/2
         nzhd = size(ffb,1)
         call PEMFIELDM32(fxyz,exyz,ffb,isign,nx,ny,nz,kstrt,nzv,nzvh,kx&
     &yp2,kyzp2,j2blok,m2blok,nzhd)
         end subroutine ipemfieldm32
!
         subroutine ipcpfieldm32(fxyz,exyz,nx,ny,nz,kstrt,j2blok)
! combines 3d electric fields, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: exyz
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok, nzvh
         nzv = size(fxyz,2); kxyp2 = size(fxyz,3); kyzp2 = size(fxyz,4)
         m2blok = size(fxyz,5)/j2blok; nzvh = size(exyz,2)/2
         call PCPFIELDM32(fxyz,exyz,nx,ny,nz,kstrt,nzv,nzvh,kxyp2,kyzp2,&
     &j2blok,m2blok)
         end subroutine ipcpfieldm32
!
         subroutine ipavpotm332(bxyz,axyz,nx,ny,nz,kstrt,j2blok)
! calculates 3d vector potential from magnetic field
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real, dimension(:,:,:,:,:), pointer :: bxyz, axyz
! local data
         integer :: nzvh, kxyp2, kyzp2, m2blok
         nzvh = size(bxyz,2)/2; kxyp2 = size(bxyz,3)
         kyzp2 = size(bxyz,4); m2blok = size(bxyz,5)/j2blok
         call PAVPOTM332(bxyz,axyz,nx,ny,nz,kstrt,nzvh,kxyp2,kyzp2,j2blo&
     &k,m2blok)
         end subroutine ipavpotm332
!
      end module pbfield32d
