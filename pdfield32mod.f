!-----------------------------------------------------------------------
!
      module pdfield32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library
! pdfield32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: october 27, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PLCGUARD32X, PLDGUARD32X, PLSCGUARD32X, PLSGUARD32X
      public :: PLACGUARD32X, PLAGUARD32X, PLSCGUARD32XL, PLSGUARD32XL
      public :: laguard, lcguard
      public :: poisd_init, poisd, cuperpd, bpoisd, ibpoisd
      public :: maxweld, cmfieldd, emfieldd, cpfieldd, avpotd
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PLCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD32X(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,zj&
     &0,nx,ny,nz,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz
         integer :: nxe, nypmx, nzpmx, idds, mblok, nblok
         real, dimension(3,nxe,nypmx,nzpmx,mblok*nblok) :: cu
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PLSGUARD32X(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,nz&
     &,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz
         integer :: nxe, nypmx, nzpmx, idds, mblok, nblok
         real, dimension(nxe,nypmx,nzpmx,mblok*nblok) :: q
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD32X(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
!        real, dimension(*) :: cu
         real :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLAGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
!        real, dimension(*) :: q
         real :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD32XL(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,z&
     &j0,nx,ny,nz,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz
         integer :: nxe, nypmx, nzpmx, idds, mblok, nblok
         real, dimension(3,nxe,nypmx,nzpmx,mblok*nblok) :: cu
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PLSGUARD32XL(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,n&
     &z,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz
         integer :: nxe, nypmx, nzpmx, idds, mblok, nblok
         real, dimension(nxe,nypmx,nzpmx,mblok*nblok) :: q
         integer, dimension(idds,mblok*nblok) :: noff, nyzp
         end subroutine
      end interface
      interface
         subroutine PPOISDX32(q,fx,fy,fz,isign,ffd,ax,ay,az,affp,we,nx,n&
     &y,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: q, fx, fy&
     &, fz
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISDX332(q,fxyz,isign,ffd,ax,ay,az,affp,we,nx,ny,n&
     &z,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: q
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISD332(q,fxyz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz&
     &,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         real, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: q
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PCUPERPDX32(cu,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blo&
     &k,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPD32(cu,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,&
     &m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISDX332(cu,bxyz,isign,ffd,ax,ay,az,affp,ci,wm,nx&
     &,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ax, ay, az, affp, ci, wm
         integer :: isign, nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxy&
     &z
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PBPOISD332(cu,bxyz,isign,ffd,ax,ay,az,affp,ci,wm,nx,&
     &ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ax, ay, az, affp, ci, wm
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxyz
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine IPBPOISDX332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,nz2d,k&
     &xyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxy&
     &z
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine IPBPOISD332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,nzv,kxy&
     &p2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         integer :: nzd
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu, bxyz
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMAXWELDX32(exyz,bxyz,cu,ffd,affp,ci,dt,wf,wm,nx,ny,&
     &nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: exyz, b&
     &xyz, cu
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMAXWELD32(exyz,bxyz,cu,ffd,affp,ci,dt,wf,wm,nx,ny,n&
     &z,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         integer :: nzd
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: exyz, bxyz,&
     & cu
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PDMFIELDD32(q3,q,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp2&
     &,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: q3
         real, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: q
         end subroutine
      end interface
      interface
         subroutine PCMFIELDD32(cu3,cu,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyz&
     &p2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: cu3
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PEMFIELDD32(fxyz,exyz,ffd,isign,nx,ny,nz,kstrt,nz2d,&
     &nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
         implicit none
         integer :: isign, nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok, nzd
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: exyz
         complex, dimension(nzd,kxyp2,kyzp2,j2blok*m2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPMFIELDD32(pot3,pot,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,k&
     &yzp2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: pot3
         real, dimension(nzv,kxyp2,kyzp2,j2blok*m2blok) :: pot
         end subroutine
      end interface
      interface
         subroutine PCPFIELDD32(fxyz,exyz,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,&
     &kyzp2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: exyz
         end subroutine
      end interface
      interface
         subroutine PAVPOTDX332(bxyz,axyz,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp&
     &2,j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: bxyz, a&
     &xyz
         end subroutine
      end interface
      interface
         subroutine PAVPOTD332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,&
     &j2blok,m2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp2, kyzp2, j2blok, m2blok
         real, dimension(3,nzv,kxyp2,kyzp2,j2blok*m2blok) :: bxyz, axyz
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!      
      interface laguard
         module procedure iplacguard32x
         module procedure iplaguard32x
      end interface
!
      interface lcguard
         module procedure iplcguard32x
         module procedure ipldguard32x
      end interface
!
      interface poisd_init
         module procedure ippoisd332init
      end interface
!
      interface poisd
         module procedure ippoisd32
         module procedure ipspoisd32
         module procedure ippoisd332
      end interface
!
      interface cuperpd
         module procedure ipcuperpdx32
      end interface
!
      interface bpoisd
         module procedure jpbpoisd332
      end interface
!
      interface ibpoisd
         module procedure jipbpoisd332
      end interface
!
      interface maxweld
         module procedure ipmaxweld32
      end interface
!
      interface cmfieldd
         module procedure ipcmfieldd32
         module procedure ipdmfieldd32
      end interface
!
      interface emfieldd
         module procedure ipemfieldd32
      end interface
!
      interface cpfieldd
         module procedure ipcpfieldd32
      end interface
!
      interface avpotd
         module procedure ipavpotd332
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine iplacguard32x(cu,nyzp,nx,inorder)
! add guard cells in x,for non-uniform 3d vector data
! copy fields, disable quadratic interpolation at edge in x direction
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:,:) :: cu
         integer, dimension(:,:) :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         mnblok = size(cu,5); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call PLACGUARD32X(cu(1,2,1,1,1),nyzp,nx-2,nxe,nypmx,nzpmx,id&
     &ds,mnblok)
         endif
         end subroutine iplacguard32x
!
         subroutine iplaguard32x(q,nyzp,nx,inorder)
! add guard cells in x for non-uniform 3d scalar data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:) :: q
         integer, dimension(:,:) :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call PLAGUARD32X(q(2,1,1,1),nyzp,nx-2,nxe,nypmx,nzpmx,idds,m&
     &nblok)
         endif
         end subroutine iplaguard32x
!
         subroutine iplcguard32x(fxyz,nyzp,nx,inorder)
! copy guard cells in x for non-uniform 3d vector data
! disable quadratic interpolation at edge in x direction
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:,:) :: fxyz
         integer, dimension(:,:) :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         mnblok = size(fxyz,5)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call PLCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine iplcguard32x
!
         subroutine ipldguard32x(q,nyzp,nx,inorder)
! copy guard cells in x for non-uniform 3d scalar data
! disable quadratic interpolation at edge in x direction
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:) :: q
         integer, dimension(:,:) :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call PLDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipldguard32x
!
         subroutine ippoisd32(q,fx,ffd,we,nx,ny,nz,kstrt,j2blok)
! poisson solver for 3d potential, conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, fx, ffd
! local data
         integer :: isign = 1, nz2d, kxyp2, kyzp2, m2blok, nzd
         real :: ax, ay, az, affp
         complex, dimension(1,1,1,1) :: fy, fz
         nz2d = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok; nzd = size(ffd,1)
         call PPOISDX32(q,fx,fy,fz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz,k&
     &strt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine ippoisd32
!
         subroutine ipspoisd32(q,fy,ffd,nx,ny,nz,kstrt,j2blok)
! smoother for 3d scalar field, conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:), pointer :: q, fy, ffd
! local data
         integer :: isign = 2, nz2d, kxyp2, kyzp2, m2blok, nzd
         real :: ax, ay, az, affp, we
         complex, dimension(1,1,1,1) :: fx, fz
         nz2d = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok; nzd = size(ffd,1)
         call PPOISDX32(q,fx,fy,fz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz,k&
     &strt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine ipspoisd32
!
         subroutine ippoisd332init(ffd,ax,ay,az,affp,nx,ny,nz,kstrt,j2bl&
     &ok)
! initialize 3d electric field solver, conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: ax, ay, az, affp
         complex, dimension(:,:,:,:), pointer :: ffd
! local data
         integer :: isign = 0, nz2d, kxyp2, kyzp2, m2blok, nzd
         real :: we
         complex, dimension(1,1,1,1) :: q
         complex, dimension(3,1,1,1,1) :: fxyz
         nzd = size(ffd,1); kxyp2 = size(ffd,2); kyzp2 = size(ffd,3)
         m2blok = size(ffd,4)/j2blok; nz2d = size(q,1)
         call PPOISDX332(q,fxyz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz,kstr&
     &t,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine ippoisd332init
!
         subroutine ippoisd332(q,fxyz,ffd,we,nx,ny,nz,kstrt,j2blok)
! poisson solver for 3d electric field, conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, ffd
         complex, dimension(:,:,:,:,:), pointer :: fxyz
! local data
         integer :: isign = -1, nz2d, kxyp2, kyzp2, m2blok, nzd
         real :: ax, ay, az, affp
         nz2d = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok; nzd = size(ffd,1)
         call PPOISDX332(q,fxyz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz,kstr&
     &t,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine ippoisd332
!
         subroutine ipcuperpdx32(cu,nx,ny,nz,kstrt,j2blok)
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nz2d, kxyp2, kyzp2, m2blok
         nz2d = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok
         call PCUPERPDX32(cu,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2bl&
     &ok)
         end subroutine ipcuperpdx32
!
         subroutine ipcuperpd32(cu,nx,ny,nz,kstrt,j2blok)
! calculates transverse part of 3d vector field, conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok
         nzv = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok
         call PCUPERPD32(cu,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok&
     &)
         end subroutine ipcuperpd32
!
         subroutine jpbpoisd332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,j2blok)
! caculates static vector potential for 3d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: ffd
         complex, dimension(:,:,:,:,:), pointer :: cu, bxyz
! local data
         integer :: isign = 1, nz2d, kxyp2, kyzp2, m2blok, nzd
         real :: ax, ay, az, affp
         nz2d = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok; nzd = size(ffd,1)
         call PBPOISDX332(cu,bxyz,isign,ffd,ax,ay,az,affp,ci,wm,nx,ny,nz&
     &,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine jpbpoisd332
!
         subroutine jipbpoisd332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,j2blok&
     &)
! calculates static magnetic field for periodic 3d vector field
! conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: ffd
         real, dimension(:,:,:,:,:), pointer :: cu, bxyz
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok, nzd
         nzv = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok; nzd = size(ffd,1)
         call IPBPOISD332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,nzv,kxyp2,kyz&
     &p2,j2blok,m2blok,nzd)
         end subroutine jipbpoisd332
!
         subroutine ipmaxweld32(exyz,bxyz,cu,ffd,affp,ci,dt,wf,wm,nx,ny,&
     &nz,kstrt,j2blok)
! calculates maxwell's equation for 3d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: affp, ci, dt, wf, wm
         complex, dimension(:,:,:,:), pointer :: ffd
         real, dimension(:,:,:,:,:), pointer :: exyz, bxyz, cu
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok, nzd
         nzv = size(cu,2); kxyp2 = size(cu,3); kyzp2 = size(cu,4)
         m2blok = size(cu,5)/j2blok; nzd = size(ffd,1)
         call PMAXWELD32(exyz,bxyz,cu,ffd,affp,ci,dt,wf,wm,nx,ny,nz,kstr&
     &t,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine ipmaxweld32
!
         subroutine ipcmfieldd32(cu3,cu,nx,ny,nz,kstrt,j2blok)
! copies from double to normal array in y,z dimension for 3d vector data
! conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: cu3
         real, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nz2d, kxyp2, kyzp2, m2blok, nzv
         nz2d = size(cu3,2); kxyp2 = size(cu3,3); kyzp2 = size(cu3,4)
         m2blok = size(cu3,5)/j2blok; nzv = size(cu,2)
         call PCMFIELDD32(cu3,cu,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp2,j2b&
     &lok,m2blok)
         end subroutine ipcmfieldd32
!
         subroutine ipdmfieldd32(q3,q,nx,ny,nz,kstrt,j2blok)
! copies from double to normal array in y,z dimension for 3d scalar data
! conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:), pointer :: q3
         real, dimension(:,:,:,:), pointer :: q
! local data
         integer :: nz2d, kxyp2, kyzp2, m2blok, nzv
         nz2d = size(q3,1); kxyp2 = size(q3,2); kyzp2 = size(q3,3)
         m2blok = size(q3,4)/j2blok; nzv = size(q,1)
         call PDMFIELDD32(q3,q,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp2,j2blo&
     &k,m2blok)
         end subroutine ipdmfieldd32
!
         subroutine ipemfieldd32(fxyz,exyz,ffd,isign,nx,ny,nz,kstrt,j2bl&
     &ok)
! combines and smooths 3d vector fields, conducting boundaries
         implicit none
         integer :: isign, nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: exyz
         complex, dimension(:,:,:,:), pointer :: ffd
! local data
         integer :: nz2d, kxyp2, kyzp2, m2blok, nzv, nzd
         nz2d = size(fxyz,2); kxyp2 = size(fxyz,3); kyzp2 = size(fxyz,4)
         m2blok = size(fxyz,5)/j2blok; nzv = size(exyz,2)
         nzd = size(ffd,1)
         call PEMFIELDD32(fxyz,exyz,ffd,isign,nx,ny,nz,kstrt,nz2d,nzv,kx&
     &yp2,kyzp2,j2blok,m2blok,nzd)
         end subroutine ipemfieldd32
!
         subroutine ipcpfieldd32(fxyz,exyz,nx,ny,nz,kstrt,j2blok)
! combines 3d electric fields, conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: exyz
! local data
         integer :: nz2d, kxyp2, kyzp2, m2blok, nzv
         nz2d = size(fxyz,2); kxyp2 = size(fxyz,3); kyzp2 = size(fxyz,4)
         m2blok = size(fxyz,5)/j2blok; nzv = size(exyz,2)
         call PCPFIELDD32(fxyz,exyz,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp2,&
     &j2blok,m2blok)
         end subroutine ipcpfieldd32
!
         subroutine ipavpotd332(bxyz,axyz,nx,ny,nz,kstrt,j2blok)
! calculates 3d vector potential from magnetic field
! conducting boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real, dimension(:,:,:,:,:), pointer :: bxyz, axyz
! local data
         integer :: nzv, kxyp2, kyzp2, m2blok
         nzv = size(bxyz,2); kxyp2 = size(bxyz,3); kyzp2 = size(bxyz,4)
         m2blok = size(bxyz,5)/j2blok
         call PAVPOTD332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok&
     &,m2blok)
         end subroutine ipavpotd332
!
      end module pdfield32d
