!-----------------------------------------------------------------------
!
      module pfield32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library pfield32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: april 2, 2008
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PSCGUARD32, PSGUARD32, PSCGUARD32L, PSGUARD32L
      public :: PACGUARD32X, PAGUARD32X, PACGUARD32XL, PAGUARD32XL
      public :: cguard, sguard, aguard, zguard
      public :: pois_init, pois, cuperp, bpois
      public :: ibpois, maxwel, emfield, avpot, gtmodes, ptmodes
      public :: addqei, imoment, ipdivf32, ipgradf32, ipcurlf32
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD32(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,id&
     &ds,mnblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PSGUARD32(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         real :: qi0
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PACGUARD32X(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PAGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PZCGUARD32(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PZGUARD32(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PCGUARD32XL(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD32L(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,i&
     &dds,mnblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PSGUARD32L(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
         implicit none
         real :: qi0
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PACGUARD32XL(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PZCGUARD32L(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(3,nxe,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PZGUARD32L(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real, dimension(nxe,nypmx,nzpmx,mnblok) :: q
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PPOIS32INIT(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,ky&
     &zp,jblok,mblok,nzhd)
         implicit none
         integer :: nx, ny, nz, kstrt, kxyp, kyzp, jblok, mblok, nzhd
         real :: ax, ay, az, affp
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PPOISP32(q,fx,fy,fz,isign,ffc,ax,ay,az,affp,we,nx,ny&
     &,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp, kyzp
         integer :: jblok, mblok, nzhd
         complex, dimension(nzv,kxyp,kyzp,jblok*mblok) :: q, fx, fy, fz
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PPOISP332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz&
     &,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         implicit none
         real :: ax, ay, az, affp, we
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp, kyzp
         integer :: jblok, mblok, nzhd
         complex, dimension(nzv,kxyp,kyzp,jblok*mblok) :: q
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: fxyz
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PDIVF32(f,df,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblo&
     &k)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: f
         complex, dimension(nzv,kxyp,kyzp,jblok*mblok) :: df
         end subroutine
      end interface
      interface
         subroutine PGRADF32(df,f,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mbl&
     &ok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         complex, dimension(nzv,kxyp,kyzp,jblok*mblok) :: df
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: f
         end subroutine
      end interface
      interface
         subroutine PCURLF32(f,g,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblo&
     &k)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PCUPERP32(cu,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblo&
     &k)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISP332(cu,bxyz,isign,ffc,ax,ay,az,affp,ci,wm,nx,&
     &ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         implicit none
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp, kyzp
         integer :: jblok, mblok, nzhd
         real :: ax, ay, az, affp, ci, wm
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: cu, bxyz
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine IPBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nzv,kxy&
     &p,kyzp,jblok,mblok,nzhd)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         integer :: nzhd
         real :: ci, wm
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: cu, bxyz
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PMAXWEL32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,nz&
     &,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         integer :: nzhd
         real :: affp, ci, dt, wf, wm
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: exyz, bxyz
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: cu
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PEMFIELD32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,nzv,kx&
     &yp,kyzp,jblok,mblok,nzhd)
         implicit none
         integer :: isign, nx, ny, nz, kstrt, nzv, kxyp, kyzp
         integer :: jblok, mblok, nzhd
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: fxyz, exyz
         complex, dimension(nzhd,kxyp,kyzp,jblok*mblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jbl&
     &ok,mblok)
         implicit none
         integer :: nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
         complex, dimension(3,nzv,kxyp,kyzp,jblok*mblok) :: bxyz, axyz
         end subroutine
      end interface
      interface
         subroutine PGTMODES32(pot,pott,nx,ny,nz,it,kstrt,modesx,modesy,&
     &modesz,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd)
         implicit none
         integer :: nx, ny, nz, it, kstrt, modesx, modesy, modesz
         integer :: nzv, kxyp, kyzp, jblok, mblok, nt
         integer :: modesxpd, modesypd, modeszd
         complex, dimension(nzv,kxyp,kyzp,jblok*mblok) :: pot
         complex, dimension(nt,modeszd,modesxpd,modesypd,jblok*mblok) ::&
     & pott
         end subroutine
      end interface
      interface
         subroutine PPTMODES32(pot,pott,nx,ny,nz,it,kstrt,modesx,modesy,&
     &modesz,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd)
         implicit none
         integer :: nx, ny, nz, it, kstrt, modesx, modesy, modesz
         integer :: nzv, kxyp, kyzp, jblok, mblok, nt
         integer :: modesxpd, modesypd, modeszd
         complex, dimension(nzv,kxyp,kyzp,jblok*mblok) :: pot
         complex, dimension(nt,modeszd,modesxpd,modesypd,jblok*mblok) ::&
     & pott
         end subroutine
      end interface
      interface
         subroutine PADDQEI32(qe,qi,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
!        real, dimension(*) :: qe, qi
         real :: qe, qi
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
      interface
         subroutine PIMOMENT32(qi,fxyz,nyzp,pxi,pyi,pzi,dt,nx,nxe,nypmx,&
     &nzpmx,idds,mnblok)
         implicit none
         integer :: nx, nxe, nypmx, nzpmx, idds, mnblok
         real :: pxi, pyi, pzi, dt
!        real, dimension(*) :: qi
!        real, dimension(*) :: fxyz
         real :: qi, fxyz
         integer, dimension(idds,mnblok) :: nyzp
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface cguard
         module procedure ipcguard32x
         module procedure ipdguard32x
      end interface
!
      interface sguard
         module procedure ipscguard32
         module procedure ipsguard32
      end interface
!      
      interface aguard
         module procedure ipacguard32x
         module procedure ipaguard32x
      end interface
!
      interface zguard
         module procedure ipzguard32
      end interface
!
       interface pois_init
!        module procedure ippois332init
         module procedure ippois32init
      end interface
! 
      interface pois
         module procedure ippois32
         module procedure ipspois32
         module procedure ippois332
      end interface
!
      interface cuperp
         module procedure ipcuperp32
      end interface
!
      interface bpois
         module procedure jpbpois332
      end interface
!
      interface ibpois
         module procedure jipbpois332
      end interface
!  
      interface maxwel
         module procedure ipmaxwel32
      end interface
!
      interface emfield
         module procedure ipemfield32
      end interface
!
      interface avpot
         module procedure ipavpot332
      end interface
!
      interface gtmodes
         module procedure ipgtmodes32
      end interface
!
      interface ptmodes
         module procedure ipptmodes32
      end interface
!
      interface addqei
         module procedure ipaddqei32
      end interface
!
      interface imoment
         module procedure ipimoment32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipcguard32x(fxyz,nyzp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 3d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: fxyz
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         mnblok = size(fxyz,5); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PCGUARD32XL(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipcguard32x
!
         subroutine ipdguard32x(q,nyzp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 3d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipdguard32x
!
         subroutine ipscguard32(cu,nyzp,xj0,yj0,zj0,nx,inorder)
! initialize periodic 3d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:,:,:), pointer :: cu
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         mnblok = size(cu,5); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PSCGUARD32L(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,idds&
     &,mnblok)
         else
            call PSCGUARD32(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,idds,&
     &mnblok)
         endif
         end subroutine ipscguard32
!
         subroutine ipsguard32(q,nyzp,qi0,nx,inorder)
! initialize periodic 3d scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PSGUARD32L(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PSGUARD32(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipsguard32
!
         subroutine ipacguard32x(cu,nyzp,nx,inorder)
! add guard cells in x for non-uniform, periodic 3d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: cu
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         mnblok = size(cu,5); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PACGUARD32XL(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PACGUARD32X(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipacguard32x
!
         subroutine ipaguard32x(q,nyzp,nx,inorder)
! add guard cells in x for non-uniform, periodic 3d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PAGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipaguard32x
!
         subroutine ipzguard32(q,nyzp,nx,inorder)
! zero out guard cells in 3d scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PZGUARD32L(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PZGUARD32(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
         end subroutine ipzguard32
!
!-----------------------------------------------------------------------
         subroutine ippois32init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,jblok)
! initialize 3d periodic electric field solver
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: ax, ay, az, affp
         complex, dimension(:,:,:,:), pointer :: ffc
! local data
         integer :: kxyp, kyzp, mblok, nzhd
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3)
         mblok = size(ffc,4)/jblok
         call PPOIS32INIT(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,kyzp,jbl&
     &ok,mblok,nzhd)
         end subroutine ippois32init
!
         subroutine ippois32(q,fx,ffc,we,nx,ny,nz,kstrt,jblok)
! poisson solver for periodic 3d potential
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, fx, ffc
! local data
         integer :: isign = 1, nzv, kxyp, kyzp, mblok, nzhd
         real :: ax, ay, az, affp
         complex, dimension(1,1,1,1) :: fy, fz
         nzv = size(q,1)
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3)
         mblok = size(ffc,4)/jblok
         call PPOISP32(q,fx,fy,fz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,ks&
     &trt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         end subroutine ippois32
!
         subroutine ipspois32(q,fy,ffc,nx,ny,nz,kstrt,jblok)
! smoother for periodic 3d scalar field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:), pointer :: q, fy, ffc
! local data
         integer :: isign = 2, nzv, kxyp, kyzp, mblok, nzhd
         real :: ax, ay, az, affp, we
         complex, dimension(1,1,1,1) :: fx, fz
         nzv = size(q,1)
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3)
         mblok = size(ffc,4)/jblok
         call PPOISP32(q,fx,fy,fz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,ks&
     &trt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         end subroutine ipspois32
!
         subroutine ippois332init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,jblok&
     &)
! initialize 3d periodic electric field solver
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: ax, ay, az, affp
         complex, dimension(:,:,:,:), pointer :: ffc
! local data
         integer :: isign = 0, nzv, kxyp, kyzp, mblok, nzhd
         real :: we
         complex, dimension(1,1,1,1) :: q
         complex, dimension(3,1,1,1,1) :: fxyz
         nzv = size(q,1)
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3)
         mblok = size(ffc,4)/jblok
         call PPOISP332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt&
     &,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         end subroutine ippois332init
!
         subroutine ippois332(q,fxyz,ffc,we,nx,ny,nz,kstrt,jblok)
! poisson solver for periodic 3d electric field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, ffc
         complex, dimension(:,:,:,:,:), pointer :: fxyz
! local data
         integer :: isign = -1, nzv, kxyp, kyzp, mblok, nzhd
         real :: ax, ay, az, affp
         nzv = size(q,1)
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3)
         mblok = size(ffc,4)/jblok
         call PPOISP332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt&
     &,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         end subroutine ippois332
!
         subroutine ipdivf32(f,df,nx,ny,nz,kstrt,jblok)
! calculates the divergence of periodic 3d vector field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: df
! local data
         integer :: nzv, kxyp, kyzp, mblok
         nzv = size(f,2); kxyp = size(f,3); kyzp = size(f,4)
         mblok = size(f,5)/jblok
         call PDIVF32(f,df,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
         end subroutine ipdivf32
!
         subroutine ipgradf32(df,f,nx,ny,nz,kstrt,jblok)
! calculates the gradient of periodic 3d scalar field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:), pointer :: df
         complex, dimension(:,:,:,:,:), pointer :: f
! local data
         integer :: nzv, kxyp, kyzp, mblok
         nzv = size(f,2); kxyp = size(f,3); kyzp = size(f,4)
         mblok = size(f,5)/jblok
         call PGRADF32(df,f,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
         end subroutine ipgradf32
!
         subroutine ipcurlf32(f,g,nx,ny,nz,kstrt,jblok)
! calculates the curl of periodic 3d vector field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:,:), pointer :: f, g
! local data
         integer :: nzv, kxyp, kyzp, mblok
         nzv = size(f,2); kxyp = size(f,3); kyzp = size(f,4)
         mblok = size(f,5)/jblok
         call PCURLF32(f,g,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
         end subroutine ipcurlf32
!
         subroutine ipcuperp32(cu,nx,ny,nz,kstrt,jblok)
! calculates the transverse part of periodic 3d vector field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:,:), pointer :: cu
! local data
         integer :: nzv, kxyp, kyzp, mblok
         nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
         mblok = size(cu,5)/jblok
         call PCUPERP32(cu,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
         end subroutine ipcuperp32
!
         subroutine jpbpois332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,jblok)
! calculates static vector potential for periodic 3d vector field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: ffc
         complex, dimension(:,:,:,:,:), pointer :: cu, bxyz
! local data
         integer :: isign = 1, nzv, kxyp, kyzp, nzhd, mblok
         real :: ax, ay, az, affp
         nzv = size(cu,2)
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3) 
         mblok = size(ffc,4)/jblok
         call PBPOISP332(cu,bxyz,isign,ffc,ax,ay,az,affp,ci,wm,nx,ny,nz,&
     &kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         end subroutine jpbpois332
!
         subroutine jipbpois332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,jblok)
! calculates static magnetic for periodic 3d vector field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: ffc
         complex, dimension(:,:,:,:,:), pointer :: cu, bxyz
! local data
         integer :: nzv, kxyp, kyzp, nzhd, mblok
         nzv = size(cu,2)
         nzhd = size(ffc,1); kxyp = size(ffc,2); kyzp = size(ffc,3) 
         mblok = size(ffc,4)/jblok
         call IPBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nzv,kxyp,kyzp&
     &,jblok,mblok,nzhd)
         end subroutine jipbpois332
!
         subroutine ipmaxwel32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,n&
     &z,kstrt,jblok)
! calculates maxwell's equation for periodic 3d vector field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         real :: affp, ci, dt, wf, wm
         complex, dimension(:,:,:,:,:), pointer :: exyz, bxyz, cu
         complex, dimension(:,:,:,:), pointer :: ffc
! local data
         integer :: nzv, kxyp, kyzp, mblok, nzhd
         nzv = size(cu,2); kxyp = size(cu,3); kyzp = size(cu,4)
         nzhd = size(ffc,1); mblok = size(ffc,4)/jblok
         call PMAXWEL32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,nz,kstrt&
     &,nzv,kxyp,kyzp,jblok,mblok,nzhd)
         end subroutine ipmaxwel32
!
         subroutine ipemfield32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,jblok&
     &)
! combines and smooths periodic 3d vector fields
         implicit none
         integer :: isign, nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:,:), pointer :: fxyz, exyz
         complex, dimension(:,:,:,:), pointer :: ffc
! local data
         integer :: nzv, kxyp, kyzp, mblok, nzhd
         nzv = size(fxyz,2); kxyp = size(fxyz,3); kyzp = size(fxyz,4)
         mblok = size(fxyz,5)/jblok; nzhd = size(ffc,1)
         call PEMFIELD32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyz&
     &p,jblok,mblok,nzhd)
         end subroutine ipemfield32
!
         subroutine ipavpot332(bxyz,axyz,nx,ny,nz,kstrt,jblok)
! calculates periodic 3d vector potential from magnetic field
         implicit none
         integer :: nx, ny, nz, kstrt, jblok
         complex, dimension(:,:,:,:,:), pointer :: bxyz, axyz
! local data
         integer :: nzv, kxyp, kyzp, mblok
         nzv = size(bxyz,2); kxyp = size(bxyz,3); kyzp = size(bxyz,4) 
         mblok = size(bxyz,5)/jblok
         call PAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mbl&
     &ok)
         end subroutine ipavpot332
!
         subroutine ipgtmodes32(pot,pott,nx,ny,nz,modesx,modesy,modesz,k&
     &strt,jblok)
! extracts lowest order modes from 3d scalar field
         implicit none
         integer :: nx, ny, nz, modesx, modesy, modesz, kstrt, jblok
         complex, dimension(:,:,:,:), pointer :: pot
         complex, dimension(:,:,:,:), pointer :: pott
! local data
         integer :: it, nzv, kxyp, kyzp, mblok, nt
         integer :: modesxpd, modesypd, modeszd
         it = 1; nt = 1
         nzv = size(pot,1);  kxyp = size(pot,2); kyzp = size(pot,3)
         mblok = size(pot,4)/jblok
         modeszd = size(pott,1); modesxpd = size(pott,2)
         modesypd = size(pott,3)
         call PGTMODES32(pot,pott,nx,ny,nz,it,kstrt,modesx,modesy,modesz&
     &,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd)
         end subroutine ipgtmodes32
!
         subroutine ipptmodes32(pot,pott,nx,ny,nz,modesx,modesy,modesz,k&
     &strt,jblok)
! extracts lowest order modes from 3d scalar field
         implicit none
         integer :: nx, ny, nz, modesx, modesy, modesz, kstrt, jblok
         complex, dimension(:,:,:,:), pointer :: pot
         complex, dimension(:,:,:,:), pointer :: pott
! local data
         integer :: it, nzv, kxyp, kyzp, mblok, nt
         integer :: modesxpd, modesypd, modeszd
         it = 1; nt = 1
         nzv = size(pot,1);  kxyp = size(pot,2); kyzp = size(pot,3)
         mblok = size(pot,4)/jblok
         modeszd = size(pott,1); modesxpd = size(pott,2)
         modesypd = size(pott,3)
         call PPTMODES32(pot,pott,nx,ny,nz,it,kstrt,modesx,modesy,modesz&
     &,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd)
         end subroutine ipptmodes32
!
         subroutine ipaddqei32(qe,qi,nyzp,nx,inorder)
! adds electron and ion densities
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: qe, qi
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         nxe = size(qe,1); nypmx = size(qe,2); nzpmx = size(qe,3)
         mnblok = size(qe,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PADDQEI32(qe(1,1,1,1),qi(1,1,1,1),nyzp,nx,nxe,nypmx,nzp&
     &mx,idds,mnblok)
         else
            call PADDQEI32(qe(2,2,2,1),qi(2,2,2,1),nyzp,nx,nxe,nypmx,nzp&
     &mx,idds,mnblok)
         endif
         end subroutine ipaddqei32
!
         subroutine ipimoment32(qi,fxyz,nyzp,id0,iunit,px,py,pz,dt,wx,wy&
     &,wz,nx,inorder)
! calculate ion momentum from integral of qi*fxyz,
! and prints it, and adds it to total momentum, for 3d code
         implicit none
         integer :: nx, id0, iunit
         real :: px, py, pz, dt, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: qi
         real, dimension(:,:,:,:,:), pointer :: fxyz
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, idds, mnblok, order
         double precision, dimension(3) :: sum3, work3
  995    format (' ion momentum = ',3e14.7)
         nxe = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         mnblok = size(fxyz,5); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! calculate and print ion momentum
         if (order==LINEAR) then
            call PIMOMENT32(qi(1,1,1,1),fxyz(1,1,1,1,1),nyzp,px,py,pz,dt&
     &,nx,nxe,nypmx,nzpmx,idds,mnblok)
         else
            call PIMOMENT32(qi(2,2,2,1),fxyz(1,2,2,2,1),nyzp,px,py,pz,dt&
     &,nx,nxe,nypmx,nzpmx,idds,mnblok)
         endif
! sum over the y and z directions
         sum3(1) = px
         sum3(2) = py
         sum3(3) = pz
         call PDSUM(sum3,work3,3,1)
         px = sum3(1)
         py = sum3(2)
         pz = sum3(3)
         if (id0==0) write (iunit,995) px, py, pz
! add to total momentum
         wx = wx + px
         wy = wy + py
         wz = wz + pz
         end subroutine ipimoment32
!
      end module pfield32d

