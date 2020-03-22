!-----------------------------------------------------------------------
!
      module pnpfield32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library
! pfield32lib.f, pdfield32lib.f, pbfield32lib.f, pcfield32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: october 27, 2005
!
      use globals, only: LINEAR, QUADRATIC
      use pfield32d
      use pdfield32d
      use pbfield32d
      use pcfield32d
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: cguard, sguard, aguard, zguard
      public :: pois_init, pois, cuperp, bpois
      public :: ibpois, maxwel, emfield, avpot, gtmodes, ptmodes
      public :: ipdivf32, ipgradf32, ipcurlf32
      public :: poisd_init, poisd, cuperpd, bpoisd, ibpoisd
      public :: maxweld, cmfieldd, emfieldd, cpfieldd, avpotd
      public :: poism_init, poism, cuperpm, bpoism, ibpoism
      public :: maxwelm, cmfieldm, emfieldm, cpfieldm, avpotm
      public :: sguardp, aguardp, cguardp
      public :: poisc3_init, poisc
!
! define generic interfaces to Fortran90 library
!
      interface sguardp
         module procedure ipscguard32xp
         module procedure ipsguard32xp
      end interface
!      
      interface aguardp
         module procedure ipacguard32xp
         module procedure ipaguard32xp
      end interface
!
      interface cguardp
         module procedure ipcguard32xp
         module procedure ipdguard32xp
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipscguard32xp(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,z&
     &j0,nx,ny,nz,mblok,ipbc,inorder)
! initialize non-uniform 3d vector field
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, ny, nz, mblok, ipbc
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:,:,:), pointer :: cu
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         integer :: ngx = 1, ngy = 1, ngz = 1, nxe, nypmx, nzpmx, idds
         integer :: nblok, mnblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         mnblok = size(cu,5); nblok = mnblok/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PSCGUARD32L(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,i&
     &dds,mnblok)
            else
               call PSCGUARD32(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,id&
     &ds,mnblok)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PLSCGUARD32XL(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,z&
     &j0,nx,ny,nz,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
            else
               call PLSCGUARD32X(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,zj&
     &0,nx,ny,nz,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PMSCGUARD32XL(cu,kstrt,nvpy,noff,nyzp,xj0,yj0,zj0,nx&
     &,ny,ngx,ngy,nxe,nypmx,nzpmx,idds,mblok,nblok)
            else
               call PMSCGUARD32X(cu,kstrt,nvpy,noff,nyzp,xj0,yj0,zj0,nx,&
     &ny,ngx,ngy,nxe,nypmx,nzpmx,idds,mblok,nblok)
            endif
         endif
         end subroutine ipscguard32xp
!
         subroutine ipsguard32xp(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,n&
     &z,mblok,ipbc,inorder)
! initialize non-uniform 3d scalar field
         implicit none
         integer :: kstrt, nvpy, nvpz, nx, ny, nz, mblok, ipbc
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         integer :: ngx = 1, ngy = 1, ngz = 1, nxe, nypmx, nzpmx, idds
         integer :: nblok, mnblok, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); nblok = mnblok/mblok; idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PSGUARD32L(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
            else
               call PSGUARD32(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PLSGUARD32XL(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,n&
     &z,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
            else
               call PLSGUARD32X(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,nz&
     &,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PMSGUARD32XL(q,kstrt,nvpy,noff,nyzp,qi0,nx,ny,ngx,ng&
     &y,nxe,nypmx,nzpmx,idds,mblok,nblok)
            else
               call PMSGUARD32X(q,kstrt,nvpy,noff,nyzp,qi0,nx,ny,ngx,ngy&
     &,nxe,nypmx,nzpmx,idds,mblok,nblok)
            endif
         endif
         end subroutine ipsguard32xp
!
         subroutine ipacguard32xp(cu,nyzp,nx,ipbc,inorder)
! add guard cells in x,for non-uniform 3d vector data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: cu
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         mnblok = size(cu,5); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PACGUARD32XL(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            else
               call PACGUARD32X(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call PLACGUARD32X(cu(1,2,1,1,1),nyzp,nx-2,nxe,nypmx,nzpmx&
     &,idds,mnblok)
            endif
         else if (ipbc==3) then
            if (order==QUADRATIC) then
               call PLACGUARD32X(cu(1,2,1,1,1),nyzp,nx-2,nxe,nypmx,nzpmx&
     &,idds,mnblok)
            endif
         endif
         end subroutine ipacguard32xp
!
         subroutine ipaguard32xp(q,nyzp,nx,ipbc,inorder)
! add guard cells in x for non-uniform 3d scalar data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4); idds = size(nyzp,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            else
               call PAGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call PLAGUARD32X(q(2,1,1,1),nyzp,nx-2,nxe,nypmx,nzpmx,idd&
     &s,mnblok)
            endif
         else if (ipbc==3) then
            if (order==QUADRATIC) then
               call PLAGUARD32X(q(2,1,1,1),nyzp,nx-2,nxe,nypmx,nzpmx,idd&
     &s,mnblok)
            endif
         endif
         end subroutine ipaguard32xp
!
         subroutine ipcguard32xp(fxyz,nyzp,nx,ipbc,inorder)
! copy guard cells in x for non-uniform 3d vector data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:,:), pointer :: fxyz
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         mnblok = size(fxyz,5)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PCGUARD32XL(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
            else
               call PCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call PLCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
            endif
         else if (ipbc==3) then
            if (order==QUADRATIC) then
               call PLCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok&
     &)
            endif
         endif
         end subroutine ipcguard32xp
!
         subroutine ipdguard32xp(q,nyzp,nx,ipbc,inorder)
! copy guard cells in x for non-uniform 3d scalar data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: q
         integer, dimension(:,:), pointer :: nyzp
! local data
         integer :: nxe, nypmx, nzpmx, mnblok, idds, order
         nxe = size(q,1); nypmx = size(q,2); nzpmx = size(q,3)
         mnblok = size(q,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            else
               call PDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call PLDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         else if (ipbc==3) then
            if (order==QUADRATIC) then
               call PLDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
            endif
         endif
         end subroutine ipdguard32xp
!
      end module pnpfield32d
