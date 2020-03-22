!-----------------------------------------------------------------------
!
      module pcfield32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library 
! pcfield32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: october 7, 2005
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: poisc3_init, poisc
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PFORMC32(ffg,f,fs,ft,bs,br,g,h,fpotc,mixup3,sct3,aff&
     &p,ar,indx1,indy1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,k&
     &yp2d,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,kxyp2d,kyzp&
     &2d,nz1d,nxhyz2,nxyzh2)
         implicit none
         integer :: indx1, indy1, indz1, kstrt, nxv, ny2d, nz2d
         integer :: kxyp2, kyp2, kyzp2, kzp2, kyp2d, kzyp2
         integer :: j2blok, k2blok, jk2blok, l2blok, m2blok, ml2blok
         integer :: kxyp2d, kyzp2d, nz1d, nxhyz2, nxyzh2
         real :: ar, affp
         real, dimension(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok) :: ffg
         real, dimension(2*nxv,kyp2d,kzp2,k2blok*l2blok) :: f
         complex, dimension(ny2d,kxyp2,kzp2,j2blok*l2blok) :: fs
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: ft
         real, dimension(5,nz1d,kyzp2,j2blok*m2blok) :: g, h
         complex, dimension(kxyp2,kzyp2,kzp2,jk2blok*l2blok) :: bs
         complex, dimension(kxyp2,kzyp2,kzp2,j2blok*ml2blok) :: br
         integer, dimension(nxhyz2) :: mixup3
         complex, dimension(nxyzh2) :: sct3
         real, external :: fpotc
         end subroutine
      end interface
      interface
         subroutine PFORMC32X(ffg,f,fs,ft,g,h,fpotc,mixup3,sct3,affp,ar,&
     &indx1,indy1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kyp2d,&
     &j2blok,k2blok,l2blok,m2blok,kxyp2d,kyzp2d,nz1d,nxhyz2,nxyzh2)
         implicit none
         integer :: indx1, indy1, indz1, kstrt, nxv, ny2d, nz2d
         integer :: kxyp2, kyp2, kyzp2, kzp2, kyp2d, j2blok, k2blok
         integer :: l2blok, m2blok, kxyp2d, kyzp2d, nz1d, nxhyz2, nxyzh2
         real :: ar, affp
         real, dimension(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok) :: ffg
         real, dimension(2*nxv,kyp2d,kzp2,k2blok*l2blok) :: f
         complex, dimension(ny2d,kxyp2,kzp2,j2blok*l2blok) :: fs
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: ft
         real, dimension(5,nz1d,kyzp2,j2blok*m2blok) :: g, h
         integer, dimension(nxhyz2) :: mixup3
         complex, dimension(nxyzh2) :: sct3
         real, external :: fpotc
         end subroutine
      end interface
      interface
         subroutine PPOISC32(q,fx,fy,fz,isign,ffg,we,nx,ny,nz,kstrt,nz2d&
     &,kxyp2,kyzp2,j2blok,m2blok,kxyp2d,kyzp2d,nz1d)
         implicit none
         real :: we
         integer :: isign, nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, kxyp2d, kyzp2d, nz1d
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: q
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: fx, fy,fz
         real, dimension(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok) :: ffg
         end subroutine
      end interface
      interface
         subroutine PPOISC332(q,fxyz,ffg,we,nx,ny,nz,kstrt,nz2d,kxyp2,ky&
     &zp2,j2blok,m2blok,kxyp2d,kyzp2d,nz1d)
         implicit none
         real :: we
         integer :: nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
         integer :: j2blok, m2blok, kxyp2d, kyzp2d, nz1d
         complex, dimension(nz2d,kxyp2,kyzp2,j2blok*m2blok) :: q
         complex, dimension(3,nz2d,kxyp2,kyzp2,j2blok*m2blok) :: fxyz
         real, dimension(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok) :: ffg
         end subroutine
      end interface
      interface
         function POTC3(r,affp,ari,ifun)
         implicit none
         integer :: ifun
         real :: POTC3, r, affp, ari
         end function
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface poisc3_init
!        module procedure ippoisc32init
         module procedure ippoisc32xinit
      end interface
!
      interface poisc
         module procedure ippoisc32
         module procedure ipspoisc32
         module procedure ippoisc332
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ippoisc32(q,fx,ffg,we,nx,ny,nz,kstrt,j2blok)
! poisson solver for 3d potential, open boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q, fx
         real, dimension(:,:,:,:,:), pointer :: ffg
! local data
         integer :: isign = 1, nz2d, kxyp2, kyzp2, m2blok, nz1d
         integer :: kxyp2d, kyzp2d
         complex, dimension(1,1,1,1) :: fy, fz
         nz2d = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok
         nz1d = size(ffg,2); kxyp2d = size(ffg,3); kyzp2d = size(ffg,4)
         call PPOISC32(q,fx,fy,fz,isign,ffg,we,nx,ny,nz,kstrt,nz2d,kxyp2&
     &,kyzp2,j2blok,m2blok,kxyp2d,kyzp2d,nz1d)
         end subroutine ippoisc32
!
         subroutine ipspoisc32(q,fy,ffg,nx,ny,nz,kstrt,j2blok)
! smoother for 3d scalar field, open boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         complex, dimension(:,:,:,:), pointer :: q, fy
         real, dimension(:,:,:,:,:), pointer :: ffg
! local data
         integer :: isign = 2, nz2d, kxyp2, kyzp2, m2blok, nz1d
         integer :: kxyp2d, kyzp2d
         real :: we
         complex, dimension(1,1,1,1) :: fx, fz
         nz2d = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok
         nz1d = size(ffg,2); kxyp2d = size(ffg,3); kyzp2d = size(ffg,4)
         call PPOISC32(q,fx,fy,fz,isign,ffg,we,nx,ny,nz,kstrt,nz2d,kxyp2&
     &,kyzp2,j2blok,m2blok,kxyp2d,kyzp2d,nz1d)
         end subroutine ipspoisc32
!
         subroutine ippoisc32init(ffg,f,fs,ft,mixup3,sct3,ar,affp,indx,i&
     &ndy,indz,kstrt,kyp2,j2blok,k2blok)
! initialize 3d poisson solver with 3d fields, open boundary conditions
         implicit none
         integer :: indx, indy, indz, kstrt, kyp2, j2blok, k2blok
         real :: ar, affp
         real, dimension(:,:,:,:,:), pointer :: ffg
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: fs
         complex, dimension(:,:,:,:), pointer :: ft
         integer, dimension(:), pointer :: mixup3
         complex, dimension(:), pointer :: sct3
! local data
         integer :: indx1, indy1, indz1, nxv, ny2d, nz2d, kxyp2, kyzp2
         integer :: kzp2, kyp2d, kzyp2, jk2blok, l2blok, m2blok, ml2blok
         integer :: kxyp2d, kyzp2d, nz1d, nxhyz2, nxyzh2
         complex, dimension(size(ft,2),max(size(ft,3),kyp2),size(f,3),ma&
     &x(size(f,4),size(fs,4))) :: bs
         complex, dimension(size(ft,2),max(size(ft,3),kyp2),size(f,3),ma&
     &x(size(fs,4),size(ft,4))) :: br
         real, dimension(5,size(ffg,2),size(ft,3),size(ft,4)) :: g, h
         real, external :: POTC3
         indx1 = indx + 1; indy1 = indy + 1; indz1 = indz + 1
         nz1d = size(ffg,2); kxyp2d = size(ffg,3); kyzp2d = size(ffg,4)
         nxv = size(f,1)/2; kyp2d = size(f,2); kzp2 = size(f,3)
         l2blok = size(f,4)/k2blok
         nz2d = size(ft,1); kxyp2 = size(ft,2); kyzp2 = size(ft,3)
         m2blok = size(ft,4)/j2blok; ny2d = size(fs,1)
         kzyp2 = max(kyzp2,kyp2)
         jk2blok = max(j2blok,k2blok); ml2blok = max(m2blok,l2blok)
         nxhyz2 = size(mixup3); nxyzh2 = size(sct3)
         call PFORMC32(ffg,f,fs,ft,bs,br,g,h,POTC3,mixup3,sct3,affp,ar,i&
     &ndx1,indy1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kyp2d,k&
     &zyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,kxyp2d,kyzp2d,nz1&
     &d,nxhyz2,nxyzh2)
         end subroutine ippoisc32init
!
         subroutine ippoisc32xinit(ffg,f,fs,ft,mixup3,sct3,ar,affp,indx,&
     &indy,indz,kstrt,kyp2,j2blok,k2blok)
! initialize 3d poisson solver with 3d fields, open boundary conditions
         implicit none
         integer :: indx, indy, indz, kstrt, kyp2, j2blok, k2blok
         real :: ar, affp
         real, dimension(:,:,:,:,:), pointer :: ffg
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: fs
         complex, dimension(:,:,:,:), pointer :: ft
         integer, dimension(:), pointer :: mixup3
         complex, dimension(:), pointer :: sct3
! local data
         integer :: indx1, indy1, indz1, nxv, ny2d, nz2d, kxyp2, kyzp2
         integer :: kzp2, kyp2d, l2blok, m2blok, kxyp2d, kyzp2d, nz1d
         integer :: nxhyz2, nxyzh2
         real, dimension(5,size(ffg,2),size(ft,3),size(ft,4)) :: g, h
         real, external :: POTC3
         indx1 = indx + 1; indy1 = indy + 1; indz1 = indz + 1
         nz1d = size(ffg,2); kxyp2d = size(ffg,3); kyzp2d = size(ffg,4)
         nxv = size(f,1)/2; kyp2d = size(f,2); kzp2 = size(f,3)
         l2blok = size(f,4)/k2blok
         nz2d = size(ft,1); kxyp2 = size(ft,2); kyzp2 = size(ft,3)
         m2blok = size(ft,4)/j2blok; ny2d = size(fs,1)
         nxhyz2 = size(mixup3); nxyzh2 = size(sct3)
         call PFORMC32X(ffg,f,fs,ft,g,h,POTC3,mixup3,sct3,affp,ar,indx1,&
     &indy1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kyp2d,j2blok&
     &,k2blok,l2blok,m2blok,kxyp2d,kyzp2d,nz1d,nxhyz2,nxyzh2)
         end subroutine ippoisc32xinit
!
         subroutine ippoisc332(q,fxyz,ffg,we,nx,ny,nz,kstrt,j2blok)
! poisson solver for 3d electric field, open boundaries
         implicit none
         integer :: nx, ny, nz, kstrt, j2blok
         real :: we
         complex, dimension(:,:,:,:), pointer :: q
         complex, dimension(:,:,:,:,:), pointer :: fxyz
         real, dimension(:,:,:,:,:), pointer :: ffg
! local data
         integer :: nz2d, kxyp2, kyzp2, m2blok, nz1d, kxyp2d, kyzp2d
         nz2d = size(q,1); kxyp2 = size(q,2); kyzp2 = size(q,3)
         m2blok = size(q,4)/j2blok
         nz1d = size(ffg,2); kxyp2d = size(ffg,3); kyzp2d = size(ffg,4)
         call PPOISC332(q,fxyz,ffg,we,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2&
     &blok,m2blok,kxyp2d,kyzp2d,nz1d)
         end subroutine ippoisc332
!
      end module pcfield32d
