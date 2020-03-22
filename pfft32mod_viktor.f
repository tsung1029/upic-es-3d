!-----------------------------------------------------------------------
!
      module pfft32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library pfft32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: june 16, 2008
!
      use globals, only: LINEAR, QUADRATIC
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: wtimer
      public :: fft_init, fft, fftn, fftc_init
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy&
     &,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,k&
     &zyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd
         real :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,ind&
     &y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,&
     &kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd
         real :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(3,kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(3,kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT32RX(f,g,h,isign,ntpose,mixup,sct,indx,indy,indz&
     &,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblok,&
     &kblok,lblok,mblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok
         integer :: nxhyzd, nxyzhd
         real :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT32RX3(f,g,h,isign,ntpose,mixup,sct,indx,indy,ind&
     &z,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblok&
     &,kblok,lblok,mblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok
         integer :: nxhyzd, nxyzhd
         real :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT32C(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy&
     &,indz,kstrt,nxv,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kz&
     &yp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxv, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxyzd, nxyzhd
         complex :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer :: indx, indy, indz
         integer :: nxhyzd, nxyzhd
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd
         real :: ttp
         real :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT32RX(f,g,h,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,j&
     &blok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok
         integer :: nxhyzd, nxyzhd
         real :: ttp
         real :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd
         real :: ttp
         real :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(3,kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(3,kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT32RX3(f,g,h,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,&
     &jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok
         integer :: nxhyzd, nxyzhd
         real :: ttp
         real :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT32RN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd,nxy&
     &zhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: ndim, nxhyzd, nxyzhd
         real :: ttp
         real :: f
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(ndim,kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(ndim,kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         complex, dimension(ndim,nxvh) :: ss
         end subroutine
      end interface
      interface
         subroutine WPFFT32RXN(f,g,h,ss,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kz&
     &pd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok
         integer :: ndim, nxhyzd, nxyzhd
         real :: ttp
         real :: f
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         complex, dimension(ndim,nxvh) :: ss
         end subroutine
      end interface
      interface
         subroutine WP2FFT32RN(f1,f2,g1,g2,h1,h2,bs,br,ss,isign,ntpose,m&
     &ixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,k&
     &xypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nd&
     &im1,ndim2,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: ndim1, ndim2, nxhyzd, nxyzhd
         real :: ttp
!        real, dimension(*) :: f1, f2
         real :: f1, f2
         complex, dimension(ndim1,nyv,kxypd,kzpd,jblok*lblok) :: g1
         complex, dimension(ndim2,nyv,kxypd,kzpd,jblok*lblok) :: g2
         complex, dimension(ndim1,nzv,kxypd,kyzpd,jblok*mblok) :: h1
         complex, dimension(ndim2,nzv,kxypd,kyzpd,jblok*mblok) :: h2
         complex, dimension(ndim1+ndim2,kxyp,kzyp,kzp,jkblok*lblok) ::  &
     &bs
         complex, dimension(ndim1+ndim2,kxyp,kzyp,kzp,jblok*mlblok) ::  &
     &br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         complex, dimension(ndim1+ndim2,nxvh) :: ss
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft_init
!        module procedure ipfft32rinit
!        module procedure ipfft32rxinit
         module procedure iwpfft32rinit
      end interface
!
      interface fft
!        module procedure ipfft32r
!        module procedure ipfft32r3
!        module procedure ipfft32rx
!        module procedure ipfft32rx3
         module procedure iwpfft32r
         module procedure iwpfft32r3
!        module procedure iwpfft32rx
!        module procedure iwpfft32rx3
         module procedure ipfft32c
      end interface
!
      interface fftn
         module procedure iwpfft32rn
!        module procedure iwpfft32rxn
         module procedure iwp1fft32rn
         module procedure iwp2fft32rn
      end interface
!
      interface fftc_init
         module procedure ipfft32cinit
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipfft32rinit(mixup,sct,indx,indy,indz)
! initialize 3d real to complex fft
         implicit none
         integer :: indx, indy, indz
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, ntpose = 1, kstrt = 1
         integer :: nxvh = 1, nyv = 1, nzv = 1
         integer :: kxyp = 1, kyp = 1, kyzp = 1, kzp = 1
         integer :: kxypd = 1, kypd = 1, kyzpd = 1, kzpd = 1, kzyp = 1
         integer :: jblok = 1, kblok = 1, lblok = 1, mblok = 1
         integer :: jkblok = 1, mlblok = 1
         integer :: nxhyzd, nxyzhd
         real :: f
         complex, dimension(1,1,1,1) :: g, h, bs, br
         nxhyzd = size(mixup); nxyzhd = size(sct)
         call PFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy,indz,&
     &kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jb&
     &lok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         end subroutine ipfft32rinit
!
         subroutine ipfft32r(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,k&
     &strt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(f,4),size(g,&
     &4))) :: bs
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(g,4),size(h,&
     &4))) :: br
         real :: tf
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT32R(f(1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         else
            call PFFT32R(f(2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft32r
!
         subroutine ipfft32r3(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,&
     &kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         real :: tf
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT32R3(f(1,1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         else
            call PFFT32R3(f(1,2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft32r3
!
         subroutine ipfft32rxinit(mixup,sct,indx,indy,indz)
! initialize optimized 3d real to complex fft
         implicit none
         integer :: indx, indy, indz
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, ntpose = 1, kstrt = 1
         integer :: nxvh = 1, nyv = 1, nzv = 1
         integer :: kxyp = 1, kyp = 1, kyzp = 1, kzp = 1
         integer :: kxypd = 1, kypd = 1, kyzpd = 1, kzpd = 1
         integer :: jblok = 1, kblok = 1, lblok = 1, mblok = 1
         integer :: nxhyzd, nxyzhd
         real :: f
         complex, dimension(1,1,1,1) :: g, h
         nxhyzd = size(mixup); nxyzhd = size(sct)
         call PFFT32RX(f,g,h,isign,ntpose,mixup,sct,indx,indy,indz,kstrt&
     &,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblok,kblok,&
     &lblok,mblok,nxhyzd,nxyzhd)
         end subroutine ipfft32rxinit
!
         subroutine ipfft32rx(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,&
     &kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform optimized 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order
         real :: tf
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT32RX(f(1,1,1,1),g,h,isign,ntpose,mixup,sct,indx,ind&
     &y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,&
     &jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         else
            call PFFT32RX(f(2,2,2,1),g,h,isign,ntpose,mixup,sct,indx,ind&
     &y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,&
     &jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft32rx
!
         subroutine ipfft32rx3(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform optimized 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order
         real :: tf
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT32RX3(f(1,1,1,1,1),g,h,isign,ntpose,mixup,sct,indx,&
     &indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kz&
     &pd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         else
            call PFFT32RX3(f(1,2,2,2,1),g,h,isign,ntpose,mixup,sct,indx,&
     &indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kz&
     &pd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft32rx3
!
         subroutine iwpfft32rinit(mixup,sct,indx,indy,indz)
! initialize 3d real to complex fft
         implicit none
         integer :: indx, indy, indz
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhyzd, nxyzhd
         nxhyzd = size(mixup); nxyzhd = size(sct)
         call WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         end subroutine iwpfft32rinit
!
         subroutine iwpfft32r(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,&
     &kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(f,4),size(g,&
     &4))) :: bs
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(g,4),size(h,&
     &4))) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT32R(f(1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,k&
     &yzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd&
     &)
         else
            call WPFFT32R(f(2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,k&
     &yzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd&
     &)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft32r
!
         subroutine iwpfft32r3(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT32R3(f(1,1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct&
     &,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kyp&
     &d,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxy&
     &zhd)
         else
            call WPFFT32R3(f(1,2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct&
     &,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kyp&
     &d,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxy&
     &zhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft32r3
!
         subroutine iwpfft32rx(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform optimized 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT32RX(f(1,1,1,1),g,h,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         else
            call WPFFT32RX(f(2,2,2,1),g,h,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft32rx
!
         subroutine iwpfft32rx3(f,g,h,isign,mixup,sct,tfft,indx,indy,ind&
     &z,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform optimized 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT32RX3(f(1,1,1,1,1),g,h,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)
         else
            call WPFFT32RX3(f(1,2,2,2,1),g,h,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd)

         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft32rx3
!
         subroutine iwpfft32rn(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d vector real to complex fft for n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         complex, dimension(size(f,1),size(f,2)/2) :: ss
         integer :: ntpose = 1, ndim, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kzpd = size(f,4); lblok = size(f,5)/kblok
         nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT32RN(f(1,1,1,1,1),g,h,bs,br,ss,isign,ntpose,mixup,&
     &sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,&
     &kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nx&
     &hyzd,nxyzhd)
         else
            call WPFFT32RN(f(1,2,2,2,1),g,h,bs,br,ss,isign,ntpose,mixup,&
     &sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,&
     &kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nx&
     &hyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft32rn
!
         subroutine iwpfft32rxn(f,g,h,isign,mixup,sct,tfft,indx,indy,ind&
     &z,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform optimized 3d vector real to complex fft for n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g
         complex, dimension(:,:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         complex, dimension(size(f,1),size(f,2)/2) :: ss
         integer :: ntpose = 1, ndim, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kzpd = size(f,4); lblok = size(f,5)/kblok
         nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT32RXN(f(1,1,1,1,1),g,h,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,&
     &kyzpd,kzpd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd)
         else
            call WPFFT32RXN(f(1,2,2,2,1),g,h,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,&
     &kyzpd,kzpd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft32rxn
!
         subroutine iwp1fft32rn(f1,f2,g1,g2,h1,h2,isign,mixup,sct,tfft,i&
     &ndx,indy,indz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d real to complex fft for a scalar and n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f1
         real, dimension(:,:,:,:,:), pointer :: f2
         complex, dimension(:,:,:,:), pointer :: g1, h1
         complex, dimension(:,:,:,:,:), pointer :: g2, h2
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim1, ndim2, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(1+size(f2,1),kxyp,max(kyzp,kyp),kzp,max(size&
     &(f1,4),size(g1,4))) :: bs
         complex, dimension(1+size(f2,1),kxyp,max(kyzp,kyp),kzp,max(size&
     &(g1,4),size(h1,4))) :: br
         complex, dimension(1+size(f2,1),size(f1,1)/2) :: ss
         real :: tf, ttp
         double precision :: dtime
         ndim1 = 1; ndim2 = size(f2,1)
         nxvh = size(f1,1)/2; kypd = size(f1,2)
         kzpd = size(f1,3); lblok = size(f1,4)/kblok
         nyv = size(g1,1); kxypd = size(g1,2)
         nzv = size(h1,1); kyzpd = size(h1,3); jblok = size(h1,4)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WP2FFT32RN(f1(1,1,1,1),f2(1,1,1,1,1),g1,g2,h1,h2,bs,br,&
     &ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kx&
     &yp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblo&
     &k,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd)
         else
            call WP2FFT32RN(f1(2,2,2,1),f2(1,2,2,2,1),g1,g2,h1,h2,bs,br,&
     &ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kx&
     &yp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblo&
     &k,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwp1fft32rn
!
         subroutine iwp2fft32rn(f1,f2,g1,g2,h1,h2,isign,mixup,sct,tfft,i&
     &ndx,indy,indz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d real to complex fft for two n component vectors
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f1, f2
         complex, dimension(:,:,:,:,:), pointer :: g1, g2, h1, h2
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim1, ndim2, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order
         complex, dimension(size(f1,1)+size(f2,1),kxyp,max(kyzp,kyp),kzp&
     &,max(size(f1,5),size(g1,5))) :: bs
         complex, dimension(size(f1,1)+size(f2,1),kxyp,max(kyzp,kyp),kzp&
     &,max(size(g1,5),size(h1,5))) :: br
         complex, dimension(size(f1,1)+size(f2,1),size(f1,2)/2) :: ss
         real :: tf, ttp
         double precision :: dtime
         ndim1 = size(f1,1); ndim2 = size(f2,1)
         nxvh = size(f1,2)/2; kypd = size(f1,3)
         kzpd = size(f1,4); lblok = size(f1,5)/kblok
         nyv = size(g1,2); kxypd = size(g1,3)
         nzv = size(h1,2); kyzpd = size(h1,4); jblok = size(h1,5)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WP2FFT32RN(f1(1,1,1,1,1),f2(1,1,1,1,1),g1,g2,h1,h2,bs,b&
     &r,ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,&
     &kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lb&
     &lok,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd)
         else
            call WP2FFT32RN(f1(1,2,2,2,1),f2(1,2,2,2,1),g1,g2,h1,h2,bs,b&
     &r,ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,&
     &kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lb&
     &lok,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwp2fft32rn
!
         subroutine ipfft32c(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,k&
     &strt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform 3d scalar complex to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real :: tfft
         complex, dimension(:,:,:,:), pointer :: f, g
         complex, dimension(:,:,:,:), pointer :: h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxv, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxyzd, nxyzhd, order
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(f,4),size(g,&
     &4))) :: bs
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(g,4),size(h,&
     &4))) :: br
         real :: tf
         double precision :: dtime
         nxv = size(f,1); kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT32C(f(1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct,ind&
     &x,indy,indz,kstrt,nxv,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxyzd,nxyzhd)
         else
            call PFFT32C(f(2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct,ind&
     &x,indy,indz,kstrt,nxv,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxyzd,nxyzhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft32c
!
         subroutine ipfft32cinit(mixup,sct,indx,indy,indz)
! initialize 3d complex to complex fft
         implicit none
         integer :: indx, indy, indz
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, ntpose = 1, kstrt = 1
         integer :: nxv = 1, nyv = 1, nzv = 1
         integer :: kxyp = 1, kyp = 1, kyzp = 1, kzp = 1
         integer :: kxypd = 1, kypd = 1, kyzpd = 1, kzpd = 1, kzyp = 1
         integer :: jblok = 1, kblok = 1, lblok = 1, mblok = 1
         integer :: jkblok = 1, mlblok = 1
         integer :: nxyzd, nxyzhd
         complex :: f
         complex, dimension(1,1,1,1) :: g, h, bs, br
         nxyzd = size(mixup); nxyzhd = size(sct)
         call PFFT32C(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy,indz,&
     &kstrt,nxv,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jbl&
     &ok,kblok,jkblok,lblok,mblok,mlblok,nxyzd,nxyzhd)
         end subroutine ipfft32cinit
!
      end module pfft32d

