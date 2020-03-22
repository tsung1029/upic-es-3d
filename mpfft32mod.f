!-----------------------------------------------------------------------
!
      module mpfft32d
!
! Fortran90 interface to 3d parallel PIC F77 library mpfft32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: june 16, 2008
!
      use globals, only: LINEAR, QUADRATIC
      use pfft32d, only: wtimer, fft_init, fftc_init
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: fft_init, fft, fftn, fftc_init
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyzip&
     &,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT32RX(f,g,h,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,j&
     &blok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
         integer ::  nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyzi&
     &p,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(3,kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(3,kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT32RX3(f,g,h,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,&
     &jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(3,nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT32RN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd,nxy&
     &zhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: ndim, nxhyzd, nxyzhd, nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd,jblok*mblok) :: h
         complex, dimension(ndim,kxyp,kzyp,kzp,jkblok*lblok) :: bs
         complex, dimension(ndim,kxyp,kzyp,kzp,jblok*mlblok) :: br
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         complex, dimension(ndim,nxvh,nmt+1) :: ss
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT32RXN(f,g,h,ss,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kz&
     &pd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask,nmt,ier&
     &r)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok, ndim, nxhyzd, nxyzhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(ndim,nyv,kxypd,kzpd,jblok*lblok) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd,jblok*mblok) :: h
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         complex, dimension(ndim,nxvh,nmt+1) :: ss
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
      interface
         subroutine MP2FFT32RN(f1,f2,g1,g2,h1,h2,bs,br,ss,isign,ntpose,m&
     &ixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,k&
     &xypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nd&
     &im1,ndim2,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, indz, kstrt
         integer :: nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
         integer :: kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: ndim1, ndim2, nxhyzd, nxyzhd
         integer :: nmt, ierr
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
         complex, dimension(ndim1+ndim2,nxvh,nmt+1) :: ss
         integer, dimension(nmt) :: kyzip, iftask
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft
         module procedure impfft32r
         module procedure impfft32r3
!        module procedure impfft32rx
!        module procedure impfft32rx3
      end interface
!
      interface fftn
         module procedure impfft32rn
!        module procedure impfft32rxn
         module procedure imp1fft32rn
         module procedure imp2fft32rn
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impfft32r(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,&
     &kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(f,4),size(g,&
     &4))) :: bs
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(g,4),size(h,&
     &4))) :: br
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32R(f(1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,k&
     &yzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd&
     &,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32R(f(2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,k&
     &yzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd&
     &,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft32r
!
         subroutine impfft32r3(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32R3(f(1,1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct&
     &,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kyp&
     &d,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxy&
     &zhd,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32R3(f(1,2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct&
     &,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kyp&
     &d,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxy&
     &zhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft32r3
!
         subroutine impfft32rx(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking optimized 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32RX(f(1,1,1,1),g,h,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32RX(f(2,2,2,1),g,h,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,&
     &kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft32rx
!
         subroutine impfft32rx3(f,g,h,isign,mixup,sct,tfft,indx,indy,ind&
     &z,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking optimized 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32RX3(f(1,1,1,1,1),g,h,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ier&
     &r)
         else
            call MPFFT32RX3(f(1,2,2,2,1),g,h,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyz&
     &pd,kzpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ier&
     &r)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft32rx3
!
         subroutine impfft32rn(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking 3d vector real to complex fft
! for n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         complex, dimension(size(f,1),size(f,2)/2,ntasks+1) :: ss
         integer :: ntpose = 1, ndim, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kzpd = size(f,4); lblok = size(f,5)/kblok; nyv = size(g,2)
         kxypd = size(g,3); nzv = size(h,2); kyzpd = size(h,4)
         jblok = size(h,5)/mblok; kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32RN(f(1,1,1,1,1),g,h,bs,br,ss,isign,ntpose,mixup,&
     &sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,&
     &kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nx&
     &hyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32RN(f(1,2,2,2,1),g,h,bs,br,ss,isign,ntpose,mixup,&
     &sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,&
     &kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nx&
     &hyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft32rn
!
         subroutine impfft32rxn(f,g,h,isign,mixup,sct,tfft,indx,indy,ind&
     &z,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking  3d vector real to complex fft
! for n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         complex, dimension(size(f,1),size(f,2)/2,ntasks+1) :: ss
         integer :: ntpose = 1, ndim, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kzpd = size(f,4); lblok = size(f,5)/kblok; nyv = size(g,2)
         kxypd = size(g,3); nzv = size(h,2); kyzpd = size(h,4)
         jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32RXN(f(1,1,1,1,1),g,h,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,&
     &kyzpd,kzpd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask&
     &,nmt,ierr)
         else
            call MPFFT32RXN(f(1,2,2,2,1),g,h,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,&
     &kyzpd,kzpd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask&
     &,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft32rxn
!
         subroutine imp1fft32rn(f1,f2,g1,g2,h1,h2,isign,mixup,sct,tfft,i&
     &ndx,indy,indz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking 3d real to complex fft for a scalar and n
! component vector
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
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         complex, dimension(1+size(f2,1),kxyp,max(kyzp,kyp),kzp,max(size&
     &(f1,4),size(g1,4))) :: bs
         complex, dimension(1+size(f2,1),kxyp,max(kyzp,kyp),kzp,max(size&
     &(g1,4),size(h1,4))) :: br
         complex, dimension(1+size(f2,1),size(f1,1)/2,ntasks+1) :: ss
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim1 = 1; ndim2 = size(f2,1)
         nxvh = size(f1,1)/2; kypd = size(f1,2)
         kzpd = size(f1,3); lblok = size(f1,4)/kblok; nyv = size(g1,1)
         kxypd = size(g1,2); nzv = size(h1,1); kyzpd = size(h1,3)
         jblok = size(h1,4)/mblok; kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MP2FFT32RN(f1(1,1,1,1),f2(1,1,1,1,1),g1,g2,h1,h2,bs,br,&
     &ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kx&
     &yp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblo&
     &k,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MP2FFT32RN(f1(2,2,2,1),f2(1,2,2,2,1),g1,g2,h1,h2,bs,br,&
     &ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kx&
     &yp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblo&
     &k,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine imp1fft32rn
!
         subroutine imp2fft32rn(f1,f2,g1,g2,h1,h2,isign,mixup,sct,tfft,i&
     &ndx,indy,indz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! perform multi-tasking 3d real to complex fft for two n component
! vectors
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
         integer :: nxhyzd, nxyzhd, nmt, order, ierr
         complex, dimension(size(f1,1)+size(f2,1),kxyp,max(kyzp,kyp),kzp&
     &,max(size(f1,5),size(g1,5))) :: bs
         complex, dimension(size(f1,1)+size(f2,1),kxyp,max(kyzp,kyp),kzp&
     &,max(size(g1,5),size(h1,5))) :: br
         complex, dimension(size(f1,1)+size(f2,1),size(f1,2)/2,ntasks+1)&
     &:: ss
         integer, dimension(ntasks) :: kyzip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim1 = size(f1,1); ndim2 = size(f2,1)
         nxvh = size(f1,2)/2; kypd = size(f1,3)
         kzpd = size(f1,4); lblok = size(f1,5)/kblok; nyv = size(g1,2)
         kxypd = size(g1,3); nzv = size(h1,2); kyzpd = size(h1,4)
         jblok = size(h1,5)/mblok; kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MP2FFT32RN(f1(1,1,1,1,1),f2(1,1,1,1,1),g1,g2,h1,h2,bs,b&
     &r,ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,&
     &kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lb&
     &lok,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MP2FFT32RN(f1(1,2,2,2,1),f2(1,2,2,2,1),g1,g2,h1,h2,bs,b&
     &r,ss,isign,ntpose,mixup,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,&
     &kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lb&
     &lok,mblok,mlblok,ndim1,ndim2,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine imp2fft32rn
!
      end module mpfft32d

