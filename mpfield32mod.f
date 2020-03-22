!-----------------------------------------------------------------------
!
      module mpfield32d
!
! Fortran90 interface to 3d parallel PIC F77 library mpfield32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: june 10, 2005
!
      use globals, only: LINEAR, QUADRATIC
      use pfield32d
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: cguard, sguard, aguard, zguard
      public :: pois_init, pois, cuperp, bpois
      public :: ibpois, maxwel, emfield, avpot, gtmodes, ptmodes
      public :: fft_init, fft, fftn, fftc_init
      public :: ipdivf32, ipgradf32, ipcurlf32
      public :: mft, mftn
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,ind&
     &y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,&
     &kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyzip,ift&
     &ask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, ierr
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
         subroutine MPFFT32RX(f,g,h,isign,ntpose,mixup,sct,indx,indy,ind&
     &z,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblok&
     &,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
         integer ::  nmt, ierr
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
         subroutine MPFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,in&
     &dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd&
     &,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyzip,if&
     &task,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: nxhyzd, nxyzhd, nmt, ierr
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
         subroutine MPFFT32RX3(f,g,h,isign,ntpose,mixup,sct,indx,indy,in&
     &dz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblo&
     &k,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
         integer :: nmt, ierr
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
         subroutine MPFFT32RN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,indx&
     &,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd,nxyzhd,&
     &kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
         integer :: jblok, kblok, jkblok, lblok, mblok, mlblok
         integer :: ndim, nxhyzd, nxyzhd, nmt, ierr
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
         subroutine MPFFT32RXN(f,g,h,ss,isign,ntpose,mixup,sct,indx,indy&
     &,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,j&
     &blok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         implicit none
         integer :: isign
         integer :: ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
         integer :: kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
         integer :: jblok, kblok, lblok, mblok, ndim, nxhyzd, nxyzhd
         integer :: nmt, ierr
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
!
! define generic interfaces to Fortran90 library
!
      interface mft
!        module procedure impfft32r
!        module procedure impfft32r3
         module procedure impfft32rx
         module procedure impfft32rx3
      end interface
!
      interface mftn
!        module procedure impfft32rn
         module procedure impfft32rxn
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impfft32r(f,g,h,isign,mixup,sct,tfft,indx,indy,indz,&
     &kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,nmt,inorder)
! perform multi-tasking 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok, nmt
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order, ierr
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(f,4),size(g,&
     &4))) :: bs
         complex, dimension(kxyp,max(kyzp,kyp),kzp,max(size(g,4),size(h,&
     &4))) :: br
         integer, dimension(nmt) :: kyzip, iftask
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
         call wtimer(tfft,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32R(f(1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct,in&
     &dx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd&
     &,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyz&
     &ip,iftask,nmt,ierr)
         else
            call MPFFT32R(f(2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct,in&
     &dx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd&
     &,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyz&
     &ip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tfft,dtime)
         end subroutine impfft32r
!
         subroutine impfft32r3(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,nmt,inorder)
! perform multi-tasking 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok, nmt
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order, ierr
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         integer, dimension(nmt) :: kyzip, iftask
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
         call wtimer(tfft,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32R3(f(1,1,1,1,1),g,h,bs,br,isign,ntpose,mixup,sct&
     &,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,ky&
     &zpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,&
     &kyzip,iftask,nmt,ierr)
         else
            call MPFFT32R3(f(1,2,2,2,1),g,h,bs,br,isign,ntpose,mixup,sct&
     &,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,ky&
     &zpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,&
     &kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tfft,dtime)
         end subroutine impfft32r3
!
         subroutine impfft32rx(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,nmt,inorder)
! perform multi-tasking optimized 3d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok, nmt
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv, kxypd
         integer :: kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order, ierr
         integer, dimension(nmt) :: kyzip, iftask
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kzpd = size(f,3)
         lblok = size(f,4)/kblok; nyv = size(g,1); kxypd = size(g,2)
         nzv = size(h,1); kyzpd = size(h,3); jblok = size(h,4)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfft,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32RX(f(1,1,1,1),g,h,isign,ntpose,mixup,sct,indx,in&
     &dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd&
     &,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32RX(f(2,2,2,1),g,h,isign,ntpose,mixup,sct,indx,in&
     &dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd&
     &,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tfft,dtime)
         end subroutine impfft32rx
!
         subroutine impfft32rx3(f,g,h,isign,mixup,sct,tfft,indx,indy,ind&
     &z,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,nmt,inorder)
! perform multi-tasking optimized 3d vector real to complex fft
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok, nmt
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         integer :: ntpose = 1, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order, ierr
         integer, dimension(nmt) :: kyzip, iftask
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kzpd = size(f,4)
         lblok = size(f,5)/kblok; nyv = size(g,2); kxypd = size(g,3)
         nzv = size(h,2); kyzpd = size(h,4); jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call MPFFT32RX3(f(1,1,1,1,1),g,h,isign,ntpose,mixup,sct,indx&
     &,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32RX3(f(1,2,2,2,1),g,h,isign,ntpose,mixup,sct,indx&
     &,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,k&
     &zpd,jblok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tfft,dtime)
         end subroutine impfft32rx3
!
         subroutine impfft32rn(f,g,h,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,nmt,inorder)
! perform multi-tasking 3d vector real to complex fft
! for n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok, nmt
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         complex, dimension(size(f,1),size(f,2)/2,nmt+1) :: ss
         integer :: ntpose = 1, ndim, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok, kzyp, jkblok, mlblok
         integer :: nxhyzd, nxyzhd, order, ierr
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(f,5),size(&
     &g,5))) :: bs
         complex, dimension(3,kxyp,max(kyzp,kyp),kzp,max(size(g,5),size(&
     &h,5))) :: br
         integer, dimension(nmt) :: kyzip, iftask
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kzpd = size(f,4); lblok = size(f,5)/kblok; nyv = size(g,2)
         kxypd = size(g,3); nzv = size(h,2); kyzpd = size(h,4)
         jblok = size(h,5)/mblok; kzyp = max(kyzp,kyp)
         jkblok = max(jblok,kblok); mlblok = max(mblok,lblok)
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfft,dtime,-1)
         if (order==LINEAR) then
            call MPFFT32RN(f(1,1,1,1,1),g,h,bs,br,ss,isign,ntpose,mixup,&
     &sct,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd&
     &,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd&
     &,nxyzhd,kyzip,iftask,nmt,ierr)
         else
            call MPFFT32RN(f(1,2,2,2,1),g,h,bs,br,ss,isign,ntpose,mixup,&
     &sct,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd&
     &,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd&
     &,nxyzhd,kyzip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tfft,dtime)
         end subroutine impfft32rn
!
         subroutine impfft32rxn(f,g,h,isign,mixup,sct,tfft,indx,indy,ind&
     &z,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,nmt,inorder)
! perform multi-tasking  3d vector real to complex fft
! for n component vector
         implicit none
         integer :: isign, indx, indy, indz, kstrt, kxyp, kyp, kyzp, kzp
         integer :: kblok, mblok, nmt
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:,:), pointer :: g, h
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         complex, dimension(size(f,1),size(f,2)/2,nmt+1) :: ss
         integer :: ntpose = 1, ndim, nxvh, kypd, kzpd, nzv, nyv
         integer :: kxypd, kyzpd, jblok, lblok
         integer :: nxhyzd, nxyzhd, order, ierr
         integer, dimension(nmt) :: kyzip, iftask
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kzpd = size(f,4); lblok = size(f,5)/kblok; nyv = size(g,2)
         kxypd = size(g,3); nzv = size(h,2); kyzpd = size(h,4)
         jblok = size(h,5)/mblok
         nxhyzd = size(mixup); nxyzhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call MPFFT32RXN(f(1,1,1,1,1),g,h,ss,isign,ntpose,mixup,sct,i&
     &ndx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzp&
     &d,kzpd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask,nmt&
     &,ierr)
         else
            call MPFFT32RXN(f(1,2,2,2,1),g,h,ss,isign,ntpose,mixup,sct,i&
     &ndx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzp&
     &d,kzpd,jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask,nmt&
     &,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tfft,dtime)
         end subroutine impfft32rxn
!
      end module mpfield32d

