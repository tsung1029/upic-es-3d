!-----------------------------------------------------------------------
!
      module mpbpush32d
!
! Fortran90 interface to parallel 3d PIC F77 library mpbpush32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: january 16, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use pbpush32d, only: wtimer, retard
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: djpost, push, retard
!
! define interface to Fortran77 procedures
!
      interface
         subroutine MPJDOST32(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idi&
     &mp,npmax,mnblok,nxv,nypmx,nzpmx,idds,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: cux, cuy, cuz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxv,nypmx,nzpmx,3,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGJPOST32(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxv,nypmx,nzpmx,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJPOST32(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,id&
     &imp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,cup,idtask,nmt,ierr&
     &)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, n81, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n81,npd,mnblok,nmt+1) :: nn
         real, dimension(n81,npd,mnblok,nmt+1) :: amxyz
         real, dimension(3*nxvyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc,cup,idt&
     &ask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n81, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n81,npd,mnblok,nmt+1) :: nn
         real, dimension(n81,npd,mnblok,nmt+1) :: amxyz
         real, dimension(3*nxvyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPJDOST32L(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,id&
     &imp,npmax,mnblok,nxv,nypmx,nzpmx,idds,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: cux, cuy, cuz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nxv,nypmx,nzpmx,3,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGJPOST32L(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxv,nypmx,nzpmx,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJPOST32L(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idi&
     &mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,i&
     &dimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,cup,idtask,nmt,ier&
     &r)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
         integer :: npd, n24, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n24,npd,mnblok,nmt+1) :: nn
         real, dimension(n24,npd,mnblok,nmt+1) :: amxyz
         real, dimension(3*nxvyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc,cup,id&
     &task,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n24, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n24,npd,mnblok,nmt+1) :: nn
         real, dimension(n24,npd,mnblok,nmt+1) :: amxyz
         real, dimension(3*nxvyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH32(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,dt&
     &,ek,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr&
     &)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: bx, by, bz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ek&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,&
     &nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,e&
     &k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask&
     &,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH32L(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,d&
     &t,ek,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ier&
     &r)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: bx, by, bz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,e&
     &k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask&
     &,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,&
     &ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH32C(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,d&
     &t,ek,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ier&
     &r)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: bx, by, bz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH32C(part,fxyz,bxyz,npp,nps,noff,qbm,dt,ek,nx&
     &,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt&
     &,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH32CL(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,&
     &dt,ek,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
         integer :: nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: fx, fy, fz
         real, dimension(nxv,nypmx,nzpmx,mnblok) :: bx, by, bz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH32CL(part,fxyz,bxyz,npp,nps,noff,qbm,dt,ek,n&
     &x,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nm&
     &t,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure impgjpost32
      end interface
!
      interface push
         module procedure impgbpush32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgjpost32(part,cu,npp,noff,qm,dt,tdjpost,nx,ny,nz,&
     &ipbc,inorder,djopt)
! multi-tasking current deposit, 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, n24 = 24, n81 = 81
         integer, dimension(ntasks) :: idtask
         real, dimension(3,size(cu,2),size(cu,3),size(cu,4),size(cu,5),n&
     &tasks) :: cup
         integer, dimension(n81,npd,size(part,3),ntasks+1) :: nn
         real, dimension(n81,npd,size(part,3),ntasks+1) :: amxyz
         real :: tj
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSJPOST32L(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idi&
     &mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n24,ipbc,cup,idt&
     &ask,nmt,ierr)
            else
               call MPGJPOST32L(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSJPOST32(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n81,ipbc,cup,idta&
     &sk,nmt,ierr)
            else
               call MPGJPOST32(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tj,dtime)
         tdjpost = tdjpost + tj
         end subroutine impgjpost32
!
         subroutine impgbpush32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,tp&
     &ush,nx,ny,nz,ipbc,inorder,popt)
! multi-tasking particle push with 3d electromagnetic fields,
! 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: fxyz, bxyz
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,&
     &ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtas&
     &k,nmt,ierr)
            else
               call MPGBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,e&
     &k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask&
     &,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,e&
     &k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask&
     &,nmt,ierr)
            else
               call MPGBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ek&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,&
     &nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine impgbpush32
!
      end module mpbpush32d
