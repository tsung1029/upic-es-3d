!-----------------------------------------------------------------------
!
      module mprbpush32d
!
! Fortran90 interface to parallel 3d PIC F77 library mprbpush32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: january 16, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use prbpush32d, only: wtimer, retard
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: rdjpost, rpush, retard
!
! define interface to Fortran77 procedures
!
      interface
         subroutine MPGRJPOST32(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,i&
     &dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxv,nypmx,nzpmx,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJPOST32(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,ci,&
     &nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc,cup&
     &,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n81, ipbc, nmt, ierr
         real :: qm, dt, ci
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
         subroutine MPGRJPOST32L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxv,nypmx,nzpmx,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJPOST32L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(3,nxyzp,mnblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,ci&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc,cu&
     &p,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n24, ipbc, nmt, ierr
         real :: qm, dt, ci
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
         subroutine MPGRPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,ny&
     &,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,i&
     &err)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,i&
     &err)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,&
     &ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,c&
     &i,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idt&
     &ask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,&
     &ci,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,id&
     &task,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,&
     &ci,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,id&
     &task,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc&
     &,ci,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,i&
     &dtask,nmt,ierr)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok):: npp, nps
         integer, dimension(idds,mnblok) :: noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdjpost
         module procedure impgrjpost32
      end interface
!
      interface rpush
         module procedure impgrpush32
         module procedure impgrbpush32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgrpush32(part,fxyz,npp,noff,qbm,dt,ci,ek,tpush,nx&
     &,ny,nz,ipbc,inorder,popt)
! multi-tasking relativistic particle push with 3d electrostatic fields,
! 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: fxyz
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
               call MPGSRPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,&
     &ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,&
     &ierr)
            else
               call MPGRPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,i&
     &err)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSRPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,i&
     &err)
            else
               call MPGRPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,ny&
     &,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ie&
     &rr)
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
         end subroutine impgrpush32
!
         subroutine impgrjpost32(part,cu,npp,noff,qm,dt,ci,tdjpost,nx,ny&
     &,nz,ipbc,inorder,djopt)
! multi-tasking, relativisitc current deposit, 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
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
               call MPGSRJPOST32L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSRJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,ci&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n24,ipbc,cup&
     &,idtask,nmt,ierr)
            else
               call MPGRJPOST32L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSRJPOST32(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSRJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,ci,&
     &nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n81,ipbc,cup,&
     &idtask,nmt,ierr)
            else
               call MPGRJPOST32(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,i&
     &dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
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
         end subroutine impgrjpost32
!
         subroutine impgrbpush32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,e&
     &k,tpush,nx,ny,nz,ipbc,inorder,popt)
! multi-tasking, relativistic particle push with 3d electromagnetic
! fields, 2d partition
         implicit none
         integer :: nx, ny, nz, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
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
               call MPGSRBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc&
     &,ci,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,i&
     &dtask,nmt,ierr)
            else
               call MPGRBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,&
     &ci,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,id&
     &task,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSRBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,&
     &ci,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,id&
     &task,nmt,ierr)
            else
               call MPGRBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,c&
     &i,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idt&
     &ask,nmt,ierr)
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
         end subroutine impgrbpush32
!
      end module mprbpush32d
