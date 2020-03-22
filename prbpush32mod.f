!-----------------------------------------------------------------------
!
      module prbpush32d
!
! Fortran90 interface to 3d PIC Fortran77 library prbpush32lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: january 16, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: rdjpost, rpush, retard
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PGRJPOST32(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,&
     &npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRJPOST32(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,ci,nx,ny&
     &,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n81, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n81,npd,mnblok) :: nn
         real, dimension(n81,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PGRJPOST32L(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRJPOST32L(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,ci,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp
         integer :: idds, npd, n24, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3*nxvyzp,mnblok) :: cu
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         integer, dimension(n24,npd,mnblok) :: nn
         real, dimension(n24,npd,mnblok) :: amxyz
         end subroutine
      end interface
      interface
         subroutine PGRPUSH32(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,i&
     &dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRPUSH32(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGRPUSH32L(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRPUSH32L(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,&
     &nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxv,nypmx,nzpmx,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp
         integer :: idds, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,mnblok) :: part
         real, dimension(3,nxyzp,mnblok) :: fxyz, bxyz
         integer, dimension(mnblok) :: npp
         integer, dimension(idds,mnblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PRRETARD32(part,npp,dtc,ci,nx,ny,nz,idimp,npmax,mnbl&
     &ok,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mnblok, ipbc
         real :: dtc, ci
         real, dimension(idimp,npmax,mnblok) :: part
         integer, dimension(mnblok) :: npp
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdjpost
         module procedure ipgrjpost32
      end interface
!
      interface rpush
         module procedure ipgrpush32
         module procedure ipgrbpush32
      end interface
!
      interface retard
         module procedure iprretard32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipgrpush32(part,fxyz,npp,noff,qbm,dt,ci,ek,tpush,nx,&
     &ny,nz,ipbc,inorder,popt)
! push relativistic particles with 3d electrostatic fields, 2d partition
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
         integer :: order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSRPUSH32L(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz&
     &,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGRPUSH32L(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSRPUSH32(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,&
     &idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGRPUSH32(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,i&
     &dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgrpush32
!
         subroutine ipgrjpost32(part,cu,npp,noff,qm,dt,ci,tdjpost,nx,ny,&
     &nz,ipbc,inorder,djopt)
! deposit relativisitc current, 2d partition
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
         integer :: order, opt
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, n24 = 24, n81 = 81
         integer, dimension(n81,npd,size(part,3)) :: nn
         real, dimension(n81,npd,size(part,3)) :: amxyz
         real :: tj
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         nxv = size(cu,2); nypmx = size(cu,3); nzpmx = size(cu,4)
         idds = size(noff,1)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSRJPOST32L(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idim&
     &p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else if (opt==VECTOR) then
               call PGSRJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,ci,nx,n&
     &y,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n24,ipbc)
            else
               call PGRJPOST32L(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSRJPOST32(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp&
     &,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else if (opt==VECTOR) then
               call PGSRJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,ci,nx,ny&
     &,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,npd,n81,ipbc)
            else
               call PGRJPOST32(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,&
     &npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         endif
! record time
         call wtimer(tj,dtime)
         tdjpost = tdjpost + tj
         end subroutine ipgrjpost32
!
         subroutine ipgrbpush32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek&
     &,tpush,nx,ny,nz,ipbc,inorder,popt)
! push relativistic particles with 3d electromagnetic fields,
! 2d partition
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
         integer :: order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); idds = size(noff,1)
         nxv = size(fxyz,2); nypmx = size(fxyz,3); nzpmx = size(fxyz,4)
         nxyzp = nxv*nypmx*nzpmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSRBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGRBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSRBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek&
     &,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
            else
               call PGRBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,&
     &nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgrbpush32
!
         subroutine iprretard32(part,npp,dtc,ci,nx,ny,nz,ipbc)
! retards relativistic particle positions half time-step
         implicit none
         integer :: nx, ny, nz, ipbc
         real :: dtc, ci
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, mnblok
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3)
         call PRRETARD32(part,npp,dtc,ci,nx,ny,nz,idimp,npmax,mnblok,ipb&
     &c)
         end subroutine iprretard32
!
      end module prbpush32d
