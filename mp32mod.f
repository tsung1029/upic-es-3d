!-----------------------------------------------------------------------
!
      module mp32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library mp32lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: april 11, 2008
!
      use p0d
      use mp0d, only: ntasks
      use p32d, only: get_funit, pwtimer, plsum, plmax, plbcast, writebf&
     &, readbf, wtimer, wrdata, rddata, fcomp, dcomp, pcguard, pcguardp,&
     &pnlcguard, paguard, paguardp, pncguardp, pnlaguard, pnaguardp,    &
     &pfmove, repart, fnoff, trpsin, haftrp, dblsin, hafdbl, ztrp
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PPINIT, PPID, PPEXIT, HARTBEAT
      public :: get_funit, pwtimer, plsum, plmax, plbcast
      public :: writebf, readbf, wtimer, wrdata, rddata
      public :: fcomp, dcomp, pmove, pcguard, pcguardp, pnlcguard
      public :: paguard, paguardp, pncguardp, pnlaguard, pnaguardp
      public :: pfmove, repart
      public :: fnoff, trpsin, haftrp, dblsin, hafdbl, ztrp, pmoves
!
! buffer data for particle managers
      real, dimension(:,:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(:,:), allocatable :: ihole
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MPMOVE32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,&
     &ihole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok&
     &,idps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp, npq
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPXMOV32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,&
     &ihole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok&
     &,idps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(npmax,mblok*nblok) :: maskp
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp, npq
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPMOVES32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl&
     &,ihole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblo&
     &k,idps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp, npq
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPXMOVS32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl&
     &,ihole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblo&
     &k,idps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer :: mblok, nblok, idps, nbmax, idds, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(npmax,mblok*nblok) :: maskp
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp, npq
         real, dimension(idimp,nbmax,mblok*nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,mblok*nblok) :: rbufl, rbufr
         integer, dimension(idds,mblok*nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,mblok*nblok) :: ihole
         integer, dimension(9) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface pmove
         module procedure impmove32
         module procedure impdmove32
      end interface
!
      interface pmoves
         module procedure impmoves32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine impmove32(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz&
     &,nbmax,idds,mblok,vt,ierr)
! multi-tasking particle manager
! non-uniform 2d partition boundaries in 3d code
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, idps, ntmax, nmt
         integer, dimension(ntasks) :: idtask
         integer, dimension(size(edges,1),size(edges,2),ntasks) :: jssp
         integer, dimension(size(npp)) :: npq
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(9) :: info
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         nmt = ntasks
         th = 0.0
! check if size of buffers has changed
         if (szbuf.ne.nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call MPXMOV32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
         else
            call MPMOVE32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
         endif
         if ((ierr /= 0) .or. (info(1) /= 0)) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine impmove32
!
         subroutine impdmove32(part,edges,npp,anpav,pibal,tmove,ny,nz,ks&
     &trt,nvpy,nvpz,nbmax,idds,mblok,vt,ierr)
! multi-tasking particle manager
! non-uniform 2d partition boundaries in 3d code
! returns load imbalance
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real :: anpav, pibal
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, idps, ntmax, nmt
         integer, dimension(ntasks) :: idtask
         integer, dimension(size(edges,1),size(edges,2),ntasks) :: jssp
         integer, dimension(size(npp)) :: npq
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(9) :: info
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         nmt = ntasks
         th = 0.0
! check if size of buffers has changed
         if (szbuf.ne.nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call MPXMOV32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
         else
            call MPMOVE32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id&
     &ps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
         endif
         if ((ierr /= 0) .or. (info(1) /= 0)) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
! calculate percent imbalance
         anpav = real(info(9))/real(nvpy*nvpz)
         if (anpav > 0.0) then
            pibal = max(real(info(2))-anpav,anpav-real(info(3)))/anpav
         endif
         ierr = info(1)
         end subroutine impdmove32
!
         subroutine impmoves32(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvp&
     &z,nbmax,idds,mblok,vt,ierr)
! multi-tasking particle manager
! non-uniform 2d partition boundaries in 3d code
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nbmax, idds, mblok, vt
         integer :: ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, idps, ntmax, nmt
         integer, dimension(ntasks) :: idtask
         integer, dimension(size(edges,1),size(edges,2),ntasks) :: jssp
         integer, dimension(size(npp)) :: npq
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(9) :: info
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         nmt = ntasks
! set maximum number of passes
         info(6) = 2; info(7) = 2
         th = 0.0
! check if size of buffers has changed
         if (szbuf.ne.nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call MPXMOVS32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ih&
     &ole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,i&
     &dps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
         else
            call MPMOVES32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ih&
     &ole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,i&
     &dps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
         endif
         if ((ierr /= 0) .or. (info(1) /= 0)) then
            call MP_END
            call PPABORT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine impmoves32
!
      end module mp32d
