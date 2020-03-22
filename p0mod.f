!-----------------------------------------------------------------------
!
      module p0d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library p0lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: may 2, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PPINIT, PPID, PPEXIT, HARTBEAT, PWRITE1, PREAD1
      public :: get_funit, pwtimer, wtimer
      public :: plsum, plmax, plscan, plbcast
      public :: writef, readf, wrdata, rddata
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PPINIT(idproc,id0,nvp)
         implicit none
         integer :: idproc, id0, nvp
         end subroutine
      end interface
      interface
         subroutine PPID(idproc,id0,nvp)
         implicit none
         integer :: idproc, id0, nvp
         end subroutine
      end interface
      interface
         subroutine PPEXIT
         end subroutine
      end interface
      interface
         subroutine HARTBEAT(msg,n)
         implicit none
         integer :: n
         double precision, dimension(n) :: msg
         end subroutine
      end interface
      interface
         subroutine PTIMERA(icntrl,time,dtime)
         implicit none
         integer :: icntrl
         real, dimension(2) :: time
         double precision :: dtime
         end subroutine
      end interface
      interface
         subroutine PWTIMERA(icntrl,time,dtime)
         implicit none
         integer :: icntrl
         real :: time
         double precision :: dtime
         end subroutine
      end interface
      interface
         subroutine PSUM(f,g,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         real, dimension(nxp,nblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PISUM(if,ig,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         integer, dimension(nxp,nblok) :: if, ig
         end subroutine
      end interface
      interface
         subroutine PCSUM(f,g,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         complex, dimension(nxp,nblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PDSUM(f,g,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         double precision, dimension(nxp,nblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PMAX(f,g,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         real, dimension(nxp,nblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PIMAX(if,ig,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         integer, dimension(nxp,nblok) :: if, ig
         end subroutine
      end interface
      interface
         subroutine PDMAX(f,g,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         double precision, dimension(nxp,nblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PSCAN(f,g,s,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         real, dimension(nxp,nblok) :: f, g, s
         end subroutine
      end interface
      interface
         subroutine PISCAN(if,ig,is,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         integer, dimension(nxp,nblok) :: if, ig, is
         end subroutine
      end interface
      interface
         subroutine PDSCAN(f,g,s,nxp,nblok)
         implicit none
         integer :: nxp, nblok
         double precision, dimension(nxp,nblok) :: f, g, s
         end subroutine
      end interface
      interface
         subroutine PBCAST(f,nxp)
         implicit none
         integer :: nxp
         real, dimension(nxp) :: f
         end subroutine
      end interface
      interface
         subroutine PBICAST(f,nxp)
         implicit none
         integer :: nxp
         integer, dimension(nxp) :: f
         end subroutine
      end interface
      interface
         subroutine PBDCAST(f,nxp)
         implicit none
         integer :: nxp
         double precision, dimension(nxp) :: f
         end subroutine
      end interface
      interface
         subroutine PWRITE1(f,nxp,iunit,nrec,name)
         implicit none
         integer :: nxp, iunit, nrec
         character(len=*) :: name
         real, dimension(nxp) :: f
         end subroutine
      end interface
      interface
         subroutine PREAD1(f,nxp,iunit,nrec,name)
         implicit none
         integer :: nxp, iunit, nrec
         character(len=*) :: name
         real, dimension(nxp) :: f
         end subroutine
      end interface
      interface
         subroutine PWRDATA(f,nvp,nxp,nblok,iunit)
         implicit none
         integer :: nvp, nxp, nblok, iunit
         real, dimension(nxp,nblok) :: f
         end subroutine
      end interface
      interface
         subroutine PRDDATA(f,nvp,nxp,nblok,iunit,ierror)
         implicit none
         integer :: nvp, nxp, nblok, iunit, ierror
         real, dimension(nxp,nblok) :: f
         end subroutine
      end interface
      interface
         subroutine PWRPART(part,npp,idimp,npmax,nblok,iunit)
         implicit none
         integer :: idimp, npmax, nblok, iunit
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PRDPART(part,npp,idimp,npmax,nblok,iunit,ierror)
         implicit none
         integer :: idimp, npmax, nblok, iunit, ierror
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface plsum
         module procedure ipsum
         module procedure ipisum
         module procedure ipcsum
      end interface
!
      interface plmax
         module procedure ipmax
         module procedure ipimax
      end interface
!
      interface plscan
         module procedure ipscan
         module procedure ipiscan
      end interface
!
      interface plbcast
         module procedure ipbcast
         module procedure ipbicast
      end interface
!
      interface writef
         module procedure ipwrite1
      end interface
!
      interface readf
         module procedure ipread1
      end interface
!
      interface wrdata
         module procedure ipwrdata0
         module procedure ipwrpart
      end interface
!
      interface rddata
         module procedure iprddata0
         module procedure iprdpart
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         function get_funit(start) result(funit)
! this function returns an unconnected fortran unit number,
! starting with unit = start.  returns -1 if none found
         integer, intent(in) :: start
         integer :: funit
! local data
         integer :: i
         logical :: connected
         funit = -1
! check connection status
         do i = start, 98
            inquire(unit=i,opened=connected)
            if (.not.connected) then
               funit = i
               exit
            endif
         enddo
         end function get_funit
!
         subroutine pwtimer(time,dtime,icntrl)
! parallel wall clock timer
         implicit none
         real, dimension(2), intent(out) :: time
         double precision, intent(inout) :: dtime
         integer, intent(in), optional :: icntrl
         integer :: ltime
         ltime = 1
         if (present(icntrl)) ltime = icntrl
         if (ltime==0) return
         call PTIMERA(ltime,time,dtime)
         end subroutine pwtimer
!
         subroutine wtimer(time,dtime,icntrl)
! local wall clock timer
         implicit none
         real, intent(out) :: time
         double precision, intent(inout) :: dtime
         integer, intent(in), optional :: icntrl
         integer :: ltime
         ltime = 1
         if (present(icntrl)) ltime = icntrl
         if (ltime==0) return
         call PWTIMERA(ltime,time,dtime)
         end subroutine wtimer
!
         subroutine ipsum(f)
! perform global sum of 1d real array
         implicit none
         real, dimension(:) :: f
         integer :: nxp, nblok
         real, dimension(size(f)) :: g
         nxp = size(f); nblok = 1
         call PSUM(f,g,nxp,nblok)
         end subroutine ipsum
!
         subroutine ipisum(if)
! perform global sum of 1d integer array
         implicit none
         integer, dimension(:) :: if
         integer :: nxp, nblok
         integer, dimension(size(if)) :: ig
         nxp = size(if); nblok = 1
         call PISUM(if,ig,nxp,nblok)
         end subroutine ipisum
!
         subroutine ipcsum(f)
! perform global sum of 1d complex array
         implicit none
         complex, dimension(:) :: f
         integer :: nxp, nblok
         complex, dimension(size(f)) :: g
         nxp = size(f); nblok = 1
         call PCSUM(f,g,nxp,nblok)
         end subroutine ipcsum
!
         subroutine ipmax(f)
! perform global maximum of 1d real array
         implicit none
         real, dimension(:) :: f
         integer :: nxp, nblok
         real, dimension(size(f)) :: g
         nxp = size(f); nblok = 1
         call PMAX(f,g,nxp,nblok)
         end subroutine ipmax
!
         subroutine ipimax(if)
! perform global maximum of 1d integer array
         implicit none
         integer, dimension(:) :: if
         integer :: nxp, nblok
         integer, dimension(size(if)) :: ig
         nxp = size(if); nblok = 1
         call PIMAX(if,ig,nxp,nblok)
         end subroutine ipimax
!
         subroutine ipscan(f)
! perform parallel prefix sum reduction of 1d real array
         implicit none
         real, dimension(:) :: f
         integer :: nxp, nblok
         real, dimension(size(f)) :: g, s
         nxp = size(f); nblok = 1
         call PSCAN(f,g,s,nxp,nblok)
         end subroutine ipscan
!
         subroutine ipiscan(if)
! perform parallel prefix sum reduction of 1d integer array
         implicit none
         integer, dimension(:) :: if
         integer :: nxp, nblok
         integer, dimension(size(if)) :: ig, is
         nxp = size(if); nblok = 1
         call PISCAN(if,ig,is,nxp,nblok)
         end subroutine ipiscan
!
         subroutine ipbcast(f)
! broadcast 1d real array
         implicit none
         real, dimension(:) :: f
         integer :: nxp
         nxp = size(f)
         call PBCAST(f,nxp)
         end subroutine ipbcast
!
         subroutine ipbicast(if)
! broadcast 1d integer array
         implicit none
         integer, dimension(:) :: if
         integer :: nxp
         nxp = size(if)
         call PBICAST(if,nxp)
         end subroutine ipbicast
!
         subroutine ipwrite1(f,iunit,nrec,name)
! collects distributed real 2d scalar array and writes it to a
! direct access binary file, for uniform 1d partitions
         implicit none
         integer :: iunit, nrec
         real, dimension(:,:,:), pointer :: f
         character(len=*) :: name
         integer :: nxp
         nxp = size(f)
         call PWRITE1(f,nxp,iunit,nrec,name)
         end subroutine ipwrite1
!
         subroutine ipread1(f,iunit,nrec,name)
! reads distributed real 2d scalar array from a direct access
! binary file and distributes it, for uniform 1d partitions
         implicit none
         integer :: iunit, nrec
         real, dimension(:,:,:), pointer :: f
         character(len=*) :: name
         integer :: nxp
         nxp = size(f)
         call PREAD1(f,nxp,iunit,nrec,name)
         end subroutine ipread1
!
         subroutine ipwrdata0(f,nvp,iunit)
! collects distributed real 1d scalar data and writes it to a
! fortran unformatted sequential file, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit
         real, dimension(:,:), pointer :: f
         integer :: nxp, nblok
         nxp = size(f); nblok = 1
         call PWRDATA(f,nvp,nxp,nblok,iunit)
         end subroutine ipwrdata0
!
         subroutine iprddata0(f,nvp,iunit,ierror)
! reads real 1d scalar data from a fortran unformatted sequential file
! and distributes it, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit, ierror
         real, dimension(:,:), pointer :: f
         integer :: nxp, nblok
         nxp = size(f); nblok = 1
         call PRDDATA(f,nvp,nxp,nblok,iunit,ierror)
         end subroutine iprddata0
!
         subroutine ipwrpart(part,npp,iunit)
! reads particle data from a fortran unformatted sequential file and
! distributes it, for 1d partitions
         implicit none
         integer :: iunit
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         integer :: idimp, npmax, nblok
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PWRPART(part,npp,idimp,npmax,nblok,iunit)
         end subroutine ipwrpart
!
         subroutine iprdpart(part,npp,iunit,ierror)
! collects distributed particle data and writes it to a fortran
! unformatted sequential file, for 1d partitions
         implicit none
         integer :: iunit, ierror
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         integer :: idimp, npmax, nblok
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PRDPART(part,npp,idimp,npmax,nblok,iunit,ierror)
         end subroutine iprdpart
!
      end module p0d
