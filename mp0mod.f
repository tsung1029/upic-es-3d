!-----------------------------------------------------------------------
!
      module mp0d
!
! Fortran90 interface for initializing Multi-tasking library MacMP.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 12, 2007
!
      implicit none
      private
      public :: mpinit, ncpus, ntasks
!
! ncpus = number of cpus found
! ntasks = number of additional tasks for threaded programming
      integer, save :: ncpus = 1, ntasks = 0
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MP_INIT(nproc)
         implicit none
         integer nproc
         end subroutine
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         function mpinit(sntasks)
! initialize for multiprocessing
         implicit none
         integer mpinit
         integer, optional :: sntasks
         call MP_INIT(ncpus)
         if (ncpus==0) then
            !write (2,*) 'MPLibrary not installed'
            !print*, 'MPLibrary not installed'
            ! flush(2) ! required to see anything on IBMs
            ntasks= 0
! return the number of processors on host computer
         else
            ntasks = ncpus - 1
! set specific number of tasks if requested
            if (present(sntasks)) then
               if (sntasks >= 0) ntasks = sntasks
            endif
         endif
         mpinit = ntasks
         end function
!
      end module mp0d
