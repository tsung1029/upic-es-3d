!-----------------------------------------------------------------------
!
      module globals
!
! Fortran90 interface to PIC library constants
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: april 26, 2005
!
      implicit none
      public
!
      integer, parameter :: LINEAR = 1, QUADRATIC = 2
      integer, parameter :: STANDARD = 1, LOOKAHEAD = 2, VECTOR = 3
      integer, parameter :: PERIODIC_3D = 1, DIRICHLET_3D = 2
      integer, parameter :: MIXED_3D = 3, VACUUM_3D = 4
!
      end module globals
