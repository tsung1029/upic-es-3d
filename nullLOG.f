c-----------------------------------------------------------------------
c Null library for a few miscellaneous routines from MacMPI
c update: april 9, 2003
c-----------------------------------------------------------------------
      subroutine LOGNAME(name)
c this subroutine records and displays user-defined label
c name = label to display
      implicit none
      character*(*) name
      return
      end
c-----------------------------------------------------------------------
      subroutine SET_MON(monval)
c this subroutine sets new monitor value and corresponding window
c monval = new monitor value
      implicit none
      integer monval
      return
      end
c-----------------------------------------------------------------------
      function GET_MON()
c this function gets current monitor value
      implicit none
      integer GET_MON
      GET_MON = 0
      return
      end
