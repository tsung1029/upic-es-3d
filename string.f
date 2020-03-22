! general module for string operation
! Written by Chengkun Huang
! Dated: 09/02/05
      module m_string
      
      interface int2str 
        module procedure i2s
      end interface
      
      contains
!---------------------------------------------------

      subroutine i2s(int_in,string,ndigits)

      implicit none
      integer, intent(in) :: int_in, ndigits
      character(len=20), intent(inout) :: string
      
! local variables      
      integer  ::  izero, i, nd, m
      character(len=20)  :: chindx 
      
      m = 1 
      izero =  ichar('0')
      if (ndigits > 20) then 
         nd = 20
      else
         nd = ndigits
      endif
      chindx = ''
      do i = nd, 1, -1 
         m = 10**(i-1)
         chindx = trim(chindx) // char(  izero + mod( int_in/m , 10 ) ) 
      enddo 
      string = trim(chindx)

      end subroutine i2s     
      
      end module m_string
