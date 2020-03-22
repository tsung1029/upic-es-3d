!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     system module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/source/os-sys-macosx.f90 $
! $Id: os-sys-macosx.f90 57 2006-02-02 19:29:00Z zamb $
!

module m_system

implicit none


! give general access to all information in this module
public


!==============================================================
! parameters essential at compile time that depend on the system
! the program is running on
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Constants defining datatypes
integer, public, parameter :: p_single = kind(1.0e0)
integer, public, parameter :: p_double = kind(1.0d0)
integer, public, parameter :: p_byte   = selected_int_kind(2)
integer, public, parameter :: p_int64  = selected_int_kind(10)

! directory seperator for this system
character, parameter :: p_dir_sep  = '/'      ! unix like systems


! log file id
integer, parameter, public :: file_id_tem  =  10
integer, parameter, private :: file_id_msg  =  20


integer, parameter :: p_stderr = 0
integer, parameter :: p_stdin  = 5
integer, parameter :: p_stdout = 6

! max size of filename
integer, parameter :: p_max_filename_len = 256


integer, parameter :: p_err_assert          = -21 ! assertion 

!=================================================================================================
! High resolution timer interfaces
! these routines are defined in os-sys-multi-c.c
!-------------------------------------------------------------------------------------------------

! number of seconds elapsed since a fixed time in the past
real(p_double), external :: timer_cpu_seconds 

! minimal difference between successive high resolution timer
! calls (note: this is usual much less than the actual resolution
! due to the conversion to seconds)
real(p_double), external :: timer_resolution

! number of ticks since a fixed time in the past
integer(p_int64), external :: timer_ticks

! Convert tick interval to seconds
real(p_double), external :: timer_interval_seconds

!unix time in seconds
integer(p_int64), external :: unix_time

! buffer variables for debug/log routines
! used by the DEBUG, LOG, ERROR and WARNING macros
character(len = 1024) :: dbg_buf__, log_buf__, err_buf__, wrn_buf__
 
! Structure to hold command line / environment options
type :: t_options
  logical :: test    = .false.
  logical :: restart = .false.
	integer(p_int64) :: maxtime = 0
  character( len = p_max_filename_len ) :: input_file = '' 
  character( len = p_max_filename_len ) :: work_dir  = '' 
  
  integer :: random_seed = 0
  integer :: algorithm   = 0
  real( p_double ) :: omega_p0
  real( p_double ) ::  n0
  real( p_double ) ::  gamma  
end type t_options
 
 
interface isnan
  module procedure isnan_single
end interface

interface isinf
  module procedure isinf_single
end interface

interface ffloor
  module procedure fastfloor_single
end interface

interface signbit
  module procedure signbit_single
end interface

interface file_exists
   module procedure file_exists
end interface

interface getopt
  module procedure getopt
end interface getopt

contains

!---------------------------------------------------------------------------------------------------
function file_exists( filename )
!---------------------------------------------------------------------------------------------------
! Check is supplied filename exists and is readable
! - This is implemented in pure fortran, a posix version using fstat() or access() is also 
!   possible
!---------------------------------------------------------------------------------------------------

   implicit none
   
   logical :: file_exists
   character(len = *), intent(in) :: filename
   
   integer :: ierr
   
   open( file_id_tem, file = trim(filename), status = 'old', iostat = ierr, action = 'read' )
   if ( ierr == 0 ) then
	 close( file_id_tem )
	 file_exists = .true.
   else
	 file_exists = .false.
   endif

end function file_exists
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine link( target, name, ierr_out )
!---------------------------------------------------------------------------------------------------
! create a link "name" to the file "target"
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  character(len=*), intent(in) :: target, name
  integer, intent(out), optional :: ierr_out

  ! local variables
  integer :: ierr

  ! remove old link
  call remove( name, ierr )

  ! create new one
  call symlink( target, name, ierr )
  
  ! return result if necessary
  if (present(ierr_out)) then
	ierr_out = ierr
  endif

end subroutine link
!---------------------------------------------------------------------------------------------------


!************************************************************************
! Fortran wrappers for os-sys-multi-c.c
!************************************************************************

!---------------------------------------------------------------------------------------------------
subroutine getenv(name, value, ierr)
!---------------------------------------------------------------------------------------------------
! gets the value of an environment variable
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: value
  integer, intent(out) :: ierr

  character(len=len_trim(name)+1) :: lname

  ! Make shure this is a proper C string
  lname = trim(name)//char(0)

  value = ""
  call getenv_f(lname, value, ierr)
  call cleanup_cstring(value)

end subroutine getenv
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine gethostname(hostname, ierr)
!---------------------------------------------------------------------------------------------------
!        gets the hostname
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(out) :: hostname
  integer, intent(out) :: ierr

  character(len=256) :: lhostname

  lhostname = " "
  call gethostname_f(lhostname, ierr)
  call cleanup_cstring(lhostname)
  hostname = trim(lhostname)

end subroutine gethostname
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function hostname()
!---------------------------------------------------------------------------------------------------
!        gets the hostname
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=256) :: hostname

  integer:: ierr
  character(len=256) :: lhostname

  lhostname = " "
  ierr = 0
  call gethostname_f(lhostname, ierr)
  call cleanup_cstring(lhostname)

  hostname = trim(lhostname)

end function hostname
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------- 
subroutine umask( numask )
!--------------------------------------------------- 
! set the current umask
!--------------------------------------------------- 
   implicit none
   
   integer, intent(in) :: numask

   call umask_f( numask )

end subroutine umask
!--------------------------------------------------- 

!--------------------------------------------------- 
subroutine mkdir(path, ierr)
!---------------------------------------------------------------------------------------------------
!        create a directory
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr
  
  character(len=len_trim(path)+1) :: lpath

  ! Make shure this is a proper C string
  lpath = trim(path)//char(0)
  
  call mkdir_f(lpath, ierr)

end subroutine mkdir
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------- 
subroutine chdir(path, ierr)
!---------------------------------------------------------------------------------------------------
!        change current working directory
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr
  
  character(len=len_trim(path)+1) :: lpath

  ! Make shure this is a proper C string
  lpath = trim(path)//char(0)
  
  call chdir_f(lpath, ierr)


end subroutine chdir
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------- 
subroutine getcwd(path, ierr)
!---------------------------------------------------------------------------------------------------
!        get current working directory
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(out) :: path
  integer, intent(out) :: ierr
  
  
  path = " "
  call getcwd_f(path, len(path), ierr)
  call cleanup_cstring(path)


end subroutine getcwd
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function strerror(ierr)
!---------------------------------------------------------------------------------------------------
!        get the name of a system error
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: ierr
  character(len = 256) :: strerror

  strerror = ""
  call strerror_f(ierr, strerror)
  call cleanup_cstring(strerror)

end function strerror
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine symlink(targetf, linkf, ierr)
!---------------------------------------------------------------------------------------------------
!        create a symbolic link
!---------------------------------------------------------------------------------------------------

implicit none

  character(len=*), intent(in) :: targetf, linkf
  integer, intent(out) :: ierr

  character(len=len_trim(targetf)+1) :: ltargetf
  character(len=len_trim(linkf)+1) :: llinkf
  
  ! Make shure this are proper C strings
  
  ltargetf = trim(targetf)//char(0)
  llinkf = trim(linkf)//char(0)

  call symlink_f(ltargetf, llinkf, ierr) 

 end subroutine symlink
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine remove(path, ierr)
!---------------------------------------------------------------------------------------------------
!        remove a file
!---------------------------------------------------------------------------------------------------

implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr
  character(len=len_trim(path)+1) :: lpath

  ! Make shure this is a proper C string
  lpath = trim(path)//char(0)

  call remove_f(trim(lpath), ierr)

 end subroutine remove
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------

! Use Fortran 2003 extensions to access command line arguments

!---------------------------------------------------------------------------------------------------
function getargc( )
!---------------------------------------------------------------------------------------------------
! Get number of command line arguments
!---------------------------------------------------------------------------------------------------

  implicit none

  integer :: getargc

  getargc = command_argument_count()

end function getargc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function getargv( i )
!---------------------------------------------------------------------------------------------------
! Get requested command line argument
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len = p_max_filename_len) :: getargv
  integer, intent(in) :: i
  
  integer :: argc  
  argc = getargc()
  
  if ( i < 1 .or. i > argc ) then
    print *, '(*error*) Trying to read argument ', i, ' of ', argc
    getargv = ""
    return
  endif
  
  call get_command_argument( i, getargv )
  
end function getargv
!---------------------------------------------------------------------------------------------------


! optstring	- string specifying options
! optidx    - command line argument index to process (start with 1)
! opt		- next known option character in optstring
! optarg	- option argument, if it is anticipated

! opt is set to '\0' when finished processing arguments

subroutine getopt( optstring, optidx, opt, optarg  )

  implicit none
  
  character(len=*), intent(in) :: optstring
  integer, intent(inout) :: optidx
  character, intent(out) :: opt
  character(len = *), intent(out) :: optarg

  integer :: argc, idx
  character(len = p_max_filename_len) :: argv
  logical :: valid

  argc = getargc()
  
  ! If over the last command line argument return '\0'
  if ( optidx > argc ) then
    opt = achar(0)
    return
  endif

  ! get argument value
  argv = getargv( optidx )

  if ( argv(1:1) == '-' ) then
    
    ! termination with '--'
    if ( len_trim(argv) == 2 .and. argv(2:2) == '-' ) then
	  opt = achar(0)
	  return
    endif
    
    ! look for options in optstring
    idx = 1
    valid = .false.
    do while ( idx <= len_trim(optstring) ) 
      ! found option
      if ( optstring(idx:idx) == argv(2:2) ) then
        ! check for option argument
        if ( idx + 1 <= len_trim(optstring) ) then
		  if ( optstring(idx+1:idx+1) == ':' ) then
			optidx = optidx + 1
			if ( optidx > argc ) then
			  exit
			else
			  optarg = getargv( optidx )
			endif
		  endif
        endif
        ! set opt to option and finish search
        opt = argv(2:2)
        valid = .true.
        exit
      endif
      idx = idx + 1
    enddo
    
    if ( .not. valid ) then
      opt = '?'
	  return
    endif
  
    optidx = optidx + 1
    
  else
    opt = achar(0)
  endif
  
  
end subroutine getopt




!---------------------------------------------------------------------------------------------------
function isnan_single( x )
!---------------------------------------------------------------------------------------------------
! Check if supplied value is a NAN
!---------------------------------------------------------------------------------------------------

  implicit none

  real( p_single ), intent(in) :: x
  logical :: isnan_single

  integer :: isnan_res

  call isnan_single_f( x , isnan_res )
  
  isnan_single = ( isnan_res == 1 )

end function isnan_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isinf_single( x )
!---------------------------------------------------------------------------------------------------
!  Check if number is infinity
!---------------------------------------------------------------------------------------------------

  implicit none

  real( p_single ), intent(in) :: x
  logical :: isinf_single

  integer :: isinf_res

  call isinf_f( x , isinf_res )
  
  isinf_single = ( isinf_res .eq. 1 )

end function isinf_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! For their usability range these are actually (MUCH) faster than the native floor implementation
! (tested in the range [-2,2])
!---------------------------------------------------------------------------------------------------
function fastfloor_single(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_single), intent(in) :: x
  
  integer :: fastfloor_single

  fastfloor_single = int(x)
  if ( x < 0 ) fastfloor_single = fastfloor_single - 1

end function fastfloor_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function signbit_single(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_single), intent(in) :: x
  
  integer :: signbit_single
  
  if ( x < 0 ) then
    signbit_single = 1
  else
    signbit_single = 0
  endif
  
end function signbit_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine cleanup_cstring( str )
!---------------------------------------------------------------------------------------------------
! Cleans up a C string (null terminated) so that it is a proper fortran string
!---------------------------------------------------------------------------------------------------

  implicit none

	character(len=*), intent(inout) :: str
	integer :: pos
	! a C string will be terminated by a null character
	
	pos = index( str, char(0) )
	
	if (pos > 0) then
	  str = str( 1:pos-1)
	endif


end subroutine cleanup_cstring
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine linebuf_stdio()
!---------------------------------------------------------------------------------------------------
! Force stdout and stderr to be line buffered
!---------------------------------------------------------------------------------------------------

  implicit none

  call linebuf_stdio_f()

 end subroutine linebuf_stdio

end module m_system
