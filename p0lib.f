c basic parallel PIC library for MPI communications
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: june 7, 2008
c update: 25 August 2010: removed status='scratch' from open
c statements because BGL doesn't like it
c-----------------------------------------------------------------------
      function NDIAN()
c this function determines whether number format is big or little endian
c assumes ieee format for numbers
c ndian = (0,1) = architecture is (little-endian,big-endian)
      implicit none
      integer NDIAN
      integer*4 i
      real*8 d
      dimension i(2)
      equivalence (i,d)
      data d /1.0d0/
      save
c  big endian
      if (i(2).eq.0) then
         NDIAN = 1
c little endian
      else if (i(1).eq.0) then
         NDIAN = 0
c error
      else
         NDIAN = -1
      endif
      return
      end
c-----------------------------------------------------------------------
      function NDPREC()
c this subroutine determines whether default reals are double precision
c ndprec = (0,1) = default reals are (normal,double-precision)
      implicit none
      integer NDPREC
      real small, prec, vresult
      data small /1.0e-12/
      save
      prec = 1.0 + small
      if (vresult(prec).gt.1.0) then
         NDPREC = 1
      else
         NDPREC = 0
      endif
      return
      end
c-----------------------------------------------------------------------
      function IDPREC()
c this subroutine determines whether default integers are double
c precision
c idprec = (0,1) = default integers are (normal,double-precision)
      implicit none
      integer IDPREC
      integer ibig, iprec, iresult
      data ibig /2147483647/
      save
      iprec = ibig + 1
      if (iresult(iprec).gt.0.0) then
         IDPREC = 1
      else
         IDPREC = 0
      endif
      return
      end
c-----------------------------------------------------------------------
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
      return
      end
c-----------------------------------------------------------------------
      function iresult(iprec)
      implicit none
      integer iprec, iresult
      iresult = iprec
      return
      end
c-----------------------------------------------------------------------
      subroutine PPINIT(idproc,id0,nvp)
c this subroutine initializes parallel processing
c creates a communicator with number of processors equal to a power of 2
c output: idproc, id0, nvp
c idproc = processor id in lgrp communicator
c id0 = processor id in MPI_COMM_WORLD
c nvp = number of real or virtual processors obtained
      implicit none
      integer idproc, id0, nvp
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
c mint = default datatype for integers
c mcplx = default datatype for complex type
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror, ndprec, i, k
      logical flag
      real small, prec, vresult
      save /PPARMS/
      data small /1.0e-12/
      prec = 1.0 + small
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (vresult(prec).gt.1.0) then
         ndprec = 1
      else
         ndprec = 0
      endif
c this segment is used for shared memory computers
c     nproc = nvp
c     id0 = 0
c     idproc = 0
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         stop
      endif
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
c initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror.ne.0) stop
c already initialized
      else
         call PPID(idproc,id0,nvp)
         return
      endif
      lworld = MPI_COMM_WORLD
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,id0,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nproc,ierror)
c set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
c single precision
      if (ndprec.eq.0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
c double precision
      else
c        mint = MPI_INTEGER8
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
c set nvp to largest power of 2 contained in nproc
      nvp = 1
      idproc = id0
   10 i = nvp + nvp
      if (i.le.nproc) then
         nvp = i
         go to 10
      endif
c check if nproc is a power of 2
      if (nproc.gt.nvp) then
         idproc = idproc - 1
      endif
c create communicator which contains power of 2 number of processors
c exclude processor 0 in MPI_COMM_WORLD if nproc is not a power of 2
c this processor can be used as a diagnostic node
      if ((idproc.ge.0).and.(idproc.lt.nvp)) then
         k = 2
      else
         idproc = MPI_PROC_NULL
         k = 1
      endif
      call MPI_COMM_SPLIT(lworld,k,1,lgrp,ierror)
c requested number of processors not obtained
      if (ierror.ne.0) then
         write (2,*) ' MPI_COMM_SPLIT error: nvp, nproc=', nvp, nproc
         call PPEXIT
         stop
      endif
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,k,ierror)
      if (idproc.ge.0) then
         idproc = k
      else
         idproc = -(k + 1)
      endif
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine PPID(idproc,id0,nvp)
c this subroutine gets processor ids and number of processors
c output: idproc, id0, nvp
c idproc = processor id in lgrp communicator
c id0 = processor id in MPI_COMM_WORLD
c nvp = number of real or virtual processors obtained
c output: all
      implicit none
      integer idproc, id0, nvp
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx,  mdouble, lworld
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,id0,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierror)
c special case if diagnostic nodes are present
      if (nvp.gt.nproc) then
         idproc = id0 - 1
         if (idproc.ge.nproc) idproc = nproc - idproc - 2
      else
         idproc = id0
      endif
      nvp = nproc
      return
      end
c-----------------------------------------------------------------------
      subroutine PPINIT0(idproc,nvp)
c this subroutine initializes parallel processing
c input: nvp, output: idproc
c idproc = processor id in lgrp communicator
c nvp = number of real or virtual processors requested
      implicit none
      integer idproc, nvp
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
c mint = default datatype for integers
c mcplx = default datatype for complex type
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror, ndprec
      logical flag
      real small, prec, vresult
      save /PPARMS/
      data small /1.0e-12/
      prec = 1.0 + small
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (vresult(prec).gt.1.0) then
         ndprec = 1
      else
         ndprec = 0
      endif
c this segment is used for shared memory computers
c     nproc = nvp
c     idproc = 0
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         stop
      endif
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
c initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror.ne.0) stop
      endif
      lworld = MPI_COMM_WORLD
      lgrp = lworld
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
c set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
c single precision
      if (ndprec.eq.0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
c double precision
      else
c        mint = MPI_INTEGER8
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
c requested number of processors not obtained
      if (nproc.ne.nvp) then
         write (2,*) ' processor number error: nvp, nproc=', nvp, nproc
         call ppexit
         stop
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPINIT1(idproc,id0,nvp)
c this subroutine initializes parallel processing
c output: idproc, id0, nvp
c idproc = processor id in lgrp communicator
c id0 = processor id in MPI_COMM_WORLD
c nvp = number of real or virtual processors obtained
      implicit none
      integer idproc, id0, nvp
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
c mint = default datatype for integers
c mcplx = default datatype for complex type
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror, ndprec, idprec
      integer ibig, iprec, iresult
      logical flag
      real small, prec, vresult
      save /PPARMS/
      data small /1.0e-12/
      data ibig /2147483647/
      prec = 1.0 + small
      iprec = ibig + 1
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (vresult(prec).gt.1.0) then
         ndprec = 1
      else
         ndprec = 0
      endif
c idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
      if (iresult(iprec).gt.0) then
         idprec = 1
      else
         idprec = 0
      endif
c this segment is used for shared memory computers
c     nproc = nvp
c     id0 = 0
c     idproc = 0
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         stop
      endif
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
c initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror.ne.0) stop
      endif
      lworld = MPI_COMM_WORLD
      lgrp = lworld
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
c set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
c single precision real
      if (ndprec.eq.0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
c double precision real
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
c single precision integer
c     if (idprec.eq.0) then
c        mint = MPI_INTEGER
c double precision integer
c     else
c        mint = MPI_INTEGER8
c     endif
      nvp = nproc
      id0 = idproc
      return
      end
c-----------------------------------------------------------------------
      subroutine PPEXIT
c this subroutine terminates parallel processing
      implicit none
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      integer ierror
      logical flag
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
c synchronize processes
         call MPI_BARRIER(lworld,ierror)
c terminate MPI execution environment
         call MPI_FINALIZE(ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPABORT
c this subroutine aborts parallel processing
      implicit none
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      integer errorcode, ierror
      logical flag
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
         errorcode = 1
c terminate MPI execution environment
         call MPI_ABORT(lworld,errorcode,ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine HARTBEAT(msg,n)
c participating nodes send a message msg to non-participating nodes
c msg = message to send
c n = length of message
c input: msg, n
c output: msg
      implicit none
      integer n
      double precision msg
      dimension msg(n)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, ierror
      dimension istatus(lstat)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierror)
c quit if communicators are the same size
      if (nvp.eq.nproc) return
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierror)
c send message
      if (idproc.eq.0) then
         call MPI_RECV(msg,n,mdouble,1,999,lworld,istatus,ierror)
      elseif (idproc.eq.1) then
         call MPI_SEND(msg,n,mdouble,0,999,lworld,ierror)
      endif
c broadcast message to non-participating nodes
      if (nvp.gt.(2*nproc)) then
         call MPI_BCAST(msg,n,mdouble,0,lgrp,ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PTIMERA(icntrl,time,dtime)
c this subroutine performs parallel wall clock timing
c input: icntrl, dtime
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c time = maximum/minimum elapsed time in seconds
c dtime = current time
c written for mpi
      implicit none
      integer icntrl
      real time
      double precision dtime
      dimension time(2)
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr
      real nclock, ttime
      double precision jclock
      dimension ttime(2)
c initialize clock
      if (icntrl.eq.(-1)) then
         call MPI_BARRIER(lgrp,ierr)
         dtime = MPI_WTIME()
c read clock and write time difference from last clock initialization
      else if (icntrl.eq.1) then
         jclock = dtime
         dtime = MPI_WTIME()
         nclock = real(dtime - jclock)
         ttime(1) = nclock
         ttime(2) = -nclock
         call MPI_ALLREDUCE(ttime,time,2,mreal,MPI_MAX,lgrp,ierr)
         time(2) = -time(2)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
c this subroutine performs local wall clock timing
c input: icntrl, dtime
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c time = elapsed time in seconds
c dtime = current time
c written for mpi
      implicit none
      integer icntrl
      real time
      double precision dtime
c local data
      double precision jclock
      double precision MPI_WTIME
      external MPI_WTIME
c initialize clock
      if (icntrl.eq.(-1)) then
         dtime = MPI_WTIME()
c read clock and write time difference from last clock initialization
      else if (icntrl.eq.1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PSUM(f,g,nxp,nblok)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output data
c g = scratch array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      real f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c perform sum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = f(j,k) + g(j,k)
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PISUM(if,ig,nxp,nblok)
c this subroutine performs a parallel sum of a vector, that is:
c if(j,k) = sum over k of if(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c if = input and output integer data
c ig = scratch integer array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      integer if, ig, nxp, nblok
      dimension if(nxp,nblok), ig(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mint = default datatype for integers
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, nsid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        ig(j,k) = if(j,kb+kxs)
c     else
c        ig(j,k) = if(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_ISEND(if,nxp,mint,kb+kxs-1,l+nxp,lgrp,nsid,ierr)
         call MPI_RECV(ig,nxp,mint,kb+kxs-1,l+nxp,lgrp,istatus,ierr)
      else
         call MPI_ISEND(if,nxp,mint,kb-kxs-1,l+nxp,lgrp,nsid,ierr)
         call MPI_RECV(ig,nxp,mint,kb-kxs-1,l+nxp,lgrp,istatus,ierr)
      endif
      call MPI_WAIT(nsid,istatus,ierr)
   30 continue
c perform sum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      if(j,k) = if(j,k) + ig(j,k)
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PCSUM(f,g,nxp,nblok)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output complex data
c g = scratch complex array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      complex f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mcplx,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mcplx,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mcplx,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mcplx,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c perform sum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = f(j,k) + g(j,k)
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PDSUM(f,g,nxp,nblok)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output double precision data
c g = scratch double precision array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      double precision f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mdouble,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mdouble,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mdouble,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mdouble,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c perform sum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = f(j,k) + g(j,k)
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PMAX(f,g,nxp,nblok)
c this subroutine finds parallel maximum for each element of a vector
c that is, f(j,k) = maximum as a function of k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output data
c g = scratch array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      real f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c find maximum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = amax1(f(j,k),g(j,k))
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PIMAX(if,ig,nxp,nblok)
c this subroutine finds parallel maximum for each element of a vector
c that is, if(j,k) = maximum as a function of k of if(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c if = input and output integer data
c ig = scratch integer array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      integer if, ig
      integer nxp, nblok
      dimension if(nxp,nblok), ig(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mint = default datatype for integers
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        ig(j,k) = if(j,kb+kxs)
c     else
c        ig(j,k) = if(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(ig,nxp,mint,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(if,nxp,mint,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(ig,nxp,mint,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(if,nxp,mint,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c find maximum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      if(j,k) = max0(if(j,k),ig(j,k))
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PDMAX(f,g,nxp,nblok)
c this subroutine finds parallel maximum for each element of a vector
c that is, f(j,k) = maximum as a function of k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output double precision data
c g = scratch double precision array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      double precision f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mdouble,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mdouble,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mdouble,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mdouble,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c find maximum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = dmax1(f(j,k),g(j,k))
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PSCAN(f,g,s,nxp,nblok)
c this subroutine performs a parallel prefix reduction of a vector,
c that is: f(j,k) = sum over k of f(j,k), where the sum is over k values
c less than idproc.
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output data
c g, s = scratch arrays
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      real f, g, s
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok), s(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c initialize global sum
      do 20 k = 1, nblok
      do 10 j = 1, nxp
      s(j,k) = f(j,k)
   10 continue
   20 continue
c main iteration loop
   30 if (kxs.ge.nproc) go to 90
c shift data
      do 60 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 40 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = s(j,kb+kxs)
c     else
c        g(j,k) = s(j,kb-kxs)
c     endif
c  40 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(s,nxp,mreal,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(s,nxp,mreal,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
c perform prefix scan
      if (lb.ne.0) then
         do 50 j = 1, nxp
         f(j,k) = f(j,k) + g(j,k)
   50    continue
      endif
   60 continue
c perform sum
      do 80 k = 1, nblok
      do 70 j = 1, nxp
      s(j,k) = s(j,k) + g(j,k)
   70 continue
   80 continue
      l = l + 1
      kxs = kxs + kxs
      go to 30
   90 return
      end
c-----------------------------------------------------------------------
      subroutine PISCAN(if,ig,is,nxp,nblok)
c this subroutine performs a parallel prefix reduction of a vector,
c that is: if(j,k) = sum over k of if(j,k), where the sum is over k
c values less than idproc.
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c if = input and output integer data
c ig, is = scratch integer arrays
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      integer if, ig, is
      integer nxp, nblok
      dimension if(nxp,nblok), ig(nxp,nblok), is(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mint = default datatype for integers
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c initialize global sum
      do 20 k = 1, nblok
      do 10 j = 1, nxp
      is(j,k) = if(j,k)
   10 continue
   20 continue
c main iteration loop
   30 if (kxs.ge.nproc) go to 90
c shift data
      do 60 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 40 j = 1, nxp
c     if (lb.eq.0) then
c        ig(j,k) = is(j,kb+kxs)
c     else
c        ig(j,k) = is(j,kb-kxs)
c     endif
c  40 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(ig,nxp,mint,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(is,nxp,mint,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(ig,nxp,mint,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(is,nxp,mint,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
c perform prefix scan
      if (lb.ne.0) then
         do 50 j = 1, nxp
         if(j,k) = if(j,k) + ig(j,k)
   50    continue
      endif
   60 continue
c perform sum
      do 80 k = 1, nblok
      do 70 j = 1, nxp
      is(j,k) = is(j,k) + ig(j,k)
   70 continue
   80 continue
      l = l + 1
      kxs = kxs + kxs
      go to 30
   90 return
      end
c-----------------------------------------------------------------------
      subroutine PDSCAN(f,g,s,nxp,nblok)
c this subroutine performs a parallel prefix reduction of a vector,
c that is: f(j,k) = sum over k of f(j,k), where the sum is over k values
c less than idproc.
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output double precision data
c g, s = scratch double precision arrays
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      double precision f, g, s
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok), s(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c initialize global sum
      do 20 k = 1, nblok
      do 10 j = 1, nxp
      s(j,k) = f(j,k)
   10 continue
   20 continue
c main iteration loop
   30 if (kxs.ge.nproc) go to 90
c shift data
      do 60 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 40 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = s(j,kb+kxs)
c     else
c        g(j,k) = s(j,kb-kxs)
c     endif
c  40 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mdouble,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(s,nxp,mdouble,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mdouble,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(s,nxp,mdouble,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
c perform prefix scan
      if (lb.ne.0) then
         do 50 j = 1, nxp
         f(j,k) = f(j,k) + g(j,k)
   50    continue
      endif
   60 continue
c perform sum
      do 80 k = 1, nblok
      do 70 j = 1, nxp
      s(j,k) = s(j,k) + g(j,k)
   70 continue
   80 continue
      l = l + 1
      kxs = kxs + kxs
      go to 30
   90 return
      end
c-----------------------------------------------------------------------
      subroutine PBCAST(f,nxp)
c this subroutine broadcasts real data f
c f = data to be broadcast
c nxp = size of data f
c input: f, nxp
c output: f
      implicit none
      integer nxp
      real f
      dimension f(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, np, ioff, i, id, ierr
      dimension istatus(lstat)
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
         endif
c first send data to remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_SEND(f,nxp,mreal,id,97,lworld,ierr)
   10    continue
c then send data to node 0
         if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nxp,mreal,id,97,lworld,ierr)
         endif
c other nodes receive data from node 0
      elseif (idproc.le.nproc) then
         call MPI_RECV(f,nxp,mreal,0,97,lworld,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PBICAST(f,nxp)
c this subroutine broadcasts integer data f
c f = data to be broadcast
c nxp = size of data f
c input: f, nxp
c output: f
      implicit none
      integer nxp
      integer f
      dimension f(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mint = default datatype for integers
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, np, ioff, i, id, ierr
      dimension istatus(lstat)
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
         endif
c first send data to remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_SEND(f,nxp,mint,id,96,lworld,ierr)
   10    continue
c then send data to node 0
         if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nxp,mint,id,96,lworld,ierr)
         endif
c other nodes receive data from node 0
      elseif (idproc.le.nproc) then
         call MPI_RECV(f,nxp,mint,0,96,lworld,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PBDCAST(f,nxp)
c this subroutine broadcasts double precision data f
c f = data to be broadcast
c nxp = size of data f
c input: f, nxp
c output: f
      implicit none
      integer nxp
      double precision f
      dimension f(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, np, ioff, i, id, ierr
      dimension istatus(lstat)
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
         endif
c first send data to remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_SEND(f,nxp,mdouble,id,97,lworld,ierr)
   10    continue
c then send data to node 0
         if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nxp,mdouble,id,97,lworld,ierr)
         endif
c other nodes receive data from node 0
      elseif (idproc.le.nproc) then
         call MPI_RECV(f,nxp,mdouble,0,97,lworld,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine P0COPY(f,g,nxp)
c this subroutine copies real data f on node 0 to g on other nodes
c f/g = data to be sent/received
c nxp = size of data f
c input: f, nxp
c output: g
      implicit none
      integer nxp
      real f, g
      dimension f(nxp), g(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, idproc, i, id, ierr
      dimension istatus(lstat)
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
         do 10 i = 2, nproc
            id = i - 1
            call MPI_SEND(f,nxp,mreal,id,95,lgrp,ierr)
   10    continue
c other nodes receive data from node 0
      else
         call MPI_RECV(g,nxp,mreal,0,95,lgrp,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PWRITE1(f,nxp,iunit,nrec,name)
c this subroutine collects distributed real data f and writes to a 
c direct access binary file
c f = input data to be written, modified on node 0
c nxp = size of data f
c iunit = fortran unit number
c nrec = current record number for write (if negative, open file with
c recl=-nren)
c name = file name (used only if nren < 0)
c input: f, nxp, iunit, nrec, fname
c output: nrec
      implicit none
      integer nxp, iunit, nrec
      real f
      character*(*) name
      dimension f(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, ierr
      dimension istatus(lstat)
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        lrec = -nrec
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='unknown')
c        nrec = 1
c     endif
c     write (unit=iunit,rec=nrec) f
c     nrec = nrec + 1
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            lrec = -nrec
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='unknown')
            nrec = 1
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
            id = 1
            call MPI_RECV(f,nxp,mreal,id,99,lworld,istatus,ierr)
         endif
c first write data for node 0
         nrec0 = nrec
         write (unit=iunit,rec=nrec) f
         nrec = nrec + 1
c then write data from remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_RECV(f,nxp,mreal,id,99,lworld,istatus,ierr)
            write (unit=iunit,rec=nrec) f
            nrec = nrec + 1
   10    continue
c read data back for node 0
         read (unit=iunit,rec=nrec0) f
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         call MPI_SEND(f,nxp,mreal,0,99,lworld,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PREAD1(f,nxp,iunit,nrec,name)
c this subroutine reads real data f from a direct access binary file
c and distributes it
c f = output data to be read
c nxp = size of data f
c iunit = fortran unit number
c nrec = current record number for read (if negative, open file with
c recl=-nren)
c name = file name (used only if nren < 0)
c input: nxp, iunit, nrec, fname
c output: f, nrec
      implicit none
      integer nxp, iunit, nrec
      real f
      character*(*) name
      dimension f(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lrec, nvp, idproc, np, ioff, nrec0, i, id, ierr
      dimension istatus(lstat)
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        lrec = -nrec
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c        nrec = 1
c     endif
c     read (unit=iunit,rec=nrec) f
c     nrec = nrec + 1
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            lrec = -nrec
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
            nrec = 1
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
         endif
c first read data for remaining nodes
         nrec0 = nrec
         nrec = nrec + 1
         do 10 i = 2, np
            read (unit=iunit,rec=nrec) f
            nrec = nrec + 1
            id = i - ioff
            call MPI_SEND(f,nxp,mreal,id,98,lworld,ierr)
   10    continue
c then read data from node 0
         read (unit=iunit,rec=nrec0) f
         if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nxp,mreal,id,98,lworld,ierr)
         endif
c other nodes receive data from node 0
      elseif (idproc.le.(nproc+1)) then
         call MPI_RECV(f,nxp,mreal,0,98,lworld,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PWRDATA(f,nvp,nxp,nblok,iunit)
c this subroutine collects distributed real data f
c and writes it to a fortran unformatted sequential file
c f = data to be written per processor
c nvp = number of real or virtual processors requested
c nxp = size of data f
c nblok = number of particle partitions
c iunit = fortran unit number
c input: all
c output: none
      implicit none
      integer nvp, nxp, nblok, iunit
      real f
      dimension f(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, j, l, idproc, lvp, joff, id, ierr
      integer istatus
      dimension istatus(lstat)
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nproc
c     do 10 l = 1, nblok
c     write (iunit) (f(j,l),j=1,nxp)
c  10 continue
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,lvp,ierr)
      do 20 l = 1, nblok
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
c special diagnostic node present
         if (joff.eq.0) then
            id = 1
            call MPI_RECV(f,nxp,mreal,id,107,lworld,istatus,ierr)
         endif
c first write data for node 0
         write (iunit) (f(j,l),j=1,nxp)
c save copy of current data
c          open(unit=99,form='unformatted',status='scratch')
         open(unit=99,form='unformatted')
         write (99) (f(j,l),j=1,nxp)
c then write data from remaining nodes
         do 10 i = 2, nvp
         id = i - joff
         call MPI_RECV(f,nxp,mreal,id,107,lworld,istatus,ierr)
         write (iunit) (f(j,l),j=1,nxp)
   10    continue
c read data back for node 0
         rewind 99
         read (99) (f(j,l),j=1,nxp)
         close (99)
c other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         call MPI_SEND(f,nxp,mreal,0,107,lworld,ierr)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PRDDATA(f,nvp,nxp,nblok,iunit,ierror)
c this subroutine reads data f from a fortran unformatted sequential
c file and distributes it
c f = data to be read
c nvp = number of real or virtual processors requested
c nxp = size of data f per processor
c nblok = number of particle partitions
c iunit = fortran unit number
c ierror = (0,1) = (no,yes) error condition exists
c input: all
c output: part
      implicit none
      integer nvp, nxp, nblok, iunit, ierror
      real f
      dimension f(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, j, l, idproc, lvp, joff, id
      integer ios, ierr
      integer istatus
      dimension istatus(lstat)
      ierror = 0
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nproc
c     do 10 l = 1, nblok
c     read (iunit) (f(j,l),j=1,nxp)
c  10 continue
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,lvp,ierr)
      do 20 l = 1, nblok
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
c first read data for node 0
         read (iunit,iostat=ios) (f(j,l),j=1,nxp)
         if (ios.ne.0) ierror = ierror + 1
c special diagnostic node present
         if (joff.eq.0) then
            id = 1
            call MPI_SEND(f,nxp,mreal,id,106,lworld,ierr)
         endif
c save copy of current data
c         open(unit=99,form='unformatted',status='scratch')
         open(unit=99,form='unformatted')
         write (99) (f(j,l),j=1,nxp)
c then send data to remaining nodes
         do 10 i = 2, nvp
         read (iunit,iostat=ios) (f(j,l),j=1,nxp)
         if (ios.ne.0) ierror = ierror + 1
         id = i - joff
         call MPI_SEND(f,nxp,mreal,id,106,lworld,ierr)
   10    continue
c read data back for node 0
         rewind 99
         read (99,iostat=ios) (f(j,l),j=1,nxp)
         if (ios.ne.0) ierror = ierror + 1
         close (99)
c other nodes receive data from node 0
      else if (idproc.le.(nvp-joff)) then
         call MPI_RECV(f,nxp,mreal,0,106,lworld,istatus,ierr)
      endif
   20 continue
c check for error condition
      call PISUM(ierror,ierr,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine PWRPART(part,npp,idimp,npmax,nblok,iunit)
c this subroutine collects distributed particle data part
c and writes it to a fortran unformatted sequential file
c part = particle data to be written
c npp(l) = number of particles in partition l
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c iunit = fortran unit number
c input: all
c output: none
      implicit none
      integer idimp, npmax, nblok, iunit
      real part
      integer npp
      dimension part(idimp,npmax,nblok), npp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, j, k, l, ndt, idproc, lvp, nvp, joff, id, nd, nd0, ierr
      integer istatus
      dimension istatus(lstat)
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nproc
c     do 10 l = 1, nblok
c     write (iunit) npp(l)
c     write (iunit) (part(:,j,l),j=1,nd)
c  10 continue
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,lvp,ierr)
      do 20 l = 1, nblok
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         ndt = idimp*npmax
c no special diagnostic node
         if (joff.eq.1) then
            nd = npp(l)
c special diagnostic node present
         else
            nvp = lvp - nproc
            id = 1
            call MPI_RECV(part,ndt,mreal,id,105,lworld,istatus,ierr)
c determine how many particles to write
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/idimp
         endif
c first write data for node 0
         write (iunit) nd
         write (iunit) ((part(j,k,l),j=1,idimp),k=1,nd)
c save copy of current data
         nd0 = nd
c        open(unit=99,form='unformatted',status='scratch')
         open(unit=99,form='unformatted')
         write (99) ((part(j,k,l),j=1,idimp),k=1,nd0)
c then write data from remaining nodes
         do 10 i = 2, nvp
         id = i - joff
         call MPI_RECV(part,ndt,mreal,id,105,lworld,istatus,ierr)
c determine how many particles to write
         call MPI_GET_COUNT(istatus,mreal,nd,ierr)
         nd = nd/idimp
         write (iunit) nd
         write (iunit) ((part(j,k,l),j=1,idimp),k=1,nd)
   10    continue
c read data back for node 0
         rewind 99
         read (99) ((part(j,k,l),j=1,idimp),k=1,nd0)
         close (99)
c other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         ndt = idimp*npp(l)
         call MPI_SEND(part,ndt,mreal,0,105,lworld,ierr)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PRDPART(part,npp,idimp,npmax,nblok,iunit,ierror)
c this subroutine reads particle data from a fortran unformatted
c sequential file and distributes it
c part = particle data to be read
c npp(l) = number of particles in partition l
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c iunit = fortran unit number
c ierror = (0,1) = (no,yes) error condition exists
c input: all
c output: part
      implicit none
      integer idimp, npmax, nblok, iunit, ierror
      real part
      integer npp
      dimension part(idimp,npmax,nblok), npp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, j, k, l, ndt, idproc, lvp, nvp, joff, id, nd, nd0
      integer ios, ierr
      integer istatus
      dimension istatus(lstat)
      ierror = 0
      ndt = idimp*npmax
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nproc
c     do 10 l = 1, nblok
c     read (iunit,iostat=ios) npp(l)
c     read (iunit,iostat=ios) (part(:,j,l),j=1,nd)
c  10 continue
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,lvp,ierr)
      do 20 l = 1, nblok
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
c first read data for node 0
         read (iunit,iostat=ios) nd
         if (ios.ne.0) ierror = ierror + 1
         read (iunit,iostat=ios) ((part(j,k,l),j=1,idimp),k=1,nd)
         if (ios.ne.0) ierror = ierror + 1
c no special diagnostic node
         if (joff.eq.1) then
            npp(l) = nd
c special diagnostic node present
         else
            nvp = lvp - nproc
            id = 1
            ndt = idimp*nd
            call MPI_SEND(part,ndt,mreal,id,104,lworld,ierr)
         endif
c save copy of current data
         nd0 = nd
c         open(unit=99,form='unformatted',status='scratch')
         open(unit=99,form='unformatted')
         write (99) ((part(j,k,l),j=1,idimp),k=1,nd0)
c then send data to remaining nodes
         do 10 i = 2, nvp
         read (iunit,iostat=ios) nd
         if (ios.ne.0) ierror = ierror + 1
         read (iunit,iostat=ios) ((part(j,k,l),j=1,idimp),k=1,nd)
         if (ios.ne.0) ierror = ierror + 1
         ndt = idimp*nd
         id = i - joff
         call MPI_SEND(part,ndt,mreal,id,104,lworld,ierr)
   10    continue
c read data back for node 0
         rewind 99
         read (99,iostat=ios) ((part(j,k,l),j=1,idimp),k=1,nd0)
         if (ios.ne.0) ierror = ierror + 1
         close (99)
c other nodes receive data from node 0
      else if (idproc.le.(nvp-joff)) then
         ndt = idimp*npmax
         call MPI_RECV(part,ndt,mreal,0,104,lworld,istatus,ierr)
c determine how many particles to write
         call MPI_GET_COUNT(istatus,mreal,nd,ierr)
         npp(l) = nd/idimp
      endif
   20 continue
c check for error condition
      call PISUM(ierror,ierr,1,1)
      return
      end
