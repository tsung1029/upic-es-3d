c 3d parallel PIC library for MPI communications
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: february 3, 2003
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
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
      integer*4 id0s, ierror, one, k
      integer ndprec, i
      logical*4 flag
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
      call MPI_COMM_RANK(lworld,id0s,ierror)
      id0 = id0s
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
      one = 1
      call MPI_COMM_SPLIT(lworld,k,one,lgrp,ierror)
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
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
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
      integer*4 nproc, lgrp, mreal, mint, mcplx,  mdouble, lworld
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer*4 id0s, ierror, npv
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,id0s,ierror)
      id0 = id0s
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,npv,ierror)
c special case if diagnostic nodes are present
      if (npv.gt.nproc) then
         idproc = id0 - 1
         if (idproc.ge.nproc) idproc = nproc - idproc - 2
      else
         idproc = id0
      endif
      nvp = nproc
      return
      end
c-----------------------------------------------------------------------
      subroutine PPEXIT
c this subroutine terminates parallel processing
      implicit none
c common block for parallel processing
      integer*4 nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lgrp = current communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      integer ierror
      logical flag
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
c synchronize processes
         call MPI_BARRIER(lgrp,ierror)
c terminate MPI execution environment
         call MPI_FINALIZE(ierror)
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer*4 istatus, npv, idproc, ierror, ns, zero, one, tag
      dimension istatus(lstat)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,npv,ierror)
c quit if communicators are the same size
      if (npv.eq.nproc) return
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierror)
c send message
      ns = n
      zero = 0
      one = 1
      tag = 999
      if (idproc.eq.0) then
         call MPI_RECV(msg,ns,mdouble,one,tag,lworld,istatus,ierror)
      elseif (idproc.eq.1) then
         call MPI_SEND(msg,ns,mdouble,zero,tag,lworld,ierror)
      endif
c broadcast message to non-participating nodes
      if (npv.gt.(2*nproc)) then
         call MPI_BCAST(msg,ns,mdouble,zero,lgrp,ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine DCOMP32(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,idds
     1,mblok,nblok)
c this subroutine determines spatial boundaries for particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c for 2D spatial decomposition
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c where m = my + mblok*(mz - 1)
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idps = number of particle partition boundaries
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
      implicit none
      real edges
      integer nyzp, noff, ny, nz, kstrt, nvpy, nvpz, idps, idds
      integer mblok, nblok
      dimension edges(idps,mblok*nblok)
      dimension nyzp(idds,mblok*nblok), noff(idds,mblok*nblok)
c local data
      integer js, ks, jb, kb, kr, m, moff, my, mz
      real at1, at2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      at1 = float(ny)/float(nvpy)
      at2 = float(nz)/float(nvpz)
      do 20 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb = mz + ks
      do 10 my = 1, mblok
      m = my + moff
      jb = my + js
      edges(1,m) = at1*float(jb)
      noff(1,m) = edges(1,m) + .5
      edges(2,m) = at1*float(jb + 1)
      kr = edges(2,m) + .5
      nyzp(1,m) = kr - noff(1,m)
      edges(3,m) = at2*float(kb)
      noff(2,m) = edges(3,m) + .5
      edges(4,m) = at2*float(kb + 1)
      kr = edges(4,m) + .5
      nyzp(2,m) = kr - noff(2,m)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCOMP32L(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,idd
     1s,mblok,nblok)
c this subroutine determines spatial boundaries for particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c for 2D spatial decomposition
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c where m = my + mblok*(mz - 1)
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idps = number of particle partition boundaries
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
      implicit none
      real edges
      integer nyzp, noff, ny, nz, kstrt, nvpy, nvpz, idps, idds
      integer mblok, nblok
      dimension edges(idps,mblok*nblok)
      dimension nyzp(idds,mblok*nblok), noff(idds,mblok*nblok)
c local data
      integer js, ks, jb, kb, kr, m, moff, my, mz
      real at1, at2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      at1 = float(ny)/float(nvpy)
      at2 = float(nz)/float(nvpz)
      do 20 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb = mz + ks
      do 10 my = 1, mblok
      m = my + moff
      jb = my + js
      edges(1,m) = at1*float(jb)
      noff(1,m) = edges(1,m)
      edges(2,m) = at1*float(jb + 1)
      kr = edges(2,m)
      nyzp(1,m) = kr - noff(1,m)
      edges(3,m) = at2*float(kb)
      noff(2,m) = edges(3,m)
      edges(4,m) = at2*float(kb + 1)
      kr = edges(4,m)
      nyzp(2,m) = kr - noff(2,m)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCGUARD32(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,n
     1blok,kyp,kzp,ngds)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs(j,k,ngds,m) = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer kyp, kzp
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer*4 istatus, msid, ierr, lenc, kls, krs, tag
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      ngc = 2
      if (kyp.eq.1) ngc = 1
      do 30 m = 1, mnblok
      do 20 k = 1, nzpmx
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+1,k,m)
      scs(j,k,2,m) = f(j,2,k,m)
      scs(j,k,3,m) = f(j,ngc+1,k,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 90 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 80 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      krr = kr + kz
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kll = kl + kz
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = kr + 1
         if (krr.ge.nvpy) krr = krr - nvpy
         krr = krr + kz
         kll = kl - 1
         if (kll.lt.0) kll = kll + nvpy
         kll = kll + kz
         ngc = 1
      endif
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 50 k = 1, nzpmx
c     do 40 j = 1, nxv
c     scs(j,k,4,m) = scs(j,k,1,kl)
c     scs(j,k,5,m) = scs(j,k,2,kr)
c     scs(j,k,6,m) = scs(j,k,3,krr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      lenc = nxvz
      kls = kl - 1
      krs = kr - 1
      tag = noff + 3
      call MPI_IRECV(scs(1,1,4,m),lenc,mreal,kls,tag,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),lenc,mreal,krs,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      lenc = ngc*nxvz
      tag = noff + 4
      call MPI_IRECV(scs(1,1,5,m),lenc,mreal,krs,tag,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,2,m),lenc,mreal,kls,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         kls = kll - 1
         krs = krr - 1
         tag = noff + 6
         call MPI_IRECV(scs(1,1,6,m),lenc,mreal,krs,tag,lgrp,msid,ierr)
         call MPI_SEND(scs(1,1,3,m),lenc,mreal,kls,tag,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      f(j,1,k,m) = scs(j,k,4,m)
      f(j,kyp+2,k,m) = scs(j,k,5,m)
      f(j,kyp+3,k,m) = scs(j,k,6,m)
   60 continue
   70 continue
   80 continue
   90 continue
c copy to guard cells in z
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = ky + nvpy*kl
      ngc = 2
c special case of only one grid per processor
      if (kzp.eq.1) then
         krr = kr + 1
         if (krr.ge.nvpz) krr = krr - nvpz
         krr = ky + nvpy*krr
         kll = kl - 1
         if (kll.lt.0) kll = kll + nvpz
         kll = ky + nvpy*kll
         ngc = 1
      endif
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 110 k = 1, nypmx
c     do 100 j = 1, nxv
c     f(j,k,1,m) = f(j,k,kzp+1,kl)
c     f(j,k,kzp+2,m) = f(j,k,2,kr)
c     f(j,k,kzp+3,m) = f(j,k,ngc+1,krr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
      lenc = nxvy
      kls = kl - 1
      krs = kr - 1
      tag = noff + 7
      call MPI_IRECV(f(1,1,1,m),lenc,mreal,kls,tag,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+1,m),lenc,mreal,krs,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      lenc = ngc*nxvy
      tag = noff + 8
      call MPI_IRECV(f(1,1,kzp+2,m),lenc,mreal,krs,tag,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,2,m),lenc,mreal,kls,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         kls = kll - 1
         krs = krr - 1
         tag = noff + 10
         call MPI_IRECV(f(1,1,kzp+3,m),lenc,mreal,krs,tag,lgrp,msid,ierr
     1)
         call MPI_SEND(f(1,1,2,m),lenc,mreal,kls,tag,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,
     1nblok,kyp,kzp,ngds)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs(j,k,m) = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer kyp, kzp
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer*4 istatus, msid, ierr, lenc, kls, krs, tag
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 30 m = 1, mnblok
      do 20 k = 1, nzpmx
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,1,k,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 90 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 80 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 50 k = 1, nzpmx
c     do 40 j = 1, nxv
c     scs(j,k,2,m) = scs(j,k,1,kr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      lenc = nxvz
      kls = kl - 1
      krs = kr - 1
      tag = noff + 3
      call MPI_IRECV(scs(1,1,2,m),lenc,mreal,krs,tag,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),lenc,mreal,kls,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      f(j,kyp+1,k,m) = scs(j,k,2,m)
   60 continue
   70 continue
   80 continue
   90 continue
c copy to guard cells in z
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 110 k = 1, nypmx
c     do 100 j = 1, nxv
c     f(j,k,kzp+1,m) = f(j,k,1,kr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
      lenc = nxvy
      kls = kl - 1
      krs = kr - 1
      tag = noff + 4
      call MPI_IRECV(f(1,1,kzp+1,m),lenc,mreal,krs,tag,lgrp,msid,ier
     1r)
      call MPI_SEND(f(1,1,1,m),lenc,mreal,kls,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mbl
     1ok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer kyp, kzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer*4 istatus, msid, ierr, lenc, kls, krs, tag
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      if (kyp.eq.1) ngc = 1
      do 30 m = 1, mnblok
      do 20 k = 1, nzpmx
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+2,k,m)
      scs(j,k,2,m) = f(j,kyp+3,k,m)
      scs(j,k,3,m) = f(j,1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 90 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 80 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      krr = kr + kz
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kll = kl + kz
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = kr + 1
         if (krr.ge.nvpy) krr = krr - nvpy
         krr = krr + kz
         kll = kl - 1
         if (kll.lt.0) kll = kll + nvpy
         kll = kll + kz
         ngc = 1
      endif
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 50 k = 1, nzpmx
c     do 40 j = 1, nxv
c     scs(j,k,4,m) = scs(j,k,1,kl)
c     scs(j,k,5,m) = scs(j,k,2,kll)
c     scs(j,k,6,m) = scs(j,k,3,kr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      lenc = ngc*nxvz
      kls = kl - 1
      krs = kr - 1
      tag = noff + 1
      call MPI_IRECV(scs(1,1,4,m),lenc,mreal,kls,tag,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),lenc,mreal,krs,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      lenc = nxvz
      kls = kl - 1
      krs = kr - 1
      tag = noff + 2
      call MPI_IRECV(scs(1,1,6,m),lenc,mreal,krs,tag,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,3,m),lenc,mreal,kls,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         lenc = ngc*nxvz
         kls = kll - 1
         krs = krr - 1
         tag = noff + 5
         call MPI_IRECV(scs(1,1,5,m),lenc,mreal,kls,tag,lgrp,msid,ierr)
         call MPI_SEND(scs(1,1,2,m),lenc,mreal,krs,tag,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,4,m)
      f(j,ngc+1,k,m) = f(j,ngc+1,k,m) + scs(j,k,5,m)
      f(j,kyp+1,k,m) = f(j,kyp+1,k,m) + scs(j,k,6,m)
   60 continue
   70 continue
   80 continue
   90 continue
c add guard cells in z
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = ky + nvpy*kl
      ngc = 2
c special case of only one grid per processor
      if (kzp.eq.1) then
         krr = kr + 1
         if (krr.ge.nvpz) krr = krr - nvpz
         krr = ky + nvpy*krr
         kll = kl - 1
         if (kll.lt.0) kll = kll + nvpz
         kll = ky + nvpy*kll
         ngc = 1
      endif
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 110 k = 1, nypmx
c     do 100 j = 1, nxv
c     scr(j,k,1,m) = f(j,k,kzp+2,kl)
c     scr(j,k,2,m) = f(j,k,kzp+3,kll)
c     scr(j,k,3,m) = f(j,k,1,kr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
      lenc = ngc*nxvy
      kls = kl - 1
      krs = kr - 1
      tag = noff + 3
      call MPI_IRECV(scr,lenc,mreal,kls,tag,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+2,m),lenc,mreal,krs,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      lenc = nxvy
      tag = noff + 4
      call MPI_IRECV(scr(1,1,3,m),lenc,mreal,krs,tag,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,m),lenc,mreal,kls,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         lenc = ngc*nxvy
         kls = kll - 1
         krs = krr - 1
         tag = noff + 6
         call MPI_IRECV(scr(1,1,2,m),lenc,mreal,kls,tag,lgrp,msid,ierr)
         call MPI_SEND(f(1,1,kzp+3,m),lenc,mreal,krs,tag,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 130 k = 1, nypmx
      do 120 j = 1, nxv
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,1,m)
      f(j,k,ngc+1,m) = f(j,k,ngc+1,m) + scr(j,k,2,m)
      f(j,k,kzp+1,m) = f(j,k,kzp+1,m) + scr(j,k,3,m)
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mb
     1lok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer kyp, kzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer*4 istatus, msid, ierr, lenc, kls, krs, tag
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 30 m = 1, mnblok
      do 20 k = 1, nzpmx
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 90 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 80 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 50 k = 1, nzpmx
c     do 40 j = 1, nxv
c     scs(j,k,2,m) = scs(j,k,1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      lenc = nxvz
      kls = kl - 1
      krs = kr - 1
      tag = noff + 1
      call MPI_IRECV(scs(1,1,2,m),lenc,mreal,kls,tag,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),lenc,mreal,krs,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   60 continue
   70 continue
   80 continue
   90 continue
c add guard cells in z
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 110 k = 1, nypmx
c     do 100 j = 1, nxv
c     scr(j,k,m) = f(j,k,kzp+1,kl)
c 100 continue
c 110 continue
c this segment is used for mpi computers
      lenc = nxvy
      kls = kl - 1
      krs = kr - 1
      tag = noff + 2
      call MPI_IRECV(scr,lenc,mreal,kls,tag,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+1,m),lenc,mreal,krs,tag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 130 k = 1, nypmx
      do 120 j = 1, nxv
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax
     2,idds,ntmax,ierr)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions with 2D spatial decomposition
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c npp(m) = number of particles in partition m
c sbufl = buffer for particles being sent to back processor
c sbufr = buffer for particles being sent to front processor
c rbufl = buffer for particles being received from back processor
c rbufr = buffer for particles being received from front processor
c ihole = location of holes left in particle arrays
c jsl(idds,m) = number of particles going back in particle partition m
c jsr(idds,m) = number of particles going front in particle partition m
c jss(idds,m) = scratch array for particle partition m
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of particle partition boundaries
c nbmax =  size of buffers for passing particles between processors
c idds = dimensionality of domain decomposition
c ntmax =  size of hole array for particles leaving processors
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, ierr
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer js, ks, mnblok, n, m, my, mz, moff, nvp
      integer iter, npr, nps, ibflg, iwork, kb, kl, kr, j, j1, j2
      integer nbsize, nter
      integer*4 istatus, msid, ierror, lenc, kls, krs, tag1, tag2
      real an, xt
      dimension kb(2), msid(4), istatus(lstat), ibflg(3), iwork(3)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      iter = 2
      nter = 0
c debugging section: count total number of particles before move
      npr = 0
      do 10 m = 1, mnblok
      npr = npr + npp(m)
   10 continue
c buffer outgoing particles, first in y then in z direction
      do 300 n = 1, 2
      if (n.eq.1) then
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         nvp = nvpz
         an = float(nz)
      endif
   20 do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 50 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 30 j = 1, npp(m)
c particles going up or forward
      xt = part(n+1,j,m)
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            part(n+1,j,m) = xt
            sbufr(1,jsr(1,m),m) = part(1,j,m)
            sbufr(2,jsr(1,m),m) = part(2,j,m)
            sbufr(3,jsr(1,m),m) = part(3,j,m)
            sbufr(4,jsr(1,m),m) = part(4,j,m)
            sbufr(5,jsr(1,m),m) = part(5,j,m)
            sbufr(6,jsr(1,m),m) = part(6,j,m)
            ihole(jsl(1,m)+jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 40
         endif
c particles going down or backward
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            part(n+1,j,m) = xt
            sbufl(1,jsl(1,m),m) = part(1,j,m)
            sbufl(2,jsl(1,m),m) = part(2,j,m)
            sbufl(3,jsl(1,m),m) = part(3,j,m)
            sbufl(4,jsl(1,m),m) = part(4,j,m)
            sbufl(5,jsl(1,m),m) = part(5,j,m)
            sbufl(6,jsl(1,m),m) = part(6,j,m)
            ihole(jsl(1,m)+jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 40
         endif
      endif
   30 continue
   40 jss(1,m) = jsl(1,m) + jsr(1,m)
   50 continue
   60 continue
c check for full buffer condition
      nps = 0
      do 100 m = 1, mnblok
      nps = nps + jss(2,m)
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
c get particles from below and above or back and front
      kl = kb(n)
      kb(n) = kl + 1
      if (kb(n).ge.nvp) kb(n) = kb(n) - nvp
      kr = kb(1) + nvpy*kb(2) + 1
      kb(n) = kl - 1
      if (kb(n).lt.0) kb(n) = kb(n) + nvp
      kl = kb(1) + nvpy*kb(2) + 1
c this segment is used for shared memory computers
c     jsl(2,m) = jsr(1,kl)
c     do 120 j = 1, jsl(2,m)
c     rbufl(1,j,m) = sbufr(1,j,kl)
c     rbufl(2,j,m) = sbufr(2,j,kl)
c     rbufl(3,j,m) = sbufr(3,j,kl)
c     rbufl(4,j,m) = sbufr(4,j,kl)
c     rbufl(5,j,m) = sbufr(5,j,kl)
c     rbufl(6,j,m) = sbufr(6,j,kl)
c 120 continue
c     jsr(2,m) = jsl(1,kr)
c     do 130 j = 1, jsr(2,m)
c     rbufr(1,j,m) = sbufl(1,j,kr)
c     rbufr(2,j,m) = sbufl(2,j,kr)
c     rbufr(3,j,m) = sbufl(3,j,kr)
c     rbufr(4,j,m) = sbufl(4,j,kr)
c     rbufr(5,j,m) = sbufl(5,j,kr)
c     rbufr(6,j,m) = sbufl(6,j,kr)
c 130 continue
c this segment is used for mpi computers
c post receive
      lenc = nbsize
      kls = kl - 1
      krs = kr - 1
      tag1 = iter - 1
      tag2 = iter
      call MPI_IRECV(rbufl,lenc,mreal,kls,tag1,lgrp,msid(1),ierror)
      call MPI_IRECV(rbufr,lenc,mreal,krs,tag2,lgrp,msid(2),ierror)
c send particles
      lenc = idimp*jsr(1,m)
      call MPI_ISEND(sbufr,lenc,mreal,krs,tag1,lgrp,msid(3),ierror)
      lenc = idimp*jsl(1,m)
      call MPI_ISEND(sbufl,lenc,mreal,kls,tag2,lgrp,msid(4),ierror)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierror)
      call MPI_GET_COUNT(istatus,mreal,lenc,ierror)
      jsl(2,m) = lenc/idimp
      call MPI_WAIT(msid(2),istatus,ierror)
      call MPI_GET_COUNT(istatus,mreal,lenc,ierror)
      jsr(2,m) = lenc/idimp
  140 continue
  150 continue
c check if particles must be passed further
      nps = 0
      do 180 m = 1, mnblok
c check if any particles coming from above or front belong here
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 160 j = 1, jsr(2,m)
      if (rbufr(n+1,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(n+1,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
  160 continue
      if (jsr(1,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jsr(1,m),' particles returning above'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jsr(1,m),' particles returning front'
         endif
      endif
c check if any particles coming from below or back belong here
      do 170 j = 1, jsl(2,m)
      if (rbufl(n+1,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(n+1,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  170 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = nps + (jsl(1,m) + jsr(1,m))
  180 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierror)
      call MPI_WAIT(msid(4),istatus,ierror)
      if (nps.eq.0) go to 240
c remove particles which do not belong here
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 220 my = 1, mblok
      m = my + moff
      kb(1) = my + js
c first check particles coming from above or front
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 190 j = 1, jsr(2,m)
      xt = rbufr(n+1,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(n+1,j,m) = xt
         sbufl(1,jsl(1,m),m) = rbufr(1,j,m)
         sbufl(2,jsl(1,m),m) = rbufr(2,j,m)
         sbufl(3,jsl(1,m),m) = rbufr(3,j,m)
         sbufl(4,jsl(1,m),m) = rbufr(4,j,m)
         sbufl(5,jsl(1,m),m) = rbufr(5,j,m)
         sbufl(6,jsl(1,m),m) = rbufr(6,j,m)
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(n+1,j,m) = xt
         sbufr(1,jsr(1,m),m) = rbufr(1,j,m)
         sbufr(2,jsr(1,m),m) = rbufr(2,j,m)
         sbufr(3,jsr(1,m),m) = rbufr(3,j,m)
         sbufr(4,jsr(1,m),m) = rbufr(4,j,m)
         sbufr(5,jsr(1,m),m) = rbufr(5,j,m)
         sbufr(6,jsr(1,m),m) = rbufr(6,j,m)
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         rbufr(1,jss(2,m),m) = rbufr(1,j,m)
         rbufr(2,jss(2,m),m) = rbufr(2,j,m)
         rbufr(3,jss(2,m),m) = rbufr(3,j,m)
         rbufr(4,jss(2,m),m) = rbufr(4,j,m)
         rbufr(5,jss(2,m),m) = rbufr(5,j,m)
         rbufr(6,jss(2,m),m) = rbufr(6,j,m)
      endif
  190 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 200 j = 1, jsl(2,m)
      xt = rbufl(n+1,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(n+1,j,m) = xt
            sbufr(1,jsr(1,m),m) = rbufl(1,j,m)
            sbufr(2,jsr(1,m),m) = rbufl(2,j,m)
            sbufr(3,jsr(1,m),m) = rbufl(3,j,m)
            sbufr(4,jsr(1,m),m) = rbufl(4,j,m)
            sbufr(5,jsr(1,m),m) = rbufl(5,j,m)
            sbufr(6,jsr(1,m),m) = rbufl(6,j,m)
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(n+1,j,m) = xt
            sbufl(1,jsl(1,m),m) = rbufl(1,j,m)
            sbufl(2,jsl(1,m),m) = rbufl(2,j,m)
            sbufl(3,jsl(1,m),m) = rbufl(3,j,m)
            sbufl(4,jsl(1,m),m) = rbufl(4,j,m)
            sbufl(5,jsl(1,m),m) = rbufl(5,j,m)
            sbufl(6,jsl(1,m),m) = rbufl(6,j,m)
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         rbufl(1,jss(2,m),m) = rbufl(1,j,m)
         rbufl(2,jss(2,m),m) = rbufl(2,j,m)
         rbufl(3,jss(2,m),m) = rbufl(3,j,m)
         rbufl(4,jss(2,m),m) = rbufl(4,j,m)
         rbufl(5,jss(2,m),m) = rbufl(5,j,m)
         rbufl(6,jss(2,m),m) = rbufl(6,j,m)
      endif
  200 continue
  210 jsl(2,m) = jss(2,m)
  220 continue
  230 continue
c check if move would overflow particle array
  240 nps = 0
      do 250 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m) - npmax
      if (jss(2,m).le.0) jss(2,m) = 0
      nps = nps + jss(2,m)
  250 continue
      ibflg(1) = nps
      call PISUM(ibflg,iwork,3,1)
      ierr = ibflg(1)
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         return
      endif
      do 290 m = 1, mnblok
c distribute incoming particles from buffers
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 260 j = 1, jss(2,m)
      part(1,ihole(j,m),m) = rbufl(1,j,m)
      part(2,ihole(j,m),m) = rbufl(2,j,m)
      part(3,ihole(j,m),m) = rbufl(3,j,m)
      part(4,ihole(j,m),m) = rbufl(4,j,m)
      part(5,ihole(j,m),m) = rbufl(5,j,m)
      part(6,ihole(j,m),m) = rbufl(6,j,m)
  260 continue
      if (jss(1,m).gt.jsl(2,m)) then
         jss(2,m) = min0(jss(1,m)-jsl(2,m),jsr(2,m))
      else
         jss(2,m) = jsl(2,m) - jss(1,m)
      endif
      do 270 j = 1, jss(2,m)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1,m).gt.jsl(2,m)) then
         part(1,ihole(j+jsl(2,m),m),m) = rbufr(1,j,m)
         part(2,ihole(j+jsl(2,m),m),m) = rbufr(2,j,m)
         part(3,ihole(j+jsl(2,m),m),m) = rbufr(3,j,m)
         part(4,ihole(j+jsl(2,m),m),m) = rbufr(4,j,m)
         part(5,ihole(j+jsl(2,m),m),m) = rbufr(5,j,m)
         part(6,ihole(j+jsl(2,m),m),m) = rbufr(6,j,m)
      else
c no more holes
c distribute remaining particles from below or back into bottom
         part(1,j+npp(m),m) = rbufl(1,j+jss(1,m),m)
         part(2,j+npp(m),m) = rbufl(2,j+jss(1,m),m)
         part(3,j+npp(m),m) = rbufl(3,j+jss(1,m),m)
         part(4,j+npp(m),m) = rbufl(4,j+jss(1,m),m)
         part(5,j+npp(m),m) = rbufl(5,j+jss(1,m),m)
         part(6,j+npp(m),m) = rbufl(6,j+jss(1,m),m)
      endif
  270 continue
      if (jss(1,m).le.jsl(2,m)) then
         npp(m) = npp(m) + (jsl(2,m) - jss(1,m))
         jss(1,m) = jsl(2,m)
      endif
      jss(2,m) = jss(1,m) - (jsl(2,m) + jsr(2,m))
      if (jss(2,m).gt.0) then
         jss(1,m) = (jsl(2,m) + jsr(2,m))
         jsr(2,m) = jss(2,m)
      else
         jss(1,m) = jss(1,m) - jsl(2,m)
         jsr(2,m) = -jss(2,m)
      endif
      do 280 j = 1, jsr(2,m)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,m).gt.0) then
         j1 = npp(m) - j + 1
         j2 = jss(1,m) + jss(2,m) - j + 1
         if (j1.gt.ihole(j2,m)) then
c move particle only if it is below current hole
            part(1,ihole(j2,m),m) = part(1,j1,m)
            part(2,ihole(j2,m),m) = part(2,j1,m)
            part(3,ihole(j2,m),m) = part(3,j1,m)
            part(4,ihole(j2,m),m) = part(4,j1,m)
            part(5,ihole(j2,m),m) = part(5,j1,m)
            part(6,ihole(j2,m),m) = part(6,j1,m)
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         part(1,j+npp(m),m) = rbufr(1,j+jss(1,m),m)
         part(2,j+npp(m),m) = rbufr(2,j+jss(1,m),m)
         part(3,j+npp(m),m) = rbufr(3,j+jss(1,m),m)
         part(4,j+npp(m),m) = rbufr(4,j+jss(1,m),m)
         part(5,j+npp(m),m) = rbufr(5,j+jss(1,m),m)
         part(6,j+npp(m),m) = rbufr(6,j+jss(1,m),m)
      endif
  280 continue
      if (jss(2,m).gt.0) then
         npp(m) = npp(m) - jsr(2,m)
      else
         npp(m) = npp(m) + jsr(2,m)
      endif
      jss(1,m) = 0
  290 continue
c check if any particles have to be passed further
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         go to 110
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         go to 20
      endif
  300 continue
c debugging section: count total number of particles after move
      nps = 0
      do 310 m = 1, mnblok
      nps = nps + npp(m)
  310 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         ierr = 1
      endif
c information
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax
     2,idds,ntmax,maskp,ierr)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions with 2D spatial decomposition
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c npp(m) = number of particles in partition m
c sbufl = buffer for particles being sent to back processor
c sbufr = buffer for particles being sent to front processor
c rbufl = buffer for particles being received from back processor
c rbufr = buffer for particles being received from front processor
c ihole = location of holes left in particle arrays
c jsl(idds,m) = number of particles going back in particle partition m
c jsr(idds,m) = number of particles going front in particle partition m
c jss(idds,m) = scratch array for particle partition m
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of particle partition boundaries
c nbmax =  size of buffers for passing particles between processors
c idds = dimensionality of domain decomposition
c ntmax =  size of hole array for particles leaving processors
c maskp = scratch array for particle addresses
c ierr = (0,1) = (no,yes) error condition exists
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, maskp, ierr
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax
      dimension part(idimp,npmax,mblok*nblok), maskp(npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer js, ks, mnblok, n, m, my, mz, moff, nvp
      integer iter, npr, nps, ibflg, iwork, kb, kl, kr, j, j1, j2
      integer nbsize, nter
      integer*4 istatus, msid, ierror, lenc, kls, krs, tag1, tag2
      real an, xt
      dimension kb(2), msid(4), istatus(lstat), ibflg(3), iwork(3)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      iter = 2
      nter = 0
c debugging section: count total number of particles before move
      npr = 0
      do 10 m = 1, mnblok
      npr = npr + npp(m)
   10 continue
c buffer outgoing particles, first in y then in z direction
      do 300 n = 1, 2
      if (n.eq.1) then
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         nvp = nvpz
         an = float(nz)
      endif
   20 do 90 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 80 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jss(1,m) = 0
      jss(2,m) = 0
c find mask function for particles out of bounds
      do 30 j = 1, npp(m)
      xt = part(n+1,j,m)
      if ((xt.ge.edges(2*n,m)).or.(xt.lt.edges(2*n-1,m))) then
         jss(1,m) = jss(1,m) + 1
         maskp(j,m) = 1
      else
         maskp(j,m) = 0
      endif
   30 continue
c set flag if hole buffer would overflow
      if (jss(1,m).gt.ntmax) then
         jss(1,m) = ntmax
         jss(2,m) = 1
      endif
c accumulate location of holes
      do 40 j = 2, npp(m)
      maskp(j,m) = maskp(j,m) + maskp(j-1,m)
   40 continue
c store addresses of particles out of bounds
      do 50 j = 2, npp(m)
      if ((maskp(j,m).gt.maskp(j-1,m)).and.(maskp(j,m).le.ntmax)) then
         ihole(maskp(j,m),m) = j
      endif
   50 continue
      if (maskp(1,m).gt.0) ihole(1,m) = 1
      jsl(1,m) = 0
      jsr(1,m) = 0
c load particle buffers
      do 60 j = 1, jss(1,m)
      xt = part(n+1,ihole(j,m),m)
c particles going forward
      if (xt.ge.edges(2*n,m)) then
         if ((kb(n)+1).eq.nvp) xt = xt - an
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            part(n+1,ihole(j,m),m) = xt
            sbufr(1,jsr(1,m),m) = part(1,ihole(j,m),m)
            sbufr(2,jsr(1,m),m) = part(2,ihole(j,m),m)
            sbufr(3,jsr(1,m),m) = part(3,ihole(j,m),m)
            sbufr(4,jsr(1,m),m) = part(4,ihole(j,m),m)
            sbufr(5,jsr(1,m),m) = part(5,ihole(j,m),m)
            sbufr(6,jsr(1,m),m) = part(6,ihole(j,m),m)
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
         else
            jss(2,m) = 1
c           go to 70
         endif
c particles going backward
      else
         if (kb(n).eq.0) xt = xt + an
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            part(n+1,ihole(j,m),m) = xt
            sbufl(1,jsl(1,m),m) = part(1,ihole(j,m),m)
            sbufl(2,jsl(1,m),m) = part(2,ihole(j,m),m)
            sbufl(3,jsl(1,m),m) = part(3,ihole(j,m),m)
            sbufl(4,jsl(1,m),m) = part(4,ihole(j,m),m)
            sbufl(5,jsl(1,m),m) = part(5,ihole(j,m),m)
            sbufl(6,jsl(1,m),m) = part(6,ihole(j,m),m)
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
         else
            jss(2,m) = 1
c           go to 70
         endif
      endif
   60 continue
   70 jss(1,m) = jsl(1,m) + jsr(1,m)
   80 continue
   90 continue
c check for full buffer condition
      nps = 0
      do 100 m = 1, mnblok
      nps = nps + jss(2,m)
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
c get particles from below and above or back and front
      kl = kb(n)
      kb(n) = kl + 1
      if (kb(n).ge.nvp) kb(n) = kb(n) - nvp
      kr = kb(1) + nvpy*kb(2) + 1
      kb(n) = kl - 1
      if (kb(n).lt.0) kb(n) = kb(n) + nvp
      kl = kb(1) + nvpy*kb(2) + 1
c this segment is used for shared memory computers
c     jsl(2,m) = jsr(1,kl)
c     do 120 j = 1, jsl(2,m)
c     rbufl(1,j,m) = sbufr(1,j,kl)
c     rbufl(2,j,m) = sbufr(2,j,kl)
c     rbufl(3,j,m) = sbufr(3,j,kl)
c     rbufl(4,j,m) = sbufr(4,j,kl)
c     rbufl(5,j,m) = sbufr(5,j,kl)
c     rbufl(6,j,m) = sbufr(6,j,kl)
c 120 continue
c     jsr(2,m) = jsl(1,kr)
c     do 130 j = 1, jsr(2,m)
c     rbufr(1,j,m) = sbufl(1,j,kr)
c     rbufr(2,j,m) = sbufl(2,j,kr)
c     rbufr(3,j,m) = sbufl(3,j,kr)
c     rbufr(4,j,m) = sbufl(4,j,kr)
c     rbufr(5,j,m) = sbufl(5,j,kr)
c     rbufr(6,j,m) = sbufl(6,j,kr)
c 130 continue
c this segment is used for mpi computers
c post receive
      lenc = nbsize
      kls = kl - 1
      krs = kr - 1
      tag1 = iter - 1
      tag2 = iter
      call MPI_IRECV(rbufl,lenc,mreal,kls,tag1,lgrp,msid(1),ierror)
      call MPI_IRECV(rbufr,lenc,mreal,krs,tag2,lgrp,msid(2),ierror)
c send particles
      lenc = idimp*jsr(1,m)
      call MPI_ISEND(sbufr,lenc,mreal,krs,tag1,lgrp,msid(3),ierror)
      lenc = idimp*jsl(1,m)
      call MPI_ISEND(sbufl,lenc,mreal,kls,tag2,lgrp,msid(4),ierror)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierror)
      call MPI_GET_COUNT(istatus,mreal,lenc,ierror)
      jsl(2,m) = lenc/idimp
      call MPI_WAIT(msid(2),istatus,ierror)
      call MPI_GET_COUNT(istatus,mreal,lenc,ierror)
      jsr(2,m) = lenc/idimp
  140 continue
  150 continue
c check if particles must be passed further
      nps = 0
      do 180 m = 1, mnblok
c check if any particles coming from above or front belong here
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 160 j = 1, jsr(2,m)
      if (rbufr(n+1,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(n+1,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
  160 continue
      if (jsr(1,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jsr(1,m),' particles returning above'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jsr(1,m),' particles returning front'
         endif
      endif
c check if any particles coming from below or back belong here
      do 170 j = 1, jsl(2,m)
      if (rbufl(n+1,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(n+1,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  170 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = nps + (jsl(1,m) + jsr(1,m))
  180 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierror)
      call MPI_WAIT(msid(4),istatus,ierror)
      if (nps.eq.0) go to 240
c remove particles which do not belong here
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 220 my = 1, mblok
      m = my + moff
      kb(1) = my + js
c first check particles coming from above or front
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 190 j = 1, jsr(2,m)
      xt = rbufr(n+1,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(n+1,j,m) = xt
         sbufl(1,jsl(1,m),m) = rbufr(1,j,m)
         sbufl(2,jsl(1,m),m) = rbufr(2,j,m)
         sbufl(3,jsl(1,m),m) = rbufr(3,j,m)
         sbufl(4,jsl(1,m),m) = rbufr(4,j,m)
         sbufl(5,jsl(1,m),m) = rbufr(5,j,m)
         sbufl(6,jsl(1,m),m) = rbufr(6,j,m)
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(n+1,j,m) = xt
         sbufr(1,jsr(1,m),m) = rbufr(1,j,m)
         sbufr(2,jsr(1,m),m) = rbufr(2,j,m)
         sbufr(3,jsr(1,m),m) = rbufr(3,j,m)
         sbufr(4,jsr(1,m),m) = rbufr(4,j,m)
         sbufr(5,jsr(1,m),m) = rbufr(5,j,m)
         sbufr(6,jsr(1,m),m) = rbufr(6,j,m)
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         rbufr(1,jss(2,m),m) = rbufr(1,j,m)
         rbufr(2,jss(2,m),m) = rbufr(2,j,m)
         rbufr(3,jss(2,m),m) = rbufr(3,j,m)
         rbufr(4,jss(2,m),m) = rbufr(4,j,m)
         rbufr(5,jss(2,m),m) = rbufr(5,j,m)
         rbufr(6,jss(2,m),m) = rbufr(6,j,m)
      endif
  190 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 200 j = 1, jsl(2,m)
      xt = rbufl(n+1,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(n+1,j,m) = xt
            sbufr(1,jsr(1,m),m) = rbufl(1,j,m)
            sbufr(2,jsr(1,m),m) = rbufl(2,j,m)
            sbufr(3,jsr(1,m),m) = rbufl(3,j,m)
            sbufr(4,jsr(1,m),m) = rbufl(4,j,m)
            sbufr(5,jsr(1,m),m) = rbufl(5,j,m)
            sbufr(6,jsr(1,m),m) = rbufl(6,j,m)
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(n+1,j,m) = xt
            sbufl(1,jsl(1,m),m) = rbufl(1,j,m)
            sbufl(2,jsl(1,m),m) = rbufl(2,j,m)
            sbufl(3,jsl(1,m),m) = rbufl(3,j,m)
            sbufl(4,jsl(1,m),m) = rbufl(4,j,m)
            sbufl(5,jsl(1,m),m) = rbufl(5,j,m)
            sbufl(6,jsl(1,m),m) = rbufl(6,j,m)
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         rbufl(1,jss(2,m),m) = rbufl(1,j,m)
         rbufl(2,jss(2,m),m) = rbufl(2,j,m)
         rbufl(3,jss(2,m),m) = rbufl(3,j,m)
         rbufl(4,jss(2,m),m) = rbufl(4,j,m)
         rbufl(5,jss(2,m),m) = rbufl(5,j,m)
         rbufl(6,jss(2,m),m) = rbufl(6,j,m)
      endif
  200 continue
  210 jsl(2,m) = jss(2,m)
  220 continue
  230 continue
c check if move would overflow particle array
  240 nps = 0
      do 250 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m) - npmax
      if (jss(2,m).le.0) jss(2,m) = 0
      nps = nps + jss(2,m)
  250 continue
      ibflg(1) = nps
      call PISUM(ibflg,iwork,3,1)
      ierr = ibflg(1)
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         return
      endif
      do 290 m = 1, mnblok
c distribute incoming particles from buffers
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 260 j = 1, jss(2,m)
      part(1,ihole(j,m),m) = rbufl(1,j,m)
      part(2,ihole(j,m),m) = rbufl(2,j,m)
      part(3,ihole(j,m),m) = rbufl(3,j,m)
      part(4,ihole(j,m),m) = rbufl(4,j,m)
      part(5,ihole(j,m),m) = rbufl(5,j,m)
      part(6,ihole(j,m),m) = rbufl(6,j,m)
  260 continue
      if (jss(1,m).gt.jsl(2,m)) then
         jss(2,m) = min0(jss(1,m)-jsl(2,m),jsr(2,m))
      else
         jss(2,m) = jsl(2,m) - jss(1,m)
      endif
      do 270 j = 1, jss(2,m)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1,m).gt.jsl(2,m)) then
         part(1,ihole(j+jsl(2,m),m),m) = rbufr(1,j,m)
         part(2,ihole(j+jsl(2,m),m),m) = rbufr(2,j,m)
         part(3,ihole(j+jsl(2,m),m),m) = rbufr(3,j,m)
         part(4,ihole(j+jsl(2,m),m),m) = rbufr(4,j,m)
         part(5,ihole(j+jsl(2,m),m),m) = rbufr(5,j,m)
         part(6,ihole(j+jsl(2,m),m),m) = rbufr(6,j,m)
      else
c no more holes
c distribute remaining particles from below or back into bottom
         part(1,j+npp(m),m) = rbufl(1,j+jss(1,m),m)
         part(2,j+npp(m),m) = rbufl(2,j+jss(1,m),m)
         part(3,j+npp(m),m) = rbufl(3,j+jss(1,m),m)
         part(4,j+npp(m),m) = rbufl(4,j+jss(1,m),m)
         part(5,j+npp(m),m) = rbufl(5,j+jss(1,m),m)
         part(6,j+npp(m),m) = rbufl(6,j+jss(1,m),m)
      endif
  270 continue
      if (jss(1,m).le.jsl(2,m)) then
         npp(m) = npp(m) + (jsl(2,m) - jss(1,m))
         jss(1,m) = jsl(2,m)
      endif
      jss(2,m) = jss(1,m) - (jsl(2,m) + jsr(2,m))
      if (jss(2,m).gt.0) then
         jss(1,m) = (jsl(2,m) + jsr(2,m))
         jsr(2,m) = jss(2,m)
      else
         jss(1,m) = jss(1,m) - jsl(2,m)
         jsr(2,m) = -jss(2,m)
      endif
      do 280 j = 1, jsr(2,m)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,m).gt.0) then
         j1 = npp(m) - j + 1
         j2 = jss(1,m) + jss(2,m) - j + 1
         if (j1.gt.ihole(j2,m)) then
c move particle only if it is below current hole
            part(1,ihole(j2,m),m) = part(1,j1,m)
            part(2,ihole(j2,m),m) = part(2,j1,m)
            part(3,ihole(j2,m),m) = part(3,j1,m)
            part(4,ihole(j2,m),m) = part(4,j1,m)
            part(5,ihole(j2,m),m) = part(5,j1,m)
            part(6,ihole(j2,m),m) = part(6,j1,m)
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         part(1,j+npp(m),m) = rbufr(1,j+jss(1,m),m)
         part(2,j+npp(m),m) = rbufr(2,j+jss(1,m),m)
         part(3,j+npp(m),m) = rbufr(3,j+jss(1,m),m)
         part(4,j+npp(m),m) = rbufr(4,j+jss(1,m),m)
         part(5,j+npp(m),m) = rbufr(5,j+jss(1,m),m)
         part(6,j+npp(m),m) = rbufr(6,j+jss(1,m),m)
      endif
  280 continue
      if (jss(2,m).gt.0) then
         npp(m) = npp(m) - jsr(2,m)
      else
         npp(m) = npp(m) + jsr(2,m)
      endif
      jss(1,m) = 0
  290 continue
c check if any particles have to be passed further
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         go to 110
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         go to 20
      endif
  300 continue
c debugging section: count total number of particles after move
      nps = 0
      do 310 m = 1, mnblok
      nps = nps + npp(m)
  310 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         ierr = 1
      endif
c information
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxy
     1pd,kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l,mx,mz) = f(j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g, s, t
      dimension f(nxv,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension s(kxyp,kyp,kzp,kblok*lblok), t(kxyp,kyp,kzp,jblok*lblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1)
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(k+koff,j,l,m) = f(j+joff,k,l,i+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      lenc = kxyp*kyp*kzp
      do 100 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 90 mx = 1, jkblok
      m = mx + moff
      do 80 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = ir + kyb*(mz + ls) + 1
         irs = ir - 1
         tag = ir + kxym + 1
         call MPI_IRECV(t(1,1,1,m),lenc,mcplx,irs,tag,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = is + kxb*(mz + ks) + 1
         do 30 l = 1, kzp
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         s(j,k,l,m) = f(j+joff,k,l,m)
   10    continue
   20    continue
   30    continue
         iss = is - 1
         tag = m + kstrt + kxym
         call MPI_SEND(s(1,1,1,m),lenc,mcplx,iss,tag,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 k = 1, kyp
         do 40 j = 1, kxyp
         g(k+koff,j,l,m) = t(j,k,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kx
     1ypd,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(j+kzp*(l-1),k,l,my,mz) = g(k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g = complex input array
c h = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h, s, t
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension s(kyzp,kxyp,kzp,jblok*lblok)
      dimension t(kyzp,kxyp,kzp,jblok*mblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(l+loff,j,k,m) = g(k+koff,j,l,mx+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      lenc = kxyp*kyzp*kzp
      do 100 my = 1, mlblok
      moff = jblok*(my - 1)
      do 90 mx = 1, jblok
      m = mx + moff
      do 80 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         ir = ir0 + kyzm*(ii - 1)
         loff = kzp*ir
         ir = mx + js + kxb*ir + 1
         irs = ir - 1
         tag = ir + kyzm + 1
         call MPI_IRECV(t(1,1,1,m),lenc,mcplx,irs,tag,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = mx + js + kxb*is + 1
         do 30 l = 1, kzp
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         s(k,j,l,m) = g(k+koff,j,l,m)
   10    continue
   20    continue
   30    continue
         iss = is - 1
         tag = m + kstrt + kyzm
         call MPI_SEND(s(1,1,1,m),lenc,mcplx,iss,tag,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 j = 1, kxyp
         do 40 k = 1, kyzp
         h(l+loff,j,k,m) = t(k,j,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kx
     1ypd,kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l,mx,mz) = f(j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g, s, t
      dimension f(3,nxv,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension s(3,kxyp,kyp,kzp,kblok*lblok)
      dimension t(3,kxyp,kyp,kzp,jblok*lblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1)
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(1,k+koff,j,l,m) = f(1,j+joff,k,l,i+ioff)
c     g(2,k+koff,j,l,m) = f(2,j+joff,k,l,i+ioff)
c     g(3,k+koff,j,l,m) = f(3,j+joff,k,l,i+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      lenc = 3*kxyp*kyp*kzp
      do 100 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 90 mx = 1, jkblok
      m = mx + moff
      do 80 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = ir + kyb*(mz + ls) + 1
         irs = ir - 1
         tag = ir + kxym + 1
         call MPI_IRECV(t(1,1,1,1,m),lenc,mcplx,irs,tag,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = is + kxb*(mz + ks) + 1
         do 30 l = 1, kzp
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         s(1,j,k,l,m) = f(1,j+joff,k,l,m)
         s(2,j,k,l,m) = f(2,j+joff,k,l,m)
         s(3,j,k,l,m) = f(3,j+joff,k,l,m)
   10    continue
   20    continue
   30    continue
         iss = is - 1
         tag = m + kstrt + kxym
         call MPI_SEND(s(1,1,1,1,m),lenc,mcplx,iss,tag,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 k = 1, kyp
         do 40 j = 1, kxyp
         g(1,k+koff,j,l,m) = t(1,j,k,l,m)
         g(2,k+koff,j,l,m) = t(2,j,k,l,m)
         g(3,k+koff,j,l,m) = t(3,j,k,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(j+kzp*(l-1),k,l,my,mz) = g(k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g = complex input array
c h = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h, s, t
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension s(3,kyzp,kxyp,kzp,jblok*lblok)
      dimension t(3,kyzp,kxyp,kzp,jblok*mblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(1,l+loff,j,k,m) = g(1,k+koff,j,l,mx+ioff)
c     h(2,l+loff,j,k,m) = g(2,k+koff,j,l,mx+ioff)
c     h(3,l+loff,j,k,m) = g(3,k+koff,j,l,mx+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      lenc = 3*kxyp*kyzp*kzp
      do 100 my = 1, mlblok
      moff = jblok*(my - 1)
      do 90 mx = 1, jblok
      m = mx + moff
      do 80 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         ir = ir0 + kyzm*(ii - 1)
         loff = kzp*ir
         ir = mx + js + kxb*ir + 1
         irs = ir - 1
         tag = ir + kyzm + 1
         call MPI_IRECV(t(1,1,1,1,m),lenc,mcplx,irs,tag,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = mx + js + kxb*is + 1
         do 30 l = 1, kzp
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         s(1,k,j,l,m) = g(1,k+koff,j,l,m)
         s(2,k,j,l,m) = g(2,k+koff,j,l,m)
         s(3,k,j,l,m) = g(3,k+koff,j,l,m)
   10    continue
   20    continue
   30    continue
         iss = is - 1
         tag = m + kstrt + kyzm
         call MPI_SEND(s(1,1,1,1,m),lenc,mcplx,iss,tag,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 j = 1, kxyp
         do 40 k = 1, kyzp
         h(1,l+loff,j,k,m) = t(1,k,j,l,m)
         h(2,l+loff,j,k,m) = t(2,k,j,l,m)
         h(3,l+loff,j,k,m) = t(3,k,j,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l,mx,mz) = f(j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c f = complex input array
c g = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g
      dimension f(nxv*kypd*kzpd*kblok*lblok)
      dimension g(nyv*kxypd*kzpd*jblok*lblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1) - 1
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1)))) = f(j+joff+nxv*(k-1+kyp
c    1d*(l-1+kzpd*(i+ioff))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      lenc = kxyp*kyp*kzp
c transpose local data
      do 70 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 60 mx = 1, kblok
      m = mx + moff
      ioff = kxb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kxym
      is0 = iand(kxym-1,ieor(mx+js,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kyb*kzb)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kyp*(l + is) - 1
         koff = kypd*(l + loff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         g(j+kxyp*(k+ir)) = f(j+joff+nxv*(k+koff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 100 mx = 1, jkblok
      m = mx + moff
      do 90 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         is = ir0 + kxym*(ii - 1)
         ir = is + kyb*(mz + ls) + 1
         irs = ir - 1
         tag = ir + kxym + 1
         call MPI_IRECV(f(1+kxyp*kyp*kzp*is),lenc,mcplx,irs,tag,lgrp,msi
     1d,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         iss = is - 1
         tag = m + kstrt + kxym
         call MPI_SEND(g(1+kxyp*kyp*kzp*ir),lenc,mcplx,iss,tag,lgrp,ierr
     1)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      do 180 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kyb*(m - 1) 
      loff = kzpd*(m - 1) - 1
      do 160 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kzb)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kyp*(l + ir) - 1
         joff = kxypd*(l + loff) - 1
         do 130 k = 1, kyp
         do 120 j = 1, kxyp
         g(k+koff+nyv*(j+joff)) = f(j+kxyp*(k+is))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd
     1,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(j+kzp*(l-1),k,l,my,mz) = g(k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c g = complex input array
c h = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h
      dimension g(nyv*kxypd*kzpd*jblok*lblok)
      dimension h(nzv*kxypd*kyzpd*jblok*mblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1) - 1
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1)))) = g(k+koff+nyv*(j-1+kx
c    1ypd*(l-1+kzpd*(mx+ioff))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      lenc = kxyp*kyzp*kzp
c transpose local data
      do 70 my = 1, lblok
      moff = jblok*(my - 1)
      do 60 mx = 1, jblok
      m = mx + moff
      ioff = kyzb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kyzm
      is0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kxb*kzb)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kxyp*(l + is) - 1
         joff = kxypd*(l + loff) - 1
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         h(k+kyzp*(j+ir)) = g(k+koff+nyv*(j+joff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 my = 1, mlblok
      moff = jblok*(my - 1)
      do 100 mx = 1, jblok
      m = mx + moff
      do 90 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         is = ir0 + kyzm*(ii - 1)
         ir = mx + js + kxb*is + 1
         irs = ir - 1
         tag = ir + kyzm + 1
         call MPI_IRECV(g(1+kxyp*kyzp*kzp*is),lenc,mcplx,irs,tag,lgrp,ms
     1id,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         iss = is - 1
         tag = m + kstrt + kyzm
         call MPI_SEND(h(1+kxyp*kyzp*kzp*ir),lenc,mcplx,iss,tag,lgrp,ier
     1r)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      nzyv = nzv*kxypd
      do 180 my = 1, mblok
      moff = jblok*(my - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kzb*(m - 1)
      koff = kyzpd*(m - 1) - 1
      do 160 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kyzb)) then
         ir = ir0 + kyzm*(ii - 1)
         joff = kzp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kxyp*(l + ir) - 1
         do 130 j = 1, kxyp
         loff = nzv*(j - 1) + joff
         do 120 k = 1, kyzp
         h(l+loff+nzyv*(k+koff)) = g(k+kyzp*(j+is))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxypd
     1,kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(1:3,k+kyp*(m-1),j,l,mx,mz) = f(1:3,j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c f = complex input array
c g = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
c optimized version
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g
      dimension f(3*nxv*kypd*kzpd*kblok*lblok)
      dimension g(3*nyv*kxypd*kzpd*jblok*lblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1) - 1
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js) - 1
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1) - 1
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(1+3*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(1+3*(j+joff+n
c    1xv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c     g(2+3*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(2+3*(j+joff+n
c    1xv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c     g(3+3*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(3+3*(j+joff+n
c    1xv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      lenc = 3*kxyp*kyp*kzp
c transpose local data
      do 70 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 60 mx = 1, kblok
      m = mx + moff
      ioff = kxb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kxym
      is0 = iand(kxym-1,ieor(mx+js,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kyb*kzb)) then
         is = is0 + kxym*(ii - 1)
         joff = 3*kxyp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kyp*(l + is) - 1
         koff = kypd*(l + loff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, 3*kxyp
         g(j+3*kxyp*(k+ir)) = f(j+joff+3*nxv*(k+koff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 100 mx = 1, jkblok
      m = mx + moff
      do 90 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         is = ir0 + kxym*(ii - 1)
         ir = is + kyb*(mz + ls) + 1
         irs = ir - 1
         tag = ir + kxym + 1
         call MPI_IRECV(f(1+3*kxyp*kyp*kzp*is),lenc,mcplx,irs,tag,lgrp,m
     1sid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         iss = is - 1
         tag = m + kstrt + kxym
         call MPI_SEND(g(1+3*kxyp*kyp*kzp*ir),lenc,mcplx,iss,tag,lgrp,ie
     1rr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      do 180 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kyb*(m - 1) 
      loff = kzpd*(m - 1) - 1
      do 160 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kzb)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kyp*(l + ir) - 1
         joff = kxypd*(l + loff) - 1
         do 130 k = 1, kyp
         do 120 j = 1, kxyp
         g(3*(k+koff+nyv*(j+joff))-2) = f(3*(j+kxyp*(k+is))-2)
         g(3*(k+koff+nyv*(j+joff))-1) = f(3*(j+kxyp*(k+is))-1)
         g(3*(k+koff+nyv*(j+joff))) = f(3*(j+kxyp*(k+is)))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxyp
     1d,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(1:3,j+kzp*(l-1),k,l,my,mz) = g(1:3,k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c g = complex input array
c h = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
c optimized version
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h
      dimension g(3*nyv*kxypd*kzpd*jblok*lblok)
      dimension h(3*nzv*kxypd*kyzpd*jblok*mblok)
c common block for parallel processing
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is
      integer*4 istatus, msid, ierr, lenc, tag, irs, iss
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks) - 1
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1) - 1
c     ioff = jblok*(i - 1) - 1
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(1+3*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(1+3*(k+koff+
c    1nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c     h(2+3*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(2+3*(k+koff+
c    1nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c     h(3+3*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(3+3*(k+koff+
c    1nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      lenc = 3*kxyp*kyzp*kzp
c transpose local data
      do 70 my = 1, lblok
      moff = jblok*(my - 1)
      do 60 mx = 1, jblok
      m = mx + moff
      ioff = kyzb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kyzm
      is0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kxb*kzb)) then
         is = is0 + kyzm*(ii - 1)
         koff = 3*kyzp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kxyp*(l + is) - 1
         joff = kxypd*(l + loff) - 1
         do 20 j = 1, kxyp
         do 10 k = 1, 3*kyzp
         h(k+3*kyzp*(j+ir)) = g(k+koff+3*nyv*(j+joff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 my = 1, mlblok
      moff = jblok*(my - 1)
      do 100 mx = 1, jblok
      m = mx + moff
      do 90 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         is = ir0 + kyzm*(ii - 1)
         ir = mx + js + kxb*is + 1
         irs = ir - 1
         tag = ir + kyzm + 1
         call MPI_IRECV(g(1+3*kxyp*kyzp*kzp*is),lenc,mcplx,irs,tag,lgrp,
     1msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         iss = is - 1
         tag = m + kstrt + kyzm
         call MPI_SEND(h(1+3*kxyp*kyzp*kzp*ir),lenc,mcplx,iss,tag,lgrp,i
     1err)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      nzyv = nzv*kxypd
      do 180 my = 1, mblok
      moff = jblok*(my - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kzb*(m - 1)
      koff = kyzpd*(m - 1) - 1
      do 160 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kyzb)) then
         ir = ir0 + kyzm*(ii - 1)
         joff = kzp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kxyp*(l + ir) - 1
         do 130 j = 1, kxyp
         loff = nzv*(j - 1) + joff
         do 120 k = 1, kyzp
         h(3*(l+loff+nzyv*(k+koff))-2) = g(3*(k+kyzp*(j+is))-2)
         h(3*(l+loff+nzyv*(k+koff))-1) = g(3*(k+kyzp*(j+is))-1)
         h(3*(l+loff+nzyv*(k+koff))) = g(3*(k+kyzp*(j+is)))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
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
      integer*4 nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer*4 one, ierr
      real nclock
      double precision jclock
c initialize clock
      if (icntrl.eq.(-1)) then
         call MPI_BARRIER(lgrp,ierr)
         dtime = MPI_WTIME()
c read clock and write time difference from last clock initialization
      else if (icntrl.eq.1) then
         one = 1
         jclock = dtime
         dtime = MPI_WTIME()
         nclock = real(dtime - jclock)
         call MPI_ALLREDUCE(nclock,time(2),one,mreal,MPI_MIN,lgrp,ierr)
         call MPI_ALLREDUCE(nclock,time(1),one,mreal,MPI_MAX,lgrp,ierr)
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer*4 istatus, msid, idproc, ierr, lenc, kls, krs, tag
      integer kstrt, ks, l, kxs, k, kb, lb, j
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
      lenc = nxp
      kls = kb - kxs - 1
      krs = kb + kxs - 1
      tag = l + nxp
      if (lb.eq.0) then
         call MPI_IRECV(g,lenc,mreal,krs,tag,lgrp,msid,ierr)
         call MPI_SEND(f,lenc,mreal,krs,tag,lgrp,ierr)
      else
         call MPI_IRECV(g,lenc,mreal,kls,tag,lgrp,msid,ierr)
         call MPI_SEND(f,lenc,mreal,kls,tag,lgrp,ierr)
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mint = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer*4 istatus, nsid, idproc, ierr, lenc, kls, krs, tag
      integer kstrt, ks, l, kxs, k, kb, lb, j
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
      lenc = nxp
      kls = kb - kxs - 1
      krs = kb + kxs - 1
      tag = l + nxp
      if (lb.eq.0) then
         call MPI_ISEND(if,lenc,mint,krs,tag,lgrp,nsid,ierr)
         call MPI_RECV(ig,lenc,mint,krs,tag,lgrp,istatus,ierr)
      else
         call MPI_ISEND(if,lenc,mint,kls,tag,lgrp,nsid,ierr)
         call MPI_RECV(ig,lenc,mint,kls,tag,lgrp,istatus,ierr)
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
      subroutine PDSUM(f,g,nxp,nblok)
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
      double precision f, g
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
      subroutine PBCAST(f,nxp)
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ioff, i
      integer*4 istatus, nvp, idproc, np, id, ierr, lenc, tag
      dimension istatus(lstat)
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
      lenc = nxp
      tag = 97
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
            call MPI_SEND(f,lenc,mdouble,id,tag,lworld,ierr)
   10    continue
c then send data to node 0
         if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,lenc,mdouble,id,tag,lworld,ierr)
         endif
c other nodes receive data from node 0
      elseif (idproc.le.(nproc+1)) then
         id = 0
         call MPI_RECV(f,lenc,mdouble,id,tag,lworld,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PWRITE(f,nxp,iunit,nrec,name)
c this subroutine collects distributed real data f and writes to a file
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer lrec, ioff, nrec0, i
      integer*4 istatus, nvp, idproc, np, id, ierr, lenc, tag
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
      lenc = nxp
      tag = 99
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
            call MPI_RECV(f,lenc,mreal,id,tag,lworld,istatus,ierr)
         endif
c first write data for node 0
         nrec0 = nrec
         write (unit=iunit,rec=nrec) f
         nrec = nrec + 1
c then write data from remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_RECV(f,lenc,mreal,id,tag,lworld,istatus,ierr)
            write (unit=iunit,rec=nrec) f
            nrec = nrec + 1
   10    continue
c read data back for node 0
c        read (unit=iunit,rec=nrec0) f
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         id = 0
         call MPI_SEND(f,lenc,mreal,id,tag,lworld,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PREAD(f,nxp,iunit,nrec,name)
c this subroutine reads real data f from a file and distributes it
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
      integer*4 nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer lrec, ioff, nrec0, i
      integer*4 istatus, nvp, idproc, np, id, ierr, lenc, tag
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
      lenc = nxp
      tag = 98
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
            call MPI_SEND(f,lenc,mreal,id,tag,lworld,ierr)
   10    continue
c then read data from node 0
         read (unit=iunit,rec=nrec0) f
         if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,lenc,mreal,id,tag,lworld,ierr)
         endif
c other nodes receive data from node 0
      elseif (idproc.le.(nproc+1)) then
         id = 0
         call MPI_RECV(f,lenc,mreal,id,tag,lworld,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PHEAD3(iunit,nx,ny,nz,nvp,fname)
c this subroutine writes header file for diagnostics
c iunit = fortran unit number
c nx = size of first dimension
c ny = size of second dimension
c nz =size of third dimension
c nvp = size of fourth dimension
c name = file name
c input: all
      implicit none
      integer iunit, nx, ny, nz,  nvp
      character*(*) fname
c common block for parallel processing
      integer*4 nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer*4 idproc, ierr
c this segment is used for shared memory computers
c     idproc = 0
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
      if (idproc.ne.0) return
      open(unit=iunit,file=fname//'head',form='formatted',status='unknow
     1n')
      write (iunit,*) nx, ny, nz, nvp, fname
      close(unit=iunit)
      return
      end
